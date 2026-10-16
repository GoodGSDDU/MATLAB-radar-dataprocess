function verify_results(mapFile)
%VERIFY_RESULTS 结果验证：打印建图结果统计量，可与论文数值逐项对照
%
%   verify_results            % 自动打开最新一次 gicp_map_*.mat
%   verify_results('xxx.mat') % 指定结果文件
%
%   参考值（100 帧原始数据集，默认参数）:
%     地图总点数        3413
%     单帧点数          均值 34.1，中位数 34，最大 48
%     水平覆盖          x 0.44 ~ 18.32 m, y -4.55 ~ 7.17 m
%     点云半径          均值 6.81 m，95% 分位 17.94 m
%     强度范围          5.25e4 ~ 1.33e6
%     配准平移量        均值 0.076 m，最大 0.336 m
%     配准旋转量        均值 1.17 deg，最大 4.24 deg
%
%   不使用任何工具箱函数（quantile/prctile 需 Statistics Toolbox，故用排序法替代）

if nargin < 1 || isempty(mapFile)
    files = dir('gicp_map_*.mat');
    if isempty(files)
        error('verify_results:noFile', '未找到 gicp_map_*.mat，请先运行 run_processing。');
    end
    [~, newest] = max([files.datenum]);
    mapFile = files(newest).name;
end
fprintf('结果文件: %s\n', mapFile);

S = load(mapFile);
mapPoints = S.map_points;
mapIntensity = S.map_intensity(:);
mapFrameIds = S.map_frame_ids(:);
transforms = S.transforms;

%% 1. 点云规模
[frameList, ~, grp] = unique(mapFrameIds);
counts = accumarray(grp, 1);

fprintf('\n──── 点云规模 ────\n');
fprintf('地图总点数                    : %d\n', size(mapPoints, 1));
fprintf('参与建图帧数                  : %d\n', numel(frameList));
fprintf('单帧点数 (均值/中位/最大/最小): %.1f / %g / %d / %d\n', ...
    mean(counts), median(counts), max(counts), min(counts));

%% 2. 空间范围与强度
radius = sqrt(sum(mapPoints(:, 1:2) .^ 2, 2));
sortedR = sort(radius);

fprintf('\n──── 空间分布 ────\n');
fprintf('x 范围                        : %.2f ~ %.2f m\n', min(mapPoints(:,1)), max(mapPoints(:,1)));
fprintf('y 范围                        : %.2f ~ %.2f m\n', min(mapPoints(:,2)), max(mapPoints(:,2)));
fprintf('z 范围                        : %.2f ~ %.2f m\n', min(mapPoints(:,3)), max(mapPoints(:,3)));
fprintf('点云半径 (均值/95%%分位/最大) : %.2f / %.2f / %.2f m\n', ...
    mean(radius), sortedR(max(1, ceil(0.95 * numel(sortedR)))), max(radius));
fprintf('强度范围                      : %.2e ~ %.2e (均值 %.2e)\n', ...
    min(mapIntensity), max(mapIntensity), mean(mapIntensity));

%% 3. 配准位姿（观测时雷达固定，位姿大小即配准误差）
nT = numel(transforms);
trans = nan(nT, 1);
rot = nan(nT, 1);
for f = 1:nT
    T = transforms{f};
    if isempty(T) || ~isequal(size(T), [4 4])
        continue;
    end
    trans(f) = norm(T(1:3, 4));
    rot(f) = acos(min(1, max(-1, (trace(T(1:3, 1:3)) - 1) / 2))) * 180 / pi;
end
valid = ~isnan(trans);

fprintf('\n──── 配准位姿（相对地图坐标系）────\n');
fprintf('平移量 (均值/最大)            : %.3f / %.3f m\n', mean(trans(valid)), max(trans(valid)));
fprintf('旋转量 (均值/最大)            : %.2f / %.2f deg\n', mean(rot(valid)), max(rot(valid)));

%% 4. 场景结构：两把椅子的位置
%   场景为一人在两把椅子之间走动，椅子是观测次数最多的静态物体，
%   在点云中表现为密度最高的两个点簇，可用于核查建图是否正确。
fprintf('\n──── 最高密度点簇（0.75 m 半径邻域计数）────\n');
xy = mapPoints(:, 1:2);
d2 = zeros(size(xy, 1), 1);
for i = 1:size(xy, 1)
    d2(i) = sum((xy(:, 1) - xy(i, 1)) .^ 2 + (xy(:, 2) - xy(i, 2)) .^ 2 < 0.75 ^ 2);
end
[~, densOrder] = sort(d2, 'descend');
shown = 0;
picked = zeros(0, 2);
for k = 1:numel(densOrder)
    idx = densOrder(k);
    c = xy(idx, :);
    if ~isempty(picked) && any(sqrt(sum((picked - c) .^ 2, 2)) <= 1.5)
        continue;                        % 与已输出簇过近，跳过
    end
    near = abs(xy(:,1) - c(1)) < 0.75 & abs(xy(:,2) - c(2)) < 0.75;
    fprintf('  簇 %d: x = %5.2f m, y = %5.2f m | 邻域点数 %4d | 距离 %5.2f m | z = %6.1f ~ %5.1f m\n', ...
        shown + 1, c(1), c(2), d2(idx), norm(c), ...
        min(mapPoints(near, 3)), max(mapPoints(near, 3)));
    picked = [picked; c];                %#ok<AGROW>
    shown = shown + 1;
    if shown >= 3
        break;
    end
end

fprintf('\n验证完成。\n');
end
