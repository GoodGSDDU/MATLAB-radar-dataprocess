function result = process_radar_data(dataFile, opts)
%PROCESS_RADAR_DATA 毫米波雷达数据处理与 G-ICP 多帧配准建图（主处理函数）
%
%   result = process_radar_data(dataFile, opts)
%
%   输入
%     dataFile —— adc_data.bin 路径（TI IWR1843 + DCA1000EVM 原始 ADC 数据）
%     opts     —— 可选参数结构体，见下方默认值
%   输出
%     result   —— 含 map_points / map_intensity / map_frame_ids / transforms /
%                 stats 的结构体，同时保存为 gicp_map_<时间戳>.mat
%
%   处理流程（逐帧）:
%     读取 ADC -> 距离维 FFT -> 多普勒维 FFT -> 方位/俯仰维 FFT
%     -> 二维 CFAR 检测 + 角度峰值校验 -> 极坐标转直角坐标
%     -> 静态点筛选 -> 与累积地图作 G-ICP 配准 -> 并入地图
%
%   相对原始版本的修复（算法本身未做改动，故结果可复现）:
%     1. 补充了缺失的 radar_constants.m 依赖（原 version 调用该函数但目录中没有该文件）。
%     2. 帧读取改为 fseek 直接定位（原实现每帧都从头读入，100 帧累计读取约 8 GB）。
%     3. 修正 I/Q 解交织的循环边界（原实现漏掉每帧最后 4 个 int16）。
%     4. 明确区分 Tc 与 TcGroup(3*Tc)，速度刻度可按物理正确值计算（见 radar_constants.m）。
%     5. 去掉 global 变量与函数内的 clear，参数全部显式传递。
%     6. 局部协方差与距离矩阵改为分块计算，避免地图增长时的内存峰值。
%     7. map_frame_ids 显式存为 double，避免帧数超过 255 时溢出。
%     8. 输出含逐帧统计与位姿，便于结果验证（见 verify_results.m）。
%
%   依赖: radar_constants.m, DCA1000_Read_Data.m（同目录）

if nargin < 1
    error('process_radar_data:notEnoughInputs', ...
        ['本文件是函数文件\n' ...
         '请运行主脚本 run_processing.m（在编辑器中打开它后按 F5），\n' ...
         '或在命令行窗口输入：  run_processing\n\n' ...
         '本函数用法: result = process_radar_data(dataFile, opts)']);
end
if nargin < 2
    opts = struct();
end
if ~isfield(opts, 'maxFrames'),          opts.maxFrames = inf;   end
if ~isfield(opts, 'gicpMaxIters'),       opts.gicpMaxIters = 20; end
if ~isfield(opts, 'gicpTolerance'),      opts.gicpTolerance = 1e-4; end
if ~isfield(opts, 'gicpKNeighbors'),     opts.gicpKNeighbors = 12;  end
if ~isfield(opts, 'gicpMaxCorrDist'),    opts.gicpMaxCorrDist = 0.8; end
if ~isfield(opts, 'pointLimitPerFrame'), opts.pointLimitPerFrame = 1000; end
if ~isfield(opts, 'staticSpeedThresh'),  opts.staticSpeedThresh = 0.2;  end
if ~isfield(opts, 'minRange'),           opts.minRange = 2.5;   end
if ~isfield(opts, 'maxRange'),           opts.maxRange = 20.0;  end
if ~isfield(opts, 'intensityRatio'),     opts.intensityRatio = 0.7; end
if ~isfield(opts, 'partialSavePeriod'),  opts.partialSavePeriod = 25; end
if ~isfield(opts, 'makeFigure'),         opts.makeFigure = true; end
if ~isfield(opts, 'verbose'),            opts.verbose = true;   end

params = radar_constants();

%% 检查数据文件与帧数
info = dir(dataFile);
if isempty(info)
    error('process_radar_data:noFile', '数据文件不存在: %s', dataFile);
end
totalFrames = floor(info.bytes / params.bytesPerFrame);
if totalFrames < 1
    error('process_radar_data:badFile', ...
        '文件大小 %d 字节小于单帧长度 %d 字节。', info.bytes, params.bytesPerFrame);
end
totalFrames = min(totalFrames, opts.maxFrames);

if opts.verbose
    fprintf('=== 毫米波雷达数据处理 (G-ICP 建图) ===\n');
    fprintf('数据文件      : %s\n', dataFile);
    fprintf('文件大小      : %.1f MB\n', info.bytes / 1024 / 1024);
    fprintf('单帧字节数    : %d (%.2f MB)\n', params.bytesPerFrame, params.bytesPerFrame/1024/1024);
    fprintf('待处理帧数    : %d\n', totalFrames);
    fprintf('距离分辨率    : %.4f m (有效带宽 %.1f MHz)\n', params.rangeRes, params.BandwidthValid/1e6);
    fprintf('速度分辨率    : %.4f m/s (PRI = %.1f us)\n', params.speedRes, params.PRI*1e6);
    fprintf('方位/俯仰分辨 : %.1f deg / %.1f deg\n\n', params.angleResAz, params.angleResEl);
end

%% 逐帧处理
map_points = [];
map_intensity = [];
map_frame_ids = [];
transforms = cell(totalFrames, 1);
stats = repmat(struct('frame', 0, 'rawPoints', 0, 'staticPoints', 0, 'mapPoints', 0, ...
                      'gicpIters', 0, 'gicpError', 0, 'translation', 0, 'rotation', 0), totalFrames, 1);
nStats = 0;
timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
mapFile = sprintf('gicp_map_%s.mat', timestamp);

gicpOpts = struct('max_iters', opts.gicpMaxIters, 'tolerance', opts.gicpTolerance, ...
                  'k_neighbors', opts.gicpKNeighbors, 'max_corr_dist', opts.gicpMaxCorrDist);

tStart = tic;
for frame = 1:totalFrames
    dataRadar = DCA1000_Read_Data(dataFile, params, frame);

    [pointsFrame, velocitiesFrame, nDetect] = process_single_frame(dataRadar, params, frame, opts);
    if isempty(pointsFrame)
        if opts.verbose, fprintf('  帧 %3d: 无有效点\n', frame); end
        continue;
    end
    [filteredPoints, ~] = filter_static_points(pointsFrame, velocitiesFrame, opts);
    if isempty(filteredPoints)
        if opts.verbose, fprintf('  帧 %3d: 静态筛选后无点\n', frame); end
        continue;
    end

    currentXyz = filteredPoints(:, 1:3);
    currentIntensity = filteredPoints(:, 4);
    gicpInfo = struct('iters', 0, 'mean_error', 0);

    if isempty(map_points)
        map_points = currentXyz;
        map_intensity = currentIntensity;
        map_frame_ids = double(frame) * ones(size(currentXyz, 1), 1);
        transforms{frame} = eye(4);
        if opts.verbose
            fprintf('  帧 %3d: 检测 %4d -> 静态 %3d 点 (地图初始帧)\n', ...
                frame, nDetect, size(currentXyz, 1));
        end
    else
        [T, gicpInfo] = gicp_register(currentXyz, map_points, gicpOpts);
        currentXyz = apply_transform(currentXyz, T);

        map_points = [map_points; currentXyz];
        map_intensity = [map_intensity; currentIntensity];
        map_frame_ids = [map_frame_ids; double(frame) * ones(size(currentXyz, 1), 1)];
        transforms{frame} = T;

        if opts.verbose
            fprintf('  帧 %3d: 检测 %4d -> 静态 %3d 点 | G-ICP 迭代 %2d 残差 %.4f m | 平移 %.3f m 旋转 %.2f deg\n', ...
                frame, nDetect, size(currentXyz, 1), gicpInfo.iters, gicpInfo.mean_error, ...
                norm(T(1:3,4)), rotation_angle_deg(T));
        end
    end

    nStats = nStats + 1;
    stats(nStats).frame = frame;
    stats(nStats).rawPoints = nDetect;
    stats(nStats).staticPoints = size(filteredPoints, 1);
    stats(nStats).mapPoints = size(map_points, 1);
    stats(nStats).gicpIters = gicpInfo.iters;
    stats(nStats).gicpError = gicpInfo.mean_error;
    stats(nStats).translation = norm(transforms{frame}(1:3, 4));
    stats(nStats).rotation = rotation_angle_deg(transforms{frame});

    if opts.partialSavePeriod > 0 && mod(frame, opts.partialSavePeriod) == 0 && frame < totalFrames
        save(mapFile, 'map_points', 'map_intensity', 'map_frame_ids', 'transforms', ...
             'params', 'opts', 'stats', 'nStats');
    end
end

elapsed = toc(tStart);

%% 汇总
stats = stats(1:nStats);
result = struct();
result.map_points = map_points;
result.map_intensity = map_intensity;
result.map_frame_ids = map_frame_ids;
result.transforms = transforms;
result.params = params;
result.opts = opts;
result.stats = stats;
result.elapsed_seconds = elapsed;
result.file = mapFile;

if opts.verbose
    fprintf('\n=== 处理完成 ===\n');
    fprintf('地图总点数    : %d\n', size(map_points, 1));
    if ~isempty(stats)
        fprintf('单帧平均点数  : %.1f (中位数 %g, 最大 %d)\n', ...
            mean([stats.staticPoints]), median([stats.staticPoints]), max([stats.staticPoints]));
    end
    fprintf('水平覆盖范围  : x %.2f~%.2f m, y %.2f~%.2f m\n', ...
        min(map_points(:,1)), max(map_points(:,1)), min(map_points(:,2)), max(map_points(:,2)));
    fprintf('耗时          : %.1f s (%.2f s/帧)\n', elapsed, elapsed/totalFrames);
end

save(mapFile, 'map_points', 'map_intensity', 'map_frame_ids', 'transforms', 'params', 'opts', 'stats');
fprintf('结果已保存    : %s\n', mapFile);

%% 结果图
if opts.makeFigure && ~isempty(map_points)
    fig = figure('Name', 'G-ICP Point Cloud Map', 'Position', [100 100 1200 520]);
    subplot(1, 2, 1);
    scatter(map_points(:,1), map_points(:,2), 6, map_intensity, 'filled');
    colormap(fig, 'jet'); colorbar;
    axis equal; grid on; box on;
    xlabel('x (m)'); ylabel('y (m)'); title('(a) Top view');
    subplot(1, 2, 2);
    scatter3(map_points(:,1), map_points(:,2), map_points(:,3), 6, map_intensity, 'filled');
    colormap(fig, 'jet'); colorbar;
    axis equal; grid on; box on;
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); title('(b) 3-D view');
    view(-30, 20);
    print(fig, sprintf('gicp_map_%s.png', timestamp), '-dpng', '-r200');
end

end

%% ======================= 单帧处理链 =======================
function [points3d, velocities, nDetect] = process_single_frame(dataRadar, params, frameId, opts)
%PROCESS_SINGLE_FRAME 单帧数据处理：距离/多普勒/角度 FFT -> CFAR -> 点云

%% 1. 距离维 FFT（加汉明窗）
rangeWin = hamming(params.Samples);
rangeProfile = zeros(params.Samples, params.Chirps, params.virtualAntenna);
for a = 1:params.virtualAntenna
    for m = 1:params.Chirps
        rangeProfile(:, m, a) = fft(dataRadar(:, m, a) .* rangeWin, params.rangeBin);
    end
end

%% 2. 多普勒维 FFT（沿慢时间，加汉明窗并 fftshift）
dopplerWin = hamming(params.Chirps);
speedProfile = zeros(params.Samples, params.Chirps, params.virtualAntenna);
for a = 1:params.virtualAntenna
    for n = 1:params.rangeBin
        speedProfile(n, :, a) = fftshift(fft(rangeProfile(n, :, a) .* dopplerWin.', params.dopplerBin));
    end
end

%% 3. 角度维 FFT（零填充到 angleBin 点）
% 方位: 同一发射通道下的 4 个接收阵元 (1:4)
% 俯仰: 3 个发射通道在同一接收阵元上的虚拟阵元 (1, 5, 9)
angleProfileAz = fftshift(fft(speedProfile(:, :, 1:params.rxNum), params.angleBin, 3), 3);
angleProfileEl = fftshift(fft(speedProfile(:, :, 1:params.rxNum:params.virtualAntenna), ...
                            params.angleBin, 3), 3);

%% 4. 检测
cfarPoints = detect_points_3d_cfar(speedProfile, angleProfileAz, angleProfileEl, params, opts);

%% 5. 极坐标 -> 直角坐标
points3d = [];
velocities = [];
nDetect = 0;
if ~isempty(cfarPoints)
    nDetect = size(cfarPoints, 1);
    [points3d, velocities] = convert_to_3d_coordinates_3d(cfarPoints, frameId);
end

%% 6. 帧内有效性筛选
if ~isempty(points3d)
    valid = (points3d(:,1) > 0.05) & (points3d(:,1) < 20) & ...
            (abs(points3d(:,2)) < 20) & (abs(points3d(:,3)) < 20);
    points3d = points3d(valid, :);
    velocities = velocities(valid, :);
end
end

%% ======================= CFAR 检测 =======================
function points = detect_points_3d_cfar(speedProfile, angleProfileAz, angleProfileEl, params, opts)
%DETECT_POINTS_3D_CFAR 距离-多普勒二维 CFAR 检测 + 角度谱峰值校验
%   门限 T = mu + K*sigma（参考窗 13x13，扣除中心 4x4 保护单元，共 153 个参考单元）

rdMap = abs(squeeze(speedProfile(:, :, 1)));      % 第一个虚拟通道的距离-多普勒图
[numRange, numDoppler] = size(rdMap);

guardCells = [2, 2];
trainingCells = [4, 4];
thresholdFactor = 0.8;
angleRatio = 1.2;
maxPoints = opts.pointLimitPerFrame;

rStart = 1 + guardCells(1) + trainingCells(1);
rEnd   = numRange - guardCells(1) - trainingCells(1);
dStart = 1 + guardCells(2) + trainingCells(2);
dEnd   = numDoppler - guardCells(2) - trainingCells(2);

mask = true(2 * (guardCells(1) + trainingCells(1)) + 1, ...
            2 * (guardCells(2) + trainingCells(2)) + 1);
mask(trainingCells(1) + 1 : trainingCells(1) + 2*guardCells(1), ...
     trainingCells(2) + 1 : trainingCells(2) + 2*guardCells(2)) = false;

buf = zeros(numel(rStart:rEnd) * numel(dStart:dEnd), 9);
nOut = 0;

for rangeIdx = rStart:rEnd
    for dopplerIdx = dStart:dEnd
        currentPower = rdMap(rangeIdx, dopplerIdx);

        trainRegion = rdMap(rangeIdx - trainingCells(1) - guardCells(1) : rangeIdx + trainingCells(1) + guardCells(1), ...
                            dopplerIdx - trainingCells(2) - guardCells(2) : dopplerIdx + trainingCells(2) + guardCells(2));
        trainVals = trainRegion(mask);

        if currentPower <= mean(trainVals) + thresholdFactor * std(trainVals)
            continue;
        end

        azimuthSpectrum = abs(squeeze(angleProfileAz(rangeIdx, dopplerIdx, :)));
        elevationSpectrum = abs(squeeze(angleProfileEl(rangeIdx, dopplerIdx, :)));
        [azimuthPeak, azimuthIdx] = max(azimuthSpectrum);
        [elevationPeak, elevationIdx] = max(elevationSpectrum);

        if azimuthPeak > angleRatio * mean(azimuthSpectrum) && ...
           elevationPeak > angleRatio * mean(elevationSpectrum)
            range = (rangeIdx - 1) * params.c / (2 * params.BandwidthValid);
            speed = (dopplerIdx - params.dopplerBin/2 - 1) * params.speedRes;
            azimuth = azimuthIdx - params.angleBin/2;
            elevation = elevationIdx - params.angleBin/2;

            nOut = nOut + 1;
            buf(nOut, :) = [range, azimuth, speed, currentPower, elevation, ...
                            rangeIdx, dopplerIdx, azimuthIdx, elevationIdx];
        end
    end
end

points = buf(1:nOut, :);
if size(points, 1) > maxPoints
    [~, order] = sort(points(:, 4), 'descend');
    points = points(order(1:maxPoints), :);
end
end

%% ======================= 坐标转换 =======================
function [points3d, velocities] = convert_to_3d_coordinates_3d(points, frameId)
%CONVERT_TO_3D_COORDINATES_3D 极坐标(距离/方位/俯仰) -> 直角坐标，并分解速度矢量

range = points(:, 1);
azimuth = deg2rad(points(:, 2));
speed = points(:, 3);
intensity = points(:, 4);
elevation = deg2rad(points(:, 5));

x = range .* cos(elevation) .* cos(azimuth);
y = range .* cos(elevation) .* sin(azimuth);
z = range .* sin(elevation);

vx = speed .* cos(elevation) .* cos(azimuth);
vy = speed .* cos(elevation) .* sin(azimuth);
vz = speed .* sin(elevation);

points3d = [x, y, z, intensity, frameId * ones(size(x)), range, points(:, 2), points(:, 5)];
velocities = [vx, vy, vz, speed];
end

%% ======================= 静态点筛选 =======================
function [filteredPoints, filteredVel] = filter_static_points(points, velocities, opts)
%FILTER_STATIC_POINTS 静态点筛选：速度门限 -> 距离门控 -> 强度门控
%   观测场景中有一人绕两把椅子走动，运动目标需在配准建图前剔除；
%   保留静态点可使多次观测到同一物体的点云在配准后相互叠加、逐渐加密。
%   points     —— [x, y, z, intensity, frameId, range, azimuth, elevation]
%   velocities —— [vx, vy, vz, 径向速度]

filteredPoints = [];
filteredVel = [];
if isempty(points)
    return;
end

speedMagnitude = sqrt(sum(velocities(:, 1:2) .^ 2, 2));       % 水平速度分量模值
staticIdx = speedMagnitude <= opts.staticSpeedThresh;
filteredPoints = points(staticIdx, :);
filteredVel = velocities(staticIdx, :);
if isempty(filteredPoints)
    return;
end

distances = sqrt(filteredPoints(:, 1).^2 + filteredPoints(:, 2).^2);
rangeIdx = (distances >= opts.minRange) & (distances <= opts.maxRange);
filteredPoints = filteredPoints(rangeIdx, :);
filteredVel = filteredVel(rangeIdx, :);
if isempty(filteredPoints)
    return;
end

intensityThreshold = mean(filteredPoints(:, 4)) * opts.intensityRatio;
strongIdx = filteredPoints(:, 4) >= intensityThreshold;
filteredPoints = filteredPoints(strongIdx, :);
filteredVel = filteredVel(strongIdx, :);
end

%% ======================= G-ICP 配准 =======================
function [T, info] = gicp_register(source, target, opts)
%GICP_REGISTER 广义 ICP：以局部协方差构造权重，加权 SVD 求刚体变换
%   源点云 source 向目标点云 target 配准，返回 4x4 变换 T（p_target = T * p_source）

if nargin < 3
    opts = struct('max_iters', 20, 'tolerance', 1e-4, 'k_neighbors', 12, 'max_corr_dist', 0.8);
end

T = eye(4);
prevError = inf;

sourceCov = point_covariances(source, opts.k_neighbors);
targetCov = point_covariances(target, opts.k_neighbors);

iter = 0;                                    %#ok<NASGU> 保证 max_iters = 0 时输出仍有效
for iter = 1:opts.max_iters
    sourceT = apply_transform(source, T);

    dmat = pairwise_dist_sq(sourceT, target);
    [minDistSq, nnIdx] = min(dmat, [], 2);
    validIdx = minDistSq <= opts.max_corr_dist ^ 2;
    if nnz(validIdx) < 6
        break;
    end

    srcCorr = sourceT(validIdx, :);
    tgtCorr = target(nnIdx(validIdx), :);

    R = T(1:3, 1:3);
    srcCovValid = sourceCov(:, :, validIdx);
    tgtCovValid = targetCov(:, :, nnIdx(validIdx));

    weights = zeros(size(srcCorr, 1), 1);
    for i = 1:size(srcCorr, 1)
        C = tgtCovValid(:, :, i) + R * srcCovValid(:, :, i) * R.';
        weights(i) = 1 / (trace(C) + 1e-6);
    end

    dT = estimate_rigid_transform(srcCorr, tgtCorr, weights);
    T = dT * T;

    diff = tgtCorr - apply_transform(srcCorr, dT);
    meanError = mean(sqrt(sum(diff .^ 2, 2)));
    if abs(prevError - meanError) < opts.tolerance
        prevError = meanError;
        break;
    end
    prevError = meanError;
end

info.iters = iter;
info.mean_error = prevError;
end

function covs = point_covariances(points, k)
%POINT_COVARIANCES 每个点 k 近邻（含自身）的 3x3 局部协方差矩阵，分块计算以控制内存
n = size(points, 1);
covs = zeros(3, 3, n);
if n == 0
    return;
end

k = min(k, n);
if k < 2
    for i = 1:n
        covs(:, :, i) = 1e-6 * eye(3);
    end
    return;
end

chunk = 512;
for s = 1:chunk:n
    e = min(s + chunk - 1, n);
    blk = points(s:e, :);
    dmat = pairwise_dist_sq(blk, points);
    [~, idx] = sort(dmat, 2, 'ascend');
    for r = 1:size(blk, 1)
        neighbors = points(idx(r, 1:k), :);
        c = cov(neighbors);
        if any(isnan(c(:)))
            c = eye(3);
        end
        covs(:, :, s + r - 1) = c + 1e-6 * eye(3);
    end
end
end

function T = estimate_rigid_transform(source, target, weights)
%ESTIMATE_RIGID_TRANSFORM 加权 SVD 求解最优刚体变换（Umeyama）
if nargin < 3 || isempty(weights)
    weights = ones(size(source, 1), 1);
end

weights = weights(:) / (sum(weights) + 1e-12);
cs = sum(source .* weights, 1);
ct = sum(target .* weights, 1);

X = source - cs;
Y = target - ct;
W = diag(weights);

S = X.' * W * Y;
[U, ~, V] = svd(S);
R = V * U.';

if det(R) < 0
    V(:, 3) = -V(:, 3);
    R = V * U.';
end

t = ct.' - R * cs.';

T = eye(4);
T(1:3, 1:3) = R;
T(1:3, 4) = t;
end

function pointsT = apply_transform(points, T)
%APPLY_TRANSFORM 对点云施加 4x4 刚体变换
pointsT = (T(1:3, 1:3) * points.' + T(1:3, 4)).';
end

function d2 = pairwise_dist_sq(A, B)
%PAIRWISE_DIST_SQ 两组点之间的平方距离矩阵（不依赖 Statistics Toolbox）
if isempty(A) || isempty(B)
    d2 = zeros(size(A, 1), size(B, 1));
    return;
end
aa = sum(A .^ 2, 2);
bb = sum(B .^ 2, 2).';
d2 = aa + bb - 2 * (A * B.');
d2(d2 < 0) = 0;
end

function ang = rotation_angle_deg(T)
%ROTATION_ANGLE_DEG 由旋转矩阵提取旋转角（度）
ang = acos(min(1, max(-1, (trace(T(1:3, 1:3)) - 1) / 2))) * 180 / pi;
end
