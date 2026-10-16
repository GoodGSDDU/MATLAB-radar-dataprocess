%% RUN_PROCESSING  毫米波雷达数据处理与建图 —— 主运行脚本
%
%   数据集 : adc_data.bin（TI IWR1843 + DCA1000EVM 采集，100 帧，约 150 MB）
%            场景为一人绕两把椅子走动
%   运行   : 在本文件所在目录下直接运行本脚本即可（按 F5 或输入 run_processing）
%   输出   : gicp_map_<时间戳>.mat   地图点云 + 各帧位姿 + 统计
%            gicp_map_<时间戳>.png   结果图
%
%   环境要求: MATLAB R2018b 及以上；全部代码不使用任何工具箱函数
%             （不使用 Statistics / Image Processing / Computer Vision Toolbox）
%
%   ★ 注意: 必须运行"本文件"（run_processing.m）。直接按 F5 时 MATLAB 只运行编辑器
%           当前打开的那个文件，如果在 DCA1000_Read_Data.m 等函数文件上按 F5，
%           会因为没有输入参数而报错。
%
%   文件清单:
%     run_processing.m        —— 本脚本（主入口）
%     process_radar_data.m    —— 处理链主体：单帧处理 + G-ICP 配准建图
%     DCA1000_Read_Data.m     —— 原始 ADC 数据读取与数据立方体重排
%     radar_constants.m       —— 雷达参数与派生指标
%     verify_results.m        —— 结果验证脚本（跑完后可运行，打印各项统计量）

clear; clc; close all;

%% 0. 路径与数据集定位
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;                       % 直接粘贴到命令行运行时
end
addpath(scriptDir);

parentDir = fileparts(scriptDir);
candidates = { ...
    fullfile(scriptDir, 'adc_data.bin'), ...
    fullfile(parentDir, 'SLAM2', 'ADC_SLAM', 'ADC_VSCode', 'adc_data.bin'), ...
    fullfile(parentDir, 'SLAM2', '01_data_generation', 'adc_data.bin'), ...
    fullfile(parentDir, 'SLAM资料', 'DataAll', 'example', 'radar', 'adc_data.bin')};

DATA_FILE = '';                            % 如需指定数据文件，在此填入完整路径
dataFile = '';
if ~isempty(DATA_FILE)
    dataFile = DATA_FILE;
else
    for k = 1:numel(candidates)
        if exist(candidates{k}, 'file') == 2
            dataFile = candidates{k};
            break;
        end
    end
end

if isempty(dataFile)
    error(['未找到 adc_data.bin。请将数据文件放到本目录下，' ...
           '或在 run_processing.m 中把 DATA_FILE 设为数据文件的完整路径。']);
end
fprintf('数据文件: %s\n', dataFile);

%% 1. 处理参数
opts = struct();
opts.maxFrames        = inf;     % 最多处理多少帧（inf 表示处理文件中的全部帧）
opts.pointLimitPerFrame = 1000;  % 单帧 CFAR 检测点数上限
opts.staticSpeedThresh  = 0.2;   % 静态点速度门限 (m/s)
opts.minRange         = 2.5;     % 距离门下限 (m)
opts.maxRange         = 20.0;    % 距离门上限 (m)
opts.intensityRatio   = 0.7;     % 强度门限系数（相对当帧均值）
opts.gicpMaxIters     = 20;      % G-ICP 最大迭代次数
opts.gicpTolerance    = 1e-4;    % G-ICP 收敛容差
opts.gicpKNeighbors   = 12;      % 局部协方差估计的近邻数
opts.gicpMaxCorrDist  = 0.8;     % 最大对应点距离 (m)
opts.partialSavePeriod = 25;     % 每多少帧保存一次中间结果（0 表示不保存）
opts.makeFigure       = true;
opts.verbose          = true;

%% 2. 运行
result = process_radar_data(dataFile, opts);

%% 3. 结果概览
fprintf('\n=== 结果概览 ===\n');
fprintf('地图总点数    : %d\n', size(result.map_points, 1));
fprintf('参与建图帧数  : %d\n', numel(result.stats));
fprintf('单帧点数均值  : %.1f\n', mean([result.stats.staticPoints]));
fprintf('结果文件      : %s\n', result.file);
fprintf('耗时          : %.1f s\n', result.elapsed_seconds);
fprintf('\n如需完整统计量，请运行 verify_results\n');
