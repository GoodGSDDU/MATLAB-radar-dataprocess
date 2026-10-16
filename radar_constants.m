function params = radar_constants()
%RADAR_CONSTANTS 雷达与波形参数（TI IWR1843 + DCA1000EVM）
%
%   本文件集中定义全部雷达参数，所有派生指标（分辨率、刻度因子、帧字节数等）
%   均由参数自动计算，避免在其他文件中硬编码。
%
%   数据立方体结构: data_radar(Samples, Chirps, txNum*rxNum)
%     Samples  —— 每个 chirp 的快时间采样点数
%     Chirps   —— 慢时间分组数（每个分组对应 txNum 个 chirp，见下）
%     12       —— 虚拟通道数 (TX1RX1 ... TX3RX4)

%% +++物理常数
params.c = 3e8;                                     % 光速 (m/s)

%% +++雷达波形参数
params.stratFreq = 77e9;                            % 起始频率 (Hz)
params.Tr        = 60e-6;                           % 发射扫频时间 (s)
params.Idle_time = 100e-6;                          % chirp 空闲时间 (s)
params.Tc        = params.Tr + params.Idle_time;    % chirp 间隔 (s)
params.Slope     = 29.982e12;                       % 调频斜率 (Hz/s)
params.Samples   = 256;                             % 每 chirp 采样点数
params.Fs        = 10e6;                            % 采样率 (Hz)
params.Tframe_set = 80e-3;                          % 帧周期 (s)

%% +++处理维度
params.rangeBin   = params.Samples;                 % 距离维 FFT 点数
params.Chirps     = 128;                            % 慢时间分组数
params.dopplerBin = params.Chirps;                  % 多普勒维 FFT 点数
params.angleBin   = 180;                            % 角度维 FFT 点数（零填充）

%% +++阵列参数 (3T4R -> 12 虚拟通道)
params.txNum = 3;                                   % 发射通道数
params.rxNum = 4;                                   % 接收通道数
params.virtualAntenna = params.txNum * params.rxNum;
params.dx = 0.5 * (params.c / (params.stratFreq + params.Slope * params.Tr / 2));  % 水平阵元间距
params.dz = params.dx;                              % 垂直阵元间距

%% +++派生量：带宽与波长
params.Bandwidth      = params.Slope * params.Tr;                        % 完整扫频带宽 (Hz)
params.BandwidthValid = params.Samples / params.Fs * params.Slope;       % 有效带宽（受采样窗限制）
params.centerFreq     = params.stratFreq + params.Bandwidth / 2;         % 中心频率 (Hz)
params.lambda         = params.c / params.centerFreq;                    % 波长 (m)

%% +++派生量：数据文件结构
% 每个 chirp 在原始文件中占用的 16bit 字数 = 每 chirp 采样点数 × 接收通道数
params.samplesPerChirp = params.Samples * params.rxNum;
% 一个慢时间分组（共 12 个虚拟通道）所覆盖的 chirp 数
params.chirpsPerGroup  = params.virtualAntenna / params.rxNum;
% 同一虚拟通道相邻两次慢时间采样的实际间隔（TDM 分时发射）
params.TcGroup         = params.Tc * params.chirpsPerGroup;
% 单帧原始数据量
params.wordsPerFrame = params.Samples * params.Chirps * params.virtualAntenna * 2;  % int16 个数
params.bytesPerFrame = params.wordsPerFrame * 2;                                    % 字节数

%% +++速度刻度
% 12 个虚拟通道由 3 个发射通道分时获得，因此同一虚拟通道的有效脉冲重复间隔
% 为 TcGroup = 3*Tc = 480us，而原实现按 Tc = 160us 计算，速度值会偏大 3 倍。
%   false —— 与原实现一致（论文数据即由此得到）
%   true  —— 按物理正确的 TcGroup 计算
params.dopplerUseGroupPRI = false;
if params.dopplerUseGroupPRI
    params.PRI = params.TcGroup;
else
    params.PRI = params.Tc;
end

%% +++派生量：性能指标（供打印与论文引用）
params.rangeRes   = params.c / (2 * params.BandwidthValid);   % 距离分辨率 (m)
params.rangeMax   = params.Samples * params.rangeRes;         % 最大不模糊距离 (m)
params.speedRes   = params.lambda / (2 * params.dopplerBin * params.PRI);  % 速度分辨率 (m/s)
params.speedMax   = params.lambda / (4 * params.PRI);                       % 最大不模糊速度 (m/s)
params.angleResAz = 2 / params.rxNum * 180 / pi;              % 方位角分辨率 (deg, 近似)
params.angleResEl = 2 / params.txNum * 180 / pi;              % 俯仰角分辨率 (deg, 近似)

end
