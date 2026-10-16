function data_radar = DCA1000_Read_Data(fname, params, frame)
%DCA1000_READ_DATA 读取 DCA1000EVM 采集的第 frame 帧原始 ADC 数据并重排
%
%   data_radar = DCA1000_Read_Data(fname, params, frame)
%
%   输入
%     fname  —— adc_data.bin 的完整路径
%     params —— 由 radar_constants() 返回的参数结构体
%     frame  —— 帧序号 (1 起)
%   输出
%     data_radar(Samples, Chirps, txNum*rxNum) —— 复数数据立方体
%
%   数据格式：16bit 有符号整数，复数 I/Q 以 4 个字为一组存放
%             [I0 I1 Q0 Q1] -> 2 个复数采样 (I0+jQ0, I1+jQ1)
%
%   与原实现的区别（均为修复，不改变算法）
%     1. 用 fseek 直接定位目标帧，不再读入前面所有帧后丢弃。
%        原实现读取量为 O(N_frame^2)，100 帧需读约 8 GB，此处仅需 150 MB。
%     2. 原实现以 4 个字为步长的循环条件为 i <= fileSize-5，会漏掉最后 4 个字
%        （即每帧最后 2 个复数采样被置零），此处处理全部分组。
%     3. 去掉了 global 变量，参数与帧序号改为显式传参。

if nargin < 3
    error('DCA1000_Read_Data:notEnoughInputs', ...
        ['本文件是函数文件\n' ...
         '请运行主脚本 run_processing.m（在编辑器中打开它后按 F5），\n' ...
         '或在命令行窗口输入：  run_processing\n\n' ...
         '本函数用法: data_radar = DCA1000_Read_Data(fname, params, frame)']);
end

nWords = params.wordsPerFrame;
nBytes = params.bytesPerFrame;

fid = fopen(fname, 'rb');
if fid < 0
    error('DCA1000_Read_Data:openFailed', '无法打开数据文件: %s', fname);
end

try
    if fseek(fid, (frame - 1) * nBytes, 'bof') ~= 0
        error('DCA1000_Read_Data:seekFailed', '定位第 %d 帧失败: %s', frame, fname);
    end
    sdata = fread(fid, nWords, 'int16=>double');
catch err
    fclose(fid);
    rethrow(err);
end
fclose(fid);

if numel(sdata) < nWords
    error('DCA1000_Read_Data:shortRead', ...
        '第 %d 帧数据不完整：期望 %d 个 int16，实际读到 %d 个。', ...
        frame, nWords, numel(sdata));
end

%% I/Q 解交织
q = reshape(sdata, 4, []);
lvds = zeros(1, 2 * size(q, 2));
lvds(1:2:end) = q(1, :) + 1i * q(3, :);
lvds(2:2:end) = q(2, :) + 1i * q(4, :);

%% 重排为数据立方体 [Samples, Chirps, virtualAntenna]
% 每 3072 个字为一个慢时间分组，组内按虚拟通道优先排列：通道1的256点、通道2的256点...
nSamples = params.Samples;
lvds = reshape(lvds, nSamples * params.virtualAntenna, params.Chirps).';   % Chirps x (Samples*12)
cdata = zeros(nSamples, params.Chirps, params.virtualAntenna);
for a = 1:params.virtualAntenna
    blk = lvds(:, (a - 1) * nSamples + (1:nSamples));       % Chirps x Samples
    cdata(:, :, a) = blk.';
end

data_radar = cdata;

end
