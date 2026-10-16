# 毫米波雷达数据处理与 G-ICP 建图（MATLAB）

基于 **TI IWR1843 + DCA1000EVM** 采集的原始 ADC 数据，完成从原始信号到点云地图的完整处理链，并用 **G-ICP** 做多帧配准建图。

采集场景为**一人在两把椅子之间走动**，共 **100 帧**（约 150 MB）。

全部代码**不使用任何工具箱函数**（不依赖 Statistics / Image Processing / Computer Vision Toolbox），MATLAB **R2018b 及以上**可直接运行。

## 快速开始

```matlab
run_processing        % 主入口：处理全部帧，输出 gicp_map_<时间戳>.mat / .png
verify_results        % 打印各项统计量，可与论文数值逐项对照
```

`run_processing.m` 会自动在若干候选路径查找 `adc_data.bin`，找不到时按提示把数据放到本目录下即可；也可以直接在该文件的 `DATA_FILE` 中填入完整路径。

运行后会生成两个输出文件：`gicp_map_<时间戳>.mat`（地图点云 + 各帧位姿 + 逐帧统计）与同名 `.png`（结果图）。

> **注意**：这两个输出未纳入仓库，因此 `verify_results` 需要先成功跑过一次 `run_processing` 才能使用——它读取的是 `.mat` 里的建图结果，而不是原始数据。

## 文件清单

| 文件 | 作用 |
| --- | --- |
| [`run_processing.m`](run_processing.m) | 主入口脚本：定位数据、设置参数、调用处理链、打印概览 |
| [`process_radar_data.m`](process_radar_data.m) | 处理链主体：逐帧处理 + G-ICP 配准建图 + 结果图 |
| [`DCA1000_Read_Data.m`](DCA1000_Read_Data.m) | 原始 ADC 读取、I/Q 解交织、数据立方体重排 |
| [`radar_constants.m`](radar_constants.m) | 雷达参数与全部派生指标（单一真源） |
| [`verify_results.m`](verify_results.m) | 结果验证：点云规模、空间分布、配准位姿、场景结构 |
| 输出 `gicp_map_<时间戳>.mat` / `.png` | 运行后生成：地图点云 + 各帧位姿 + 统计 / 结果图 |

## 处理流程

```
读取 ADC
  └─ DCA1000_Read_Data.m        fseek 定位目标帧 -> I/Q 解交织 -> 立方体 (Samples, Chirps, 12)
      ↓
距离维 FFT（汉明窗）
      ↓
多普勒维 FFT（沿慢时间，汉明窗 + fftshift）
      ↓
角度维 FFT（方位 / 俯仰，零填充至 180 点）
      ↓
三维 CFAR 检测 + 角度峰值校验        detect_points_3d_cfar
      ↓
极坐标 -> 直角坐标                   convert_to_3d_coordinates_3d
      ↓
静态点筛选（速度门限）                filter_static_points
      ↓
与累积地图作 G-ICP 配准               gicp_register
      └─ point_covariances / estimate_rigid_transform / apply_transform
      ↓
并入地图，进入下一帧
```

## 雷达参数

集中定义在 [`radar_constants.m`](radar_constants.m)，所有派生指标由参数自动计算，不在别处硬编码。

| 参数 | 值 |
| --- | --- |
| 起始频率 | 77 GHz |
| 调频斜率 | 29.982 THz/s |
| 扫频 / 空闲时间 | 60 / 100 µs |
| chirp 间隔 Tc | 160 µs |
| 采样点数 / 采样率 | 256 / 10 MHz |
| 慢时间分组数 Chirps | 128 |
| 阵列 | 3 TX × 4 RX = 12 虚拟通道 |
| 帧周期 | 80 ms |
| 距离维 / 多普勒维 / 角度维点数 | 256 / 128 / 180 |
| 单帧原始数据量 | 1.5 MiB |

数据立方体统一为 `data_radar(Samples, Chirps, txNum*rxNum)`。

### 关于速度刻度的一处口径问题

12 个虚拟通道由 3 个发射通道**分时（TDM-MIMO）**获得，因此同一虚拟通道相邻两次慢时间采样的实际间隔为

```
TcGroup = 3 × Tc = 480 µs
```

而原实现按 `Tc = 160 µs` 计算，**速度值会偏大 3 倍**。`radar_constants.m` 中的开关控制这一选择：

- `params.dopplerUseGroupPRI = false`（默认）—— 与原实现一致，**论文数据即由此得到**
- `params.dopplerUseGroupPRI = true` —— 按物理正确的 `TcGroup` 计算，速度值相应变为 1/3

改动此开关会改变速度门限的判据，请留意其影响。

## 结果参考值

[`verify_results.m`](verify_results.m) 内置的参考值（100 帧原始数据集、默认参数）：

| 指标 | 数值 |
| --- | --- |
| 地图总点数 | 3413 |
| 参与建图帧数 | 100 |
| 单帧点数（均值 / 中位 / 最大） | 34.1 / 34 / 48 |
| 水平覆盖 | x 0.44 ~ 18.32 m，y −4.55 ~ 7.17 m |
| 点云半径（均值 / 95% 分位） | 6.81 / 17.94 m |
| 强度范围 | 5.25×10⁴ ~ 1.33×10⁶ |
| 配准平移（均值 / 最大） | 0.076 / 0.336 m |
| 配准旋转（均值 / 最大） | 1.17 / 4.24° |

> 雷达在采集中固定不动，因此上表的配准位姿量级即**配准误差**——平均平移 7.6 cm、平均旋转 1.17°，可用于判断建图是否收敛。

## 相对原始版本的修复

算法本身未做改动，因此结果与论文一致、可复现。以下修复记录在源码注释中：

| # | 修复 | 影响 |
| --- | --- | --- |
| 1 | 补充缺失的 `radar_constants.m` 依赖 | 原版本调用该函数但目录中没有该文件，无法运行 |
| 2 | 帧读取改为 `fseek` 直接定位 | 原实现每帧都从头读入，100 帧累计读取约 8 GB；现仅需 150 MB |
| 3 | 修正 I/Q 解交织的循环边界 | 原实现漏掉每帧最后 4 个 `int16`（即最后 2 个复数采样被置零） |
| 4 | 明确区分 `Tc` 与 `TcGroup`(3·Tc) | 见上文「关于速度刻度的一处口径问题」 |
| 5 | 去掉 `global` 变量与函数内的 `clear` | 参数全部显式传递，函数可重入 |
| 6 | 局部协方差与距离矩阵改为分块计算 | 避免地图增长时的内存峰值 |
| 7 | `map_frame_ids` 显式存为 `double` | 原实现帧数超过 255 时溢出 |
| 8 | 输出含逐帧统计与位姿 | 便于结果验证（配合 `verify_results.m`） |

## 数据

`adc_data.bin`（TI IWR1843 + DCA1000EVM 原始采集，100 帧，约 150 MB）**未纳入仓库**——超出 GitHub 单文件 100 MB 上限。使用时把它放到本目录下即可，`run_processing.m` 会自动找到。

建图结果（`gicp_map_*.mat` / `.png`）为程序输出，同样未纳入仓库。

## 致谢

- 雷达平台：Texas Instruments IWR1843 + DCA1000EVM
