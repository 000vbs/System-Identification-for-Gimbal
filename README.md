# System Identification for Gimbal

面向 RoboMaster 云台伺服关节的系统辨识与 PID 整定仓库。Yaw 和 Pitch 共用同一条主线：

```text
安全限幅与工作点选择
  -> 开环力矩正弦扫频
  -> 数据质量审查与复频响
  -> 灰箱最小二乘（J、阻尼、摩擦、延迟）
  -> 局部/全行程刚度标定
  -> 按相位裕度设计串级 PI
  -> 摩擦/位置前馈
  -> 1 ms 非线性仿真
  -> 实机阶跃/扫频验证
  -> 参数持久化与报告
```

两轴不能直接套用同一组参数，主要差异来自传动比、行程、重力和静摩擦。

## 目录

- [`gimbal-joint-tuning/`](gimbal-joint-tuning/)：通用辨识、设计方法和固件接入说明。
- [`gimbal-joint-tuning/references/identification_and_design.md`](gimbal-joint-tuning/references/identification_and_design.md)：扫频处理、灰箱 LS、摩擦描述函数、相位裕度设计。
- [`gimbal-joint-tuning/references/pitch_axis_sweep.md`](gimbal-joint-tuning/references/pitch_axis_sweep.md)：Pitch 开环扫频的行程约束。
- [`pitch-axis/`](pitch-axis/)：本仓库 Pitch 轴的完整案例、原始 CSV、MATLAB 脚本和结果图。

## 1. 两轴差异

| 项目 | Yaw | Pitch |
| --- | --- | --- |
| 传动 | 2:1，关节角 = 电机角 / 2；关节侧力矩 = 电机侧 × 2 | 1:1 |
| 行程 | 连续旋转，无机械限位，但要注意线缆/滑环 | 约 0.98 rad，必须做行程保护 |
| 重力 | 基本不受重力 | 受重力和机械重力补偿弹簧影响 |
| 低频扫频 | 可做到 0.3 Hz | 低频摆幅大，通常从 1 Hz 左右开始 |
| 刚度项 | 通常可作为 `k≈0` 的惯性+摩擦对象 | 有位置相关残余刚度，需要静态/全行程标定 |
| 外环对象 | 含 0.5 减速比系数 | 不含减速比，系数为 1 |
| 前馈 | 摩擦前馈为主 | 摩擦前馈 + 位置/弹簧前馈 |
| 主要非线性 | 库仑摩擦、线缆阻尼 | 库仑摩擦、静摩擦死区、弹簧滞回、机械下限 |

## 2. 共同辨识流程

### 2.1 安全配置

开始前先确认：

- 机械限位、线缆余量、最大允许力矩；
- 电机在线和使能状态；
- 角度/速度/力矩软限幅；
- 故障码和饱和计数可观测；
- 记录率足够，不能用几 Hz 的实时波形代替内部高频记录。

### 2.2 选工作点

工作点应避开机械端点，并覆盖实际常用姿态：

- 小信号阶跃用于检查线性和超调；
- 大信号阶跃/斜坡用于检查速度和力矩；
- 开环扫频围绕一个中间位置进行，避免低频摆动撞限位；
- 有位置相关刚度时，至少选两个工作点，分别拟合局部刚度。

本仓库 Pitch 轴使用：

```text
0.15 rad：0.1 rad 瞬时阶跃 + 0.3 rad / 300 ms 斜坡
0.30 rad：开环扫频 + 0.07 rad 角度阶跃
全行程：-0.3244 ~ +0.6256 rad，0.05 rad 静态标定点
```

### 2.3 开环力矩扫频

辨识对象时必须绕过 PID，直接给力矩：

```text
J*theta'' + b_eff*theta' + k*theta = tau(t - Td)
b_eff = b0 + c/|Omega|
c = 4*tau_c/pi
```

逐频点驻留正弦，再做：

1. `[1, t, sin, cos]` 最小二乘，频率在记录值 ±2% 内搜索；
2. 计算复频响 `G = theta / tau`；
3. 剔除 R² 低、三次谐波高、正负不对称或净漂移大的点；
4. 位置通道和速度通道分别做灰箱 LS；
5. 对纯延迟 `Td=0~20 ms` 做一维扫描；
6. 用描述函数估计库仑摩擦 `tau_c`，并与等效阻尼回归交叉验证。

数据判据：

| 指标 | 建议阈值 |
| --- | --- |
| 位置通道 R² | > 0.85 |
| 每周期采样点数 | ≥ 4，最好 ≥ 10 |
| 三次谐波 / 基波 | ≤ 12% |
| 正负摆幅不对称 | < 10% |
| 窗口净漂移 / 振荡速度 | < 5% |

### 2.4 静态/局部刚度

动态扫频给出 `J`、阻尼、摩擦和延迟；静态保持实验给出位置相关保持力矩：

```text
tau_hold(theta)
```

Pitch 的全行程表在两端比二次拟合可靠，二次拟合只适合中间段。局部阶跃可在每个工作点拟合 `kLow/kHigh`，但不要把静摩擦死区强行解释成库仑摩擦。

### 2.5 串级 PID 设计

先用最坏情况 `b=0` 设计内环，保证快速摆动时仍有裕度。速度环：

```text
L = (Kp + Ki/s) / (J*s)
theta_req = PM + omega_c*Td
Ki = J*omega_c^2*cos(theta_req)
Kp = J*omega_c*sin(theta_req)
```

外环使用速度闭环的一阶近似：

```text
Yaw:   |G| = 0.5 / (omega_c*sqrt(1+(omega_c/omega_ci)^2))
Pitch: |G| = 1.0 / (omega_c*sqrt(1+(omega_c/omega_ci)^2))
```

建议：

- 内环相位裕度约 65°；
- 外环相位裕度约 60°；
- 内/外带宽比 ≥ 2.5；
- 机械纹波处 `|C*S| < 0.8`；
- 外环用 PI，不加 D。

### 2.6 前馈

摩擦前馈两轴共用：

```text
tau_friction_ff = tau_c * sat(w_ref / w_lin)
```

必须使用参考角速度，避免实测速度噪声导致符号抖动。Pitch 还需要位置前馈：

```text
tau_cmd = PID_speed(w_ref - w)
        + tau_friction_ff(w_ref)
        + Pitch_Spring_FF(target, center)
```

### 2.7 实机验证

至少记录：

- 参考和目标；
- 实际角度/角速度；
- 总力矩；
- 模式状态和故障码；
- 力矩/速度/积分饱和计数。

被限幅削顶的数据不能用于线性模型验证。

## 3. Yaw 轴流程

Yaw 无机械限位，重点在线缆和传动安全。典型流程：

1. 开环力矩扫频覆盖约 0.3~7 Hz，每点 10~12 个周期；幅值按频段分档。
2. 灰箱拟合 `J、b0、tau_c、Td`，确认 `Re(Z)*|Omega|` 近似常数，区分库仑摩擦与粘性阻尼。
3. 用 2:1 传动换算关节侧力矩和角度环对象。
4. 速度环按 65° 相位裕度设计；外环按 60° 裕度，并检查内外带宽比。
5. 摩擦前馈 `tau_c*sat(w_ref/w_lin)`。
6. 速度阶跃、双环角度阶跃和闭环扫频验证。

通用文档中记录的 Yaw 参考结果：

- 内环穿越频率约 3.5 Hz；
- 总闭环带宽随阻尼工况约 1.55~2.79 Hz；
- 库仑摩擦三路交叉验证约 0.145 / 0.157 / 0.11~0.15 N·m。

Yaw 的具体对象参数、PID 和实测 CSV 应在对应主分支或实机工程中维护，不能直接套用到 Pitch。

## 4. Pitch 轴流程

Pitch 为 1:1 传动，行程只有约 0.98 rad，并且有机械重力补偿、残余弹簧刚度和明显静摩擦。

### 4.1 开环扫频结果

本仓库 Pitch 数据：

- `pitch_rec_6pt.csv`：1.5、2.2、3.2、4.7、6.0、7.5 Hz；
- `pitch_rec_low.csv`：1.0、1.3 Hz 加密点；
- 内部记录约 143 Hz；
- 8 个频点 IMU 拟合 R² 为 0.939~0.996。

模型 B+延迟：

```text
J       = 0.07008 kg*m^2
b0      = 0.7784  N*m*s/rad
tau_c   = 0.9391  N*m
k       = 15.783  N*m/rad
Td      = 16.25 ms
逐点中位相对误差 = 10.4%
留一中位误差     = 18.6%
```

这组数据直接得到固件中的：

```text
pitch_ff_tau_c = 0.94 N*m
pitch_ff_w_lin = 0.30 rad/s
pitch_ff_k     = 15.78 N*m/rad
```

### 4.2 全行程静态标定

9 N·m 标定结果：

```text
20 个目标点，19 个有效点
覆盖 -0.2508 ~ +0.6177 rad
负端 -0.3244 rad 未达到
峰值力矩 9.0 N*m
```

负端不再继续加力矩，优先检查弹簧预紧、线束和机械下限。

二次拟合：

```text
a = -6.9334
b = +8.4749
c = +0.3031
```

两端实际保持力矩明显高于二次曲线，因此前馈优先使用原始表，中间段再用平滑曲线。

### 4.3 角度阶跃联合拟合

4 段独立数据：

- 修复前 0.07 rad；
- 修复后 0.07 rad；
- 0.1 rad 瞬时阶跃；
- 0.3 rad / 300 ms 斜坡。

联合模型：

```text
J       = 0.040535 kg*m^2
b0      = 0.16552  N*m*s/rad
tau_c   = 2.24448  N*m
kLow    = 9.930    N*m/rad  (0.15 rad 工作区)
kHigh   = 42.635   N*m/rad  (0.30 rad 工作区)
tau_act = 5.45 ms
Td      = 16.25 ms
```

注意：阶跃拟合出的 `tau_c=2.24 N*m` 不能覆盖开环扫频得到的 `0.939 N*m`。阶跃中的静摩擦、粘滑释放和位置相关负载会把等效摩擦抬高；动态参数优先相信开环扫频。

拟合误差：

| 数据段 | NRMSE | R² |
| --- | ---: | ---: |
| 0.1 rad 瞬时阶跃 | 6.39% | 0.971 |
| 0.3 rad / 300 ms 斜坡 | 2.94% | 0.995 |
| 修复前 0.07 rad | 18.59% | 0.694 |
| 修复后 0.07 rad | 9.34% | 0.952 |

### 4.4 实机参数与结果

```text
角度环 P        Kpa = 15
速度环 PI       Kpw = 3, Ki_cont = 20 /s
速度环输出限幅  3.5 N*m
速度环积分限幅  2.0 N*m
角度环速度限幅  3.0 rad/s
```

0.07 rad 阶跃：

```text
上升时间   约 51 ms（22 Hz 采样重复估计；45 ms 采样约 80 ms）
超调       2.6%
2% 稳定    209 ms
峰值力矩   6.02 N*m
峰值速度   1.00 rad/s
饱和 tick  0
```

0.1 rad 瞬时阶跃实测超调 18.37%，仍未满足 10% 目标；0.3 rad / 300 ms 斜坡更适合大角度运动。

闭环角度扫频 −3 dB 带宽约 5.5~6 Hz。

完整 Pitch 数据、脚本和图片见 [`pitch-axis/README.md`](pitch-axis/README.md)。

## 5. 固件模式

| `sys_id_method` | 模式 | 用途 |
| ---: | --- | --- |
| 1 | 开环力矩正弦 | 辨识植物，绕过 PID |
| 2 | 速度闭环正弦 | 检查速度环跟踪 |
| 3 | 双环角度正弦 | 检查角度闭环带宽 |
| 4 | Yaw 速度阶跃 | 速度环验证 |
| 5 | 双环角度阶跃 | 角度环验证 |
| 6 | Pitch 速度阶跃 | Pitch 速度环验证 |
| 7 | Pitch 角度阶跃 | 单次角度阶跃 |
| 9 | Pitch 闭环角度扫频 | 角度 Bode |
| 10 | Pitch 多幅值阶跃 | 0.1/0.3 rad |
| 11 | Pitch 全行程静态标定 | 保持力矩表 |

启动钩子：

```text
force_sys_id_axis = 0
sys_id_method     = <模式>
force_sys_id_axis = 1    # Yaw
# 或
force_sys_id_axis = 2    # Pitch
```

## 6. MATLAB 与数据

Pitch 脚本位于 [`pitch-axis/matlab/`](pitch-axis/matlab/)：

| 脚本 | 用途 |
| --- | --- |
| `pitch_ident_261018.m` | 开环扫频质量审查、灰箱 LS、延迟扫描、留一验证 |
| `analyze_pitch_angle_sweep_261019.m` | 闭环角度扫频 Bode |
| `pitch_step_all_csv_joint_fit_261019.m` | 自动去重并联合拟合全部阶跃 |
| `pitch_multi_fit_261019.m` | 0.1/0.3 rad 局部植物拟合 |
| `pitch_static_cal_261018.m` | 静态保持力矩拟合模板 |
| `pitch_closed_loop_sim_261018.m` | 1 ms 非线性串级 PID 仿真 |
| `pitch_gain_sweep_261019.m` | 角度阶跃增益族扫描 |

在仓库根目录运行：

```matlab
addpath('pitch-axis/matlab');
pitch_ident_261018
analyze_pitch_angle_sweep_261019
pitch_step_all_csv_joint_fit_261019
pitch_closed_loop_sim_261018
pitch_gain_sweep_261019
```

输出写入 `pitch-axis/generated/`。

## 7. 分支

- `main`：通用方法和 Yaw 基线。
- `pitch-axis`：Pitch 工作点、扫频、静态全行程标定、阶跃拟合、前馈、PID 和实机数据。

## 8. 已知问题

- Pitch 全行程负端在 9 N·m 下仍无法保持，应先检查机械结构；
- 角度阶跃采样只有约 22~32 Hz，速度峰值和 50~80 ms 上升沿存在量化误差；
- 回中后的静摩擦死区尚未进入普通 `tanh` 动态模型；
- 全行程表与二次系数在 RAM，掉电或重新烧录会丢失；
- `pitch_full_valid_mask` 的高位会被 `(uint16_t)` 截断，仅影响诊断显示。

## 9. 参考

- [灰箱辨识单关节动力学参数（知乎）](https://zhuanlan.zhihu.com/p/164413818)
- [`gimbal-joint-tuning/references/identification_and_design.md`](gimbal-joint-tuning/references/identification_and_design.md)
- [`gimbal-joint-tuning/references/pitch_axis_sweep.md`](gimbal-joint-tuning/references/pitch_axis_sweep.md)
