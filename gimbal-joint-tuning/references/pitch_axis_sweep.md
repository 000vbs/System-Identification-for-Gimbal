# Pitch 轴专用：正弦力矩扫频（行程有限，但已有机械重力补偿）

先分清两种 Pitch，做法完全不同：

| 机构 | 对象结构 | 开环扫频可行性 |
| --- | --- | --- |
| **有机械重力补偿（全行程）** ← 本工程 | `J·θ̈ + b·θ̇ = τ − 摩擦` —— 重力连同它的等效刚度都被机构抵消，**与 Yaw 同构** | 可以直接做，和 Yaw 一样 |
| 无补偿 / 只有弹簧 | `J·θ̈ + b·θ̇ + k·θ = τ`，k = 弹簧刚度 − 重力等效刚度；有共振点 `f_n = (1/2π)√(k/J)` | 必须先做重力补偿标定（见文末），频点还要跨过 f_n |

判断方法：手轻推离平衡位置再放开——**回位并来回振荡**=稳定（弹簧主导，按下面第二种处理）；**朝一个方向滑走**=重力负刚度占上风（务必用很小的幅值＋0.2 rad 的角度限幅）。

## 本工程的做法（全行程机械补偿 → 与 Yaw 同构）

结论：**不需要重力补偿偏置**，`sweep_ol_pitch_gA / gPhi / gC` 保持 0；也**不需要共振点分析**。剩下的唯一约束是"行程有限"。

### 唯一的硬约束：低频点必须砍掉

Pitch 总行程只有 `GIMBAL_PITCH_MIN..MAX` ≈ 0.98 rad；而低频摆幅 ∝ 1/f。用"刚能克服摩擦"的幅值（0.3 N·m 量级）外推：

| 频点 | 0.3 N·m 下预计摆幅 | 结论 |
| --- | --- | --- |
| 0.30 Hz | ≈1.2 rad | **顶限位，不可用** |
| 0.45 Hz | ≈0.66 rad | 不可行 |
| 0.70 Hz | ≈0.29 rad | 贴边可用 |
| 1.00 Hz | ≈0.14 rad | 安全 |

所以默认从 **0.70 Hz 起扫到 3.20 Hz**（4.7/7 Hz 因采样率不够也砍掉）。

### 固件里的 Pitch 专用默认值（不改任何变量即可直接扫）

| 变量 | 默认 | 含义 |
| --- | --- | --- |
| `sweep_ol_start_index_pitch` | **2** | 起始频点索引（0.70 Hz） |
| `sweep_ol_end_index_pitch` | **6** | 结束频点索引（3.20 Hz） |
| `sweep_ol_excursion_limit_pitch` | **0.30 rad** | 相对起点的角度限幅，比机械软限位更早、更干净地停 |
| `sweep_ol_amp_scale_pitch` | 1.0 | 幅值缩放；低频段摆幅 <0.02 rad 时往上加 |
| `sweep_ol_pitch_torque_limit` | **3.0 N·m** | 总力矩（正弦+补偿）限幅，0=关闭 |
| `sweep_ol_pitch_clip_ticks` | 0 | 本轮被力矩限幅的周期数，**非 0 说明咬住了** |

（Yaw 侧对应的是 `sweep_ol_start_index`/`sweep_ol_end_index`/`sweep_ol_excursion_limit`，默认 0/8/1.5 rad —— 两轴配置彼此独立，改一个不影响另一个。）

### 操作：零改动流程

```
force_sys_id_axis = 0        # 复位（确认 sweep_ol_state=0、sweep_done=0）
sys_id_method     = 1        # 开环力矩正弦扫频
（先把 Data Sampling 打开）
force_sys_id_axis = 2        # 2 = Pitch（1 = Yaw）
```

运行中看 `sweep_ol_state`(1→2)、`sweep_ol_freq_index`(2→6)、`sweep_ol_actual_hz`、`sweep_ol_fault`(必须 0)。随时 `force_sys_id_axis=0` 中止；中止即复位，下次从 `start_index` 重新来。全程约 6 个频点 × 12 周期 ≈ 78 s。

采样变量（纯符号）：`sweep_ol_freq_index`、`sweep_ol_actual_hz`、`sweep_ol_u`（正弦分量）、`sweep_ol_u_total`（实际总力矩）、`sweep_ol_unwrap_pos`（连续角度，位置拟合用这个）、`sweep_ol_state`、`sweep_ol_fault`、`sweep_ol_pitch_clip_ticks`、`BPRS_log_u_pitch`、`BPRS_log_w_pitch`、`gimbal_task_freq`。

### 辨识与设计注意

- 模型 A（`1/(J s + b)`）就够；脚本报出的 **k 应≈0** —— 若明显非零，那正好量化了"机械补偿的残余误差"，值得记下来。
- 辨出的 J、b、τ_c 是 Pitch 自己的数值，**别沿用 Yaw 的 PID**。
- Pitch 是 **1:1 传动**：外环对象里没有 Yaw 那个 0.5 的减速比系数，角度环增益要用 `references/identification_and_design.md` 的外环公式重算（`|G| = 1/(ωc·√(1+(ωc/ωci)²))`）。

## 附：无机械补偿的 Pitch 才需要的重力补偿标定

1. 让正弦不产生（`start_index > end_index`，或把 `amp_scale_pitch` 设 0），只调 `sweep_ol_pitch_gC` 直到轴能停在某个中间角；
2. 在 3~4 个角度重复，得到 (θ, τ) 点，拟合 `τ = gA·cos(θ − gPhi) + gC`；
3. 填回三个量后手动验证：轻推应"发软"但不自己滑到限位（有静摩擦死区，别只看它停在哪）。
4. 补偿项等效于负刚度（`gA·cos(θ)` 的导数是 `−gA·sin θ`）：做增量辨识用 `sweep_ol_u`（纯正弦）当输入更干净，要建"实际系统"就用 `sweep_ol_u_total`；两个量都记录了。

## 实机案例

本工程 Pitch 轴的完整结果、原始 CSV、MATLAB 辨识/阶跃拟合脚本和实机验证图见：

[`../../pitch-axis/README.md`](../../pitch-axis/README.md)

已完成的 Pitch 开环扫频得到 `J=0.07008 kg·m²`、`b0=0.7784 N·m·s/rad`、`tau_c=0.9391 N·m`、`k=15.783 N·m/rad`、`Td=16.25 ms`，逐点中位误差 10.4%。固件中的摩擦前馈 `pitch_ff_tau_c=0.94` 即来自这组数据。
