# Pitch 轴专用：正弦力矩扫频（有重力 + 有限位）

Yaw 轴无重力、无限位，所以"给力矩就能摆动"。Pitch 完全不同，**照搬 yaw 的参数会没有数据**：轴被重力/弹簧压向限位，力矩全被软限位切断。

## 本工程固件已实现的支撑（2026-10-15）

| 变量 | 含义 |
| --- | --- |
| `sweep_ol_amp_scale_pitch` | Pitch 专用幅值缩放（与 yaw 的 `sweep_ol_amp_scale` 独立）。默认 1.0 |
| `sweep_ol_pitch_gA` | 重力补偿幅值 (N·m)，默认 0 = 关闭 |
| `sweep_ol_pitch_gPhi` | 重力补偿相位 (rad)：力矩峰值对应的 IMU 俯仰角 |
| `sweep_ol_pitch_gC` | 重力补偿常量偏置 (N·m) |
| `sweep_ol_pitch_torque_limit` | **Pitch 开环扫频总力矩限幅 (N·m)，默认 3.0**；填 0 或负值关闭 |
| `sweep_ol_pitch_clip_ticks` | 本轮被限幅的周期数；**非 0 说明限幅咬住**（多半是 gC 填太大） |
| `sweep_ol_u_total` | 实际施加的总力矩 = 正弦 + 重力补偿，Ozone 直接可读 |

补偿形式：`τ_bias = gA·cos(θ_imu − gPhi) + gC`（刚性臂的重力力矩就是这个形状）。

软限位：`Gimbal_Sine_Sweep_Pitch()` 里超出 `[GIMBAL_PITCH_MIN, GIMBAL_PITCH_MAX]` 会把**正弦与补偿一起清零**——安全，但**那一段数据作废**，事后要按时间剔除。

记录量：`BPRS_log_u_pitch`（总力矩）、`BPRS_log_pos_pitch`、`BPRS_log_w_pitch`、`sweep_ol_u`（正弦分量）、`sweep_ol_u_total`、`sweep_ol_unwrap_pos`（连续角度）、`sweep_ol_freq_index`、`sweep_ol_actual_hz`、`gimbal_task_freq`。

## 上机顺序（先标定重力补偿，再扫频）

1. `manual_target_mode = 0`（方法 1 **不会**自动清它），`force_sys_id_axis = 0`，确认 `sweep_ol_state = 0`。
2. **台架标定重力补偿——不要一上来就扫频**：
   - 先把轴摆到行程中间，把 `sweep_ol_amp_scale_pitch` 设为 0 或让 `start_index > end_index`（即不产生正弦），只用 `sweep_ol_pitch_gC` 慢慢加，找到"轴能停住"的值与对应角度；
   - 在 3~4 个不同角度重复，得到 (θ, τ) 点，用 `τ = gA·cos(θ−gPhi) + gC` 拟合（3 个参数，至少 3 个点）；
   - 填回三个参数后手动验证：轻推轴应"发软"但不会自己滑到限位（有静摩擦死区，别只看它停在哪）。
3. **定幅值**：`sweep_ol_amp_scale_pitch = 1.0` 起，只看低频段（0.3~0.45 Hz）实际摆幅。判据：**摆幅 >0.02 rad** 且记录里 `sweep_ol_u_total` 没有长时间被限位切断。不够就提到 2~3（注意 DM4310 峰值约 7 N·m，且要留行程余量：行程总共 0.98 rad）。
4. **启动**：`sys_id_method = 1` → 开采样 → `force_sys_id_axis = 2`（**2 = Pitch**，不是 1）。`sweep_ol_state` 会 1→2。
5. **判读**：与 yaw 相同的五个质量判据（R²>0.85、每周期采样≥4、3 次谐波≤12%、摆幅不对称<10%、漂移<5%）。`sweep_ol_u_total − sweep_ol_u` 就是重力补偿量。

## 辨识时的两个陷阱

- **补偿项等效于负刚度**：`gA·cos(θ)` 对角度求导是 `−gA·sin(θ)`。做**增量辨识**（为了控制器设计）建议拿 `sweep_ol_u`（纯正弦）当输入；要建模"实际系统"就用 `sweep_ol_u_total`。两个都记录了，分析时选一个并写明。
- **重力让低频"变刚"**：<0.5 Hz 的响应由重力/弹簧主导，可能测不出 J、b。对策：把最低频点抬到 0.5 Hz 以上，另用**准静态法**（缓慢斜坡给力矩、看稳态角）单独定刚度/重力项。对有回中弹簧的轴，k 会被抬高，这也是"pitch 与 yaw 模型结构不同"的常见原因。

## 验证

- Pitch 目前**没有专用阶跃模式**（本工程模式 5 是 yaw 双环角度阶跃）；要验证就先用模式 1 的低频点近似方波，或照 `Gimbal_Yaw_Speed_Step` 的样子加一个 pitch 阶跃模式。
- 验收仍看：`saturated_ticks = 0`（没被削顶）、`fault = 0`、实机波形与模型趋势一致、力矩不反复大幅跳变。
