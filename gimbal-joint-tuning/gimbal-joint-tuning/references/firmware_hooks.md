# 固件侧要加的东西、限位表、编译与烧录

以下以 RoboMaster 云台工程（STM32H723 + FreeRTOS + 达妙 DM4310 + `HDL_Gimbal.*`）为例；换项目时查同名符号即可。

## 1. 辨识/测试模式（`sys_id_method`）

| 值 | 模式 | 用途 |
| --- | --- | --- |
| 0 | PRBS | 旧的开环/闭环 PRBS 辨识，噪声大，一般不用于定模型 |
| **1** | **开环力矩正弦扫频** | **辨识对象的主力**：力矩直接给，绕过 PID |
| 2 | 速度闭环正弦 | 环路内扫频，查闭环跟踪 |
| 3 | 双环角度正弦 | 角度参考→角度 PID→速度 PID，查瞄准带宽 |
| 4 | 速度阶跃 | 验证速度环：300 ms 零速基线 → 阶跃保持 N ms → 回零 |
| 5 | 双环角度阶跃 | 验证角度环：0.3 s 后跳到目标并保持（不会自动回零） |

另需 `force_sys_id_axis`（0=关 / 1=强制 Yaw / 2=强制 Pitch）作为启动钩子：遥控任务每周期都会覆盖 `Gimbal.cmd.mode`，所以必须在云台线程里每周期重写回 SYS_ID 模式。

## 2. 辨识模式必须自带的东西

- **纯符号记录量**（供 Ozone 采样）：输入力矩、实际角速度、连续角度（`unwrap_pos`，不要用 ±π 缠绕值）、当前频点索引、实测激励频率、任务频率。
- **实测频率 vs 标称频率**：相位按"每周期整数个采样点"推进，于是 `实际频率 = 标称频率 × (任务频率/1000)`。必须把任务实际频率也记录下来，否则离线拟合会找错频率（3 Hz 以上尤其明显）。
- **保护**：相对起点的行程限位、电机侧超速、反馈超时、力矩硬限幅；并暴露 `*_saturated_ticks` / `*_fault`。**限幅触发过的数据不能用于拟合**。
- **摩擦前馈**：`τ_ff = τ_c·sat(ω_ref/ω_lin)`，`τ_c` 与 `ω_lin` 都要 Ozone 可写（`τ_c=0` 关闭、负值反向用于验证符号）；用**参考**速度驱动。
- **镜像变量**：`float yaw_target_vel_log; float yaw_torque_log;`，在正常串级、模式 4、模式 5 里同步写。

## 3. 限位与力矩（本工程实测值；换项目必须重新读源码）

先读源码，不要凭记忆：`rg -n "PITCH_|YAW_|MAX_TOR|GEAR_RATIO" Users/2-HDL/HDL_Gimbal/HDL_Gimbal.hpp`

| 项 | Yaw | Pitch |
| --- | --- | --- |
| 机械限位 | **无**（可连续旋转；注意线缆缠绕/滑环） | `GIMBAL_PITCH_MAX = +0.645448804 rad`（+36.98°）<br>`GIMBAL_PITCH_MIN = −0.334366739 rad`（−19.15°）<br>行程 ≈ **0.98 rad（56°）**，IMU 空间绝对角度 |
| 力矩限幅 | `GIMBAL_SPEED_MAX_OUT = 9.0 N·m`（速度环输出）<br>`YAW_LQI_MAX_TOR = 9.0 N·m`（旧 LQI 遗留值） | `PITCH_MAX_TOR = 6.0 N·m` |
| 传动比 | `GIMBAL_YAW_GEAR_RATIO = 2.0` → 关节角 = 电机角/2，**关节侧力矩 = 电机侧×2** | 无（1:1） |
| 电机力矩刻度 | `dm_motor[...].tmp.TMAX = 10.0 N·m`（MIT 帧满刻度）；DM4310 额定约 3、峰值约 7 N·m | 同 |

要点：

- **限位是角度的，不是力矩的**。Pitch 开环扫频必须在软限位处切断激励力矩（本工程在 `Gimbal_Sine_Sweep_Pitch` 里做）；Pitch 受重力，低频开环激励会一直压向限位，**做 Pitch 开环扫频前先确认行程、缓冲、或加弹簧/配重**。
- **力矩的真实约束**：DM4310 峰值约 7 N·m；Yaw 经 2:1 减速后关节侧是电机侧的 2 倍（9 N·m 电机 ≈ 18 N·m 关节），先确认皮带/结构承受能力——"测试保护值"和"物理上限"是两回事。
- **`TMAX` 只是 CAN 打包刻度**：驱动器内部实际配置不同则"指令 ≠ 实际"。**校核方法**：同时记录 `ctrl.tor_set`（指令）与 `para.tor`（驱动器回报），稳恒力矩下两者之比就是真实刻度；辨识出的 τ_c 依赖这个刻度。
- **测试模式的力矩限幅**要在代码里做合法性校验（本工程：模式 3/4 ≤1.5 N·m、模式 5 ≤4.0 N·m），并把默认值调到"不削顶"的水平，否则每次上机都要手动改。

## 4. 编译与烧录（Windows + CMake/Ninja + OpenOCD + DAPLink）

```powershell
# 编译（工程里通常已有可用 build 目录）
ninja -C build/<预设目录>
# 或： mingw32-make -C build

# 烧录（CMSIS-DAP/DAPLink，不进入调试会话）
openocd -f config/daplink.cfg -c "gdb_port disabled" `
        -c "program build/GXU-ROBOTZ-Season-26-hero.elf verify reset exit"
```

- 探针被 Ozone 占用时会报 `CMSIS-DAP command 0x0 not implemented`：关掉 Ozone，或直接在 Ozone 里 Download。
- **烧完必须校验镜像里是新值**（源码对 ≠ 烧对）：

```powershell
arm-none-eabi-nm <elf> | Select-String "<新符号>"
arm-none-eabi-objdump -s -j .data <elf> | Select-String "<符号地址>"
# 或直接用 Python 搜 float32 字节： struct.pack("<f", 0.3052)
```

## 5. 改这份固件时的两个坑

- 目标行若含**全角标点**（`，` `。`）或特殊括号，`apply_patch` 可能匹配失败。对策：把改动落在相邻的纯 ASCII 行上；或在 `Gimbal_Init()` 里用显式赋值覆盖默认值（例如模式限幅默认值写在 `Gimbal_Init`，而不是改那行带全角逗号的初始化）。
- 多个模式共用同一组 PID 结构体时，**别让前馈污染辨识链路**：前馈只加在正常控制路径与验证模式（4/5），不要加进模式 2/3 的辨识路径，否则辨识出的对象会被前馈掩盖。
