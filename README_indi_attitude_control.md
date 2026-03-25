# 基于 INDI 的滑翔飞行器 BTT 姿态控制仿真说明

## 1. 文件说明

主脚本：

- `indi_attitude_control_sim.m`

该脚本用于验证滑翔飞行器在时变法向过载指令、侧向过载指令作用下的姿态控制能力。当前版本采用：

- BTT 转弯控制
- 角速度内环 INDI
- 4 个 X 布局尾舵模态分解
- 更完整的平动 + 转动耦合动力学

当前版本不包含：

- 对称增阻模态
- 时间约束控制
- 舵机动态与舵速限制

---

## 2. 为什么这版比上一版动力学更完整

上一版的问题在于：

- `gamma` 没有单独积分
- 直接用 `gamma = theta - alpha` 代替航迹动力学
- 导致 `theta` 可以在缺少平动约束的情况下漂移到不合理值

当前版本做了以下修正：

1. 单独积分位置：

- `x`
- `y`
- `z`

2. 单独积分平动速度与速度方向：

- `V`
- `gamma`
- `chi`

3. 单独积分姿态角：

- `phi`
- `theta`
- `psi`

4. 单独积分角速度：

- `p`
- `q`
- `r`

5. 单独积分气动状态：

- `alpha`
- `beta`

因此现在：

- `theta` 是姿态角
- `gamma` 是航迹倾角
- `alpha` 是攻角

三者不再被强行绑成一个代数关系，而是通过动力学自然耦合。

---

## 3. 当前动力学模型

### 3.1 平动动力学

采用速度坐标系下的简化点质点平动模型：

```matlab
x_dot = V * cos(gamma) * cos(chi)
y_dot = V * cos(gamma) * sin(chi)
z_dot = V * sin(gamma)
V_dot = -D/m - g*sin(gamma)
gamma_dot = (L*cos(phi)/m - g*cos(gamma)) / V
chi_dot = (L*sin(phi)/(m*V*cos(gamma)))
```

其中：

- `L` 由攻角产生
- `phi` 决定升力在法向和侧向的分解

这正是 BTT 转弯的核心物理机制。

---

### 3.2 转动动力学

刚体角运动采用：

```matlab
J * omega_dot = M_aero + M_ctrl - cross(omega, J*omega)
```

其中：

- `omega = [p; q; r]`
- `M_ctrl` 来自 4 个 X 尾舵
- `M_aero` 为简化的阻尼 + 静稳定力矩

---

### 3.3 攻角与侧滑角动力学

为了兼顾控制验证与动力学一致性，脚本采用：

```matlab
alpha_dot = q - gamma_dot - a_alpha*alpha + b_alpha*up
beta_dot  = r - chi_dot   - a_beta *beta + b_beta *uy
```

这里的关键是：

- `alpha_dot` 中显式减去了 `gamma_dot`
- `beta_dot` 中显式减去了 `chi_dot`

这使得 `alpha`、`beta` 不再脱离平动单独漂移，而是反映“机体姿态相对速度方向”的变化。

---

## 4. BTT 控制逻辑

### 4.1 控制思想

当前版本改成了更符合 BTT 的控制结构：

1. `beta_cmd = 0`
   表示尽量保持协调转弯。

2. 攻角 `alpha`
   主要决定总升力，也就是总过载大小。

3. 滚转角 `phi`
   主要决定总升力在法向和侧向上的分配方向。

因此过载输出采用：

```matlab
n_total = L / (m*g)
n_n = n_total * cos(phi)
n_l = n_total * sin(phi) + Y/(m*g)
```

其中侧力项 `Y` 只保留为由侧滑引起的弱修正项。

---

### 4.2 外环

给定时变法向、侧向过载指令：

```matlab
n_n_cmd(t)
n_l_cmd(t)
```

先构造：

```matlab
n_total_cmd = sqrt(n_n_cmd^2 + n_l_cmd^2)
phi_cmd = atan2(n_l_cmd, n_n_cmd)
alpha_cmd = n_total_cmd / nn_alpha_gain
beta_cmd = 0
```

然后得到角速度指令：

```matlab
p_cmd = k_phi_p*(phi_cmd - phi) - k_phi_d*p
q_cmd = k_alpha_q*(alpha_cmd - alpha)
r_cmd = k_beta_r*(beta_cmd - beta)
```

这样就形成了：

- `phi` 控制侧向机动方向
- `alpha` 控制总过载大小
- `beta` 保持协调转弯

---

### 4.3 内环 INDI

角速度内环采用：

```matlab
nu = omega_cmd_dot + Komega*(omega_cmd - omega)
Delta_ua = pinv(Ba) * J * (nu - omega_dot_meas_prev)
ua = ua_prev + Delta_ua
```

其中：

- `ua = [ur; up; uy]`
- `Ba` 为姿态模态控制效能矩阵
- `omega_dot_meas_prev` 由角速度差分加低通得到

这保持了 INDI 的标准增量控制逻辑。

---

## 5. 当前结果为什么更合理

在 MATLAB 实跑后，当前结果大致为：

- 法向过载峰值误差约 `0.206 g`
- 侧向过载峰值误差约 `0.633 g`
- 最大舵偏约 `25 deg`
- 末端速度约 `180 m/s`
- 末端俯仰角约 `25 deg`
- 末端航迹角约 `8 deg`

这比之前“俯仰角飘到 80 度”的结果合理得多，原因是：

1. `theta` 现在受平动与转动共同约束
2. `gamma` 不再被姿态角代替
3. `alpha` 通过 `q - gamma_dot` 动态与航迹变化耦合
4. BTT 侧向机动主要由滚转建立，而不是靠不合理的姿态漂移

---

## 6. 绘图说明

脚本按如下顺序绘图，并采用子图方式组织：

1. 弹道曲线
2. 速度曲线
3. 姿态曲线
4. 过载跟踪曲线
5. 角度跟踪曲线
6. 舵偏角曲线
7. 角速度曲线

这样可以把：

- 平动结果
- 姿态结果
- 过载结果
- 控制输入

分开观察，避免信息堆在一张图中。

---

## 7. 当前版本的简化边界

虽然动力学结构已经更完整，但当前版本仍然是控制律验证模型，不是高保真工程模型。主要简化包括：

1. 气动系数仍是低阶简化形式
2. `alpha`、`beta` 动态仍是近似模型
3. 未加入舵机动态和舵速限制
4. 未加入对称增阻模态
5. 未加入时间控制与能量管理

---

## 8. 后续扩展建议

如果继续扩展，可以按下面顺序推进：

1. 增加对称增阻模态 `ud`
   扩展到姿态通道 + 阻力通道联合控制。

2. 将常值 `Ba` 扩展为状态相关 `Ba(x)`
   提高 INDI 局部逆的准确性。

3. 增加舵机动态与舵速限制
   提高工程真实性。

4. 将当前平动模型继续扩展为更完整的气动/轨迹模型
   用于时间约束和能量管理研究。

