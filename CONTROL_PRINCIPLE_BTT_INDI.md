# 基于 INDI 的滑翔飞行器 BTT 姿态与阻力联合控制方法说明

## 1. 文档对象与方法目标

本文档对应脚本：

- `indi_attitude_drag_control_sim.m`

原始文件 `indi_attitude_control_sim.m` 保持为“仅姿态控制”的基线版本；本文档说明的是在此基础上扩展得到的“姿态 + 阻力联合控制”方法。

该方法面向无动力滑翔飞行器，目标为：

- 跟踪恒定法向过载指令 `n_n_cmd`
- 跟踪恒定速度指令 `V_cmd`
- 保持侧向过载指令 `n_l_cmd = 0`
- 保持协调转弯约束 `beta_cmd = 0`

总体结构为：

- 外环 1：BTT 姿态外环，负责“过载大小 + 过载方向”
- 外环 2：速度/阻力外环，负责“对称增阻模态”
- 内环：角速度 INDI
- 执行机构：4 个 X 布局尾舵

---

## 2. 状态与动力学模型

### 2.1 状态定义

脚本采用如下状态：

```math
x =
\begin{bmatrix}
x & y & z & V & \gamma & \chi & \phi & \theta & \psi & p & q & r & \alpha & \beta
\end{bmatrix}^T
```

其中：

- `x, y, z` 为位置
- `V` 为速度大小
- `gamma` 为航迹倾角
- `chi` 为航向角
- `phi, theta, psi` 为欧拉角
- `p, q, r` 为机体角速度
- `alpha, beta` 为攻角与侧滑角

---

### 2.2 平动动力学

在速度坐标系下，平动动力学写为：

```math
\dot{x} = V \cos\gamma \cos\chi
```

```math
\dot{y} = V \cos\gamma \sin\chi
```

```math
\dot{z} = V \sin\gamma
```

```math
\dot{V} = -\frac{D}{m} - g\sin\gamma
```

```math
\dot{\gamma} =
\frac{L\cos\phi / m - g\cos\gamma}{V}
```

```math
\dot{\chi} =
\frac{L\sin\phi}{m V \cos\gamma}
```

这里：

- `L` 为总升力
- `D` 为总阻力
- `m` 为飞行器质量
- `g` 为重力加速度

代码中对 `V` 和 `cos(gamma)` 都加入了下界保护，以避免低速或大航迹角时的数值发散。

---

### 2.3 欧拉角运动学

欧拉角运动学采用标准形式：

```math
\dot{\eta} = T(\phi,\theta)\,\omega
```

其中：

```math
\eta = [\phi,\theta,\psi]^T,\qquad
\omega = [p,q,r]^T
```

```math
T(\phi,\theta)=
\begin{bmatrix}
1 & \sin\phi\tan\theta & \cos\phi\tan\theta \\
0 & \cos\phi & -\sin\phi \\
0 & \sin\phi/\cos\theta & \cos\phi/\cos\theta
\end{bmatrix}
```

代码中同样对 `cos(theta)` 做了下界保护。

---

### 2.4 转动动力学

刚体角运动方程为：

```math
J\dot{\omega} = M_{aero} + M_{ctrl} - \omega \times (J\omega)
```

其中：

- `J = diag(Jx, Jy, Jz)` 为转动惯量矩阵
- `M_aero` 为气动力矩
- `M_ctrl` 为舵面控制力矩

---

### 2.5 `alpha / beta` 低阶耦合动态

当前脚本中，`alpha` 与 `beta` 仍采用低阶近似动态，但显式与平动通道耦合：

```math
\dot{\alpha} = q - \dot{\gamma} - a_{\alpha}\alpha + b_{\alpha}u_p
```

```math
\dot{\beta} = r - \dot{\chi} - a_{\beta}\beta + b_{\beta}u_y
```

其中：

- `u_p` 为俯仰姿态模态输入
- `u_y` 为偏航姿态模态输入

这种写法的含义是：

- `alpha` 反映机体俯仰相对速度方向的偏差，因此减去 `gamma_dot`
- `beta` 反映机体航向相对速度方向的偏差，因此减去 `chi_dot`

---

## 3. 气动与过载模型

### 3.1 气动系数模型

脚本中的气动系数采用低阶解析模型：

```math
C_L = C_{L0} + C_{L\alpha}\alpha
```

```math
C_{D,ctrl} = C_{D,u_d}\,\max(u_d,0)
```

```math
C_D = C_{D0} + K_{induced} C_L^2 + C_{D\beta}\beta^2 + C_{D,ctrl}
```

```math
C_Y = C_{Y\beta}\beta
```

其中：

- `u_d` 为对称增阻模态
- `C_{D,ctrl}` 为阻力控制通道额外引入的阻力

对应气动力为：

```math
q = \frac{1}{2}\rho V^2
```

```math
L = q S C_L
```

```math
D = q S C_D
```

```math
Y = q S C_Y
```

---

### 3.2 BTT 过载分解

总过载定义为：

```math
n_{tot} = \frac{L}{mg}
```

法向过载与侧向过载分解为：

```math
n_n = n_{tot}\cos\phi
```

```math
n_l = n_{tot}\sin\phi + \frac{Y}{mg}
```

因此：

- `alpha` 主要决定总升力大小，也即 `n_tot`
- `phi` 主要决定总升力在法向/侧向方向上的分配
- `beta` 主要作为协调转弯约束项，而不是主要机动来源

---

## 4. 气动力矩模型与舵面力矩模型

### 4.1 气动力矩模型

脚本采用线性阻尼 + 静稳定项构造气动力矩：

```math
M_{aero} =
\begin{bmatrix}
-L_p p - L_{\beta}\beta \\
-M_q q - M_{\alpha}\alpha \\
-N_r r - N_{\beta}\beta
\end{bmatrix}
```

含义是：

- 滚转通道主要由滚转角速度阻尼 `L_p p` 和侧滑耦合项 `L_beta beta` 构成
- 俯仰通道主要由俯仰角速度阻尼 `M_q q` 和攻角静稳定项 `M_alpha alpha` 构成
- 偏航通道主要由偏航角速度阻尼 `N_r r` 和侧滑静稳定项 `N_beta beta` 构成

---

### 4.2 X 尾翼姿态模态

定义姿态模态输入为：

```math
u_a =
\begin{bmatrix}
u_r \\ u_p \\ u_y
\end{bmatrix}
```

4 个舵面的姿态差动部分由矩阵 `T_a` 映射：

```math
\delta_{att} = T_a u_a
```

其中：

```math
T_a =
\begin{bmatrix}
1 & 1 & -1 \\
-1 & 1 & 1 \\
1 & -1 & 1 \\
-1 & -1 & -1
\end{bmatrix}
```

展开后即：

```math
\delta_1 = u_r + u_p - u_y
```

```math
\delta_2 = -u_r + u_p + u_y
```

```math
\delta_3 = u_r - u_p + u_y
```

```math
\delta_4 = -u_r - u_p - u_y
```

---

### 4.3 舵面力矩模型

舵面对三轴控制力矩的作用由矩阵 `B4` 给出：

```math
M_{ctrl} = B_4 \delta
```

其中：

```math
B_4 =
\begin{bmatrix}
K_r & -K_r & K_r & -K_r \\
K_p & K_p & -K_p & -K_p \\
-K_y & K_y & K_y & -K_y
\end{bmatrix}
```

姿态模态下的等效控制效能矩阵为：

```math
B_a = B_4 T_a
```

因此：

```math
M_{ctrl} = B_a u_a
```

这就是 INDI 内环里使用的局部控制效能矩阵。

---

### 4.4 对称增阻模态

在姿态差动模态之外，再引入对称增阻模态：

```math
\delta_d = T_d u_d
```

其中：

```math
T_d =
\begin{bmatrix}
1 \\ 1 \\ 1 \\ 1
\end{bmatrix}
```

即 4 个舵面同向偏转，用于增加阻力而不直接产生期望姿态力矩。

最终舵偏角为：

```math
\delta = \delta_{att} + \delta_d = T_a u_a + T_d u_d
```

由于 `T_d` 为纯对称模态，在理想线性力矩模型下，其主要作用是增阻而不是提供净姿态力矩。

---

## 5. 外环控制律

### 5.1 姿态外环：法向过载与滚转分配

测试输入为恒值：

```math
n_{n,cmd} = \text{const},\qquad
n_{l,cmd} = \text{const}
```

总过载指令构造为：

```math
n_{tot,cmd} = \sqrt{\max(n_{n,cmd},0)^2 + n_{l,cmd}^2}
```

滚转角指令为：

```math
\phi_{cmd} = \operatorname{atan2}(n_{l,cmd}, n_{n,cmd})
```

攻角指令由升力近似反求：

```math
\alpha_{cmd} = \frac{n_{tot,cmd}}{k_{n\alpha}}
```

其中：

```math
k_{n\alpha} = \frac{qSC_{L\alpha}}{mg}
```

并对 `alpha_cmd` 与 `phi_cmd` 做限幅。

---

### 5.2 协调转弯约束

当前脚本中：

```math
\beta_{cmd} = 0
```

于是角速度指令写为：

```math
p_{cmd} = k_{\phi p}(\phi_{cmd} - \phi) - k_{\phi d}p
```

```math
q_{cmd} = k_{\alpha q}(\alpha_{cmd} - \alpha)
```

```math
r_{cmd} = k_{\beta r}(\beta_{cmd} - \beta)
```

这里没有额外加入偏航率前馈，偏航通道的主要作用是压低侧滑角。

---

### 5.3 速度 / 阻力外环

阻力外环目标是通过 `u_d` 维持速度跟踪。

首先，脚本根据参考平衡点构造初始配平：

```math
q_{ref} = \frac{1}{2}\rho V_{cmd}^2
```

```math
C_{L,ref} = \frac{n_{n,cmd}mg}{q_{ref}S}
```

```math
\alpha_{ref} = \frac{C_{L,ref} - C_{L0}}{C_{L\alpha}}
```

```math
\gamma_{ref} = -\arccos(n_{n,cmd})
```

```math
C_{D,base,ref} = C_{D0} + K_{induced} C_{L,ref}^2
```

```math
C_{D,req,ref} =
\frac{mg\max(-\sin\gamma_{ref},0)}{q_{ref}S}
```

```math
u_{d,trim} =
\operatorname{sat}\left(
\frac{C_{D,req,ref} - C_{D,base,ref}}{C_{D,u_d}}
\right)
```

这一步的目的是在仿真初始时给出一个与目标速度相一致的增阻配平值。

---

### 5.4 在线阻力前馈与速度反馈

在每个控制周期，脚本根据当前状态重算阻力前馈：

```math
C_{L,now} = C_{L0} + C_{L\alpha}\alpha
```

```math
C_{D,base,now} = C_{D0} + K_{induced} C_{L,now}^2 + C_{D\beta}\beta^2
```

```math
C_{D,req,now} =
\frac{mg\max(-\sin\gamma,0)}{qS}
```

```math
u_{d,ff} =
\operatorname{sat}\left(
\frac{C_{D,req,now} - C_{D,base,now}}{C_{D,u_d}}
\right)
```

再叠加速度误差反馈：

```math
u_{d,ref} = \operatorname{sat}\left(u_{d,ff} + k_{V,u_d}(V - V_{cmd})\right)
```

最后通过一阶滤波得到平滑指令：

```math
u_{d,cmd}(k) =
u_{d}(k-1) +
\frac{\Delta t}{\tau_{u_d}}
\left(u_{d,ref}(k) - u_{d}(k-1)\right)
```

因此该通道本质上是：

- 用在线前馈提供当前滑翔工况所需的基本阻力
- 用速度误差反馈修正剩余误差

---

## 6. INDI 内环控制律

### 6.1 期望角加速度

定义角速度指令：

```math
\omega_{cmd} =
\begin{bmatrix}
p_{cmd} \\ q_{cmd} \\ r_{cmd}
\end{bmatrix}
```

期望角加速度取为：

```math
\nu = \dot{\omega}_{cmd} + K_{\omega}(\omega_{cmd} - \omega)
```

其中：

- `\dot{\omega}_{cmd}` 由指令差分得到
- `K_omega` 为角速度误差反馈增益矩阵

---

### 6.2 角加速度估计

角加速度估计采用“差分 + 一阶低通”：

```math
\dot{\omega}_{raw} \approx \frac{\omega_k - \omega_{k-1}}{\Delta t}
```

```math
\dot{\omega}_{meas} =
\dot{\omega}_{meas,prev} +
\frac{\Delta t}{\tau_{\dot{\omega}}}
(\dot{\omega}_{raw} - \dot{\omega}_{meas,prev})
```

---

### 6.3 INDI 增量控制律

在采样间隔内做局部增量近似：

```math
J(\dot{\omega}_k - \dot{\omega}_{k-1}) \approx B_a \Delta u_a
```

则姿态模态增量为：

```math
\Delta u_a \approx B_a^\dagger J(\nu - \dot{\omega}_{meas})
```

因此：

```math
u_{a,k} = u_{a,k-1} + \Delta u_a
```

这就是脚本中的 INDI 内环核心公式。

---

## 7. 舵面分配与姿态优先策略

### 7.1 姿态模态先分配

先根据姿态模态得到差动舵偏：

```math
\delta_{att} = T_a u_a
```

若姿态差动部分已超过单舵限幅 `\delta_{max}`，则整体等比例缩放：

```math
\delta_{att} \leftarrow
\delta_{att}\frac{\delta_{max}}{\max|\delta_{att}|}
```

---

### 7.2 对称增阻仅占用剩余裕度

对称增阻模态不允许破坏姿态控制优先级，因此仅使用“正向剩余裕度”：

```math
u_{d,margin} = \min(\delta_{max} - \delta_{att,i})
```

```math
u_d = \operatorname{sat}\left(u_{d,cmd}, 0, \min(u_{d,max}, u_{d,margin})\right)
```

最终舵偏角为：

```math
\delta = \delta_{att} + T_d u_d
```

等效姿态模态反馈值写为：

```math
u_{a,eff} = T_a^\dagger \delta_{att}
```

这意味着：

- 姿态通道始终优先满足
- 阻力通道只能使用姿态通道未占用的舵偏空间
- 即使阻力控制有需求，也不会主动挤占姿态控制所需控制量

---

## 8. 仿真测试项与结果指标

当前测试项为：

- 恒定法向过载指令 `n_n_cmd_const`
- 恒定侧向过载指令 `n_l_cmd_const = 0`
- 恒定速度指令 `V_cmd_const`

脚本记录并输出的关键指标包括：

- 法向过载峰值误差
- 侧向过载峰值误差
- 速度峰值误差
- 最大舵偏角
- 最大对称增阻模态
- 末端速度
- `cond(Ba)`

这些指标分别反映：

- 姿态控制精度
- 速度 / 阻力控制精度
- 舵面使用程度
- 姿态模态控制效能矩阵的数值条件

---

## 9. 方法特点与当前简化

### 9.1 当前方法的特点

当前方法的核心特点是：

- 采用 BTT 思想，将“过载大小”和“过载方向”分开处理
- 用 INDI 构造姿态角速度内环，降低对完整精确模型的依赖
- 在 X 尾翼上增加对称增阻模态，实现姿态与阻力联合控制
- 用“姿态优先、阻力占剩余裕度”的策略避免通道冲突

---

### 9.2 当前模型的主要简化

当前版本仍然是控制律验证模型，而非高保真工程模型。主要简化包括：

- `alpha / beta` 仍为低阶近似动态
- 气动系数仍为低阶解析模型
- 未引入舵机动态与舵速限制
- `B_a` 未显式建模为状态相关矩阵 `B_a(x)`
- 未引入时间约束控制与能量管理逻辑

---

## 10. 后续可扩展方向

后续若继续扩展，可沿以下方向推进：

1. 将 `B_a` 扩展为随工况变化的在线控制效能矩阵。
2. 在阻力通道中引入更真实的舵面阻力模型和非线性气动模型。
3. 加入舵机动态、舵速限制和舵面饱和恢复逻辑。
4. 将速度 / 阻力通道进一步扩展为时间约束控制或能量管理通道。
5. 将当前低阶 `alpha / beta` 模型替换为更完整的机体系速度动力学模型。
