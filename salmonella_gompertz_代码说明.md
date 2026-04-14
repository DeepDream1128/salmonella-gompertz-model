# `salmonella_gompertz.m` 代码说明

本文档说明主脚本 `salmonella_gompertz.m` 的数学模型、数据流、参数化策略与统计输出，便于复现与阅读代码。

---

## 1. 研究问题与目标

脚本在**正弦温度历程**下，用 **Gompertz 型生长**描述沙门氏菌数量（以 \(\log_{10} N\) 表示），通过实验数据估计模型参数，并进行：

- 可辨识性分析（缩放灵敏度系数 SSC）
- 非线性最小二乘拟合（`nlinfit`）
- 渐近置信/预测区间与自助法（bootstrap）区间对比
- 简单最优试验设计（\(\Delta\) 准则：\(\det(\mathbf{X}^\top\mathbf{X})\)）

---

## 2. 数学模型

### 2.1 主模型（微分形式）

代码注释给出的形式为：

\[
\frac{d(\log N)}{dt} = \mu\, K\, C\, e^{-K},
\quad K = e^{-\mu(t-M)}.
\]

状态变量为 \(y = \log_{10} N\)。\(K\) 在 `gompODE` 中**直接由当前 \(t\) 与 \(\mu, M\) 计算**，不作为单独状态积分。

### 2.2 次级模型（比生长速率 \(\mu\) 与温度 \(T\)）

\[
\mu = a\,(T-T_{\min})^2\,\bigl(1 - e^{b\,(T-T_{\max})}\bigr),
\quad T_{\min} < T < T_{\max},
\]

在 `secModel` 中：若 \(T \le T_{\min}\) 或 \(T \ge T_{\max}\)，则 \(\mu = 0\)。

温度 \(T(t)\) 来自全局变量 `tTemp`：列为 \([\)时间 (h), 温度 (°C)\(]\)，在 ODE 右端用 `interp1` 线性插值（端点外推）。

### 2.3 数值实现

- 前向问题：`ode45` 在**唯一时间点**上积分，再通过 `unique` 的索引 `ic` 映射回原始 `t` 向量（支持重复时刻、多观测点）。
- 初值：`y0 = A`，其中 \(A\) 为 \(\log_{10} N\) 在 \(t=0\) 附近的参数化初值。

---

## 3. 数据输入

| 文件 | 内容 |
|------|------|
| `Salmonella sin growth.xlsx` | 第1列时间 (h)，第2、3列为两次重复的 CFU/mL |
| `Salmonella sin growth Temps.xlsx` | 第1列时间 (min)，第2列温度 (°C) → 脚本中换算为小时 |

生长数据：对每个有效、正的 CFU，取 \(\log_{10}\)，时间与观测一一展开，得到向量 `x`（时间）与 `yobs`（\(\log_{10}N\)）。

---

## 4. 参数向量与“多轮”估计策略

完整参数集为 7 个：

\[
\boldsymbol{\beta} = [A,\, C,\, M,\, a,\, b,\, T_{\min},\, T_{\max}].
\]

### 4.1 Round 0（可选）：7 参数同时估计

用于对比；若 `nlinfit` 失败会 `catch` 并打印信息。

### 4.2 SSC 分析与固定参数

- `SSC_V3` 计算**缩放灵敏度** \(X'_j \approx \beta_j \,\partial y/\partial \beta_j\)（前向差分，相对扰动 `d=0.001`）。
- 脚本根据各参数 \(\max |SSC|\) 及参数对 SSC **比值是否近似常数**（红色背景提示近似共线），说明 **\(T_{\min}\)** 与 **\(b\)** 的 SSC 相对较小，故在后续轮次中**固定为文献/初值**。

### 4.3 Round 1：5 参数 \([A,C,M,a,T_{\max}]\)

`Tmin_fixed`、`b_fixed` 固定，`gompertzFOR_5` / `gompertzINV_5` 调用同一 `gompODE`。

### 4.4 Round 2：4 参数 \([A,C,M,T_{\max}]\)

将 Round 1 得到的 **\(a\)** 固定为 `a_fixed`，估计 **\(A,C,M,T_{\max}\)**。这是脚本中**主报告结果**所用的参数化（`beta`、`fnameFOR`、`fnameINV`）。

**设计动机**：降低参数间相关与 Jacobian 病态，使 `cond(J)` 可接受、\(T_{\max}\) 等参数更可辨识。

---

## 5. 统计输出（Round 2）

- **拟合**：`nlinfit` 得到 \(\hat{\boldsymbol{\beta}}\)、残差、`J`、协方差 `COVB`、`mse`。
- **优度**：\(R^2\)、调整 \(R^2\)、RMSE、相对 RMSE 等。
- **参数不确定度**：`nlparci` 渐近 95% 区间；`corrcov` 得相关矩阵与标准误。
- **曲线与观测**：`nlpredci` 分别计算**回归线**的置信带（`'curve'`）与**观测**的预测带（`'observation'`）。
- **残差诊断**：残差图、直方图、基于排序残差的 **runs** 统计（与最小期望 runs 比较）。
- **Bootstrap**：对残差有放回抽样构造 \(y_{\mathrm{boot}}\)，重复拟合（`optsBoot` 放宽 `TolFun`/`TolX` 以减少迭代失败）；汇总参数分位数区间；并绘制 bootstrap 意义下的 CB/PB，与渐近带对比宽度。

---

## 6. 最优试验设计（脚本末尾）

对 **等间距** 时间点数目 `npts_vec = 5:50`：

1. 在 `tk` 上算 `SSC_V3` 得到缩放灵敏度矩阵列。
2. 按注释：**未缩放 Jacobian** 满足 \(J_{:,j} = X_{\mathrm{SSC},j}/\beta_j\)。
3. 计算 \(\det(\mathbf{J}^\top\mathbf{J})\) 及 \((\mathbf{J}^\top\mathbf{J})^{-1}\) 的对角元 \(C_{ii}\)（参数方差相关），作图作为 \(\Delta\) 准则与局部精度参考。

---

## 7. 图形输出

所有图保存到工程目录下 **`figures/`**（若不存在则创建），文件名包括：

- 初值与数据、SSC 曲线、SSC 比值热图式子图
- 各轮拟合、双轴（logN + 温度）、残差、OED、bootstrap 带等

---

## 8. 依赖与运行注意

- 需要 Statistics and Machine Learning Toolbox（`nlinfit`、`nlparci`、`nlpredci` 等）。
- **全局变量** `tTemp` 在 `gompODE` 中读取；运行前须先读入温度表并赋值（脚本已做）。
- 修改 Excel 路径或列布局时需同步改 `readmatrix` 与列索引说明。

---

## 9. 文件内函数一览（文末局部函数）

| 函数 | 作用 |
|------|------|
| `gompertzFOR_7` | 7 参数前向 \(\log N(t)\) |
| `gompertzFOR_5` / `gompertzINV_5` | 5 参数（固定 \(T_{\min},b\)） |
| `gompertzFOR_4b` / `gompertzINV_4b` | 4 参数（再固定 \(a\)） |
| `gompODE` | ODE 右端项 |
| `secModel` | \(\mu(T)\) |
| `SSC_V3` | 缩放灵敏度数值近似 |

前向与“逆”函数在本脚本中**形式相同**，均用于 `nlinfit` 预测 \(\log_{10}N\)；命名上的 INV 表示在非线性回归中作为用户提供的模型函数句柄。

---

## 10. 阅读代码的建议顺序

1. 开头注释与数据读入（第 1–47 行）
2. 初值与 `gompertzFOR_7` + `SSC_V3` 分析（第 48–145 行）
3. Round 0 → Round 1 → Round 2 的 `nlinfit` 与诊断（第 146–365 行）
4. 作图与残差、runs（第 366–458 行）
5. OED 与 bootstrap（第 492–614 行）
6. 文件末尾局部函数定义（第 618–704 行）

如需与另一脚本 `inv_soln5.m` 或 `global_example.m` 对照，可关注 OED 部分注释中“与 inv_soln5.m 一致”的等间距设计与 Jacobian 换算约定。
