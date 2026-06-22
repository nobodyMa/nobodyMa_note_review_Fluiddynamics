### 基本定义

1.  **笛卡尔坐标系**：记为 $\mathbf{x} = (x^1, x^2, x^3) = (x, y, z)$。这是我们的基准坐标系。
2.  **曲线坐标系**：记为 $\mathbf{u} = (u^1, u^2, u^3)$。存在一个光滑、可逆的坐标变换：
    $$ \mathbf{x} = \mathbf{x}(\mathbf{u}) \quad \text{或} \quad x^i = x^i(u^1, u^2, u^3), \quad i=1,2,3 $$
    我们要求雅可比行列式 $J = \det(\mathbf{J}) \neq 0$ 处处成立，以保证变换的可逆性。
3.  **雅可比矩阵 (Jacobian Matrix)**：
    $$ \mathbf{J} = \frac{\partial \mathbf{x}}{\partial \mathbf{u}} = \begin{bmatrix} \frac{\partial x^1}{\partial u^1} & \frac{\partial x^1}{\partial u^2} & \frac{\partial x^1}{\partial u^3} \\ \frac{\partial x^2}{\partial u^1} & \frac{\partial x^2}{\partial u^2} & \frac{\partial x^2}{\partial u^3} \\ \frac{\partial x^3}{\partial u^1} & \frac{\partial x^3}{\partial u^2} & \frac{\partial x^3}{\partial u^3} \end{bmatrix} $$
    其逆矩阵为 $\mathbf{J}^{-1} = \frac{\partial \mathbf{u}}{\partial \mathbf{x}}$。
4.  **局部协变基向量**：
    $$ \mathbf{g}_i = \frac{\partial \mathbf{x}}{\partial u^i} = \left( \frac{\partial x^1}{\partial u^i}, \frac{\partial x^2}{\partial u^i}, \frac{\partial x^3}{\partial u^i} \right)^T, \quad i=1,2,3 $$
    这些向量是雅可比矩阵的列向量。它们构成了曲线坐标系中任意一点切空间的一组基。
5.  **度规张量**：欧氏空间的度规在曲线坐标下体现为：
    $$ g_{ij} = \mathbf{g}_i \cdot \mathbf{g}_j = \sum_{k=1}^3 \frac{\partial x^k}{\partial u^i} \frac{\partial x^k}{\partial u^j} $$
    度规张量的行列式 $g = \det(g_{ij})$ 与雅可比行列式 $J$ 的关系为 $g = J^2$。
    对于**正交曲线坐标系**，有 $g_{ij} = 0 \, (i \neq j)$，此时定义**拉梅系数**：
    $$ h_i = \sqrt{g_{ii}} = |\mathbf{g}_i| $$
    并且 $ \sqrt{g} = h_1 h_2 h_3 $。
6.  **局部标准正交基（物理基）**：
    $$ \hat{\mathbf{e}}_i = \frac{\mathbf{g}_i}{h_i} = \frac{\mathbf{g}_i}{|\mathbf{g}_i|}, \quad i=1,2,3 $$

---

### 1. 梯度 (Gradient) 的推导

**目标**：将 $\nabla f = \left( \frac{\partial f}{\partial x}, \frac{\partial f}{\partial y}, \frac{\partial f}{\partial z} \right)$ 用 $u^i$ 坐标表示。

**推导**：
标量场 $f$ 的梯度是一个向量，其分量在坐标变换下必须满足一定的变换规则。我们可以通过链式法则和基向量的对偶关系来找到它。

最直接的方法是使用链式法则：
$$ \frac{\partial f}{\partial x^j} = \sum_{k=1}^3 \frac{\partial f}{\partial u^k} \frac{\partial u^k}{\partial x^j} $$
这意味着梯度算子在坐标变换下为：
$$ \frac{\partial}{\partial x^j} = \sum_{k=1}^3 \frac{\partial u^k}{\partial x^j} \frac{\partial}{\partial u^k} $$

然而，梯度向量 $\nabla f$ 本身是一个不依赖于坐标的几何实体。为了在曲线坐标下找到它的**分量**，我们将其表示为局部标准正交基 $\{ \hat{\mathbf{e}}_i \}$ 的线性组合：
$$ \nabla f = \sum_{i=1}^3 (\nabla f)^i \hat{\mathbf{e}}_i $$
我们需要找到分量 $(\nabla f)^i$。

考虑 $f$ 沿某个方向的方向导数。沿 $\hat{\mathbf{e}}_i$ 方向的一个无穷小位移 $d\mathbf{l}$ 可以表示为：
$$ d\mathbf{l} = h_i du^i \hat{\mathbf{e}}_i $$
（因为 $ds^2 = \sum_i (h_i du^i)^2$）。

该方向的方向导数为：
$$ \frac{\partial f}{\partial l_i} = \lim_{\Delta l_i \to 0} \frac{f(\mathbf{u} + \Delta u^i \mathbf{e}_i) - f(\mathbf{u})}{\Delta l_i} = \frac{1}{h_i} \frac{\partial f}{\partial u^i} $$

而这个方向导数正好是梯度 $\nabla f$ 在 $\hat{\mathbf{e}}_i$ 方向上的投影：
$$ \frac{\partial f}{\partial l_i} = \hat{\mathbf{e}}_i \cdot \nabla f = (\nabla f)^i $$

因此，我们得到：
$$ (\nabla f)^i = \frac{1}{h_i} \frac{\partial f}{\partial u^i} $$

**最终公式（正交系）**：
$$ \boxed{\nabla f = \sum_{i=1}^3 \frac{1}{h_i} \frac{\partial f}{\partial u^i} \hat{\mathbf{e}}_i = \frac{1}{h_1} \frac{\partial f}{\partial u^1} \hat{\mathbf{e}}_1 + \frac{1}{h_2} \frac{\partial f}{\partial u^2} \hat{\mathbf{e}}_2 + \frac{1}{h_3} \frac{\partial f}{\partial u^3} \hat{\mathbf{e}}_3} $$

---

### 2. 散度 (Divergence) 的推导

**目标**：将 $\nabla \cdot \mathbf{F} = \frac{\partial F_x}{\partial x} + \frac{\partial F_y}{\partial y} + \frac{\partial F_z}{\partial z}$ 用 $u^i$ 坐标表示。

**推导**：
散度公式最严谨的推导来自于**散度定理**和**体积元的表达式**。

1.  **曲线坐标下的体积元**：
    在曲线坐标下，一个无穷小的体积元 $dV$ 为：
    $$ dV = |dx^1 dx^2 dx^3| = |\det(\mathbf{J})| \, du^1 du^2 du^3 = J \, du^1 du^2 du^3 = \sqrt{g} \, du^1 du^2 du^3 $$
    对于正交系，$\sqrt{g} = h_1 h_2 h_3$，所以：
    $$ dV = h_1 h_2 h_3 \, du^1 du^2 du^3 $$

2.  **散度定理的积分形式**：
    $$ \iiint_V (\nabla \cdot \mathbf{F}) \, dV = \oiint_S \mathbf{F} \cdot d\mathbf{S} $$
    我们考虑一个由坐标面围成的无穷小平行六面体，其体积 $dV = h_1 h_2 h_3 \, du^1 du^2 du^3$。

3.  **计算通量**：
    计算向量场 $\mathbf{F} = F_1 \hat{\mathbf{e}}_1 + F_2 \hat{\mathbf{e}}_2 + F_3 \hat{\mathbf{e}}_3$ 穿过这个六面体六个面的净通量。
    *   **垂直于 $u^1$ 方向的两个面**：
        面积： $dS_1 = h_2 h_3 \, du^2 du^3$
        通量差： $ \frac{\partial}{\partial u^1}(F_1 \, h_2 h_3) du^1 du^2 du^3 $
    *   **同理**，垂直于 $u^2$ 和 $u^3$ 方向的通量差分别为：
        $ \frac{\partial}{\partial u^2}(F_2 \, h_3 h_1) du^1 du^2 du^3 $
        $ \frac{\partial}{\partial u^3}(F_3 \, h_1 h_2) du^1 du^2 du^3 $

4.  **总通量与散度**：
    总净通量为三者之和：
    $$ \oiint_S \mathbf{F} \cdot d\mathbf{S} = \left[ \frac{\partial}{\partial u^1}(F_1 h_2 h_3) + \frac{\partial}{\partial u^2}(F_2 h_3 h_1) + \frac{\partial}{\partial u^3}(F_3 h_1 h_2) \right] du^1 du^2 du^3 $$
    根据散度定理，这等于 $(\nabla \cdot \mathbf{F}) \, dV = (\nabla \cdot \mathbf{F}) \, h_1 h_2 h_3 \, du^1 du^2 du^3$。

    令两者相等：
    $$ (\nabla \cdot \mathbf{F}) \, h_1 h_2 h_3 \, du^1 du^2 du^3 = \left[ \frac{\partial}{\partial u^1}(F_1 h_2 h_3) + \ldots \right] du^1 du^2 du^3 $$

**最终公式（正交系）**：
$$ \boxed{\nabla \cdot \mathbf{F} = \frac{1}{h_1 h_2 h_3} \left[ \frac{\partial}{\partial u^1}(F_1 h_2 h_3) + \frac{\partial}{\partial u^2}(F_2 h_3 h_1) + \frac{\partial}{\partial u^3}(F_3 h_1 h_2) \right]} $$

---

### 3. 旋度 (Curl) 的推导

**目标**：将 $\nabla \times \mathbf{F}$ 用 $u^i$ 坐标表示。

**推导**：
旋度公式最严谨的推导来自于**斯托克斯定理**和**面积元的表达式**。

1.  **曲线坐标下的面积元**：
    一个垂直于 $\hat{\mathbf{e}}_1$ 的无穷小矩形面积元为：
    $$ d\mathbf{S}_1 = (h_2 du^2) (h_3 du^3) \hat{\mathbf{e}}_1 = h_2 h_3 du^2 du^3 \hat{\mathbf{e}}_1 $$
    同理，$d\mathbf{S}_2 = h_3 h_1 du^3 du^1 \hat{\mathbf{e}}_2$, $d\mathbf{S}_3 = h_1 h_2 du^1 du^2 \hat{\mathbf{e}}_3$。

2.  **斯托克斯定理的积分形式**：
    $$ \iint_S (\nabla \times \mathbf{F}) \cdot d\mathbf{S} = \oint_C \mathbf{F} \cdot d\mathbf{l} $$
    我们选择一个垂直于 $\hat{\mathbf{e}}_1$ 的无穷小矩形回路 $C$。

3.  **计算环量**：
    计算 $\mathbf{F} \cdot d\mathbf{l}$ 沿该回路的环量。
    *   沿 $u^2$ 增加方向的边： $F_2 h_2 du^2$
    *   沿 $u^3$ 增加方向的边： $F_3 h_3 du^3$
    *   考虑方向和导数，总环量为：
        $$ \oint_C \mathbf{F} \cdot d\mathbf{l} \approx \left[ \frac{\partial}{\partial u^3}(F_2 h_2) - \frac{\partial}{\partial u^2}(F_3 h_3) \right] du^2 du^3 $$

4.  **旋度与环量**：
    根据斯托克斯定理，这等于 $(\nabla \times \mathbf{F}) \cdot d\mathbf{S}_1 = (\nabla \times \mathbf{F})_1 \, h_2 h_3 du^2 du^3$。
    令两者相等：
    $$ (\nabla \times \mathbf{F})_1 \, h_2 h_3 du^2 du^3 = \left[ \frac{\partial}{\partial u^3}(F_2 h_2) - \frac{\partial}{\partial u^2}(F_3 h_3) \right] du^2 du^3 $$
    解得：
    $$ (\nabla \times \mathbf{F})_1 = \frac{1}{h_2 h_3} \left[ \frac{\partial}{\partial u^3}(F_2 h_2) - \frac{\partial}{\partial u^2}(F_3 h_3) \right] $$
   
    同理可得 $(\nabla \times \mathbf{F})_2$ 和 $(\nabla \times \mathbf{F})_3$。

**最终公式（正交系）**：
旋度的分量由以下行列式给出：
$$ \boxed{\nabla \times \mathbf{F} = \frac{1}{h_1 h_2 h_3} \begin{vmatrix} h_1 \hat{\mathbf{e}}_1 & h_2 \hat{\mathbf{e}}_2 & h_3 \hat{\mathbf{e}}_3 \\ \frac{\partial}{\partial u^1} & \frac{\partial}{\partial u^2} & \frac{\partial}{\partial u^3} \\ h_1 F_1 & h_2 F_2 & h_3 F_3 \end{vmatrix}} $$
展开后为：
$$ \nabla \times \mathbf{F} = \frac{1}{h_2 h_3} \left[ \frac{\partial}{\partial u^2}(h_3 F_3) - \frac{\partial}{\partial u^3}(h_2 F_2) \right] \hat{\mathbf{e}}_1 + \frac{1}{h_3 h_1} \left[ \frac{\partial}{\partial u^3}(h_1 F_1) - \frac{\partial}{\partial u^1}(h_3 F_3) \right] \hat{\mathbf{e}}_2 + \frac{1}{h_1 h_2} \left[ \frac{\partial}{\partial u^1}(h_2 F_2) - \frac{\partial}{\partial u^2}(h_1 F_1) \right] \hat{\mathbf{e}}_3 $$

---

### 4. 拉普拉斯算子 (Laplacian) 的推导

**目标**：将 $\nabla^2 f = \nabla \cdot (\nabla f)$ 用 $u^i$ 坐标表示。

**推导**：
拉普拉斯算子是梯度的散度。我们直接将前面推导的梯度公式和散度公式组合起来。

1.  **梯度**： $\nabla f = \sum_{i=1}^3 \frac{1}{h_i} \frac{\partial f}{\partial u^i} \hat{\mathbf{e}}_i$。
    令 $\mathbf{G} = \nabla f$，则其物理分量为 $G_i = \frac{1}{h_i} \frac{\partial f}{\partial u^i}$。

2.  **散度**：将 $\mathbf{G}$ 代入散度公式：
    $$ \nabla \cdot \mathbf{G} = \frac{1}{h_1 h_2 h_3} \left[ \frac{\partial}{\partial u^1}(G_1 h_2 h_3) + \frac{\partial}{\partial u^2}(G_2 h_3 h_1) + \frac{\partial}{\partial u^3}(G_3 h_1 h_2) \right] $$
    $$ \nabla \cdot \mathbf{G} = \frac{1}{h_1 h_2 h_3} \left[ \frac{\partial}{\partial u^1}\left( \frac{h_2 h_3}{h_1} \frac{\partial f}{\partial u^1} \right) + \frac{\partial}{\partial u^2}\left( \frac{h_3 h_1}{h_2} \frac{\partial f}{\partial u^2} \right) + \frac{\partial}{\partial u^3}\left( \frac{h_1 h_2}{h_3} \frac{\partial f}{\partial u^3} \right) \right] $$

**最终公式（正交系）**：
$$ \boxed{\nabla^2 f = \frac{1}{h_1 h_2 h_3} \left[ \frac{\partial}{\partial u^1} \left( \frac{h_2 h_3}{h_1} \frac{\partial f}{\partial u^1} \right) + \frac{\partial}{\partial u^2} \left( \frac{h_3 h_1}{h_2} \frac{\partial f}{\partial u^2} \right) + \frac{\partial}{\partial u^3} \left( \frac{h_1 h_2}{h_3} \frac{\partial f}{\partial u^3} \right) \right]} $$
