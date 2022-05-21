## Homework #1

<div align = "center">陈珂旸 1700012118</div>

### 代码

#### 运行方式

本次作业用Julia 1.7.2完成，运行下方代码得到$H_2$ 在1.4 a.u.下的总能量和结合能。

```julia
>julia H2.jl -d 1.4
# Orbital 1: |ψ₁⟩ = 0.5489340404350302|φ₁⟩ + 0.5489340404350302|φ₂⟩, energy = -0.5782029768532935.
# Orbital 2: |ψ₂⟩ = -1.2114640729141275|φ₁⟩ + 1.2114640729141275|φ₂⟩, energy = 0.6702677605933027.
# Total Energy at 1.4 a.u. is -1.1167143251757694.
# Bond Energy at 1.4 a.u. is -0.1835506243534839.
```

可选择的参数有：

- -d: 两个原子核距离，默认为$1.4$；
- --iter, -i: 最大迭代次数，默认为$100$；
- --eps, -e: 收敛精度，默认$10^{-9}$；
- --criteria, -c: 收敛条件，可选则能量收敛“energy”或密度矩阵收敛“density”;
- --prt, -p: 是否打印迭代过程，无参数，默认为“否”。
- --help, -h: 查看参数说明。

#### 代码说明

全部代码在gaussian.jl, RHFSCF.jl, H2.jl这三个文件中。

**gaussian.jl**:：定义了计算STO-3G基积分的函数，用于计算$\langle\alpha|\beta\rangle$，$\langle\alpha|-\frac{1}{2}\nabla^2|\beta\rangle$，$\langle\alpha|\frac{1}{r_{12}}|\beta\rangle$，$\langle\alpha|-\sum_I\frac{Z_I}{|r-R_I|}|\beta\rangle$，$\langle\alpha\beta|\gamma\delta\rangle$这几个形式的积分。

**RHFSCF.jl**：执行自洽场方程迭代的主要程序。定义结构**Molecule**存储原子核、基、分子轨道、轨道能量、总能量等信息。函数**RHFSCF!(::Molecule)**运行自洽场迭代，更新分子轨道的系数和能量。

**H2.jl**：处理$H_2$的程序。给定两原子核的距离$d$，程序调用前两个文件执行RHF算法啊，最终打印输出轨道能量，总能量和结合能。 

### 算法流程

1. 初始化系统，将初始的系数矩阵设为零矩阵
   $$
   C_0 = \left[\begin{array}{lll}
   0 & 0  \\
   0 & 0 
   \end{array}\right]
   $$

   2. 计算基的重叠矩阵$S$，单体哈密顿量$h$和电子-电子相互作用张量$g$，并且对$S$特征值分解$U^{\dagger}SU=s$，令$X=Us^{-1/2}$。
      $$
      S=\langle\phi_i|\phi_j\rangle = \left[\begin{array}{lll}
      1.00000 & 0.65932  \\
      0.65932 & 1.00000 
      \end{array}\right]\\
      h = \langle\phi_i|-\frac{1}{2}\nabla^2+\sum_I\frac{Z_I}{|r-R_I|}|\phi_j\rangle = \left[\begin{array}{lll}
      -1.12041 & -0.95838  \\
      -0.95838 & -1.14041 
      \end{array}\right] \\
      g_{\mu\nu\lambda\sigma} = (\phi_\mu\phi_\nu|\phi_\sigma\phi_\lambda)-\frac{1}{2}(\phi_\mu\phi_\lambda|\phi_\sigma\phi_\nu)
      $$
      

3. 计算密度矩阵$D$和Fock矩阵$F$：
   $$
   D_{\lambda\sigma} = 2C^*_{\lambda 1}C_{\sigma 1}\\
   F_{\mu\nu} = h_{\mu\nu} + \sum_{\lambda\sigma}D_{\lambda\sigma}g_{\mu\nu\lambda\sigma}
   $$
   
4. 求解广义本征值问题$FC=\varepsilon SC$，求出$X^{\dagger}FX$的本征值$\varepsilon$和本征向量$C'$，进而求出广义本征向量$C=XC'$。

   5. 如果（能量或密度矩阵）误差$\delta$是否小于给定值，程序结束；否则回到第3步。

      ```julia
      C = zeros(2,2)
      Calculate S,h,g.
      s,U = eigen(S)
      X = U / s^(1/2)
      for i=1:iter_max 
          D = 2*C[:,1]*C[:,1]'
          @einsum F[μ,ν] := h[μ,ν] + g[μ,ν,σ,λ]*D[σ,λ]
          ε, C' = eigen(X'*F*X)
          C = X*C'
          Calculate error δ.
          if δ < eps
              break
          end
      end
      ```

      ### 结果

      #### i. Energy at equilibrium bond length (1.4 a.u.)

        - Output of each SCF step

          ![image-20220521230855279](C:\Users\Locky\AppData\Roaming\Typora\typora-user-images\image-20220521230855279.png)

        - Converged total energy: -1.116714325 a.u.  Converged criterion: $|E^{i-1}_{tot}-E^{i}_{tot}|\leq 10^{-9}$.

        - Orbital energies: 
      
          ​		Orbital 1: |ψ₁⟩ = 0.5489340404350302|φ₁⟩ + 0.5489340404350302|φ₂⟩, energy = -0.578202977 a.u.
      
          ​		Orbital 2: |ψ₂⟩ = -1.2114640729141275|φ₁⟩ + 1.2114640729141275|φ₂⟩, energy = 0.670267761 a.u.
      
      #### ii. Bonding curve of $H_2$
      
      ![](D:\多体系统的量子理论\many-body code\RHF\STO-3G_H2.png)
      
      

| d (a.u.) | Total Energy (a.u.) | Binding Energy (a.u.) |
| :------: | :-----------------: | :-------------------: |
|   0.7    |      -0.837130      |       0.096033        |
|   1.0    |      -1.065999      |       -0.132836       |
|   1.3    |      -1.116871      |       -0.183707       |
|   1.4    |      -1.116714      |       -0.183551       |
|   1.5    |      -1.111696      |       -0.178532       |
|   2.0    |      -1.049171      |       -0.116007       |
|   3.0    |      -0.885275      |       0.047889        |
|   4.0    |      -0.761082      |       0.172081        |

### 讨论

1. $H_2$轨道系数矩阵的严格解为：

$$
C^{*} = \left(\begin{array}{cc}
{\left[2\left(1+S_{12}\right)\right]^{-1 / 2}} & {\left[2\left(1-S_{12}\right)\right]^{-1 / 2}} \\
{\left[2\left(1+S_{12}\right)\right]^{-1 / 2}} & -\left[2\left(1-S_{12}\right)\right]^{-1 / 2}
\end{array}\right) = \left(\begin{array}{cc}
0.548934 & 1.211464 \\
0.548934 & -1.211464
\end{array}\right)
$$

和自洽场迭代得到的结果一致。

2. STO-3G本身存在误差，用STO-3G基计算氢原子能量结果为$-0.466582$，与理论值$-0.5$相差$6.7\%$。如果用精度更高的STO-6G的基计算氢分子得到的总能量和结合能为 $-1.125324$和$-0.192161$，和STO-3G相比相差$1\%$。
3. 求解广义本征值方程$FC=\varepsilon SC$，用Julia LinearAlgebra中的eigen函数或Python中scipy.linalg.eigh函数求解广义本征值问题都会遇到收敛到错误结果的情况，这两个函数都建立在LAPACK之上，原因不明。用书中的方法或者解$S^{-1}FC=\varepsilon C$的本征值问题都能得到正确的结果。 

