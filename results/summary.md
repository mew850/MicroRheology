# 总结

## 20260604
### 单位问题
两个文献中，粘度计算公式计算直接得到的eta1和eta2，单位都是 Pa s


关于eta2的换算问题：
```python
eta1 = miu * (1 - phi_a/phi_m)**(-2.5*phi_m)
eta2 = ((alpha * fc * hc * (k0**(-1/2)) * (phi_a**2)) / ((aa**3) * gama)) * ((phi_p / phi_a) ** ((5 - 2 * df - dl) / (3 - df))) # 单位直接是 Pa s
```
hc的单位是 m ，所以eta2的单位是$\frac{N \times m}{s^{-1}\times m^2} = Pa·s$，单位没有问题

### Ma 和 Chen计算粘度公式的差别

最关键的差别：
1. 多了个$\boldsymbol{k_0}$: eta2计算上，Ma采用的是$\eta_s = \frac{\alpha F_c h_c}{\gamma a^3} \times \phi_a^2 (\frac{\phi_p}{\phi_a})^{\frac{5-2d_f-d_l}{2-d_f}}$，而Chen采用的是$\eta_s = \frac{\alpha F_c h_c}{\gamma a^3 k_0^{1/2}} \times \phi_a^2 (\frac{\phi_p}{\phi_a})^{\frac{5-2d_f-d_l}{2-d_f}}$，
2. $\phi_{int}$中人工参数k的使用:$\phi_{int}=(\frac{Rg}{a})^{df-3}$, Ma使用人工参数k，$q=k R_g$计算得到$\phi_{int}=(\frac{q}{ka})^{df-3}$，而Chen直接使用$q = R_g$计算得到$\phi_{int}=(\frac{q}{a})^{df-3}$，没有使用参数k；但是，根据代码，Chen等人在$phi_p$前面莫名其妙多乘了个系数，系数的含义不明晰，严重怀疑就是在这里做了修正！
3. 体系引起的DLVO计算方式不同：Chen使用born电势对正极浆料的DLVO建模，而Ma是采用简单的高分子链导出的空间排斥力

决定复杂问题简单化，不然做不完了，做如下简化：
* CFD建模的时候，PVDF影响的是吸附力参数，不单独建模，假设初始时就分散均匀了
* 重点考察CB和AM浓度的影响：AM只影响偏置；由于CB浓度严重影响微观流变模型收敛性，且实际上微观流变模型也确实是对均匀浆料进行的建模，因此考虑当CB体积分数5%时，将CB和溶剂分开考虑，都保留自身的“粘度”(此时CB是连续固体)；只有体积分数小于5%时，才启动微观流变模型的本构方程，这样就能近似出从浆料从离散的牛顿流体逐步演变为非牛顿流体的过程


