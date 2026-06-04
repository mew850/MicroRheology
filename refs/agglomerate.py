"""计算的是碳黑的粘度贡献, 但还没弄懂如何计算不同PVDF浓度下的粘度值"""

import math

pi      =   3.141592653
kk      =   1.38e-23    # 玻尔兹曼常数

hc      =   2.29    #* 特征长度
uc_kt   =   32.5    #*
fc      =   41.5    #* 特征力
k0      =   12.9    #*

hc      *=  1e-9
fc      *=  1e-12

temp    =   298.15
uc      =   uc_kt*kk*temp
zz      =   5 # minimum number of bond to be broken
aa      =   75e-9   # primary particle 半径, 对应原文中的 rp
dl      =   1.55    # 分形维数
df      =   1.8
miu     =   3.6
miu     *=  1e-3

# phi_p   =   3.20345 * 0.02 #* primary particle 浓度
phi_p = 3.20345 * 0.0133
phi_m   =   0.4

tau_0   =   (6*pi*miu*(aa**3))/(kk*temp)    # 布朗运动的时间常数

# 把qq 和 aa单独拎出来, 迭代找到稳定离子间距qq
tau_sic =   (tau_0**-1)*(k0**0.5)*((uc*aa)/(kk*temp*hc))*math.exp((-zz*uc)/(kk*temp)) # tau_rs不需要迭代的部分
tau_si  =   tau_sic ** -1   # 可能与粒子间相互作用或能量势垒相关
phi_ref =   1.2 # 自己设置的检测变量?

eta_list = []
gama_list = []

print('     q/a \tgama \teta1 \teta2 \teta \tetar \tphiA \tmu')
for i in  range(13):
    gm      =   i*0.5-2
    gama    =   10**gm  # 遍历gama值
    qq_min = aa
    qq_max = aa * 30

    while True: # 二分法调整`qq`, 快速找到`gama_c` ~= `gama`的`qq`, 即方程的解
        qq = (qq_max + qq_min) / 2 # aggregate particle, 对应原文中的a
        phi_int = (qq / aa) ** (df - 3) # 公式(12)
        phi_a = phi_p / phi_int
        if  phi_a > phi_m*0.95/(phi_ref**3):
            phi_a = phi_m*0.95/(phi_ref**3)

        tau_sc = (tau_0 ** -1) * (k0 ** 0.5) * (uc_kt ** 1.5) \
                 * ((aa / hc) ** 2) * (aa / qq) * math.exp(-uc_kt)
        tau_s = tau_sc ** -1
        tau_rc = (k0 ** 0.5) * gama * (aa / hc) * ((qq / aa) ** -dl)
        tau_r = tau_rc ** -1
        tau_sr = tau_0 * qq / aa
        tau_rs = tau_si * ((aa / qq) ** dl)

        alpha = (tau_s * (tau_r + tau_rs)) / ((tau_s * tau_rs) + (tau_sr * (tau_r + tau_rs)))

        eta1 = miu * ((1 - (phi_a / phi_m)*(phi_ref**3)) ** (-2)) # hydro
        # 经过推导, phi_ref的作用只是保证 phi_a <= 0.95*phi_m
        eta2 = ((alpha * fc * hc * (k0 ** -0.5) * (phi_a ** 2))
                / ((aa ** 3) * gama)) * ((phi_p/phi_a) ** ((5 - 2 * df - dl)/(3 - df))) # struct 和文献对不上啊?
        eta = eta1 + eta2
        etar = eta / miu
        # 选用 gamma 作为检验变量, 对qq进行二分法查找
        gama_c = (2 * fc * aa) / (5 * pi * eta * (qq ** 3)) #   临界剪切速率, 
        gama_if = gama_c / gama
        if gama_if > 1.0001:
            qq_min = qq

        elif gama_if < 0.9999:
            qq_max = qq

        else:
            eta_list.append(eta)
            gama_list.append(gama)
            qa = qq / aa
            line    =   'q/a: '+str('%.3f' % qa)+'\t '+str('%.3f' % gama)+'\t '+str('%.3f' % eta1)\
                        +'\t '+str('%.3f' % eta2)+'\t '+str('%.3f' % eta)+'\t '+str('%.3f' % etar)\
                        +'\t '+str('%.3f' % phi_a) +'\t '+str('%.3f' % miu)
            print(line)
            break


# save
import pandas as pd
tmp = pd.DataFrame({
    'gama': gama_list,
    'eta': eta_list,
})
tmp.to_csv('./results/eta_gama.csv', index=False)