import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 常数
PI      =   3.141592653 # 圆周率
kb      =   1.38e-23 # 玻尔兹曼常数
TT      =   298.15 # 温度25℃
hp      =   6.626e-34 # 普朗克常数
v0      =   4e15 # 绝对吸收频率
eps0    =   8.854e-12

def cal_dlvo(aa=75e-9, eps_cb=2.78, eps_nmp=32.2, n_cb=2.0, n_nmp=1.47, sigma=6e-9, nb=12, lmd=1, zeta=10e-3, verbose=False):
    """
    基于拓展DLVO理论, 考虑粘结剂的影响, 计算颗粒间距-颗粒间吸附能的函数关系,
    输入: 颗粒参数
        - 半径: aa (mm)
        - 相对介电常数: CB: eps_cb; NMP: eps_nmp
        - 折射率: CB: n_cb; NMP: n_nmp
        - ZETA电势: zeta, 静电势能
    输出: 吸附作用力参数
    """

    list_x  =   []
    list_y  =   []
    list_vdw = []
    list_born = []
    list_ele = []

    ah1 = (3/4)*kb*TT*((eps_cb - eps_nmp) / (eps_cb + eps_nmp))**2
    ah2 = 3/(16*math.sqrt(2))*hp*v0*(n_cb**2-n_nmp**2)**2 / (n_cb**2+n_nmp**2)**(3/2)
    Ah = ah1 + ah2
    # print(Ah)

    dx = 2e-12 # 差分步长

    for i in range(5000):
        hh = dx * (i + 1)
        rc = 2 + (hh / aa)

        vdw1 = 2 * aa**2 / (hh*(4*aa + hh))
        vdw2 = 2 * aa**2 / (2*aa + hh)**2
        vdw3 = math.log((hh*(4*aa + hh)) / ((2*aa + hh)**2))
        V_vdw = -(Ah/6) * (vdw1 + vdw2 + vdw3)

        V_ele = (4 * PI * eps0 * eps_nmp * (zeta**2) * (aa**2)) / (hh + 2*aa)

        born1 = 4 * Ah * (sigma / aa)**(nb - 6)
        born2 = math.factorial(nb - 6) / math.factorial(nb - 2)
        born3 = lmd / ((lmd + 1) * (rc - 1 - lmd)**(nb - 5))
        V_born = born1 * born2 * born3
        V_total = V_vdw + V_ele + V_born

        list_x.append(hh)
        
        if verbose:
            if i % 100 == 0 or i < 10:
               print(f"hh(e-9): {hh/1e-9:.3f}\t vdw(kT): {(V_vdw/(kb*TT)):.3f}\t ele(kT): {(V_ele/(kb*TT)):.3f}\t born(kT): {(V_born/(kb*TT)):.3f}\t Total(kT): {((V_total)/(kb*TT)):.3f}")

        list_vdw.append(V_vdw)
        list_born.append(V_born)
        list_ele.append(V_ele)
        list_y.append(V_total) #% v_total
    
    return list_x, list_vdw, list_born, list_ele, list_y

def cal_param(list_x, list_y, dx=2e-12, verbose=False):
    delta1 = np.diff(list_y, n=1)
    list_y1 = list(map(lambda x: x/dx, delta1))

    delta2 = np.diff(list_y1, n=1)
    list_y2 = list(map(lambda x: x/dx, delta2))

    # 求 uc, fc, k0
    # uc 和 fc 不在同一点
    # 势能函数最小时, 为uc
    # 一阶导数最大时, 为fc
    Uc, Uci = -np.min(list_y), np.argmin(list_y)
    hc = list_x[Uci]
    Fc = np.max(list_y1)
    k0 = list_y2[Uci]*hc/Fc # k0
    if verbose:
        print(f"Uc(kbT): {Uc/(kb*TT):.2f}\t hc(e-9): {hc/1e-9:.2f}\t Fc(e-12): {Fc/1e-12:.2f}\t k0: {k0:.2f}")
    return Uc, hc, Fc, k0

def cal_rheology(hc=2.29e-9, # 颗粒临界距离
                 uc=32.5*kb*TT, # 临界吸附能
                 fc=41.5e-12,  # 临界吸附力
                 k0=12.9,  # 无量纲梯度因子
                 aa=75e-9,  # 胶体半径
                 df=1.8,    # 分形维度
                 dl=1.55,   # 化学维度
                 miu=3.6e-3, # 溶剂粘度
                 phi_p=3.20345 * 0.02,  # 原始颗粒浓度
                 phi_m=0.6,  # 颗粒最大堆积浓度
                 zz=5, # 平衡条件下最小键断裂数
                 verbose=False):
    
    """
    基于吸附力参数, 计算当前颗粒、溶剂、粘结剂体系下的理想流变曲线
    输入: 吸附力参数
    输出: 理想流变曲线

    """
    # miu = 3.6e-3 # 根据PVDF提前计算好的溶液粘度 eta_0
    phi_ref = 1.2 # ??? 没弄懂哪里来的
    zz = 5

    tau_0 = (6*PI*miu*(aa**3))/(kb*TT)
    tau_sic =   (tau_0**-1)*(k0**0.5)*((uc*aa)/(kb*TT*hc))*math.exp((-zz*uc)/(kb*TT)) # tau_rs不需要迭代的部分
    tau_si  =   tau_sic ** -1   # 可能与粒子间相互作用或能量势垒相关

    eta_list = []
    gama_list = []

    phi_a_list = []

    for i in range(13):
        gm = i * 0.5 - 2
        gama = 10**gm  # 遍历gama值
        qq_min = aa
        qq_max = aa * 30
        
        counter = 0
        while True:
            counter += 1
            qq = (qq_max + qq_min) / 2
            phi_int = (qq / aa) ** (df - 3)
            phi_a = phi_p / phi_int

            if phi_a > phi_m * 0.95 / (phi_ref**3): # 不知道哪里来的限制
                phi_a = phi_m * 0.95 / (phi_ref**3)
            
            tau_sc = (tau_0 ** -1) * (k0 ** 0.5) * (uc/(kb*TT))**1.5 \
                    * ((aa / hc) ** 2) * (aa / qq) * math.exp(-uc/(kb*TT))
            tau_s = tau_sc ** -1
            tau_rc = (k0 ** 0.5) * gama * (aa / hc) * ((qq / aa) ** -dl)
            tau_r = tau_rc ** -1
            tau_sr = tau_0 * qq / aa
            tau_rs = tau_si * ((aa / qq) ** dl)

            alpha = tau_s * (tau_r + tau_rs) / (tau_s * tau_rs + (tau_r + tau_rs) * tau_sr)
            
            eta1 = miu * ((1 - (phi_a / phi_m) * (phi_ref**3)) ** (-2))
            # eta1 = miu * (1 - phi_a*(qq/aa)**(3-df)/phi_m)**(-2)
            # eta1 = miu * (1 - phi_a/phi_m)**(-2.5*phi_m) # 单位直接是 Pa s
            eta2 = ((alpha * fc * hc * (k0**(-1/2)) * (phi_a**2)) / ((aa**3) * gama)) * ((phi_p / phi_a) ** ((5 - 2 * df - dl) / (3 - df))) # 单位直接是 Pa s
            # eta2 = ((alpha * fc * hc * (phi_a**2)) / ((aa**3) * gama)) * ((phi_p / phi_a) ** ((5 - 2 * df - dl) / (3 - df)))
            # eta2 *= 10
            eta = eta1 + eta2

            gama_c = (2 * fc * aa) / (5 * PI * eta * (qq**3))
            gama_if = gama_c / gama
            if gama_if > 1.0001 and counter < 10000:
            # if gama_if > 1.0001:
                qq_min = qq
            elif gama_if < 0.9999 and counter < 10000:
            # elif gama_if < 0.9999:
                qq_max = qq
            else:
                if counter >= 10000:
                    print("Counter Overflow!")
                eta_list.append(eta)
                gama_list.append(gama)
                phi_a_list.append(phi_a)

                qa = qq / aa
                if verbose:
                    print(f'gama: {gama:.3f}\t eta: {eta:.3f}\t qa: {qa:.3f}\t eta1 {eta1:.3f}\t eta2 {eta2:.3f}')
                break
            # print('alpha: ', alpha)
        # print(f'phi_a: {phi_a:.3f}\t phi_max: {phi_m * 0.95 / (phi_ref**3):.3f}\t phi_ref: {phi_ref:.3f}')
        # break
        
    return gama_list, eta_list, phi_a_list

def show_rheology(gama_list, eta_list):
    
    fig, ax = plt.subplots(figsize=(8, 6))  # 可以调整figsize来控制图像大小
    ax.plot(gama_list, eta_list, color='blue')
    # 添加标签和标题
    ax.set_xlabel('gama (/s)', fontsize=12)
    ax.set_ylabel('eta (Pa s)', fontsize=12)
    # ax.set_title('', fontsize=14)

    # 添加网格线
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.set_yscale('log')
    ax.set_xscale('log')

    # 调整坐标轴范围 (可选，根据你的数据范围调整)
    ax.set_xlim(0.01, 100)
    ax.set_ylim(0.01, 1000)

    # 显示图形
    plt.tight_layout()  # 自动调整子图参数, 避免标签重叠
    plt.show()

def cal_miu(w_pvdf, eta_sol=1.65e-3, mode='5130'):
    """
        根据eta_sp - w_pvdf关系, 求加入PVDF后纯液体的粘度miu.cal_miu
        Args:
            w_pvdf (float): PVDF的质量分数
            eta_sol (float): 溶剂的粘度, 单位为Pa.s, 默认为NMP, 粘度为1.65e-3 Pa.s
            mode (str): 模式, 可选值为'5130', '6020', '1015'
    """
    if mode in ['5130', 'SOLEF5130', 'Solef5130']:
        k_eta = 2.11e6
    elif mode in ['6020', 'SOLEF6020', 'Solef6020']:
        k_eta = 1.18e6
    elif mode in ['1015', 'SOLEF1015', 'Solef1015']:
        k_eta = 4.95e6
    else:
        raise ValueError(f"Unsupported mode: {mode}. Supported modes are '5130', '6020', '1015'.")  
    eta_sp = k_eta * w_pvdf**3
    miu = eta_sp * eta_sol + eta_sol
    return miu

if __name__ == "__main__":
    # 测试函数
    # cal_rheology(miu=5e-2, verbose=True)
    # Fc = 17.4e-12 # N
    # Uc = 18.42 * kb * TT
    # hc = 3.46e-9
    # k0 = 12.81

    # # gama_list, eta_list = cal_rheology(verbose=False)
    # gama_list, eta_list = cal_rheology(hc=hc, uc=Uc, fc=Fc, k0=k0, dl=1.3, phi_p=0.4, phi_m=0.61, verbose=True, miu=0.04)
    # cal_dlvo(verbose=True)
    cal_rheology(verbose=True)
