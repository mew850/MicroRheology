"""
File: utils
Abstract: 微观流变模型的复现
Author: Kuraying Pat

"""
import math
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np

PI = math.pi
kb = 1.380649e-23  # J/K
hp = 6.62607015e-34  # J·s
TT = 298.15  # K

def cal_dlvo(L, Mw, aa=275e-9, eps_cb=2.78, eps_nmp=32.2, n_cb=2.0, n_nmp=1.47, v0=3e15, 
            verbose=False, is_check=False, is_plot=False):
    """
    主要参考Ma的DLVO模型
        L: 
        Mw: PVDF分子量 g/mol
        所有物理量均取国际单位
    """
    ah1 = (3/4)*kb*TT*((eps_cb - eps_nmp) / (eps_cb + eps_nmp)) ** 2
    ah2 = 3/(16*math.sqrt(2))*hp*v0*(n_cb**2-n_nmp**2)**2 / (n_cb**2+n_nmp**2)**(3/2)
    Ah = ah1 + ah2
    # print(f"Ah:{Ah:.4g}, ah1:{ah1:.4g}, ah2:{ah2:.4g}")
    list_x = []
    list_vdw = []
    list_osm = []
    list_el = []
    list_y = []

    dx = 2e-9 # 2 nm
    for i in range(1000):
        hh = dx * (i+1)

        # VDW
        vdw_1 = 2 * aa**2 / (hh*(4*aa + hh))
        vdw_2 = 2 * aa**2 / (2*aa + hh)**2
        vdw_3 = math.log((hh*(4*aa + hh)) / ((2*aa + hh)**2))
        V_vdw = -(Ah/6) * (vdw_1 + vdw_2 + vdw_3)

        # Steric
        ## Osm
        CHI = 0.2 # TODO: 参考值, 看影响大不大， 后续可能需要替换
        # v_solvent = cal_v_solvent()*(1e-9)**3 # TODO: NMP体系, 根据公式可计算
        v_solvent = 0.16e-27 # nm^3 -> m^3
        # print(v_solvent)
        phi_poly = 0.05 # TODO: 文献表明可忽略其影响
        osm_1 = (4*PI*aa*kb*TT/v_solvent) * phi_poly**2 * (0.5-CHI)
        
        if hh > 2*L:
            V_osm = 0
        elif hh >= L and hh < 2*L: 
            V_osm = osm_1 * (L-hh/2)**2
        else:
            V_osm = osm_1 * L**2 * ((hh/(2*L)) - (1/4) - math.log(hh/L))
        ## el
        if hh > L:
            V_el = 0
        else:
            rho_poly = 1.78 * 1000 # g/cm^3 TODO: 和PVDF自身有关
            V_el = 2*PI*aa*kb*TT/(Mw/1000) * phi_poly * L**2 * rho_poly \
                    *(hh/L * math.log(hh/L * ((3-hh/L)/2)**2) - 6*math.log(((3-hh/L)/2)) + 3*(1-hh/L))
            
        V_steric = V_osm + V_el

        if verbose:
            if i % 50 == 0 or i < 10:
                print(f"hh(e-9): {hh/1e-9:.3f}\t vdw(kT): {(V_vdw/(kb*TT)):.3f}\t osm(kT): {(V_osm/(kb*TT)):.3f}\t el(kT): {(V_el/(kb*TT)):.3f}\t Total(kT): {((V_vdw+V_steric)/(kb*TT)):.3f}")

        # output
        list_x.append(hh)
        list_vdw.append(V_vdw)
        list_osm.append(V_osm)
        list_el.append(V_el)
        list_y.append(V_vdw+V_steric)

    if is_plot:
            # plot
        hh_list = [hh /aa  for hh in list_x]  # 将距离转换为纳米单位
        # hh_list = [hh / 1e-8 for hh in list_x]
        vdw_list = [i / (kb*TT) for i in list_vdw]  # 将势能转换为kBT单位
        osm_list = [i / (kb*TT) for i in list_osm]  # 将势能转换为kBT单位
        el_list = [i / (kb*TT) for i in list_el]  # 将势能转换为kBT单位
        steric_list = [osm_list[i]+el_list[i] for i in range(len(osm_list))]
        v_list = [i / (kb*TT) for i in list_y]  # 将总势能转换为kBT单位

        fig, ax = plt.subplots(figsize=(8, 6))  # 可以调整figsize来控制图像大小
        ax.plot(hh_list, vdw_list, label='Van der Waals', color='blue', linestyle='dotted')
        ax.plot(hh_list, osm_list, label='osm', color='green')
        ax.plot(hh_list, el_list, label='el ', color='red')
        ax.plot(hh_list, steric_list, label='Poly Steric', color='purple', linestyle='-.')  # 添加Steric Repulsion曲线
        ax.plot(hh_list, v_list, label='Total', color='orange')  # 添加总势能曲线

        # 添加标签和标题
        ax.set_xlabel('h/a', fontsize=12)
        ax.set_ylabel('Potential Energy (kbT)', fontsize=12)
        ax.set_title('DLVO Potential vs. Distance', fontsize=14)

        # 添加图例
        ax.legend(fontsize=10)

        # 添加网格线
        ax.grid(True, linestyle='--', alpha=0.5)

        # 调整坐标轴范围 (可选，根据你的数据范围调整)
        ax.set_xlim(0.05, 0.3)
        ax.set_ylim(-15, 10)

        # 显示图形
        plt.tight_layout()  # 自动调整子图参数, 避免标签重叠
        plt.show()
    if is_check:
        return list_x, list_y, list_vdw, list_osm, list_el
    else:
        return list_x, list_y
    
def cal_param(list_x, list_y, aa=275e-9, dx=2e-9, verbose=False, is_plot=False):
    delta1 = np.diff(list_y, n=1)
    list_y1 = list(map(lambda x: x/dx, delta1)) # -dU / dt

    delta2 = np.diff(list_y1, n=1)
    list_y2 = list(map(lambda x: x/dx, delta2))

    # 求 uc, fc, k0
    # uc 和 fc 不在同一点
    # 势能函数最小时, 为uc
    # 一阶导数最大时, 为fc
    Uc, Uci = -np.min(list_y), np.argmin(list_y)
    hc = list_x[Uci]
    Fc = np.max(list_y1)
    k0 = -list_y2[Uci]*hc/Fc # k0
    # TODO: k0这里由于差分顺序, 可能会有问题Uci-1, Uci, 和Uc对应的值分别为-141.49, 1.609, 1.946
    if verbose:
        print(f"Uc(kbT): {Uc/(kb*TT):.2f}\t hc(e-9): {hc/1e-9:.2f}\t Fc(e-12): {Fc/1e-12:.2f}\t k0: {k0:.2f}")
    if is_plot:
        fig, ax = plt.subplots(figsize=(8, 6))  # 可以调整figsize来控制图像大小
        
        # 显示图形
        plt.tight_layout()  # 自动调整子图参数, 避免标签重叠
        hh_list = [hh /aa  for hh in list_x]  # 将距离转换为纳米单位
        force_list = [i / 1e-11 for i in list_y1]
        ax.plot(hh_list[:len(force_list)], force_list, label='Net Force', color='red')
        ax.set_xlim(0.05, 0.3)
        ax.set_ylim(-1, 1)

        plt.show()

    return Uc, hc, Fc, k0

def cal_qq(qq, gamma, Fc, hc, k0, aa, phi_p, miu, alpha, k, df, dl, phi_m=0.6):
    """
    Input:
        qq: [m] 颗粒间距
        gamma: [s^-1] 应变速率
        Fc: [N] 临界吸附力
        hc: [m] 颗粒临界距离
        k0: [-] 无量纲梯度因子
        
        aa: [m] secondary aggregate radius
        phi_p: [-] primary particle浓度
        miu: [Pa.s] 溶剂粘度

        alpha: [-] capture efficiency
        k: [-] q ~ k*Rg fitting paramter
        df: [-] fractal dimension
        dl: [-] chemical dimension
        phi_m: [-] maximum packing concentration, default 0.61

    Return:
        eta_s: [Pa.s] eta_s(gamma)
    """
    phi_a = phi_p * (qq/(k*aa))**(3-df)
    eta_s = ((alpha * Fc * hc * (k0**(-1/2)) * (phi_a**2)) / ((aa**3) * gamma)) * ((phi_p / phi_a) ** ((5 - 2 * df - dl) / (3 - df)))
    # eta_s = (alpha*Fc*hc/((aa**3)*gamma)) * phi_a**2 * (phi_p/phi_a)**((5-2*df-dl)/(3-df))
    eta_h = miu * (1 - phi_a / phi_m)**(-2.5*phi_m)
    eta_CB = 2*Fc*aa / (5*PI*gamma*(qq**3))
    return eta_CB - (eta_s + eta_h)

def cal_viscosity(alpha=0.9, gamma=1.0, Fc=5.72e-13, hc=35.639e-9, k0=1.2, aa=275e-9, phi_a=0.001,
                phi_p=0.009, phi_m=0.61, df=1.6, dl=1.0, miu=0.08, verbose=False
):
    eta_struct = ((alpha * Fc * hc * (k0**(-1/2)) * (phi_a**2)) / ((aa**3) * gamma)) * ((phi_p / phi_a) ** ((5 - 2 * df - dl) / (3 - df)))
    eta_hydro = miu * (1 - (phi_a/phi_m))**(-2.5*phi_m)
    eta = eta_struct + eta_hydro
    if verbose:
        print(f"gamma:{gamma:.3f}\t eta_struct:{eta_struct:.3f}\t eta_hydro:{eta_hydro:.3f}\t eta:{eta:.3f}\t eta_r:{eta/miu:.3f}")
    return eta
