# -*- coding: utf-8 -*-
"""
シナプス後AChR, AChおよび筋弛緩薬の競合シミュレーション
(Ropen(ARA*)の追加後のモデル)
"""
import sympy as sym
import numpy as np
import matplotlib.pyplot as plt

import sys
ui_for_scipy_path = r"C:\Users\hoshino\Documents\Python Scripts\ui_for_scipy"
if not ui_for_scipy_path in sys.path:
    sys.path.append(ui_for_scipy_path)
from ui_for_scipy import InitialValueProb

#変数
t = sym.Symbol('t')
A = sym.Symbol('A')
ARO, ORA = sym.symbols('ARO, ORA')
ARA = sym.symbols('ARA')
DRO, ORD = sym.symbols('DRO, ORD')
DRD = sym.symbols('DRD')
DRA, ARD    = sym.symbols('DRA, ARD')
Ropen,Rdec  = sym.symbols('Ropen, Rdec')
D = sym.Symbol('D')

# パラメータを無次元化する際の基準量
R_base = 7.75*(10**-5) #[mol/L]
Tbase = 10**-3 #[ms]

class PostSynapticModel:
    """
    シナプス後のモデル
    :param kinet: 反応の速度定数(kinetics)を表すパラメータ
    :param pnt:  患者(patient)の特性を表すパラメータ
    """    
    def __init__(self, kinet=None, pnt=None):

        #反応の速度定数                
        if kinet:
            self.k_dissA1 = kinet['k_dissA1'] *Tbase 
            self.k_dissA2 = kinet['k_dissA2'] *Tbase 
            self.k_assocA1= kinet['k_assocA1'] *R_base*Tbase
            self.k_assocA2= kinet['k_assocA2'] *R_base*Tbase
            self.k_dissD1 = kinet['k_dissD1'] *Tbase
            self.k_dissD2 = kinet['k_dissD2'] *Tbase
            self.k_assocD1= kinet['k_assocD1'] *R_base*Tbase
            self.k_assocD2= kinet['k_assocD2'] *R_base*Tbase
            self.k_hydrol = kinet['k_hydrol'] *Tbase
            self.k_close = kinet['k_close']*Tbase
            self.k_open = kinet['k_open']*Tbase
            self.k_dpls = kinet['k_dpls']*Tbase
            self.k_dmns = kinet['k_dmns']*Tbase

        # 患者の個人差に関わるパラメータ
        if pnt:
            self.R_total = pnt['R_total'] / R_base
            self.gamma = pnt['gamma']
            self.ARA50 = pnt['ARA50'] / R_base
        
        self._formulate_ODE()

    def set_patient_dependent_params(self,pnt):
        
        self.R_total = pnt['R_total'] / R_base
        self.gamma = pnt['gamma']
        self.ARA50 = pnt['ARA50'] / R_base
        self._formulate_ODE()        
        

    def _formulate_ODE(self):
        
        # パラメータの設定
        R_total  = self.R_total
        k_dissA1 = self.k_dissA1
        k_dissA2 = self.k_dissA2
        k_assocA1= self.k_assocA1
        k_assocA2= self.k_assocA2
        k_dissD1 = self.k_dissD1
        k_dissD2 = self.k_dissD2
        k_assocD1= self.k_assocD1
        k_assocD2= self.k_assocD2
        k_hydrol = self.k_hydrol
        alpha    = self.k_close
        beta     = self.k_open
        k_dpls   = self.k_dpls
        k_dmns   = self.k_dmns

        ORO = R_total - (ARO+ORA+ARA+DRO+ORD+DRD+ARD+DRA+Ropen+Rdec)
        dxdt = {A: k_dissA1*(ARA+ARO+ARD) + k_dissA2*(ARA+ORA+DRA) - k_hydrol*A \
                    - (k_assocA1+k_assocA2)*A*ORO - k_assocA1*A*(ORA+ORD) - k_assocA2*A*(ARO+DRO),
                ARO: k_assocA1*A*ORO + k_dissA2*ARA + k_dissD2*ARD \
                    - k_dissA1*ARO - k_assocA2*A*ARO - k_assocD2*D*ARO, 
                ORA: k_assocA2*A*ORO + k_dissA1*ARA + k_dissD1*DRA \
                    - k_dissA2*ORA - k_assocA1*A*ORA - k_assocD1*D*ORA, 
                ARA: k_assocA1*A*ORA + k_assocA2*A*ARO - (k_dissA1+k_dissA2)*ARA \
                    - beta*ARA + alpha*Ropen, 
                DRO: k_assocD1*D*ORO + k_dissA2*DRA + k_dissD2*DRD \
                    - k_dissD1*DRO - k_assocA2*A*DRO - k_assocD2*D*DRO, 
                ORD: k_assocD2*D*ORO + k_dissA1*ARD + k_dissD1*DRD \
                    - k_dissD2*ORD - k_assocA1*A*ORD - k_assocD1*D*ORD, 
                DRD: k_assocD1*D*ORD + k_assocD2*D*DRO - (k_dissD1+k_dissD2)*DRD, 
                ARD: k_assocA1*A*ORD + k_assocD2*D*ARO - (k_dissA1+k_dissD2)*ARD, 
                DRA: k_assocD1*D*ORA + k_assocA2*A*DRO - (k_dissD1+k_dissA2)*DRA,
                D  : 0,
                Ropen: -alpha*Ropen + beta*ARA + k_dmns*Rdec  - k_dpls*Ropen,
                Rdec : -k_dmns*Rdec  + k_dpls*Ropen
                }                
        
        self.ivp = InitialValueProb((A, ARO, ORA, ARA, DRO, ORD, DRD, DRA, ARD, D, Ropen, Rdec), t)
        self.ivp.set_derivative(dxdt)

    def set_steady_state(self, consD):
        """
        与えられた筋弛緩薬濃度の下での定常状態を計算します. 引数consDに単位 M での数値を与えてください. 
        """
        # 与えられた筋弛緩薬濃度を単位化
        consD = consD / R_base
        
        # 定常状態における濃度を計算
        K_A1 = self.k_dissA1/self.k_assocA1
        K_A2 = self.k_dissA2/self.k_assocA2
        K_D1 = self.k_dissD1/self.k_assocD1
        K_D2 = self.k_dissD2/self.k_assocD2
        common_denom = (consD*K_A1 + K_D1*K_A1)*(consD*K_A2 + K_D2*K_A2)

        # 初期値の計算                    
        self.init = {A:  0, 
                    ARO: 0,
                    ORA: 0,
                    ARA: 0, 
                    DRO: self.R_total*consD*K_A1*K_A2*K_D2/common_denom,
                    ORD: self.R_total*consD*K_A1*K_A2*K_D1/common_denom,
                    DRD: self.R_total*(consD**2)*K_A1*K_A2/common_denom,
                    ARD: 0,
                    DRA: 0,
                    D:   consD,
                    Ropen:0,
                    Rdec:0,
                    }
        self.occ = self.init[DRO] + self.init[ORD] + self.init[DRD]
        
    def solve_competition_problem(self, time, initACh):
        """
        シナプス後受容体でのAChとNMBの競合を表すODEを解いて結果をself.solutionに格納します
        """
        # 初期値の設定(ACh放出直後の状況を再現)
        self.init[A] = initACh / R_base
        self.ivp.set_initial_value(self.init)
        # 微分方程式の計算(ACh放出後の競合の様子をシミュレーション)
        self.ivp.set_time_span(0,time)
        self.solution = self.ivp.get_solution(method='LSODA') #max_stepは最大値を正確に得るために設定
        #self.solution = self.ivp.get_solution(max_step=0.01) #max_stepは最大値を正確に得るために設定
        # シミュレーション結果の出力
        return self.solution
        
    def peak_ARA(self):
        """
        ARAの濃度のピーク値を取得します
        """
        timecourse = self.ivp.get_timecourse(ARA)        
        if isinstance(timecourse, np.ndarray):
            peak = timecourse.max()
        else:
            raise TypeError("対象の時系列データがndarrayで与えられていません")
        
        return peak

    def peak_Ropen(self):
        """
        ARAの濃度のピーク値を取得します
        """
        timecourse = self.ivp.get_timecourse(Ropen)        
        if isinstance(timecourse, np.ndarray):
            peak = timecourse.max()
        else:
            raise TypeError("対象の時系列データがndarrayで与えられていません")
        
        return peak

    def twitch_height(self):
        """
        Twitchの高さを計算します
        """
        gamma = self.gamma
        ARA50 = self.ARA50       
        peakARA = self.peak_Ropen() # ARA → Ropen
        twitch = np.power(peakARA, gamma) / (np.power(peakARA, gamma)+ np.power(ARA50, gamma))
        return twitch
    
    def final_value(self):
        """
        最終の値を取得します
        """
        last = self.ivp.get_lastvalue_bydict()
        return last


def estimate_ara50_gamma(model, occ50, Twmax, initACh):
    
    if isinstance(model, PostSynapticModel):
        R_total = model.R_total
        K_A1 = model.k_dissA1/model.k_assocA1
        K_A2 = model.k_dissA2/model.k_assocA2
        K_D1 = model.k_dissD1/model.k_assocD1
        K_D2 = model.k_dissD2/model.k_assocD2
    
    # まず, 占拠率がocc50となる濃度Dnmb50を求める
    def func(consD):
        
        common_denom = (consD*K_A1 + K_D1*K_A1)*(consD*K_A2 + K_D2*K_A2)
        DRO = R_total*consD*K_A1*K_A2*K_D2/common_denom
        ORD = R_total*consD*K_A1*K_A2*K_D1/common_denom
        DRD = R_total*(consD**2)*K_A1*K_A2/common_denom
        occ = DRO + ORD + DRD
        return occ - occ50
    
    from scipy.optimize import root    
    Dnmb50 = root(func, [0]).x[0] * R_base

    # この濃度下でのARAのピーク値を計算
    model.set_steady_state(Dnmb50)
    model.solve_competition_problem(time=5, initACh=initACh)
    ARA50 = model.peak_ARA() * R_base    
    #print(ARA50)

    #　ARAmaxを用いてgammaを計算
    model.set_steady_state(0)
    model.solve_competition_problem(time=5, initACh=initACh)
    ARAmax = model.peak_ARA() * R_base    
    gamma_A = - np.log((1-Twmax)/Twmax)/np.log(ARAmax/ARA50)
    #print(gamma_A)
    
    return ARA50, gamma_A

    

##################################################################
if __name__ == "__main__":

    # モデルを作成
    kinetics = {'k_dissA1' : 1.8*(10**4),  #[1/s]
                'k_dissA2' : 1.8*(10**4), #[1/s]
                'k_assocA1': 1.1*(10**10), #[1/Ms]
                'k_assocA2': 1.1*(10**10), #[1/Ms]
                'k_dissD1' : 4*10,  #[1/s]
                'k_dissD2' : 4*10,  #[1/s]
                'k_assocD1': 4.0*(10**6)*2*10, #[1/Ms]
                'k_assocD2': 4.0*(10**8)*2*10, #[1/Ms]
                'k_hydrol' : 12000*100, # 12000, #[1/s]
                'k_close'    : 1200, #[1/s]
                'k_open'     : 50000, #[1/s]
                'k_dpls'   : 26,
                'k_dmns'   : 0.13,
                }
    patient =  {'R_total'  : 7.75*(10**-5), #[M]
                'gamma'    : 4.167,  #4.732,
                'ARA50'    : 1.45*(10**-9),  #9.524*(10**-9), #[M]
                } 
    post_syn = PostSynapticModel(kinet=kinetics, pnt=patient)

    # 筋弛緩薬の結合
    consD  =  0 *(10**-6)  #[M]
    post_syn.set_steady_state(consD)

    # ACh放出後の計算
    initACh= 7.75*(10**-6) #[M]
    sol = post_syn.solve_competition_problem(time=15, initACh=initACh)
    twitch = post_syn.twitch_height()
    
    # 結果の表示
    print(f'initial cond = \n{post_syn.init}')    
    print(f'twitch = {twitch}')

    # グラフのプロット
    def write_graph(sol, filename):
        plt.clf()
        plt.plot(sol.t, sol.y[0]*10, label='Free ACh', ls='--', lw=3, c='b')
        plt.plot(sol.t, sol.y[1], label='ARO', ls='-', lw=3, c='r')
        plt.plot(sol.t, sol.y[2], label='ORA', ls='--', lw=3, c='r')
        plt.plot(sol.t, sol.y[3], label='ARA', ls=':', lw=3, c='k')
        plt.plot(sol.t, sol.y[4], label='ORO', ls='-', lw=3, c='c')
        plt.plot(sol.t, sol.y[5], label='DRO', ls='--', lw=3, c='c')
        plt.plot(sol.t, sol.y[6], label='DRD', ls=':', lw=3, c='c')
        plt.plot(sol.t, sol.y[7], label='ARD', ls='--', lw=3, c='m')
        plt.plot(sol.t, sol.y[8], label='DRA', ls='-', lw=3, c='m')
        plt.plot(sol.t, sol.y[10]*100, label='R*', ls='-', lw=3, c='k')

    
        plt.legend(loc='upper right', fontsize=15)
        plt.margins(0)
        plt.xlabel('Time / ms', fontsize=15)
        plt.ylabel('Fraction of AChR',fontsize=15)
        #plt.xlim(0,0.01)
        plt.ylim(0,1)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.tick_params(direction='inout', width=1, length=8)
        #plt.savefig(filename)

    write_graph(sol, "temp")
    
    ARA50, gamma = estimate_ara50_gamma(model=post_syn, 
                                        occ50=0.875, 
                                        Twmax=0.9995,
                                        initACh=7.75*10**-6)
    print(f'Assumed ARA50={ARA50:.3e}')
    print(f'Assume gamma={gamma:.3f}')
    
    
    