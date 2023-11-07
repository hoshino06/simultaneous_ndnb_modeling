# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import matplotlib.pyplot as plt

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
            self.k_dissD1 = kinet['k_dissD1'] *Tbase
            self.k_dissD2 = kinet['k_dissD2'] *Tbase
            self.k_assocD1= kinet['k_assocD1'] *R_base*Tbase
            self.k_assocD2= kinet['k_assocD2'] *R_base*Tbase

        # 患者の個人差に関わるパラメータ
        if pnt:
            self.R_total = pnt['R_total'] / R_base
            self.gamma = pnt['gamma']
            self.ARA50 = pnt['ARA50'] / R_base

    def set_steady_state(self, consD):
        """
        与えられた筋弛緩薬濃度の下での定常状態を計算します. 引数consDに単位 M での数値を与えてください. 
        """
        # 与えられた筋弛緩薬濃度を単位化
        consD = consD / R_base
        
        # 定常状態における濃度を計算
        KD1 = self.k_dissD1/self.k_assocD1
        KD2 = self.k_dissD2/self.k_assocD2

        self.occ = 1- KD1*KD2/(KD1*KD2+KD1*consD+KD2*consD+consD*consD)
        
    def solve_competition_problem(self, time, initACh):
        pass
            
    def peak_Ropen(self):
        """
        ARAの濃度のピーク値を取得します
        """        
        return 1-self.occ

    def twitch_height(self):
        """
        Twitchの高さを計算します
        """
        gamma = self.gamma
        ARA50 = self.ARA50       
        peakARA = self.peak_Ropen() * R_base
        twitch = np.power(peakARA, gamma) / (np.power(peakARA, gamma)+ np.power(ARA50, gamma))
        return twitch
    
    def final_value(self):
        """
        最終の値を取得します
        """
        last = self.ivp.get_lastvalue_bydict()
        return last


##################################################################
if __name__ == "__main__":

    # モデルを作成
    kinetics = {'k_dissA1' : 1.80*(10**4),  #[1/s]
                'k_dissA2' : 1.80*(10**4), #[1/s]
                'k_assocA1': 5.0*(10**10), #[1/Ms]
                'k_assocA2': 5.0*(10**10), #[1/Ms]
                'k_dissD1' : 4*10,  #[1/s]
                'k_dissD2' : 4*10,  #[1/s]
                'k_assocD1': 4.0*(10**6), #[1/Ms]
                'k_assocD2': 4.0*(10**8), #[1/Ms]
                'k_hydrol' : 12000, # 12000, #[1/s]
                'k_close'  : 317000,#24800, #[1/s]
                'k_open'   : 337000, #[1/s]
                'k_dpls'   : 26, #26,
                'k_dmns'   : 0.13,
                'k_dissA_ast' : 1.8*(10**4), #[1/s]
                'k_assocA_ast': 1.8*(10**10), #[1/Ms]
                }
    patient =  {'R_total'  : 7.75*(10**-5), #[M]
                'gamma'    : 4.167,  #4.732,
                'ARA50'    : 1.45*(10**-9),  #9.524*(10**-9), #[M]
                } 
    post_syn = PostSynapticModel(kinet=kinetics, pnt=patient)

    # 筋弛緩薬の結合
    consD  =  3 *(10**-7)  #[M]
    post_syn.set_steady_state(consD)

    # ACh放出後の計算
    initACh= 7.75*(10**-6) #[M]
    sol = post_syn.solve_competition_problem(time=100, initACh=initACh)
    rocc = post_syn.occ
    twitch = post_syn.twitch_height()
    
    # 結果の表示
    print(f'occupancy = {rocc}')
    print(f'twitch = {twitch}')

