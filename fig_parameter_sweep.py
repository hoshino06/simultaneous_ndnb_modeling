# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 10:59:02 2022
@author: hoshino
"""
import numpy as np
import pandas as pd
from modules.concn_effect_relationship import concentration_effect_relationship

# モデルの構造の選択
MODEL_TYPE = {'C':'Cyclic', 'R':'Reciprocal', 'B':'BindingModel'}['R'] 

# In Vivo と In Vitro の選択
InVivo_InVitro = ['InVivo', 'InVitro'][1]

# 反応の計算時間
cal_time = 5.0 

# 最適化結果のパラメータの読み込み
if MODEL_TYPE == 'Cyclic':
    from minimized_parameters_cyclic import parameters
elif MODEL_TYPE == 'Reciprocal':    
    from minimized_parameters_reciprocal import parameters
elif MODEL_TYPE == 'BindingModel':    
    from minimized_parameters_bindingmodel import parameters


###########################
# KD1とKD2の変化に対する EC50, gammaの変化特性
KD2_LIST = (10**-8)*(10**np.linspace(0,5,100)) 
KD1_LIST = (10**-8)*np.ones(100)

result = pd.DataFrame()


if InVivo_InVitro == 'InVivo':
    D = 10**np.linspace(-8.0 , -4, 101)

if InVivo_InVitro == 'InVitro':
    D = 10**np.linspace(-10 , -5, 101)

if MODEL_TYPE == 'Cyclic' or MODEL_TYPE == 'Reciprocal':
    k_diss_D_list = [1,10,60]
elif MODEL_TYPE == 'BindingModel':
    k_diss_D_list = [9999] # dummy for avoiding not implemented error


for k_dissD in k_diss_D_list:

    for KD1, KD2 in zip(KD1_LIST,KD2_LIST):
        
        parameters['k_dissD1'] = k_dissD
        parameters['k_dissD2'] = k_dissD
        parameters['k_assocD1'] = k_dissD/KD1
        parameters['k_assocD2'] = k_dissD/KD2

        (c50,gamma), cod, _ = \
            concentration_effect_relationship(InVivo_InVitro, MODEL_TYPE, parameters,
                                              free_fraction=1.0, fitting=True,
                                              d_list=D, cal_time = cal_time)
    
        print(f'C50={c50:.3e}, gamma={gamma:.3f}')    
    
        result = result.append(
                {'KD1': KD1,
                 'KD2': KD2,
                 'k_dissD': k_dissD,
                 'muD': KD1/KD2,
                 'C50': c50,
                 'gamma': gamma,
                 'cod': cod,
                  }, 
                ignore_index=True)   

result.to_csv(f'fig_{MODEL_TYPE.lower()}/parameter_sweep_{InVivo_InVitro.lower()}.csv')
