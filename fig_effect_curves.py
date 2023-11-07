# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 10:03:11 2022
@author: hoshino
"""
import numpy as np
import matplotlib.pyplot as plt
from modules.concn_effect_relationship import concentration_effect_relationship

# モデルの構造の選択
MODEL_TYPE = {'C':'Cyclic', 'R':'Reciprocal', 'B':'BindingModel'}['R'] # select 0 or 1 or 2

# In Vivo と In Vitro の選択
InVivo_InVitro = ['InVivo', 'InVitro'][1]

# 最適化結果のパラメータの読み込み
if MODEL_TYPE == 'Cyclic':    
    from minimized_parameters_cyclic import parameters, k_dissD, KD1, KD2
elif MODEL_TYPE == 'Reciprocal':    
    from minimized_parameters_reciprocal import parameters, k_dissD, KD1, KD2
elif MODEL_TYPE == 'BindingModel':    
    from minimized_parameters_bindingmodel import parameters, k_dissD, KD1, KD2

# 凡例
color = {'Cisatracurium': 'r', 
         'Vecuronium'   : 'g',
         'Rocuronium'   : 'b' }


# PKPDの実験データの読み込み
from minimization_ import EC50, GAMMA, IC50, NH

for NMBD in ['Cisatracurium', 'Vecuronium', 'Rocuronium']:
        
    # Ideal Curves
    if InVivo_InVitro == 'InVivo':
        D = 10**np.linspace(-8.0 , -5, 51)
        c50 = EC50[NMBD][0]; gamma = GAMMA[NMBD][0];   

    if InVivo_InVitro == 'InVitro':
        D = 10**np.linspace(-10 , -5, 51)
        c50   = IC50[NMBD][0]; gamma = NH[NMBD][0];   

    ideal = 1 - D**gamma / (D**gamma + c50**gamma)    
    plt.plot(D, ideal, color[NMBD]+'--', label = '')
    

    # 結果の計算
    parameters['k_dissD1'] = k_dissD[NMBD]
    parameters['k_dissD2'] = k_dissD[NMBD]
    parameters['k_assocD1'] = k_dissD[NMBD]/KD1[NMBD]
    parameters['k_assocD2'] = k_dissD[NMBD]/KD2[NMBD]
    
    (c50,gamma), cod, model = \
        concentration_effect_relationship(InVivo_InVitro, MODEL_TYPE, parameters,
                                          free_fraction=1.0, fitting=True, 
                                          d_list=D, cal_time = 5.0 )

    print(f'{c50:.2e}, {gamma:.3f}, {cod:.2f}')

    plt.plot(model['[D]'], 
             model['res'],
             color[NMBD]+'-', label = NMBD)           

    # Confidence interavl of EC50
    if InVivo_InVitro == 'InVivo':
        plt.errorbar(EC50[NMBD][0], 0.5, xerr=EC50[NMBD][1]*1.96, capsize=3, ecolor='black')


if InVivo_InVitro == 'InVivo':
    xlim = [3*10**-8,1*10**-5]
    ylabel = 'Twitch Height'
    
if InVivo_InVitro == 'InVitro':
    xlim = [10**-10, 10**-5]
    ylabel = 'Relative Current'

plt.xscale('log')
plt.xlim(xlim[0], xlim[1])
plt.xlabel('$log_{10}([D] / \mathrm{M})$',fontsize=13)
plt.ylabel(ylabel,fontsize=13)
plt.ylim(-0.05,1.05)
plt.yticks(np.arange(0.0, 1.01, step=0.2))
plt.tick_params(direction='in', top=True, right=True, 
                labelbottom=True)
#plt.margins(0)
plt.legend(loc='upper right', bbox_to_anchor=(1.0, 0.98), 
           fontsize=10 )
plt.savefig(f'fig_{MODEL_TYPE.lower()}/effect_curves_{InVivo_InVitro.lower()}.pdf', bbox_inches="tight")
plt.show()



