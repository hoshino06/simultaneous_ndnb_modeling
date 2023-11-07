# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 23:59:36 2022
@author: hoshino
"""
import numpy as np
import matplotlib.pyplot as plt
from modules.concn_effect_relationship import concentration_effect_relationship

# モデルの構造の選択
MODEL_TYPE = {'C':'Cyclic', 'R':'Reciprocal', 'B':'BindingModel'}['C'] # select 0 or 1 or 2
NMBD = {'C':'Cisatracurium', 'V':'Vecuronium', 'R':'Rocuronium'}['C']

# 最適化結果のパラメータの読み込み
if MODEL_TYPE == 'Cyclic':    
    from minimized_parameters_cyclic import parameters, k_dissD, KD1, KD2
elif MODEL_TYPE == 'Reciprocal':    
    from minimized_parameters_reciprocal import parameters, k_dissD, KD1, KD2

parameters['k_dissD1'] = k_dissD[NMBD]
parameters['k_dissD2'] = k_dissD[NMBD]
parameters['k_assocD1'] = k_dissD[NMBD]/KD1[NMBD]
parameters['k_assocD2'] = k_dissD[NMBD]/KD2[NMBD]

D = 10**np.linspace(-10 , -5, 101)
(c50,gamma), cod, model \
    = concentration_effect_relationship('InVivo', MODEL_TYPE,
                                        parameters,
                                        free_fraction=1.0,
                                        fitting=True,
                                        d_list=D)
    
print('calculation completed')
print(f'C50={c50:.3e}, gamma={gamma:.3f}')    
        

# 濃度反応関係のプロット
plt.plot(np.log10(model['[D]']), model['res'], 
         'rx', label = r'Twitch')

plt.plot(np.log10(model['[D]']), 
         1 - model['[D]']**gamma / (model['[D]']**gamma + c50**gamma), 
         'g', label = 'Fitted')
   
plt.plot(np.log10(model['[D]']), 1-model['Occ'], 
         'b', label = 'Vacancy')
   
plt.plot(np.log10(model['[D]']), model['Rop']/np.max(model['Rop']), 
         'm', label = '$[ARA*]/[ARA*]_\mathrm{max}$')
   

plt.xlabel('$log_{10}([D] / \mathrm{M})$',fontsize=13)
plt.ylabel('Twitch Height',fontsize=13)
plt.ylim(0.0,1.0)
#plt.yticks(np.arange(0.95, 1.201, step=0.05))
plt.tick_params(direction='in', top=True, right=True, 
                labelbottom=True)
plt.margins(0)
plt.legend(loc='upper right', bbox_to_anchor=(0.9, 0.9), 
           fontsize=13 )
plt.show()

