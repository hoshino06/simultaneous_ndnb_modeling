# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 11:32:25 2022
@author: hoshino
"""
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# モデルの構造の選択
MODEL_TYPE = {'C':'Cyclic', 'R':'Reciprocal', 'B':'BindingModel'}['B'] 

# In Vivo と In Vitro の選択
InVivo_InVitro = ['InVivo', 'InVitro'][1]

result = pd.read_csv(f'fig_{MODEL_TYPE.lower()}/parameter_sweep_{InVivo_InVitro.lower()}.csv')    

##############################################
# プロットの準備
k_diss_D_list = [1,10,60]


legend = {1: "k-", 
          10: "k--", 
          60: "k:"}

color = {'Cisatracurium': 'r',
         'Vecuronium'   : 'g',
         'Rocuronium'   : 'b'}

if MODEL_TYPE == 'Cyclic':
    from minimized_parameters_cyclic import KD1, KD2
    EC50 = {'Cisatracurium': 1.29e-07, 'Vecuronium': 2.49e-07, 'Rocuronium': 1.34e-06 }
    IC50 = {'Cisatracurium': 9.68e-09, 'Vecuronium': 1.49e-08, 'Rocuronium': 1.75e-08 }
    gammaE = {'Cisatracurium': 6.92, 'Vecuronium': 5.81, 'Rocuronium': 5.00 }
    gammaI = {'Cisatracurium': 0.96, 'Vecuronium': 0.98, 'Rocuronium': 0.68 }

elif MODEL_TYPE == 'Reciprocal':
    from minimized_parameters_reciprocal import KD1, KD2
    EC50 = {'Cisatracurium': 1.32e-07, 'Vecuronium': 2.52e-07, 'Rocuronium': 1.21e-06 }
    IC50 = {'Cisatracurium': 9.51e-09, 'Vecuronium': 1.55e-08, 'Rocuronium': 1.77e-08 }
    gammaE = {'Cisatracurium': 7.16, 'Vecuronium': 6.55, 'Rocuronium': 5.23 }
    gammaI = {'Cisatracurium': 0.99, 'Vecuronium': 1.00, 'Rocuronium': 0.69 }

elif MODEL_TYPE == 'BindingModel':
    from minimized_parameters_bindingmodel import KD1, KD2
    EC50 = {'Cisatracurium': 1.44e-07, 'Vecuronium': 2.87e-07, 'Rocuronium': 1.08e-06 }
    IC50 = {'Cisatracurium': 8.56e-09, 'Vecuronium': 1.45e-08, 'Rocuronium': 1.88e-08 }
    gammaE = {'Cisatracurium': 7.15, 'Vecuronium': 6.90, 'Rocuronium': 4.05 }
    gammaI = {'Cisatracurium': 1.195, 'Vecuronium': 1.14, 'Rocuronium': 1.00 }


############################################################################
# gamma のプロット
############################################################################
plt.clf()

# 設定    
if InVivo_InVitro == 'InVivo':
    index = 'gammaE'
    res_index = gammaE
    plt.ylabel('$\gamma_\mathrm{E}$',fontsize=13)
    if MODEL_TYPE == 'BindingModel': plt.ylim(3,8) 
    else: plt.ylim(0,20)    

if InVivo_InVitro == 'InVitro':
    index = 'gammaI'
    res_index = gammaI
    plt.ylabel('$\gamma_\mathrm{I}$',fontsize=13)
    if MODEL_TYPE == 'BindingModel': plt.ylim(0.9,1.3) 
    else: plt.ylim(0.4,1.4)

plt.xlabel('$log_{10}({K_\mathrm{D2}/K_\mathrm{D1}})$',fontsize=13)


# シミュレーション結果のプロット
if MODEL_TYPE == 'Cyclic' or MODEL_TYPE == 'Reciprocal':

    for k_dissD in k_diss_D_list:
    
        res = result.query(f'k_dissD == {k_dissD}')
        plt.plot(np.log10(res['KD2']/res['KD1']), 
                 res['gamma'], legend[k_dissD],
                 label = '$k_\mathrm{dissD} = '+f'{k_dissD}$'
                 )

elif MODEL_TYPE == 'BindingModel':
    plt.plot(np.log10(result['KD2']/result['KD1']), 
             result['gamma'], 'k-')


# 3種類のNDNBの値
for NDNB in ['Cisatracurium', 'Vecuronium', 'Rocuronium']:
    plt.scatter(np.log10(KD2[NDNB]/KD1[NDNB]), res_index[NDNB], c=color[NDNB], s=50)

    if MODEL_TYPE == 'BindingModel':
        plt.scatter(np.log10(KD2[NDNB]/KD1[NDNB]), res_index[NDNB]/KD1[NDNB], 
                    c=color[NDNB], s=50, label=NDNB)    

    
# if MODEL_TYPE == 'BindingModel':
#     plt.ylim(3,8)


#plt.yticks(np.arange(0.95, 1.201, step=0.05))
#plt.tick_params(direction='in', top=True, right=True, 
#                labelbottom=True)
plt.margins(0)

if MODEL_TYPE == 'Cyclic' or MODEL_TYPE == 'Reciprocal':    
    plt.legend(loc='upper right', #bbox_to_anchor=(0.05, 0.2), 
               fontsize=13 )

if MODEL_TYPE == 'BindingModel':    
    plt.legend(loc='upper right', #bbox_to_anchor=(0.05, 0.2), 
               fontsize=13 )

plt.savefig(f'fig_{MODEL_TYPE.lower()}/parameter_sweep_{index}.pdf', bbox_inches="tight")
plt.show()

#################################################################################
# C50のプロット
################################################################################
plt.clf()

#
if InVivo_InVitro == 'InVivo':
    index = 'EC50'
    res_index = EC50
    plt.ylabel(r'${\rm EC50} \, /\, K_\mathrm{D1}$',fontsize=13)
    if MODEL_TYPE == 'BindingModel': plt.ylim(0,60)
    else: plt.ylim(0,120)

if InVivo_InVitro == 'InVitro':
    index = 'IC50'
    res_index = IC50
    plt.ylabel(r'${\rm IC50} \, /\, K_\mathrm{D1}$',fontsize=13)
    if MODEL_TYPE == 'BindingModel': plt.ylim(0.2,1.2) 
    else: plt.ylim(0,4)

# シミュレーション結果のプロット
if MODEL_TYPE == 'Cyclic' or MODEL_TYPE == 'Reciprocal':

    for k_dissD in k_diss_D_list:
    
        res = result.query(f'k_dissD == {k_dissD}')   
        plt.plot(np.log10(res['KD2']/res['KD1']), 
                 res['C50']/res['KD1'], 
                 legend[k_dissD],
                 label = '$k_\mathrm{dissD} = '+f'{k_dissD}$')

elif MODEL_TYPE == 'BindingModel':
    plt.plot(np.log10(result['KD2']/result['KD1']), 
             result['C50']/result['KD1'], 'k-')


# 3種類のNDNBの値
for NDNB in ['Cisatracurium', 'Vecuronium', 'Rocuronium']:

    print(np.log10(KD2[NDNB]/KD1[NDNB]))
    plt.scatter(np.log10(KD2[NDNB]/KD1[NDNB]), res_index[NDNB]/KD1[NDNB], 
                c=color[NDNB], s=50)    
    
    if MODEL_TYPE == 'BindingModel':
        plt.scatter(np.log10(KD2[NDNB]/KD1[NDNB]), res_index[NDNB]/KD1[NDNB], 
                    c=color[NDNB], s=50, label=NDNB)    


plt.xlabel('$log_{10}({K_\mathrm{D2}/K_\mathrm{D1}})$',fontsize=13)    


#plt.ylim(5,45)
plt.margins(0)

if MODEL_TYPE == 'BindingModel':    
    plt.legend(loc='lower right', #bbox_to_anchor=(0.9, 0.7),
               fontsize=13 )


if MODEL_TYPE == 'Cyclic' or MODEL_TYPE == 'Reciprocal':    
    plt.legend(loc='upper right', #bbox_to_anchor=(0.99, 0.89),
               fontsize=13 )

# 図の出力
plt.savefig(f'fig_{MODEL_TYPE.lower()}/parameter_sweep_{index}.pdf', bbox_inches="tight")
plt.show()

