# -*- coding: utf-8 -*-
"""
Competition Modelの場合の EC50 と gamma の変化を計算するプログラム
@author: hoshino
"""
import pandas as pd
import numpy as np
from modules.model_cyclic     import PostSynapticModel as CyclicModel
from modules.model_reciprocal import PostSynapticModel as ReciprocalModel
from modules.model_bindingmdl import PostSynapticModel as BindingModel
from modules.utility_tools import sigmoid_fitting

DEBUG = False # デバッグ用のフラグ

def concentration_effect_relationship(invivo_or_invitro,
                                      mdl_type, mdl_param, 
                                      free_fraction=1,
                                      fitting=False, 
                                      d_list = [],
                                      cal_time = None,
                                      ):
    """
    Parameters
    ----------
    invivo_or_invitro: select 'InVivo' or 'InVitro'
    mdl_type:  type of the model ('Cyclic' / 'Reciprocal' / 'BindingModel' )
    mdl_param: parameters of the model
    free_frac: free_fraction of NDNB (default 1)
    fitting:   if True, calculate C50 and gamma

    Returns
    -------
    
    """
    if DEBUG == True:
        print(f'Start calculation of {invivo_or_invitro} effects by {mdl_type} model')

    # モデルとパラメータの設定
    if invivo_or_invitro == 'InVivo':

        initACh = 7.75*(10**-6)
        parameters = mdl_param.copy()
        
    elif invivo_or_invitro == 'InVitro':            

        initACh = 7.75*(10**-3)
        parameters = mdl_param.copy()
        parameters['k_hydrol'] = 0

    MODEL = \
        {'Cyclic':       CyclicModel,
         'Reciprocal':   ReciprocalModel,
         'BindingModel': BindingModel
         }  

    model = MODEL[mdl_type](kinet=parameters, pnt=parameters)


    # 濃度反応関係の計算
    result = pd.DataFrame()
    
    if not d_list == []:
        D_list = d_list
    else:
        D_list = 10**np.linspace(-9 , -4, 101) # 濃度反応関係を計算する際の濃度刻み

    if cal_time == None:
        cal_time = 300

    for d in D_list:

        model.set_steady_state( d*free_fraction ) 
        model.solve_competition_problem(time=cal_time, initACh=initACh)
        Occ    = model.occ
        Rop    = model.peak_Ropen()

        if invivo_or_invitro == 'InVivo':
            res = model.twitch_height()
        elif invivo_or_invitro == 'InVitro':            
            res = model.peak_Ropen()

        # result = result.append({'[D]': d,
        #                         'Occ':Occ,
        #                         'Rop':Rop,
        #                         'res':res}, 
        #                        ignore_index=True) 

        result = pd.concat( [result, pd.DataFrame( [{'[D]': d,
                                'Occ':Occ,
                                'Rop':Rop,
                                'res':res}] ) ], 
                               ignore_index=True ) 

    if fitting != True: 
        return result        

    # フィッティングする場合には(C50,gamma)をフィッティング
    if fitting == True:
        KD1 = mdl_param['k_dissD1']/mdl_param['k_assocD1']
        KD2 = mdl_param['k_dissD2']/mdl_param['k_assocD2']
        init = [ np.log10(np.min([KD1, KD2])), 3] 
        (c50, gamma), cod = sigmoid_fitting(result['[D]'], 
                                            result['res'],
                                            init)
        if cod <  0.95:
            print(f'coeffieicnt of determination = {cod}')
    
        return (c50,gamma), cod, result
    

############################################################################
if __name__ == "__main__":

    ####################
    invivo_or_invitro = ['InVivo', 'InVitro'][1]  # select 0 or 1
    model_type = ['Cyclic', 'Reciprocal', 'BindingModel'][2] # select 0 or 1 or 2

    ####################    
    # パラメータ
    KD1 = 1.82*10**-8
    KD2 = 1.82*10**-8
    k_dissD = 10
    
    ARA50 = 7.18*10**-8*0.01
    gammaA = 10.25
    k_dissA  = 1.8*10**4
    k_assocA = 2.99*10**10
    k_dissA_ast  = k_dissA
    k_assocA_ast = k_assocA
    k_close = 4.44*10**5
    k_open  = 3.86*10**5
    
    # Competition-based Modelの作成
    parameters = {
        # kinetic parameters
        'k_dissA1' : k_dissA, #1.8*10**4*10,  #[1/s]
        'k_dissA2' : k_dissA, #1.8*10**4*10, #[1/s]
        'k_assocA1': k_assocA, #1.1*10**8*100, #[1/Ms]
        'k_assocA2': k_assocA, #1.1*10**8*100, #[1/Ms]
        'k_dissD1' : k_dissD,  #[1/s]
        'k_dissD2' : k_dissD,  #[1/s]
        'k_assocD1': k_dissD/KD1, #[1/Ms]
        'k_assocD2': k_dissD/KD2, #[1/Ms]
        'k_close'  : k_close, #1200*30, #[1/s]
        'k_open'   : k_open,  #50000, #[1/s]
        'k_dpls'   : 26,
        'k_dmns'   : 0.13,  
        'k_dissA_ast' : k_dissA_ast, #[1/s]
        'k_assocA_ast': k_assocA_ast, #[1/Ms]
        # other parameters
        'R_total'  : 7.75*(10**-5), #[M]
        'gamma'    : gammaA,
        'ARA50'    : ARA50, #[M] 
        } 
    

    (c50,gamma), cod, model \
        = concentration_effect_relationship(invivo_or_invitro, model_type,
                                            parameters,
                                            free_fraction=1.0,
                                            fitting=True)
    print('calculation completed')
    print(f'C50={c50:.3e}, gamma={gamma:.3f}')    
            
    
    # 濃度反応関係のプロット
    import matplotlib.pyplot as plt

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
    
