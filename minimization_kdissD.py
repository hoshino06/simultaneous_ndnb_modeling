# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 15:06:37 2021
@author: hoshino
k_dissDがサイト1とサイト2で異なる場合
"""
import numpy as np
from scipy.optimize import minimize
from modules.concn_effect_relationship import concentration_effect_relationship

# モデルの構造の選択
MODEL_TYPE = {'C':'Cyclic', 'R':'Reciprocal', 'B':'BindingModel'}['R']


# 臨床データ
EC50 = {'Cisatracurium': (153/1243.5*10**-6, 33/1243.5*10**-6/np.sqrt(12)),
        'Vecuronium'   : (165/637.7*10**-6,  66/637.7*10**-6/np.sqrt(12) ),
        'Rocuronium'   : (823/609.7*10**-6, 157/609.7*10**-6/np.sqrt(8)  )}

GAMMA = {'Cisatracurium': (6.9,  1.3/np.sqrt(12) ), 
         'Vecuronium'   : (7.6,  3.8/np.sqrt(12) ),
         'Rocuronium'   : (4.79, 1.7/np.sqrt(8)  )}

IC50  = {'Cisatracurium': (10*10**-9, 1*10**-9/np.sqrt(4) ),
         'Vecuronium':    (15*10**-9, 2*10**-9/np.sqrt(4) ),
         'Rocuronium':    (17*10**-9, 2*10**-9/np.sqrt(4) )}

NH    = {'Cisatracurium': (1.02, 0.09/np.sqrt(4) ),
         'Vecuronium':    (1.03, 0.12/np.sqrt(4) ),
         'Rocuronium':    (0.67, 0.05/np.sqrt(4) )}

# FF    = {'Cisatracurium': 0.62,
#          'Vecuronium':    0.31,
#          'Rocuronium':    0.54}

# AChの速度定数
#k_dissA  = 1.0*(10**4)
#k_assocA = 1.0*(10**8)*1000
#alpha    = 4871*1000 #1.2*10**3*1000
#beta     = 106373*1000 #5.0*10**4*1000

##########################
## 評価関数
##########################
DEBUG_fun = False # 関数funのデバッグ用のフラグ

def fun(x, log_file=None, param_file=None):
    """
    評価関数
    """
    log = "Parameters:\n"

    ###########################
    # 入力 x をパラメータに変換
    ###########################
    ARA50  = 10**x[0]
    gammaA =    x[1]
    log += f'ARA50={ARA50:.2e}, gammaA={gammaA:.2f} \n'

    KD1 = {'Cisatracurium': 10**x[2], 'Vecuronium': 10**x[4], 'Rocuronium': 10**x[6]}
    KD2 = {'Cisatracurium': 10**x[3], 'Vecuronium': 10**x[5], 'Rocuronium': 10**x[7]} 

    if MODEL_TYPE == 'Cyclic' or MODEL_TYPE == 'Reciprocal' :

        k_dissD = {'Cisatracurium': [10**x[8], 10**x[15]],
                   'Vecuronium'   : [10**x[9], 10**x[16]],
                   'Rocuronium'   : [10**x[10],10**x[17]]}

        k_assocA = 10**x[11]
        k_dissA  = 10**x[12]
        k_close    = 10**x[13]
        k_open     = 10**x[14]  
        log += f'k_dissA={k_dissA:.2e},  k_assocA={k_assocA:.2e}, K_A={k_dissA/k_assocA:.2e}\n'
        log += f'k_close={k_close:.2e},  k_open={k_open:.2e} \n'
        k_hydrol = 12000
        log += f'k_hydrol={k_hydrol:.2e} \n'
        
    else: # BindingModel
        k_dissD = {'Cisatracurium': 1,
                   'Vecuronium'   : 1,
                   'Rocuronium'   : 1} #dummy
        k_assocA = 1.1*10**8 # dummy
        k_dissA  = 1.8*10**4 # dummy
        k_hydrol = None
        k_close  = None
        k_open   = None          

    if MODEL_TYPE == 'Cyclic':
        k_assocA_ast = 10**x[18] #k_assocA
        k_dissA_ast  = 10**x[19] #k_dissA
        log += f'k_dissA={k_dissA_ast:.2e},  k_assocA={k_assocA_ast:.2e}, K_A={k_dissA_ast/k_assocA_ast:.2e} \n'
    else:
        k_assocA_ast = None
        k_dissA_ast  = None   

    #############################
    ## NMBD毎の評価関数の値の計算
    #############################
    error_ec50  = []
    error_gamma = []
    error_ic50  = []
    error_nH    = []

    for NMBD in ['Cisatracurium', 'Vecuronium', 'Rocuronium']:
        
        log += f'{NMBD}: \n'
        log += f'KD1={KD1[NMBD]:.2e}, KD2={KD2[NMBD]:.2e}, kdissD1={k_dissD[NMBD][0]:.1f}, kdissD2={k_dissD[NMBD][1]:.1f}, mu={KD2[NMBD]/KD1[NMBD]:.4f} \n'

        # パラメータ
        parameters = {
            # kinetic parameters
            'k_dissD1' : k_dissD[NMBD][0],  #[1/s]
            'k_dissD2' : k_dissD[NMBD][1],  #[1/s]
            'k_assocD1': k_dissD[NMBD][0]/KD1[NMBD], #[1/Ms]
            'k_assocD2': k_dissD[NMBD][1]/KD2[NMBD], #[1/Ms]
            'k_dissA1' : k_dissA, #1.8*10**4*10,  #[1/s]
            'k_dissA2' : k_dissA, #1.8*10**4*10, #[1/s]
            'k_assocA1': k_assocA, #1.1*10**8*100, #[1/Ms]
            'k_assocA2': k_assocA, #1.1*10**8*100, #[1/Ms]
            'k_close'  : k_close, #1200*30, #[1/s]
            'k_open'   : k_open,  #50000, #[1/s]
            'k_dpls'   : 26,
            'k_dmns'   : 0.13,  
            'k_dissA_ast' : k_dissA_ast, #[1/s]
            'k_assocA_ast': k_assocA_ast, #[1/Ms]
            'k_hydrol' : k_hydrol,
            # other parameters
            'R_total'  : 7.75*(10**-5), #[M]
            'gamma'    : gammaA,
            'ARA50'    : ARA50, #[M] 
            } 

        ### In vivo simulation
        if DEBUG_fun == True:
            print('Start InVivo Simulation')

        D = 10**np.linspace(-8.0 , -4, 51)
        (ec50,gamma), cod, _ = \
            concentration_effect_relationship('InVivo', MODEL_TYPE, parameters,
                                              free_fraction=1.0, fitting=True, 
                                              d_list=D, cal_time=5.0 )
        
        if cod <  0.95 and DEBUG_fun == True:
            print(f'coeffieicnt of determination = {cod}: {NMBD}, in vivo')

        EC50_M, EC50_SE = EC50[NMBD]
        error_ec50.append( (ec50 - EC50_M)/(EC50_SE*1.96) ) 
        
        GAMA_M, GAMA_SE = GAMMA[NMBD]
        error_gamma.append( (gamma - GAMA_M)/(GAMA_SE*1.96) )

        log += f'  EC50={ec50:.2e}, Gm={gamma:.2f} \n'

        # In vitro simulation                    
        if DEBUG_fun == True:
            print('Start InVitro Simulation')

        D = 10**np.linspace(-10 , -5, 51)
        (ic50,nH), cod, _ = \
            concentration_effect_relationship('InVitro', MODEL_TYPE, parameters,
                                              free_fraction=1.0, fitting=True, 
                                              d_list=D, cal_time=5.0 )

        if cod <  0.96 and DEBUG_fun == True:
            print(f'coeffieicnt of determination = {cod}: {NMBD}, in vitro')

        if cod < 0.96:
            pnl_cod = 0
            print(f'penalty for cod added: pnl_cod={pnl_cod:.1f}')
        else:
            pnl_cod = 0
                
        IC50_M, IC50_SE = IC50[NMBD]
        error_ic50.append( (ic50 - IC50_M)/(IC50_SE*1.96) ) 
        
        NH_M, NH_SE = NH[NMBD]
        error_nH.append( (nH - NH_M)/(NH_SE*1.96) )

        log += f'  IC50={ic50:.2e}, nH={nH:.2f} \n'
    
    # 評価関数
    error_sum_ec50 = np.sum(np.power(error_ec50,2))
    error_sum_gamE = np.sum(np.power(error_gamma,2))
    error_sum_ic50 = np.sum(np.power(error_ic50,2))
    error_sum_gamI = np.sum(np.power(error_nH,2)) 
    
    val =(error_sum_ec50 + error_sum_gamE +  error_sum_ic50 + error_sum_gamI)/(3*4)

    W = 0.5  # 原稿のWと値が異なることに注意(line 211を参照) 2023/10/24
    pnl = W/2 *( np.power( np.log10(1.1*10**8)-np.log10(k_assocA), 2) \
               +np.power( np.log10(1.8*10**4)-np.log10(k_dissA), 2) ) 

    log += '---\nObjective Function:\n'
    log += f' Error = {val:.4f}, Pnl = {pnl:.4f}, Total = {val+pnl:.3f} (W={W:.2f}) \n'
    log += f' (ec50:{error_sum_ec50:.2f}, gammaE:{error_sum_gamE:.2f},'
    log += f' ic50:{error_sum_ic50:.2f}, gammaI:{error_sum_gamI:.2f})'

    print('============\n' + log)
    
    if log_file:
        print(log, file=log_file)

    if param_file: 
        NDNB_param =  'KD1 =' + repr(KD1) + '\n' 
        NDNB_param += 'KD2 =' + repr(KD2) + '\n'
        NDNB_param += 'k_dissD =' + repr(k_dissD) + '\n'
        parameters['k_hydrol'] = k_hydrol
        other_params ='parameters =' + repr(parameters) + '\n'
        print(NDNB_param + other_params, file=param_file)

    return val + W*pnl + pnl_cod



##################################################################################    
## 最適化の計算
############################################################################
if __name__ == "__main__":
    
    # 初期値
    x0 = [ np.log10(8.50*10**-8), 13.0, #ARA50 , GAMMA_A
           np.log10(1.04*10**-8), np.log10(7.33*10**-7), # KD1,KD2 (Cisatracurium)
           np.log10(1.41*10**-8), np.log10(3.83*10**-6), # KD1,KD2 (Vecuronium)
           np.log10(1.50*10**-8), np.log10(8.37*10**-8), # KD1,KD2 (Rocuronium)
           np.log10(2), np.log10(2), np.log10(40),       # k_diss 
           np.log10(3.0*10**10), np.log10(1.8*10**4),    # k_assocA, k_dissA
           np.log10(8.0*10**4), np.log10(8.0*10**4),     # k_close,  k_open
           np.log10(3.0*10**10), np.log10(1.8*10**4),    # k_assocA*, k_dissA*
           #np.log10(12000) # k_hydrol
           ]
    
    # cyclic: completed 2022-03-09 08:48:57.710029
    x0 = [ np.log10(4.50e-08), 8.89, #ARA50 , GAMMA_A
            np.log10(9.669e-09), np.log10(5.057e-06), # KD1,KD2 (Cisatracurium)
            np.log10(1.596e-08), np.log10(2.789e-06), # KD1,KD2 (Vecuronium)
            np.log10(1.547e-08), np.log10(2.306e-07), # KD1,KD2 (Rocuronium)
            np.log10(2), np.log10(2), np.log10(2),       # k_diss 
            np.log10(1.71e+10), np.log10(1.78e+04),    # k_assocA, k_dissA
            np.log10(9.78e+04), np.log10(1.39e+05),     # k_close,  k_open
            np.log10(3.47e+11), np.log10(2.72e+04),    # k_assocA*, k_dissA*
            #np.log10(12000) # k_hydrol
            ]
    
    # completed 2022-07-13 13:14:54.539597
    x0 = [ 
            np.log10(1.65e-08), 9.00, #ARA50 , GAMMA_A
            np.log10(9.77e-09), np.log10(1.08e-07), # KD1,KD2 (Cisatracurium)
            np.log10(1.63e-08), np.log10(2.81e-07), # KD1,KD2 (Vecuronium)
            np.log10(9.08e-09), np.log10(1.32e-08), # KD1,KD2 (Rocuronium)
            np.log10(2.0), np.log10(2.0), np.log10(60),       # k_diss 
            np.log10(1.13e+10), np.log10(1.49e+04),    # k_assocA, k_dissA
            np.log10(1.07e+05), np.log10(3.06e+04),     # k_close,  k_open
            np.log10(2.0), np.log10(2.0), np.log10(60),   # k_diss (site #2) 
            np.log10(5.34e+11), np.log10(6.82e+04),   # k_assocA*, k_dissA*
            ]      

    # minimization (cyclic) の結果
    x0 = [ 
            np.log10(1.82e-08), 9.04, #ARA50 , GAMMA_A
            np.log10(1.02e-08), np.log10(3.60e-05), # KD1,KD2 (Cisatracurium)
            np.log10(1.63e-08), np.log10(1.58e-06), # KD1,KD2 (Vecuronium)
            np.log10(1.22e-08), np.log10(1.76e-07), # KD1,KD2 (Rocuronium)
            np.log10(4.0), np.log10(5.4), np.log10(61.5),       # k_diss 
            np.log10(5.9e+9), np.log10(2.62e+03),    # k_assocA, k_dissA
            np.log10(4.24e+04), np.log10(1.06e+04),     # k_close,  k_open
            np.log10(4.0), np.log10(5.4), np.log10(61.5),   # k_diss (site #2) 
            np.log10(923913043478), np.log10(1.70e+04),   # k_assocA*, k_dissA*
            ]      

    if MODEL_TYPE == "BindingModel":

        x0 = [ np.log10(1.50*10**-10), 5.0, #ARA50 , GAMMA_A
               np.log10(1.04*10**-8), np.log10(7.33*10**-7), # KD1,KD2 (Cisatracurium)
               np.log10(1.41*10**-8), np.log10(3.83*10**-7), # KD1,KD2 (Vecuronium)
               np.log10(1.50*10**-8), np.log10(8.37*10**-4), # KD1,KD2 (Rocuronium)
               ]

    if MODEL_TYPE == "Reciprocal":
        
            # reciprocal completed 2022-03-06 23:23:27.692762
        x0 = [ np.log10(6.09e-07), 14.53, #ARA50 , GAMMA_A
                np.log10(1.116e-08), np.log10(9.436e-07), # KD1,KD2 (Cisatracurium)
                np.log10(1.851e-08), np.log10(8.849e-07), # KD1,KD2 (Vecuronium)
                np.log10(1.916e-08), np.log10(1.847e-07), # KD1,KD2 (Rocuronium)
                np.log10(2.6), np.log10(1.8), np.log10(50.3),       # k_diss 
                np.log10(2.79e+11 ), np.log10(1.57e+03),    # k_assocA, k_dissA
                np.log10(1.31e+05), np.log10(2.80e+06),     # k_close,  k_open
                np.log10(2.6), np.log10(1.8), np.log10(50.3),       # k_diss (site#2)
                ]

        # minimization (reciprocal) の結果 
        x0 = [ np.log10(2.08e-07), 9.14, #ARA50 , GAMMA_A
               np.log10(9.57e-09), np.log10(6.75e-06), # KD1,KD2 (Cisatracurium)
               np.log10(1.58e-08), np.log10(2.76e-06), # KD1,KD2 (Vecuronium)
               np.log10(1.23e-08), np.log10(1.37e-07), # KD1,KD2 (Rocuronium)
               np.log10(2.6), np.log10(1.9), np.log10(64.3),       # k_diss 
               np.log10(4.43e+02/(1.58e-08) ), np.log10(4.43e+02),    # k_assocA, k_dissA
               np.log10(2.22e+07), np.log10(1.48e+010),     # k_close,  k_open
               np.log10(2.6), np.log10(1.9), np.log10(64.3),       # k_diss (site#2)
               ]
    
    # 関数funのデバッグ用(funの定義の直前でフラグを定義)
    if DEBUG_fun == True:
        fun(x0)
        import sys
        sys.exit() # デバッグ時はここで終了
    
    
    # 最適化計算の条件
    Method = 'Nelder-Mead' #'BFGS' #'Powell' 
    Tolerance = 1e-2
    
    # ログの開始
    log_file = f"minimization_{MODEL_TYPE.lower()}_w_kdissD.log" 
    f = open(log_file, mode='a')
    import datetime
    print('=====================================================\n' \
          + f'Optimization start at {datetime.datetime.now()}\n' \
          + f'with {Method} method with tol={Tolerance:.2e}\n'
          +'====================================================', file=f)
    
    print('Initial conditions', file=f)
    #fun(x0, log_file=f) #初期値での評価を記録
    print(f'{x0}', file=f)
    
    # 最適化計算
    res = minimize(fun, x0, method=Method, tol=Tolerance)
            
    # ログの終了とパラメータの出力
    print('========================\n',
          f' completed {datetime.datetime.now()}')
    
    print('========================\n'\
          +f' completed {datetime.datetime.now()}',file=f)
    
        
    print('Optimization results\n---')
    
    res_file = f"minimized_parameters_{MODEL_TYPE.lower()}_w_kdissD.py"
    ff = open(res_file, mode='a')
    print(f'# completed {datetime.datetime.now()}',file=ff)
    
    fun(res.x, log_file=f, param_file=ff)
    print('\n', file=f)
    
    f.close()
    ff.close()
    
