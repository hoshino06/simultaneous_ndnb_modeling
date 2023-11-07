# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 17:16:12 2021

@author: hoshino
"""
import numpy as np
from scipy.optimize import least_squares

def sigmoid_fitting(concentration, effect, init, fitted_model=False):
    """
    Least Squeresにより濃度反応関係のデータからPDモデルパラメータを決定します
    
    Args:
        concentraion (list): フィッティングに使う濃度データ
        effect (list): フィッティングに使う反応のデータ
        init: (log10(EC50) と gamma の初期推定値

    Returns:
        (EC50, gamma): フィッティング結果
        cod: Coefficient of determination
    """        

    # Residuals を計算する関数
    # (x: パラメータ, t: 独立変数, y:縦続変数)
    def residu(x, t, y):

        D = t
        ec50  = 10**x[0]
        gamma = x[1]
        yy = 1 - ( D**gamma / ( D**gamma + ec50**gamma ) )
        
        return yy-y

    # Least Squeresによるパラメータ同定
    res_lsq = least_squares(residu, init, 
                            args=(concentration, effect)) 
    
    # 結果の取り出し
    ec50  = 10**res_lsq['x'][0]
    gamma = res_lsq['x'][1]
    residuals   = res_lsq['fun']
    cod = 1 - np.sum(residuals**2) / np.sum( (effect - np.mean(effect))**2)

    if fitted_model == False:
        return [ec50,gamma], cod

    else:
        D = concentration
        model = 1 - D**gamma/( D**gamma + ec50**gamma )
        return [ec50,gamma], cod, model
    


##################################################################
if __name__ == "__main__":
    
    D = 10**np.linspace(-10, -5, 101) # 濃度反応関係を計算する際の濃度刻み
    
    KD1 = 1.0*10**-8; KD2 = 1.0*10**-8
    model = KD1*KD2/(KD1*KD2+KD1*D+KD2*D+D*D)

    # # IC50とnHのフィッティング
    init = [ np.amin(np.log10([KD1,KD2])), 2] 
    (ec50, gamma), cod, fitted = sigmoid_fitting(D, model, init, fitted_model=True)
    print(f'EC50 = {ec50:.3e}, gamma = {gamma:.3f}')
    print(f'coeffieicnt of determination = {cod}')
    
    fitted = 1 - D**gamma/( D**gamma + ec50**gamma )
        
    def plot_graph():
    
        import matplotlib.pyplot as plt
        
        plt.plot(np.log10(D), model, 'r', label='data')
        plt.plot(np.log10(D), fitted, 'b', label='fitted')
    
        plt.xlabel('$D$',fontsize=15)
        plt.ylabel('$I$',fontsize=15)
        #plt.ylim(0.3,1.1)
        #plt.yticks(np.arange(0, 1.1, step=0.5))
        plt.tick_params(direction='in', top=True, right=True, labelbottom=True)
        plt.margins(0)
    
        plt.legend(loc='upper right', bbox_to_anchor=(0.95,1.04), #凡例の位置
                frameon=False, #凡例の囲みは不必要
                fontsize=10, #凡例のフォントサイズ
                ncol=1, labelspacing=0.2 #凡例の並び方
                )
    
        plt.show()
    
    plot_graph()
