# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 00:16:45 2022
@author: hoshino
"""
import matplotlib.pyplot as plt

###########################################################################
# Cyclic model 

# モデルの構造の選択
MODEL_TYPE = {'C':'Cyclic', 'R':'Reciprocal', 'B':'BindingModel'}['C'] # select 0 or 1 or 2
NMBD = {'C':'Cisatracurium', 'V':'Vecuronium', 'R':'Rocuronium'}['R']

# モデルとパラメータの読み込み
# if MODEL_TYPE == 'Cyclic':    
#     from modules.model_cyclic     import PostSynapticModel 
#     from minimized_parameters_cyclic import parameters, k_dissD, KD1, KD2
# elif MODEL_TYPE == 'Reciprocal':    
#     from modules.model_reciprocal import PostSynapticModel
#     from minimized_parameters_reciprocal import parameters, k_dissD, KD1, KD2

from modules.model_cyclic     import PostSynapticModel 
from minimized_parameters_cyclic import parameters, k_dissD, KD1, KD2

parameters['k_dissD1'] = k_dissD[NMBD]
parameters['k_dissD2'] = k_dissD[NMBD]
parameters['k_assocD1'] = k_dissD[NMBD]/KD1[NMBD]
parameters['k_assocD2'] = k_dissD[NMBD]/KD2[NMBD]
initACh = 7.75*(10**-6) # 7.75*(10**-6)
parameters['k_hydrol'] = 12000


model = PostSynapticModel(kinet=parameters, pnt=parameters)

d = 0*0.1e-6
model.set_steady_state( d ) 
sol    = model.solve_competition_problem(time=10, initACh=initACh)
twitch = model.twitch_height()
scale  = 1e-3

# グラフのプロット
plt.clf()
#plt.plot(sol.t, sol.y[0]/scale, label='Free ACh', ls='--', lw=3, c='b')
#plt.plot(sol.t, sol.y[1], label='ARO', ls='-', lw=3, c='r')
#plt.plot(sol.t, sol.y[2], label='ORA', ls='--', lw=3, c='r')
#plt.plot(sol.t, sol.y[4], label='ORO', ls='-', lw=3, c='c')
#plt.plot(sol.t, sol.y[5], label='DRO', ls='--', lw=3, c='c')
#plt.plot(sol.t, sol.y[6], label='DRD', ls=':', lw=3, c='c')
#plt.plot(sol.t, sol.y[7], label='ARD', ls='--', lw=3, c='m')
#plt.plot(sol.t, sol.y[8], label='DRA', ls='-', lw=3, c='m')
plt.plot(sol.t, (sol.y[10]+sol.y[11]+sol.y[12])/scale, label='R*(Cyclic)', ls='-', lw=3, c='k')
plt.plot(sol.t, sol.y[10]/scale, label='ARA*(Cyclic)', ls='--', lw=3, c='m')
plt.plot(sol.t, sol.y[11]/scale, label='ARO*(Cyclic)', ls='--', lw=3, c='c')
plt.plot(sol.t, sol.y[12]/scale, label='ORA*(Cyclic)', ls=':', lw=3, c='y')
plt.plot(sol.t, sol.y[3]/scale, label='ARA(Cyclic)', ls=':', lw=3, c='k')
plt.plot(sol.t, sol.y[13]/scale, label='RD(Cyclic)', ls='--', lw=3, c='b')

ax = plt.gca()
axins = ax.inset_axes([0.15, 0.4, 0.35, 0.45])
axins.plot(sol.t, (sol.y[10]+sol.y[11]+sol.y[12])/scale, label='R*(Cyclic)', ls='-', lw=3, c='k')
axins.plot(sol.t, sol.y[10]/scale, label='ARA*(Cyclic)', ls='--', lw=3, c='m')
axins.plot(sol.t, sol.y[11]/scale, label='ARO*(Cyclic)', ls='--', lw=3, c='c')
axins.plot(sol.t, sol.y[12]/scale, label='ORA*(Cyclic)', ls=':', lw=3, c='y')
axins.plot(sol.t, sol.y[3]/scale, label='ARA(Cyclic)', ls=':', lw=3, c='k')
x1, x2, y1, y2 = 0, 0.4, 0.4, 1.0
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
ax.indicate_inset_zoom(axins)


plt.legend(loc='upper right', fontsize=14)
plt.margins(0)
plt.xlabel('Time / ms', fontsize=15)
plt.ylabel(r'Fraction of AChR $\times 10^3$',fontsize=15)
#plt.xlim(0,0.01)
plt.ylim(0,3)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tick_params(direction='inout', width=1, length=8)
plt.savefig('fig_cyclic/time_course.pdf', bbox_inches="tight")
plt.show()

#plt.savefig(filename)


###########################################################################
# Reciprocal model 
from modules.model_reciprocal import PostSynapticModel
from minimized_parameters_reciprocal import parameters, k_dissD, KD1, KD2

parameters['k_dissD1'] = k_dissD[NMBD]
parameters['k_dissD2'] = k_dissD[NMBD]
parameters['k_assocD1'] = k_dissD[NMBD]/KD1[NMBD]
parameters['k_assocD2'] = k_dissD[NMBD]/KD2[NMBD]
initACh = 7.75*(10**-6) # 7.75*(10**-6)
parameters['k_hydrol'] = 12000

model = PostSynapticModel(kinet=parameters, pnt=parameters)

d = 0*0.1e-6
model.set_steady_state( d ) 
sol    = model.solve_competition_problem(time=100, initACh=initACh)
twitch = model.twitch_height()

scale  = 1e-2

# グラフのプロット
plt.clf()
#plt.plot(sol.t, sol.y[0], label='Free ACh', ls='--', lw=3, c='b')
#plt.plot(sol.t, sol.y[1], label='ARO', ls='-', lw=3, c='r')
#plt.plot(sol.t, sol.y[2], label='ORA', ls='--', lw=3, c='r')
#plt.plot(sol.t, sol.y[4], label='ORO', ls='-', lw=3, c='c')
#plt.plot(sol.t, sol.y[5], label='DRO', ls='--', lw=3, c='c')
#plt.plot(sol.t, sol.y[6], label='DRD', ls=':', lw=3, c='c')
#plt.plot(sol.t, sol.y[7], label='ARD', ls='--', lw=3, c='m')
#plt.plot(sol.t, sol.y[8], label='DRA', ls='-', lw=3, c='m')
#plt.plot(sol.t, sol.y[10], label='ARA*', ls='-', lw=3, c='m')
plt.plot(sol.t, sol.y[10]/scale, label='R*(Reciprocal)', ls='-', lw=3, c='k')
plt.plot(sol.t, sol.y[10]/scale, label='ARA*(Reciprocal)', ls='--', lw=3, c='m')
# plt.plot(sol.t, sol.y[11], label='AROopen', ls='-', lw=3, c='m')
# plt.plot(sol.t, sol.y[12], label='ORAopen', ls='-', lw=3, c='y')
plt.plot(sol.t, sol.y[3]/scale, label='ARA(Reciprocal)', ls=':', lw=3, c='k')
plt.plot(sol.t, sol.y[11]/scale, label='RD(Reciprocal)', ls='--', lw=3, c='b')

plt.legend(loc='upper left', fontsize=14)
plt.margins(0)
plt.xlabel('Time / ms', fontsize=15)
plt.ylabel(r'Fraction of AChR $\times 10^2$',fontsize=15)
#plt.xlim(0,0.01)
plt.ylim(0,4)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tick_params(direction='inout', width=1, length=8)
plt.savefig('fig_reciprocal/time_course.pdf', bbox_inches="tight")

#plt.savefig(filename)





