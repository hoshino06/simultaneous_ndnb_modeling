# This file is automatically generated: do not edit this. 
import numpy as np 
def func(t, x):
    dxdt = np.zeros_like(x)
    dxdt[0] = -2317.25*x[0]*(x[1] + x[4]) - 2317.25*x[0]*(x[2] + x[5]) - 4634.5*x[0]*(-x[10] - x[11] - x[1] - x[2] - x[3] - x[4] - x[5] - x[6] - x[7] - x[8] + 1.0) - 12.0*x[0] + 18.0*x[1] + 18.0*x[2] + 36.0*x[3] + 18.0*x[7] + 18.0*x[8]
    dxdt[1] = -2317.25*x[0]*x[1] + 2317.25*x[0]*(-x[10] - x[11] - x[1] - x[2] - x[3] - x[4] - x[5] - x[6] - x[7] - x[8] + 1.0) - 42.5824175824176*x[1]*x[9] - 18.0*x[1] + 18.0*x[3] + 0.01*x[8]
    dxdt[2] = -2317.25*x[0]*x[2] + 2317.25*x[0]*(-x[10] - x[11] - x[1] - x[2] - x[3] - x[4] - x[5] - x[6] - x[7] - x[8] + 1.0) - 42.5824175824176*x[2]*x[9] - 18.0*x[2] + 18.0*x[3] + 0.01*x[7]
    dxdt[3] = 2317.25*x[0]*x[1] + 2317.25*x[0]*x[2] + 444.0*x[10] - 422.0*x[3]
    dxdt[4] = -2317.25*x[0]*x[4] - 42.5824175824176*x[4]*x[9] - 0.01*x[4] + 0.01*x[6] + 18.0*x[7] + 42.5824175824176*x[9]*(-x[10] - x[11] - x[1] - x[2] - x[3] - x[4] - x[5] - x[6] - x[7] - x[8] + 1.0)
    dxdt[5] = -2317.25*x[0]*x[5] - 42.5824175824176*x[5]*x[9] - 0.01*x[5] + 0.01*x[6] + 18.0*x[8] + 42.5824175824176*x[9]*(-x[10] - x[11] - x[1] - x[2] - x[3] - x[4] - x[5] - x[6] - x[7] - x[8] + 1.0)
    dxdt[6] = 42.5824175824176*x[4]*x[9] + 42.5824175824176*x[5]*x[9] - 0.02*x[6]
    dxdt[7] = 2317.25*x[0]*x[4] + 42.5824175824176*x[2]*x[9] - 18.01*x[7]
    dxdt[8] = 2317.25*x[0]*x[5] + 42.5824175824176*x[1]*x[9] - 18.01*x[8]
    dxdt[9] = 0
    dxdt[10] = -444.026*x[10] + 0.00013*x[11] + 386.0*x[3]
    dxdt[11] = 0.026*x[10] - 0.00013*x[11]
    return dxdt