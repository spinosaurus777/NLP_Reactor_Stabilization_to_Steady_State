# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 22:35:40 2025

@author: carol
"""

from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt
    
print('Imported libraries: \nGekko \nNumpy')
#%%

print('hello world')
# MODEL
model = GEKKO()
    
# PARAMETERS


C_init = model.Const(value=0.1367, name='Initial Concentration Value')
T_init = model.Const(value=0.7293, name='Initial Temperature Values')
alpha = model.Const(value=1.95e-4, name='Energy Balance - Constant')
theta = model.Const(value=20, name='Balances - Constant')
Tf = model.Const(value=0.3947, name='Input Temperature')
Tc = model.Const(value=0.3916, name='Heat Exchanger Temperature')
k_10 = model.Const(value=10, name='Reaction Rate Constant')
eta = model.Const(value=1, name='Constant')
thao = model.Const(value=10, name='Adimensional Time')
C_des = model.Const(value=0.0194, name='Concentration on Stable State')
T_des = model.Const(value=0.8, name='Temperature on Stable State')
u_des = model.Const(value=340, name='Heat Exchanger on Stable State')
alpha_1 = model.Const(value=1e4, name='Constant for Concentration FO')
alpha_2 = model.Const(value=2e4, name='Constant for Temperature FO')
alpha_3 = model.Const(value=1e-3, name='Constant for Refrigerant FO')


model.time = np.linspace(0, 20, 100)


# VARIABLES

C = model.Var(C_init)
T = model.Var(T_init)
w = model.Var()
u = model.Var(lb=0, ub=500)
t = model.Var(0)
r = model.Var()

# EQUATIONS

model.Equations([C.dt() == ((1/theta)*(1-C)) - r,
                T.dt() == ((1/theta)*(Tf-T)) + r - (alpha*u*(T-Tc)),
                r == k_10*model.exp(w)*C,
                eta == -(w*T)])

model.options.IMODE = 6
model.Minimize((alpha_1*(C_des-C)**2) + (alpha_2*(T_des - T)**2) + (alpha_3*(u_des-u)**2))

model.solve()

#%%

print(C.VALUE)

#%%

plt.figure()
plt.plot(model.time, T.VALUE)