# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 22:35:40 2025

@author: carol
"""

from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt
    
print('Imported libraries: \nGekko \nNumpy \nMatplotlib')
#%%

# INPUT

def load_parameters()->dict:
    '''
    This fucntion created a dictionary containing the reactro parameters.
    This allows the change of the parameters values outside de model, and
    in the eventual case, it can be used for asking the user to input
    the values.
    
    Returns
    -------
    dict
        A dictionary containing the reactor parameter for the optimization
        

    '''
    
    reactor_parameters = {
        
        'C_init': 0.1367, # Initial concentration on the reactor [C(t=0) = C_init]
        'T_init': 0.7293, # Initial temperature on the reactor [T(t=0) = T_init]
        'u_init': 340, # Initial coolant flow on the heat exchanger [u(t=0) = u_init]
        'alpha': 1.95e-4, # Heat exchange constant (T - Tc) for energy balance
        'theta': 20, # Heat exchange constant (Tf - T) for energy balance
        'Tf': 0.3947, # Temperature of the entering flow
        'Tc': 0.3916, # Temperature of the coolant in the heat exchanger
        'k_10': 10, # Reaction rate constant
        'eta': 1, # Adimensional reaction rate constant
        'thao': 20, # Adimensional final time
        'C_des': 0.0194, # Concentration on steady state
        'T_des': 0.8, # Temperature on steady state
        'u_des': 340, # Heat exchanger coolant flow on steady state
        'alpha_1': 1e4, # Podering weight for concetration on the objective function
        'alpha_2': 2e4, # Pondeting weight for temperature on the objective fucntion
        'alpha_3':1e-3, # Pondering weight for the coolan flow on the objevtive function
        'adi_time': [0, 30], # Adimensional time window
        'steps': 50 # Number of steps
    
        }
    
    return reactor_parameters
    


#%%
# MODEL

def run_model(reactor:dict) -> dict:
    
    '''
    Creates the reacotr models and optimize the coolant flow through time
    in order to reach steady state in concetration, temperature and coolant
    flow in the least amount of time posible.
    

    Parameters
    ----------
    reactor : dict
        A dictionary describing the reactor

    Returns
    -------
    dict
        A dictionary saving variables solution

    '''
    
    # Create model
    
    model = GEKKO()
    
    # PARAMETERS
    
    C_init = model.Const(value=reactor['C_init'], name='Initial Concentration on the Reactor')
    T_init = model.Const(value=reactor['T_init'], name='Initial Temperature on the Reactor')
    u_init = model.Const(value=reactor['u_init'], name='Initial Coolant Flow on the Heat Exchanger')
    alpha = model.Const(value=reactor['alpha'], name='Heat Exchange Constant (T - Tc)')
    theta = model.Const(value=reactor['theta'], name='Heat Exchange Constant (Tf - T)')
    Tf = model.Const(value=reactor['Tf'], name='Temperature of the Entering Flow')
    Tc = model.Const(value=reactor['Tc'], name='Temperature of the Coolant in the Heat Exchanger')
    k_10 = model.Const(value=reactor['k_10'], name='Reaction Rate Constant')
    eta = model.Const(value=reactor['eta'], name='Adimensional Reaction Rate Constant')
    C_des = model.Const(value=reactor['C_des'], name='Concentration on Steady State')
    T_des = model.Const(value=reactor['T_des'], name='Temperature on Steady State')
    u_des = model.Const(value=reactor['u_des'], name='Coolant Flow on Stable State')
    alpha_1 = model.Const(value=reactor['alpha_1'], name='Podering Weight for Concetration')
    alpha_2 = model.Const(value=reactor['alpha_2'], name='Pondering Weight for Temperature')
    alpha_3 = model.Const(value=reactor['alpha_3'], name='Pondering Weight for Coolant Flow')


    model.time = np.linspace(reactor['adi_time'][0], reactor['adi_time'][1], reactor['steps'])
    
    
    # VARIABLES
    
    C = model.Var(C_init) # Concentration through time
    T = model.Var(T_init) # Temperature through time
    w = model.Var() # Reaction function through time
    u = model.Var(u_init, lb=0, ub=500) # Coolant flow through time
    r = model.Var() # Reaction rate through time
    
    # EQUATIONS
    
    model.Equations([C.dt() == ((1/theta)*(1-C)) - r, # Mass balance
                    T.dt() == ((1/theta)*(Tf-T)) + r - (alpha*u*(T-Tc)), # Energy balance
                    r == k_10*model.exp(w)*C, # Reaction rate
                    eta == -(w*T)]) # Reaction function
    
    
    # Model IMODE for Dynamic Trayectory
    model.options.IMODE = 6
    
    # OBJECTIVE FUNCTION
    model.Minimize((alpha_1*(C_des-C)**2) + (alpha_2*(T_des - T)**2) + (alpha_3*(u_des-u)**2))
    
    # RUN
    model.solve()
    
    # STORE SOLUTION
    
    opt_sol = {
        
        'Concetration': C.VALUE,
        'Temperature': T.VALUE,
        'Coolant flow': u.VALUE,
        'Reaction rate': r.VALUE,
        'Reaction function': w.VALUE,
        'Adimensional time': model.time
        
        }

    return opt_sol
#%%

def plot_reactor(reactor: dict, reactor_sol: dict)-> None:
    '''
    Plots the optimized trayectories of concetration, temperature
    and coolant flows for the reactors.

    Parameters
    ----------
    reactor : dict
        Dictionary containing the reactor parameters
    reactor_sol : dict
        Dictionary containing the optimize trayectory solution

    '''
    
    plt.figure(1, figsize=(8,6))
    plt.plot(reactor_sol['Adimensional time'], reactor_sol['Concetration'], label='Concentration')
    plt.grid(True)
    plt.title('Reactor Concetration', weight='bold')
    plt.xlabel('Adimensional Time')
    plt.ylabel('Adimensional Concentration')
    plt.axhline(reactor['C_des'], ls='--', label='Steady State Concetration')
    plt.legend()
    plt.show()

    plt.figure(2, figsize=(8,6))
    plt.plot(reactor_sol['Adimensional time'], reactor_sol['Temperature'], label='Temperature', color='g')
    plt.grid(True)
    plt.title('Reactor Temperature', weight='bold')
    plt.xlabel('Adimensional Time')
    plt.ylabel('Adimensional Temperature')
    plt.axhline(reactor['T_des'], ls='--', label='Steady State Temperature', color='g')
    plt.legend()
    plt.show()
    
    plt.figure(3, figsize=(8,6))
    plt.plot(reactor_sol['Adimensional time'], reactor_sol['Coolant flow'], label='Coolant flow', color='orange')
    plt.grid(True)
    plt.ylim(0, 500)
    plt.title('Heat Exchanger Reactor Flow', weight='bold')
    plt.xlabel('Adimensional Time')
    plt.ylabel('Adimensional Coolant FLow')
    plt.axhline(reactor['u_des'], ls='--', label='Steady State Coolant Flow', color='orange')
    plt.legend()
    plt.show()



#%%

# MAIN

def main():
    reactor_values = load_parameters()
    reactor_solution = run_model(reactor_values)
    plot_reactor(reactor_values, reactor_solution)
    print(reactor_solution['Coolant flow'])
    

#%%

# RUN

if __name__ == '__main__':
    main()
    
    
    
