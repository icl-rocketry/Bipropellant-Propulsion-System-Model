import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
import pandas as pd
import diffeqpy
from diffeqpy import de

System = pd.read_csv('Prop_System_Network_Setup_1.csv')

#print(System)


#Find the number of nodes and elements in the system
Nodes = System['Node Number'].max() + 1
Elements = System['Element Number'].max() + 1

# Solve a series of non-linear algebraic equations at a single timestep (basically algebraic constraints for each timestep of the differential equation solver)
def algebraic_equations(x, z, t):
    F = np.zeros(len(x)) #Define an array of residuals (1 for each equation)
    P = x[:Nodes-7]
    h = x[Nodes-8:]
    m_dot = x[]

    # ---------------------------------------------------------------------------------------------------------------------------------
    # Put all the equations in here and make them equal to the residuals (so that if residuals are all zero, you have an exact solution)
    # ---------------------------------------------------------------------------------------------------------------------------------

    return F

# Define the DAE system (Differential Algebraic Equation System)
def dae_system(t, z):
    
    m_0, m_1, m_2, m_3, m_4, U_0, U_1, U_2, U_3, U_4 = z

    #Algebraic equations to determine other state variables
    
    P_0 = System[6:Nodes-1]['IC 1 Parameter']
    h_0 = System[6:Nodes-1]['IC 2 Parameter']
    m_dot_0 = System[2:Elements-1]['IC 3 Parameter']
    x0 = np.concatenate(P_0,h_0,m_dot_0)

    solution = root(algebraic_equations, x0, args=(z, t), method='hybr', tol=1e-5)
    x = solution.x
    
    # Differential Equations
    m_dot_0 = -m_dot[System.index[System['Upstream Node'] == 0]]
    m_dot_1 = m_dot[System.index[System['Downstream Node'] == 1]]
    m_dot_2 = m_dot[System.index[System['Downstream Node'] == 2]]
    if m_3 <= 0: #Check for fuel depletion
        m_3 = 0
        m_dot_3 = 0
    else:
        m_dot_3 = -m_dot[System.index[System['Upstream Node'] == 3]]
    if m_4 <= 0: #Check for oxidiser depletion
        m_4 = 0
        m_dot_4 = 0
    else:
        m_dot_4 = -m_dot[System.index[System['Upstream Node'] == 4]]

    U_dot_0 = m_dot_0*h[System.index[System['Upstream Node'] == 0]]
    U_dot_1 = m_dot_1*h[System.index[System['Downstream Node'] == 1]]
    U_dot_2 = m_dot_2*h[System.index[System['Downstream Node'] == 2]]
    U_dot_3 = m_dot_0*h[System.index[System['Upstream Node'] == 3]]
    U_dot_4 = m_dot_0*h[System.index[System['Upstream Node'] == 4]]

    return [m_dot_0, m_dot_1, m_dot_2, m_dot_3, m_dot_4, U_dot_0, U_dot_1, U_dot_2, U_dot_3, U_dot_4]

'''
# Define initial conditions of ODEs (mass, m and internal energy U)
z0 = np.zeros(len(ODEs)*2)
for i in range(0,len(ODEs)):
    IC1_Given = System.loc[ODEs[i]]['IC 1 Parameter']
    IC2_Given = System.loc[ODEs[i]]['IC 2 Parameter']
    T = System.loc[ODEs[i]]['IC 2 Value']
    if IC1_Given == 'P':
        P = System.loc[ODEs[i]]['IC 1 Value']*1e5
        if System.loc[ODEs[i]]['Param. 1 Name'] == 'Volume':
            V = System.loc[ODEs[i]]['Param. 1 Value']*1e-3
        else:
            V = System.loc[System.index[System['Upstream Node'] == ODEs[i]]]['Param. 3 Value']*1e-3
        m = V*PropsSI('D','P',P,'T',T,'Nitrogen') #Calculate the initial mass of fluid (nitrogen) using the volume of fluid and the density 
        U = m*PropsSI('U','P',P,'T',T,'Nitrogen') #Calculate initial internal energy U using the total mass, and the pressure and temperature of the fluid (nitrogen)
    else:
        m = System.loc[ODEs[i]]['IC 1 Value']



    z0[i] = m
    z0[i + len(ODEs)] = U
'''

# Define initial conditions of ODEs (mass, m and internal energy U of fluid in tanks)

#Get initial fluid pressures (note that liquid pressures are equal to the ullage pressure above them)
P_0_0 = System.loc[0]['IC 1 Value']*1e5
P_1_0 = System.loc[1]['IC 1 Value']*1e5
P_2_0 = System.loc[2]['IC 1 Value']*1e5
P_3_0 = P_1_0
P_4_0 = P_2_0

#Get initial fluid masses (for liquid masses)
m_3_0 = System.loc[3]['IC 1 Value']
m_4_0 = System.loc[4]['IC 1 Value']

#Get initial fluid temperatures [K]
T_0_0 = System.loc[0]['IC 2 Value']
T_1_0 = System.loc[1]['IC 2 Value']
T_2_0 = System.loc[2]['IC 2 Value']
T_3_0 = System.loc[3]['IC 2 Value']
T_4_0 = System.loc[4]['IC 2 Value']

#Get tank volumes and convert to SI units [m^3]
V_0 = System.loc[0]['Param. 1 Value']*1e-3
V_1 = System.loc[0]['Param. 3 Value']*1e-3
V_2 = System.loc[1]['Param. 3 Value']*1e-3

#Calculate the initial mass of fluid (nitrogen) using the volume of fluid and the density 
m_0_0 = V_0*PropsSI('D','P',P_0_0,'T',T_0_0,'Nitrogen')
m_1_0 = V_1*PropsSI('D','P',P_1_0,'T',T_1_0,'Nitrogen')
m_2_0 = V_2*PropsSI('D','P',P_2_0,'T',T_2_0,'Nitrogen')

#Calculate initial internal energy U using the total mass, and the pressure and temperature of the fluid (nitrogen)
U_0_0 = m_0_0*PropsSI('U','P',P_0_0,'T',T_0_0,'Nitrogen')
U_1_0 = m_0_0*PropsSI('U','P',P_1_0,'T',T_1_0,'Nitrogen')
U_2_0 = m_0_0*PropsSI('U','P',P_2_0,'T',T_2_0,'Nitrogen')
U_3_0 = m_0_0*PropsSI('U','P',P_3_0,'T',T_3_0,'Nitrogen')
U_4_0 = m_0_0*PropsSI('U','P',P_4_0,'T',T_4_0,'Nitrogen')

#Create vector of initial conditions 
z0 = [m_0_0, m_1_0, m_2_0, m_3_0, m_4_0, U_0_0, U_1_0, U_2_0, U_3_0, U_4_0]


# Define the time span for integration
t_span = (0, 10)

# Generate evaluation times
t_eval = np.linspace(*t_span, 100)

# Solve the DAE
sol = solve_ivp(dae_system, t_span, z0, method='LSODA', t_eval=t_eval)

# Access the solution
t = sol.t
m_0_0, m_1_0, m_2_0, m_3_0, m_4_0, U_0_0, U_1_0, U_2_0, U_3_0, U_4_0 = sol.y

#Plot some stuff