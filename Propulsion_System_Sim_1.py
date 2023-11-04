import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
import pandas as pd
import diffeqpy
from diffeqpy import de

#-------------------------------------------------------
#------------------- Constants -------------------------
#-------------------------------------------------------

rho_3_L = 786 #Fuel density [kg/m^3]
rho_4_L = 800 #Ox density [kg/m^3]
R_N2 = 296 #Gas constant for nitrogen [J/kgK]
gamma_N2 = 1.4 #Specific Heat Ratio for nitrogen
P_atm = 0 #Atmospheric pressure (gauge) [Pa]
Tst = 273.15
Pst = 101325
rhost = Pst/(R_N2*Tst)


#-------------------------------------------------------
#-------------- Initial Conditions ---------------------
#-------------------------------------------------------

#Import the system network .csv file
System = pd.read_csv('Prop_System_Network_Setup_1.csv')

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

#Calculate the initial density of fluid (nitrogen) 
rho_0_0 = PropsSI('D','P',P_0_0,'T',T_0_0,'Nitrogen')
rho_1_0 = PropsSI('D','P',P_1_0,'T',T_1_0,'Nitrogen')
rho_2_0 = PropsSI('D','P',P_2_0,'T',T_2_0,'Nitrogen')

#Calculate the initial mass of fluid (nitrogen) using the volume of fluid and the density 
m_0_0 = V_0*rho_0_0
m_1_0 = V_1*rho_1_0
m_2_0 = V_2*rho_2_0


'''
#Calculate initial internal energy U using the total mass, and the pressure and temperature of the fluid (nitrogen)
U_0_0 = m_0_0*PropsSI('U','P',P_0_0,'T',T_0_0,'Nitrogen')
U_1_0 = m_1_0*PropsSI('U','P',P_1_0,'T',T_1_0,'Nitrogen')
U_2_0 = m_2_0*PropsSI('U','P',P_2_0,'T',T_2_0,'Nitrogen')
U_3_0 = m_3_0*PropsSI('U','P',P_3_0,'T',T_3_0,'Methanol') #Idk if this is the right thing
U_4_0 = m_4_0*PropsSI('U','P',P_4_0,'T',T_4_0,'N2O')
'''



#Find the number of nodes and elements in the system
Nodes = System['Node Number'].max() + 1
NoMassNodes = Nodes - 5
Elements = System['Element Number'].max() + 1
FlowElements = Elements - 2

def explicit_algebraic_equations(z):
    m_0, m_1, m_2, m_3, m_4 = z
    V_0 = System[0]['Param. 1 Value']*1e-3
    V_1 = System[0]['Param. 3 Value']*1e-3
    V_2 = System[1]['Param. 3 Value']*1e-3

    rho_0 = m_0/V_0
    P_0 = P_0_0*((rho_0/rho_0_0)**gamma_N2)
    T = P_0/(rho_0*R_N2)
    rho_1 = m_1/(V_1 - (m_3/rho_3_L))
    P_1 = rho_1*T*R_N2
    rho_2 = m_2/(V_2 - (m_4/rho_4_L))
    P_2 = rho_2*T*R_N2
    P_3 = P_1
    P_4 = P_2

    return P_0, P_1, P_2, P_3, P_4, T

#Function to find the nozzle mass flow rate at a given chamber pressure (approximately a linear relationship)
def Nozzle(P_chamber, D_throat, D_exit):
    gamma = 1.3 #Combustion gas SHR (ideally calculate based on OF)
    T_c = 2500 #Combustion gas temperature (ideally calculate based on OF)
    R = 296 #Combustion gas Gas Constant (ideally calculate based on OF)
    A_throat = np.pi*(D_throat**2)/4
    m_dot = (A_throat*P_chamber/np.sqrt(T_c))*np.sqrt(gamma/R)*((gamma+1)/2)**(-(gamma+1)/(2*gamma-2)) #Mass flow vs chamber pressure relationship
    P_exit = 0.8e5 #Replace with an equation for nozzle exit pressure vs chamber pressure (requires iterative solve of M_e currently)
    return m_dot

#Function to find the mass flow rate through an injector orifice for a Single Phase Incompressible (SPI) flow
def SPI_Injector(P_1,P_2,rho,A,C_d):
    m_dot = C_d*A*(2*rho*(P_1 - P_2))/np.sqrt(2*rho*(P_1 - P_2))
    return m_dot

#Function to find the dP across a pipe of a certain length, diameter and roughness at some mass flow rate
def Pipe(m_dot,rho,mu,L,D,k):
    U = m_dot/(rho*(D**2)*np.pi/4) #Find fluid velocity
    Re_D = rho*U*D/mu #Find Reynolds number of flow
    f = 0.25/(np.log10(k/(3.7*D) + 5.74/(Re_D**0.9))**2) #Find Darcy friction factor of flow
    dP = (m_dot**2)*8*f*L/(rho*(np.pi**2)*D**5) #Find pressure drop across pipe
    return dP

#Function to find mass flow rate across a component with a defined Kv in Nitrogen Gas flow
def Kv_Component_N2(P_1, P_2, T, K_v):
    P_1,P_2 = P_1*1e-5, P_2*1e-5 #Convert pressure to bar
    # determine flow regime through valve
    if (P_2 >= 0.528*P_1):
        # non-choked flow
        Q_N = K_v*514*(((P_1-P_2)*P_2)/(rhost*T))/np.sqrt(abs(((P_1-P_2)*P_2)/(rhost*T)))
    else:
        # choked flow
        Q_N = K_v*257*P_1/np.sqrt(abs(rhost*T))
    m_dot = Q_N*(Pst/(R_N2*Tst))*(1/3600) #Convert standard flow rate to mass flow rate in kg/s
    return m_dot

#Function to find mass flow rate across a component with a defined Kv in SPI Liquid flow
def Kv_Component_Liq(P_1, P_2, rho, K_v):
    P_1,P_2 = P_1*1e-5, P_2*1e-5 #Convert pressure to bar
    SG = rho_L/1000 #Calculate specific gravity of propellant (value for water is 1)
    Q = K_v*((P_1 - P_2)/SG)/np.sqrt(abs((P_1 - P_2)/SG)) #Calculate volumetric flow rate based on Kv
    m_dot = Q*rho_L*(1/3600) #Convert standard flow rate to mass flow rate ing kg/s
    return m_dot

# Solve a series of non-linear algebraic equations at a single timestep (basically algebraic constraints for each timestep of the differential equation solver)
def algebraic_equations(x, z, t, P_guess, T):
    F = np.zeros(len(x)) #Define an array of residuals (1 for each equation)
    P = x[:Nodes-1]
    #h = x[Nodes-8:]
    m_dot = x[Nodes:]

    m_0, m_1, m_2, m_3, m_4 = z

    # ---------------------------------------------------------------------------------------------------------------------------------
    # Put all the equations in here and make them equal to the residuals (so that if residuals are all zero, you have an exact solution)
    # ---------------------------------------------------------------------------------------------------------------------------------

    # Mass conservation at each massless node
    for i in range(0,len(P)): #iterate through every node
        if System[i]['Node Description'] == 'Nozzle Exit': #Apply boundary condition at the nozzle exit (don't test for continuity here)
            F[i] = P[i] - 0.8*1e5
        elif i < 5: #Apply boundary conditions for pressure at all the tank nodes (don't test of rcontinuity here)
            F[i] = P[i] - P_guess[i] 
        else: #Otherwise, it is an internal node where mass continuity equations are needed
            for j in range(2,len(m_dot)+2): #iterate through every element 
                if System.index[System['Upstream Node'] == i] == j:
                    F[i] += -m_dot[j]
                elif System.index[System['Downstream Node'] == i] == j:
                    F[i] += m_dot[j]
    
    # Momentum conservation across every element
    for i in range(2,len(m_dot)+2):
        P_1 = P[int(System[i]['Upstream Node'])] #Find upstream pressure from element
        P_2 = P[int(System[i]['Downstream Node'])] #Find downstream pressure from element

        if System[i]['Element Type'] == 'Regulator':
            #Regulator dP vs mdot equation


        if System[i]['Element Type'] == 'Kv Component':
            #Kv component dP vs mdot equation
            K_v = System[i]['Param. 3 Value']
            if System[System[i]['Upstream Node']]['Fluid'] == 'Nitrogen': #Check if the fluid is N2
                F[i + Nodes] = m_dot[i] - Kv_Component_N2(P_1, P_2, T, K_v)
            elif System[System[i]['Upstream Node']]['Fluid'] == 'Nitrous Oxide' #Check if the fluid is nitrous
                F[i + Nodes] = m_dot[i] - Kv_Component_Liq(P_1, P_2, rho_4_L, K_v)
            else: #Otherwise, fluid is methanol
                F[i + Nodes] = m_dot[i] - Kv_Component_Liq(P_1, P_2, rho_3_L, K_v)
        if System[i]['Element Type'] == 'Pipe':
            #Pipe dP vs mdot equation
            D, L, k = System[i]['Param. 3 Value']*1e-3, System[i]['Param. 4 Value'], System[i]['Param. 5 Value']*1e-3
            if System[System[i]['Upstream Node']]['Fluid'] == 'Nitrogen': #Check if the fluid is N2
                rho = PropsSI('D','P',P_1,'T',T,'Nitrogen')
                mu = PropsSI('V','P',P_1,'T',T,'Nitrogen')
            elif System[System[i]['Upstream Node']]['Fluid'] == 'Nitrous Oxide' #Check if the fluid is nitrous
                rho = rho_4_L
                mu = PropsSI('V','P',P_1,'T',270,'N2O')
            else: #Otherwise, fluid is methanol
                rho = rho_3_L
                mu = PropsSI('V','P',P_1,'T',298,'methanol')
            F[i + Nodes] = P_1 - P_2 - Pipe(m_dot[i],rho,mu,L,D,k)
        if System[i]['Element Type'] == 'Valve':
            #Valve dP vs mdot equation


        if System[i]['Element Type'] == 'Injector':
            #Injector dP vs mdot equation
            A = System[i]['Param. 3 Value']
            C_d = System[i]['Param. 4 Value']
            if System[System[i]['Upstream Node']]['Fluid'] == 'Nitrous Oxide' #Check if the fluid is nitrous
                F[i + Nodes] = m_dot[i] - SPI_Injector(P_1,P_2,rho_4_L,A,C_d)
            else: #Otherwise, fluid is methanol
                F[i + Nodes] = m_dot[i] - SPI_Injector(P_1,P_2,rho_3_L,A,C_d)
        if System[i]['Element Type'] == 'Nozzle':
            #Nozzle mdot and P_exit vs Pc equations
            F[i + Nodes] = m_dot[i] - Nozzle(P_chamber=P_1, D_throat=System[i]['Param. 3 Value']*1e-3, D_exit=System[i]['Param. 4 Value']*1e-3)
        else:
            pass
            #Some other basic equation to prevent errors



    return F

# Define the DAE system (Differential Algebraic Equation System)
def dae_system(t, z):
    
    m_0, m_1, m_2, m_3, m_4 = z

    #Algebraic equations to determine other state variables
    
    P_0, P_1, P_2, P_3, P_4, T = explicit_algebraic_equations(z)
    
    P_guess = System[:Nodes-1]['IC 1 Parameter']*1e5
    P_guess[0] = P_0 #Assign Boundary Conditions
    P_guess[1] = P_1
    P_guess[2] = P_2
    P_guess[3] = P_3
    P_guess[4] = P_4
    m_dot_guess = System[2:Elements-1]['IC 3 Parameter']
    x0 = np.concatenate(P_guess,m_dot_guess)

    solution = root(algebraic_equations, x0, args=(z, t, P_guess, T), method='hybr', tol=1e-5)
    x = solution.x

    # Get m_dot and P from the solution
    
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

    '''
    U_dot_0 = m_dot_0*h[System.index[System['Upstream Node'] == 0]]
    U_dot_1 = m_dot_1*h[System.index[System['Downstream Node'] == 1]]
    U_dot_2 = m_dot_2*h[System.index[System['Downstream Node'] == 2]]
    U_dot_3 = m_dot_0*h[System.index[System['Upstream Node'] == 3]]
    U_dot_4 = m_dot_0*h[System.index[System['Upstream Node'] == 4]]
    '''

    return [m_dot_0, m_dot_1, m_dot_2, m_dot_3, m_dot_4]


#Create vector of initial conditions 
z0 = [m_0_0, m_1_0, m_2_0, m_3_0, m_4_0]


# Define the time span for integration
t_span = (0, 10)

# Generate evaluation times
t_eval = np.linspace(*t_span, 100)

# Solve the DAE
sol = solve_ivp(dae_system, t_span, z0, method='LSODA', t_eval=t_eval)

# Access the solution
t = sol.t
m_0_0, m_1_0, m_2_0, m_3_0, m_4_0 = sol.y

#Plot some stuff