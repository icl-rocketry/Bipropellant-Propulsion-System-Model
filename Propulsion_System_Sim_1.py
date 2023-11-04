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

rho_L = 1000 #Liquid propellant density [kg/m^3]
R_N2 = 296 #Gas constant for nitrogen [J/kgK]
gamma_N2 = 1.4 #Specific Heat Ratio for nitrogen
P_atm = 0 #Atmospheric pressure (gauge) [Pa]
Tst = 273.15
Pst = 101325
rhost = Pst/(R_N2*Tst)


System = pd.read_csv('Prop_System_Network_Setup_1.csv')

#Find the number of nodes and elements in the system
Nodes = System['Node Number'].max() + 1
NoMassNodes = Nodes - 5
Elements = System['Element Number'].max() + 1
FlowElements = Elements - 2

def explicit_algebraic_equations(m_2, m_3, rho_1):
    P_1 = P_1_0*((rho_1/rho_1_0)**gamma_N2)
    T = P_1/(rho_1*R_N2)
    rho_2 = m_2/(V_2 - (m_3/rho_L))
    P_2 = rho_2*T*R_N2
    m_1 = V_1*rho_1
    return P_1, T, rho_2, P_2, m_1

#Function to find the nozzle mass flow rate at a given chamber pressure (approximately a linear relationship)
def Nozzle(P_chamber, D_throat, D_exit):
    gamma = 1.3 #Combustion gas SHR (ideally calculate based on OF)
    T_c = 2500 #Combustion gas temperature (ideally calculate based on OF)
    R = 296 #Combustion gas Gas Constant (ideally calculate based on OF)
    A_throat = np.pi*(D_throat**2)/4
    m_dot = (A_throat*P_chamber/np.sqrt(T_c))*np.sqrt(gamma/R)*((gamma+1)/2)**(-(gamma+1)/(2*gamma-2)) #Mass flow vs chamber pressure relationship
    P_exit = 0.8e5 #Replace with an equation for nozzle exit pressure vs chamber pressure (requires iterative solve of M_e currently)
    return m_dot, P_exit

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
def Kv_Component_N2(P_1, P_2, rho_1, rho_2, T, K_v):
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
def algebraic_equations(x, z, t):
    F = np.zeros(len(x)) #Define an array of residuals (1 for each equation)
    P = x[:NoMassNodes-1]
    #h = x[Nodes-8:]
    m_dot = x[NoMassNodes:]

    m_0, m_1, m_2, m_3, m_4 = z

    # ---------------------------------------------------------------------------------------------------------------------------------
    # Put all the equations in here and make them equal to the residuals (so that if residuals are all zero, you have an exact solution)
    # ---------------------------------------------------------------------------------------------------------------------------------

    # Mass conservation at each massless node
    for i in range(5,len(P)+5): #iterate through every node
        if SystemSystem[i]['Node Description'] == 'Nozzle Exit': #If 
            pass
        else:
            for j in range(2,len(m_dot)+2): #iterate through every element
                if System.index[System['Upstream Node'] == i] == j:
                    F[i] += -m_dot[j]
                elif System.index[System['Downstream Node'] == i] == j:
                    F[i] += m_dot[j]
    
    # Momentum conservation across every element
    for i in range(2,len(m_dot)+2):
        if SystemSystem[i]['Element Type'] == 'Regulator':
            #Regulator dP vs mdot equation

        if SystemSystem[i]['Element Type'] == 'Kv Component':
            #Kv component dP vs mdot equation
        
        if SystemSystem[i]['Element Type'] == 'Pipe':
            #Pipe dP vs mdot equation

        if SystemSystem[i]['Element Type'] == 'Valve':
            #Valve dP vs mdot equation

        if SystemSystem[i]['Element Type'] == 'Injector':
            #Injector dP vs mdot equation
        
        if SystemSystem[i]['Element Type'] == 'Nozzle':
            #Nozzle mdot and P_exit vs Pc equations

        else:
            #Some other basic equation to prevent errors



    return F

# Define the DAE system (Differential Algebraic Equation System)
def dae_system(t, z):
    
    m_0, m_1, m_2, m_3, m_4 = z

    #Algebraic equations to determine other state variables
    
    

    P_0 = System[5:Nodes-1]['IC 1 Parameter']
    #h_0 = System[5:Nodes-1]['IC 2 Parameter']
    m_dot_0 = System[2:Elements-1]['IC 3 Parameter']
    x0 = np.concatenate(P_0,m_dot_0)

    solution = root(algebraic_equations, x0, args=(z, t), method='hybr', tol=1e-5)
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

'''
#Calculate initial internal energy U using the total mass, and the pressure and temperature of the fluid (nitrogen)
U_0_0 = m_0_0*PropsSI('U','P',P_0_0,'T',T_0_0,'Nitrogen')
U_1_0 = m_1_0*PropsSI('U','P',P_1_0,'T',T_1_0,'Nitrogen')
U_2_0 = m_2_0*PropsSI('U','P',P_2_0,'T',T_2_0,'Nitrogen')
U_3_0 = m_3_0*PropsSI('U','P',P_3_0,'T',T_3_0,'Methanol') #Idk if this is the right thing
U_4_0 = m_4_0*PropsSI('U','P',P_4_0,'T',T_4_0,'N2O')
'''

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