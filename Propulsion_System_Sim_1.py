import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
import pandas as pd
#import diffeqpy
#from diffeqpy import de

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
#System = pd.read_csv('Prop_System_Network_Setup_1.csv')
System = pd.read_csv('Simplified_Prop_System_1.csv')

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
V_1 = System.loc[1]['Param. 3 Value']*1e-3
V_2 = System.loc[0]['Param. 3 Value']*1e-3

#Calculate the initial density of fluid (nitrogen) 
rho_0_0 = PropsSI('D','P',P_0_0,'T',T_0_0,'Nitrogen')
rho_1_0 = PropsSI('D','P',P_1_0,'T',T_1_0,'Nitrogen')
rho_2_0 = PropsSI('D','P',P_2_0,'T',T_2_0,'Nitrogen')

#Calculate initial entropy of the high pressure nitrogen
S_0_0 = PropsSI('S','P',P_0_0,'T',T_0_0,'Nitrogen')

#Calculate the initial mass of fluid (nitrogen) using the volume of fluid and the density 
m_0_0 = V_0*rho_0_0
m_1_0 = rho_1_0*(V_1 - m_3_0/rho_3_L)
m_2_0 = rho_2_0*(V_2 - m_4_0/rho_4_L)


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

#Find the number of nodes and elements in the system
Nodes = int(System['Node Number'].max()) + 1
NoMassNodes = Nodes - 5
Elements = int(System['Element Number'].max()) + 1
FlowElements = Elements - 2



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
    if abs(P_1 - P_2) < 1e-5:
        m_dot = 0
    else:
        m_dot = C_d*A*(2*rho*(P_1 - P_2))/np.sqrt(abs(2*rho*(P_1 - P_2)))
    return m_dot

#Function to find the dP across a pipe of a certain length, diameter and roughness at some mass flow rate
def Pipe(m_dot,rho,mu,L,D,k):
    if abs(m_dot) < 1e-10:
        dP = 0
    else:
        U = abs(m_dot)/(rho*(D**2)*np.pi/4) #Find fluid velocity
        Re_D = rho*U*D/mu #Find Reynolds number of flow
        f = 0.25/(np.log10(k/(3.7*D) + 5.74/(Re_D**0.9))**2) #Find Darcy friction factor of flow
        dP = (m_dot*abs(m_dot))*8*f*L/(rho*(np.pi**2)*D**5) #Find pressure drop across pipe
    return dP

#Function to find mass flow rate across a component with a defined Kv in Nitrogen Gas flow
def Kv_Component_N2(P_1, P_2, T, K_v):
    P_1,P_2 = P_1*1e-5, P_2*1e-5 #Convert pressure to bar
    # determine flow regime through valve
    #if (P_1 - P_2) <= 1e-5:
    #    m_dot = 0
    if (P_2 >= 0.528*P_1):
        # non-choked flow
        Q_N = K_v*514*(((P_1-P_2)*P_2)/(rhost*T))/np.sqrt(abs(((P_1-P_2)*P_2)/(rhost*T)))
        m_dot = Q_N*(Pst/(R_N2*Tst))*(1/3600) #Convert standard flow rate to mass flow rate in kg/s
    else:
        # choked flow
        Q_N = K_v*257*P_1/np.sqrt(abs(rhost*T))
        m_dot = Q_N*(Pst/(R_N2*Tst))*(1/3600) #Convert standard flow rate to mass flow rate in kg/s
    return m_dot

#Function to find mass flow rate across a component with a defined Kv in SPI Liquid flow
def Kv_Component_Liq(P_1, P_2, rho, K_v):
    P_1,P_2 = P_1*1e-5, P_2*1e-5 #Convert pressure to bar
    SG = rho/1000 #Calculate specific gravity of propellant (value for water is 1)
    if P_1 - P_2 < 1e-5:
        m_dot = 0
    else:
        Q = K_v*((P_1 - P_2)/SG)/np.sqrt(abs((P_1 - P_2)/SG)) #Calculate volumetric flow rate based on Kv
        m_dot = Q*rho*(1/3600) #Convert standard flow rate to mass flow rate ing kg/s
    return m_dot

def Regulator(P_1, P_2, T, K_v_max, t):
    P_1,P_2 = P_1*1e-5, P_2*1e-5
    '''
    #Equations to approximate a mechanical regulator with a certain maximum Kv and spring
    Kp = 0.5
    K_v_controlled = 0.05 + (50 - P_2)*Kp
    
    if K_v_controlled < 0:
        K_v_controlled = 0
    elif K_v_controlled > K_v_max:
        K_v_controlled = K_v_max
    
    K_v = K_v_controlled
    '''
    K_v = K_v_max
    #K_v = np.interp(t, [0, 0.6, 0.8, 10], [0, 0, K_v_controlled, K_v_controlled]) # Define the Kv of the valve based on a predefined path

    # determine flow regime through valve
    if (P_2 >= 0.528*P_1):
        # non-choked flow
        Q_N = K_v*514*(((P_1-P_2)*P_2)/(rhost*T))/np.sqrt(abs(((P_1-P_2)*P_2)/(rhost*T)))
    else:
        # choked flow
        Q_N = K_v*257*P_1/np.sqrt(abs(rhost*T))
    m_dot = Q_N*(Pst/(R_N2*Tst))*(1/3600) #Convert standard flow rate to mass flow rate
    return m_dot

def Ox_Valve(P_1, P_2, K_v_max, t):
    #K_v = np.interp(t, [0, 0.7, 0.9, 2, 2.1, 4, 4.1, 10], [0, 0, 10, 10, 1, 1, 0.3, 0.3]) # Define the Kv of the valve based on a predefined path
    #K_v = np.interp(t, [0, 0.7, 0.9, 2, 2.1, 4, 4.1, 10], [0, 0, 10, 10, 1, 1, 4, 4]) # Define the Kv of the valve based on a predefined path

    K_v = K_v_max*np.interp(t, [0, 0.1, 0.2, 2, 2.1, 4, 4.1, 10], [1, 1, 1, 1, 1, 1, 1, 1]) # Define the Kv of the valve based on a predefined path

    P_1,P_2 = P_1*1e-5, P_2*1e-5
    SG = rho_4_L/1000 #Calculate specific gravity of propellant (value for water is 1)
    Q = K_v*((P_1 - P_2)/SG)/np.sqrt(abs((P_1 - P_2)/SG)) #Calculate volumetric flow rate based on Kv
    m_dot = Q*rho_4_L*(1/3600) #Convert standard flow rate to mass flow rate
    return m_dot

def Fuel_Valve(P_1, P_2, K_v_max, t):
    #K_v = np.interp(t, [0, 0.7, 0.9, 2, 2.1, 4, 4.1, 10], [0, 0, 10, 10, 1, 1, 0.3, 0.3]) # Define the Kv of the valve based on a predefined path
    #K_v = np.interp(t, [0, 0.7, 0.9, 2, 2.1, 4, 4.1, 10], [0, 0, 10, 10, 1, 1, 4, 4]) # Define the Kv of the valve based on a predefined path

    K_v = K_v_max*np.interp(t, [0, 0.1, 0.2, 2, 2.1, 4, 4.1, 10], [1, 1, 1, 1, 1, 1, 1, 1]) # Define the Kv of the valve based on a predefined path

    P_1,P_2 = P_1*1e-5, P_2*1e-5
    SG = rho_3_L/1000 #Calculate specific gravity of propellant (value for water is 1)
    Q = K_v*((P_1 - P_2)/SG)/np.sqrt(abs((P_1 - P_2)/SG)) #Calculate volumetric flow rate based on Kv
    m_dot = Q*rho_3_L*(1/3600) #Convert standard flow rate to mass flow rate
    return m_dot

def explicit_algebraic_equations(z):
    #print('mark 4')
    print(z)
    m_0, m_1, m_2, m_3, m_4 = z
    V_0 = System.loc[0]['Param. 1 Value']*1e-3
    V_1 = System.loc[1]['Param. 3 Value']*1e-3
    V_2 = System.loc[0]['Param. 3 Value']*1e-3

    rho_0 = m_0/V_0
    P_0 = PropsSI('P','S',S_0_0,'D',rho_0,'Nitrogen')
    T = PropsSI('T','S',S_0_0,'D',rho_0,'Nitrogen')
    rho_1 = m_1/(V_1 - (m_3/rho_3_L))
    P_1 = PropsSI('P','T',T,'D',rho_1,'Nitrogen')
    rho_2 = m_2/(V_2 - (m_4/rho_4_L))
    P_2 = PropsSI('P','T',T,'D',rho_2,'Nitrogen')
    P_3 = P_1
    P_4 = P_2

    #print('mark 4.5')
    return P_0, P_1, P_2, P_3, P_4, T

# Solve a series of non-linear algebraic equations at a single timestep (basically algebraic constraints for each timestep of the differential equation solver)
def implicit_algebraic_equations(x, z, t, P_guess, T):
    #print('mark 6')
    F = np.zeros(len(x)) #Define an array of residuals (1 for each equation)
    P = x[:Nodes]
    #h = x[Nodes-8:]
    m_dot = x[Nodes:]

    m_0, m_1, m_2, m_3, m_4 = z

    # ---------------------------------------------------------------------------------------------------------------------------------
    # Put all the equations in here and make them equal to the residuals (so that if residuals are all zero, you have an exact solution)
    # ---------------------------------------------------------------------------------------------------------------------------------

    # Mass conservation at each massless node
    for i in range(0,len(P)): #iterate through every node
        if System.loc[i]['Node Description'] == 'Nozzle Exit': #Apply boundary condition at the nozzle exit (don't test for continuity here)
            F[i] = (P[i] - 0.8*1e5) * 1e-7
        elif i < 5: #Apply boundary conditions for pressure at all the tank nodes (don't test for continuity here)
            F[i] = (P[i] - P_guess[i]) * 1e-5
        else: #Otherwise, it is an internal node where mass continuity equations are needed
            for j in range(2,len(m_dot)+2): #iterate through every element 
                if System.loc[j]['Upstream Node'] == i:
                    F[i] += -m_dot[j-2] #If the node is upstream of the element, then the mass flow is flowing out of the node through that element (defined as negative)
                elif System.loc[j]['Downstream Node'] == i:
                    F[i] += m_dot[j-2] #If the node is downstream of the element, then the mass flow is flowing into the node through that element (defined as positive)
        if System.loc[i]['Fluid'] == 'Nitrogen': #Check if the fluid is N2
            F[i] *= 10 #Multiply nitrogen mass flow equations 
    
    # Momentum conservation across every element
    for i in range(2,len(m_dot)+2):
        P_1 = P[int(System.loc[i]['Upstream Node'])] #Find upstream pressure from element
        P_2 = P[int(System.loc[i]['Downstream Node'])] #Find downstream pressure from element

        if System.loc[i]['Element Type'] == 'Regulator':
            #Regulator dP vs mdot equation
            K_v_max = System.loc[i]['Param. 4 Value']
            F[i - 2 + Nodes] = m_dot[i-2] - Regulator(P_1, P_2, T, K_v_max, t)

        elif System.loc[i]['Element Type'] == 'Kv Component':
            #Kv component dP vs mdot equation
            K_v = System.loc[i]['Param. 3 Value']
            if System.loc[int(System.loc[i]['Upstream Node'])]['Fluid'] == 'Nitrogen': #Check if the fluid is N2
                F[i - 2 + Nodes] = m_dot[i-2] - Kv_Component_N2(P_1, P_2, T, K_v)
            elif System.loc[int(System.loc[i]['Upstream Node'])]['Fluid'] == 'Nitrous Oxide': #Check if the fluid is nitrous
                F[i - 2 + Nodes] = m_dot[i-2] - Kv_Component_Liq(P_1, P_2, rho_4_L, K_v)
            else: #Otherwise, fluid is methanol
                F[i - 2 + Nodes] = m_dot[i-2] - Kv_Component_Liq(P_1, P_2, rho_3_L, K_v)

        elif System.loc[i]['Element Type'] == 'Pipe':
            #Pipe dP vs mdot equation
            D, L, k = System.loc[i]['Param. 3 Value']*1e-3, System.loc[i]['Param. 4 Value'], System.loc[i]['Param. 5 Value']*1e-3
            if System.loc[int(System.loc[i]['Upstream Node'])]['Fluid'] == 'Nitrogen': #Check if the fluid is N2
                rho = PropsSI('D','P',P_1,'T',T,'Nitrogen')
                mu = PropsSI('V','P',P_1,'T',T,'Nitrogen')
            elif System.loc[int(System.loc[i]['Upstream Node'])]['Fluid'] == 'Nitrous Oxide': #Check if the fluid is nitrous
                rho = rho_4_L
                mu = PropsSI('V','P',P_1,'T',298,'Water')
            else: #Otherwise, fluid is methanol
                rho = rho_3_L
                mu = PropsSI('V','P',P_1,'T',298,'methanol')
            F[i - 2 + Nodes] = P_1 - P_2 - Pipe(m_dot[i-2],rho,mu,L,D,k)

        elif System.loc[i]['Element Type'] == 'Valve':
            #Valve dP vs mdot equation
            K_v_max = System.loc[i]['Param. 4 Value']
            if System.loc[int(System.loc[i]['Upstream Node'])]['Fluid'] == 'Nitrous Oxide': #Check if the fluid is nitrous
                F[i - 2 + Nodes] = m_dot[i-2] - Ox_Valve(P_1, P_2, K_v_max, t)
            else: #Otherwise, fluid is methanol
                F[i - 2 + Nodes] = m_dot[i-2] - Fuel_Valve(P_1, P_2, K_v_max, t)

        elif System.loc[i]['Element Type'] == 'Injector':
            #Injector dP vs mdot equation
            A = System.loc[i]['Param. 3 Value']
            C_d = System.loc[i]['Param. 4 Value']
            if System.loc[int(System.loc[i]['Upstream Node'])]['Fluid'] == 'Nitrous Oxide': #Check if the fluid is nitrous
                F[i - 2 + Nodes] = m_dot[i-2] - SPI_Injector(P_1,P_2,rho_4_L,A,C_d)
            else: #Otherwise, fluid is methanol
                F[i - 2 + Nodes] = m_dot[i-2] - SPI_Injector(P_1,P_2,rho_3_L,A,C_d)

        elif System.loc[i]['Element Type'] == 'Nozzle':
            #Nozzle mdot and P_exit vs Pc equations
            F[i - 2 + Nodes] = m_dot[i-2] - Nozzle(P_chamber=P_1, D_throat=System.loc[i]['Param. 3 Value']*1e-3, D_exit=System.loc[i]['Param. 4 Value']*1e-3)

        else:
            print('error in equation definition')
            print(System.loc[i]['Element Type'])
            #Some other basic equation to prevent errors

        #print(m_dot[6],m_dot[8],m_dot[9],m_dot[11])


    #print(F)

    return F

def algebraic_equations(t, z):
    #print('mark 2')
    P_0, P_1, P_2, P_3, P_4, T = explicit_algebraic_equations(z)
    #print(explicit_algebraic_equations(z))
    
    P_guess = np.array(System.loc[:Nodes]['IC 1 Value'])*1e5
    #print(P_0)
    P_guess[0] = P_0 #Assign Boundary Conditions
    P_guess[1] = P_1
    P_guess[2] = P_2
    P_guess[3] = P_3
    P_guess[4] = P_4
    #print(P_guess)
    #print(T)
    m_dot_guess = np.array(System.loc[2:]['IC 3 Value'])
    #print(m_dot_guess)
    x0 = np.concatenate((P_guess,m_dot_guess))
    #print(x0)
    #print('mark 5')
    solution = root(implicit_algebraic_equations, x0, args=(z, t, P_guess, T), method='lm', tol=1e-8)
    x = solution.x
    print('Residuals: ', implicit_algebraic_equations(x, z, t, P_guess, T))

    #print('mark 7')
    P = x[:Nodes]
    m_dot = x[Nodes:]

    return P, m_dot, T

# Define the DAE system (Differential Algebraic Equation System)
def dae_system(t, z):
    #print('mark 1')
    print('t = ',t)
    m_0, m_1, m_2, m_3, m_4 = z

    #Algebraic equations to determine other state variables
    P, m_dot, T = algebraic_equations(t, z)
    

    # Differential Equations
    m_dot_0 = -m_dot[System.index[System['Upstream Node'] == 0]-2]
    m_dot_1 = m_dot[System.index[System['Downstream Node'] == 1]-2]
    m_dot_2 = m_dot[System.index[System['Downstream Node'] == 2]-2]
    if m_3 <= 0: #Check for fuel depletion
        m_3 = 0
        m_dot_3 = 0
    else:
        m_dot_3 = -m_dot[System.index[System['Upstream Node'] == 3]-2]
    if m_4 <= 0: #Check for oxidiser depletion
        m_4 = 0
        m_dot_4 = 0
    else:
        m_dot_4 = -m_dot[System.index[System['Upstream Node'] == 4]-2]

    '''
    U_dot_0 = m_dot_0*h[System.index[System['Upstream Node'] == 0]]
    U_dot_1 = m_dot_1*h[System.index[System['Downstream Node'] == 1]]
    U_dot_2 = m_dot_2*h[System.index[System['Downstream Node'] == 2]]
    U_dot_3 = m_dot_0*h[System.index[System['Upstream Node'] == 3]]
    U_dot_4 = m_dot_0*h[System.index[System['Upstream Node'] == 4]]
    '''

    return [m_dot_0, m_dot_1, m_dot_2, m_dot_3, m_dot_4]



# Define the time span for integration
t_span = (0, 5)

# Generate evaluation times
t_eval = np.linspace(*t_span, 100)


'''
t = 4

P, m_dot, T = algebraic_equations(t, z0)


plt.barh(System.loc[:]['Node Description'],P[:]*1e-5)
for i in range(len(P)):
    plt.text(x= P[i]*1e-5,y= i,s= round(P[i]*1e-5,1), c='k')
plt.xlim(0,50)
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Pressure [bar]')
plt.show()


plt.barh(System.loc[2:]['Element Description'],m_dot)
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Mass Flow Rate [kg/s]')
plt.show()

#print('Pressures: ',P*1e-5)
#print('Mass Flows: ',m_dot)
'''

# Solve the DAE
sol = solve_ivp(dae_system, t_span, z0, t_eval=t_eval, method='LSODA')

# Access the solution
t = sol.t
m_0, m_1, m_2, m_3, m_4 = sol.y

#Plot some stuff


P = np.zeros((len(m_0),Nodes))
m_dot = np.zeros((len(m_0),FlowElements))
T = np.zeros(len(m_0))

for i in range(0, len(P)):
    P[i,:], m_dot[i,:], T[i] = algebraic_equations(t[i], [m_0[i],m_1[i],m_2[i],m_3[i],m_4[i]])

for i in range(0, FlowElements):
    label = 'm_dot_' + str(i+2)
    plt.plot(t,m_dot[:,i], label=label)
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Time [s]')
plt.ylabel('Mass Flow Rate [kg/s]')
plt.legend()
plt.show()

plt.plot(t,T, label='T')
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Time [s]')
plt.ylabel('Temperature [K]')
plt.legend()
plt.show()

for i in range(0, Nodes):
    label = 'P_' + str(i)
    plt.plot(t,P[:,i]*1e-5, label=label)
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Time [s]')
plt.ylabel('Pressure [bar]')
plt.ylim(0,55)
plt.legend()
plt.show()

for i in range(0, Nodes):
    label = 'P_' + str(i)
    plt.plot(t,P[:,i]*1e-5, label=label)
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Time [s]')
plt.ylabel('Pressure [bar]')
plt.ylim(0,400)
plt.legend()
plt.show()

plt.plot(t,m_0, label='m_0')
plt.plot(t,m_1, label='m_1')
plt.plot(t,m_2, label='m_2')
plt.plot(t,m_3, label='m_3')
plt.plot(t,m_4, label='m_4')
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Time [s]')
plt.ylabel('Mass [kg]')
plt.legend()
plt.show()


plt.plot(t,m_dot[:,3]/m_dot[:,4])
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Time [s]')
plt.ylabel('OF Ratio')
plt.show()
