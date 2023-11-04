import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp

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

V_1 = 2.5 #HP Tank volume [L]
V_2 = 13 #Prop tank volume [L]
A_3 = 45*1e-6 #Injector Orifice Area [m^2]
Cd_3 = 0.7 #Injector Orifice Discharge Coefficient

Kv_1_max = 1.6 #Reg Valve Maximum Kv

#-------------------------------------------------------
#-------------- Initial Conditions ---------------------
#-------------------------------------------------------

T_0_N2 = 300 #Initial temperature of nitrogen in the HP tank and ullage [K]
P_1_0 = 270 #Initial pressure of nitrogen in the HP tank [bar]
P_2_0 = 40 #Initial pressure of nitrogen ullage in prop tank [bar]
P_3_0 = 0 #Initial injector upstream pressure [bar]
m_dot_N2_0 = 0 #Initial N2 mass flow rate [kg/s]
m_dot_L_0 = 0 #Initial Propellant mass flow rate [kg/s]
m_3_0 = 11 #Initial mass of liquid propellant [kg]

#Opening_vs_norm_Kv = [[0,0.0322,0.0567,0.0967,0.1611,0.25,0.3833,0.64,0.86,1],[0,10,20,30,40,50,60,70,80,90]] taken from (scaled 8"): https://www.mydatabook.org/fluid-mechanics/flow-coefficient-opening-and-closure-curves-of-full-bore-ball-valves/
Opening_vs_norm_Kv = [[0,0.0145,0.0290,0.0580,0.0870,0.1304,0.2174,0.3188,0.5507,0.7826,1],[0,13.5,18,27,36,45,54,63,72,81,90]] # taken from (scaled 1/2"): https://www.valvesonline.co.uk/technical-information/mars-v-sector-manual-control-ball-valve-series-88.pdf

#------------------------------------------------------
#------------- Convert to SI units --------------------
#------------------------------------------------------

V_1, V_2 = V_1*1e-3, V_2*1e-3 #convert volumes to m^3
P_1_0, P_2_0 = P_1_0*1e5, P_2_0*1e5 #convert pressure to Pa

#------------------------------------------------------
#-------------- Derive calculated ICs -----------------
#------------------------------------------------------

rho_1_0 = (P_1_0)/(R_N2*T_0_N2) #Initial density of nitrogen in the HP tank [kg/m^3]
rho_2_0 = (P_2_0)/(R_N2*T_0_N2) #Initial density of nitrogen in the prop tank [kg/m^3]

m_2_0 = rho_2_0*(V_2-(m_3_0/rho_L))
m_1_0 = V_1*rho_1_0


def K1_Flow_Rate(P_1, P_2, rho_1, rho_2, T, t):
    P_1,P_2 = P_1*1e-5, P_2*1e-5
    
    Kp = 0.5
    K_v_controlled = 0.05 + (30 - P_2)*Kp
    
    if K_v_controlled < 0:
        K_v_controlled = 0
    elif K_v_controlled > Kv_1_max:
        K_v_controlled = Kv_1_max
    
    K_v = np.interp(t, [0, 0.6, 0.8, 10], [0, 0, K_v_controlled, K_v_controlled]) # Define the Kv of the valve based on a predefined path
    #K_v = 1
    # determine flow regime through valve
    if (P_2 >= 0.528*P_1):
        # non-choked flow
        Q_N = K_v*514*(((P_1-P_2)*P_2)/(rhost*T))/np.sqrt(abs(((P_1-P_2)*P_2)/(rhost*T)))
    else:
        # choked flow
        Q_N = K_v*257*P_1/np.sqrt(abs(rhost*T))
    m_dot = Q_N*(Pst/(R_N2*Tst))*(1/3600) #Convert standard flow rate to mass flow rate
    return m_dot, K_v

def K2_Flow_Rate(P_2, P_3, rho_1, rho_2, T, t):
    #K_v = np.interp(t, [0, 0.7, 0.9, 2, 2.1, 4, 4.1, 10], [0, 0, 10, 10, 1, 1, 0.3, 0.3]) # Define the Kv of the valve based on a predefined path
    #K_v = np.interp(t, [0, 0.7, 0.9, 2, 2.1, 4, 4.1, 10], [0, 0, 10, 10, 1, 1, 4, 4])
    K_v = np.interp(t, [0, 0.7, 0.9, 2, 2.1, 4, 4.1, 10], [0, 0, 10, 10, 10, 10, 10, 10])
    #K_v = 20
    P_3,P_2 = P_3*1e-5, P_2*1e-5
    SG = rho_L/1000 #Calculate specific gravity of propellant (value for water is 1)
    Q = K_v*((P_2 - P_3)/SG)/np.sqrt(abs((P_2 - P_3)/SG)) #Calculate volumetric flow rate based on Kv
    m_dot = Q*rho_L*(1/3600) #Convert standard flow rate to mass flow rate
    return m_dot

def explicit_algebraic_equations(m_2, m_3, rho_1):
    P_1 = P_1_0*((rho_1/rho_1_0)**gamma_N2)
    T = P_1/(rho_1*R_N2)
    rho_2 = m_2/(V_2 - (m_3/rho_L))
    P_2 = rho_2*T*R_N2
    m_1 = V_1*rho_1
    return P_1, T, rho_2, P_2, m_1


def algebraic_equations(x, P_1, P_2, rho_1, rho_2, m_3, T, t):
    P_3, m_dot_N2, m_dot_L = x
    F = np.zeros(3)
    if m_3 <= 0.1:
        F[0] = m_dot_L - 0
        F[1] = P_3 - P_atm
    else:
        F[0] = m_dot_L - Cd_3*A_3*2*rho_L*(P_3 - P_atm)/np.sqrt(abs(2*rho_L*(P_3 - P_atm)))
        F[1] = m_dot_L - K2_Flow_Rate(P_2, P_3, rho_1, rho_2, T, t)
    m_dot_reg, K_v_reg = K1_Flow_Rate(P_1, P_2, rho_1, rho_2, T, t)
    F[2] = m_dot_N2 - m_dot_reg
    return F


# Define the DAE system
def dae_system(t, z):
    m_2, m_3, rho_1 = z  # State variables: Mass of N2 in ullage, Mass of Liquid in tank, Density of N2 in HP tank
    
    #Algebraic equations to determine other state variables
    
    P_1, T, rho_2, P_2, m_1 = explicit_algebraic_equations(m_2, m_3, rho_1)
    
    solution = root(algebraic_equations, (1,1,1), args=(P_1, P_2, rho_1, rho_2, m_3, T, t), method='hybr', tol=1e-5)
    P_3, m_dot_N2, m_dot_L = solution.x
    
    # Differential Equations
    m_dot_2 = m_dot_N2
    
    if m_3 <= 0:
        m_3 = 0
        m_dot_3 = 0
    else:
        m_dot_3 = -m_dot_L
        
    rho_dot_1 = -(m_dot_N2/V_1)*((rho_1/rho_1_0)**(1-gamma_N2))

    return [m_dot_2, m_dot_3, rho_dot_1]

# Define initial conditions
z0 = [m_2_0, m_3_0, rho_1_0]


# Define the time span for integration
t_span = (0, 5)

# Generate evaluation times
t_eval = np.linspace(*t_span, 100)

# Solve the DAE
sol = solve_ivp(dae_system, t_span, z0, method='LSODA', t_eval=t_eval)

# Access the solution
t = sol.t
m_2, m_3, rho_1= sol.y
P_1, T, rho_2, P_2, m_1 = explicit_algebraic_equations(m_2, m_3, rho_1)


P_3 = np.zeros(len(t_eval))
m_dot_N2 = np.zeros(len(t_eval))
m_dot_L = np.zeros(len(t_eval))
K_v_reg = np.zeros(len(t_eval))
for i in range(0,len(t_eval)):
    solution = root(algebraic_equations, (1,1,1), args=(P_1[i], P_2[i], rho_1[i], rho_2[i], m_3[i], T[i], t[i]), method='hybr', tol=1e-5)
    m_dot_N2[i], K_v_reg[i] = K1_Flow_Rate(P_1[i], P_2[i], rho_1[i], rho_2[i], T[i], t[i])
    P_3[i], m_dot_N2[i], m_dot_L[i] = solution.x
    
Opening = np.interp(K_v_reg, np.transpose(Opening_vs_norm_Kv[0])*Kv_1_max, np.transpose(Opening_vs_norm_Kv[1]))

# Plot the results
import matplotlib.pyplot as plt


plt.plot(t,m_3, label='Prop Mass')
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Time [s]')
plt.ylabel('Mass [kg]')
plt.legend()
plt.show()
'''
plt.plot(t,rho_1, label='rho_1')
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Time [s]')
plt.ylabel('Density [kg/m^3]')
plt.legend()
plt.show()
'''
plt.plot(t,T, label='Nitrogen Temp.')
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Time [s]')
plt.ylabel('Temperature [K]')
plt.legend()
plt.show()


plt.plot(t,m_dot_L, label='Prop Mass Flow')
plt.plot(t,m_dot_N2, label='N2 Mass Flow')
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Time [s]')
plt.ylabel('Mass Flow Rate [kg/s]')
plt.legend()
plt.show()

plt.plot(t,K_v_reg, label='Regulator Kv')
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Time [s]')
plt.ylabel('Kv [m^3/hr]')
plt.legend()
plt.show()

plt.plot(np.transpose(Opening_vs_norm_Kv[0]), np.transpose(Opening_vs_norm_Kv[1]), label='Ball Valve')
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Kv (normalised)')
plt.ylabel('Opening Angle [deg]')
plt.legend()
plt.show()

reg_label = 'Reg Valve, Kv = ' + str(Kv_1_max)
plt.plot(t,Opening, label=reg_label)
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Time [s]')
plt.ylabel('Opening Angle [deg]')
plt.legend()
plt.show()

plt.plot(t,P_1*1e-5, label='HP N2 Tank')
plt.plot(t,P_2*1e-5, label='Prop Tank')
#plt.plot(t,P_3*1e-5, label='Injector')
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Time [s]')
plt.ylabel('Pressure [bar]')
plt.legend()
plt.show()

#plt.plot(t,P_1*1e-5, label='HP N2 Tank')
plt.plot(t,P_2*1e-5, label='Prop Tank')
plt.plot(t,P_3*1e-5, label='Injector')
plt.grid(which='major',axis='both',linewidth = 0.8)
plt.minorticks_on()
plt.grid(which='minor',axis='both',linewidth = 0.2)
plt.xlabel('Time [s]')
plt.ylabel('Pressure [bar]')
plt.legend()
plt.show()
