# Bipropellant-Propulsion-System-Model
Modelling and Simulation of the Bipropellant Feedsystem and Engine. Use for tuning the operating point of the engine, sizing components and active controller design

For initial testing of the solver, a simplified network representing the engine and plumbing system is used (see below)
![image](https://github.com/icl-rocketry/Bipropellant-Propulsion-System-Model/assets/87128082/481b7028-a13c-446a-89eb-f48c857762ae)

Node Number	Node Description
0	HP N2 Tank
1	Fuel Tank Ullage
2	Ox Tank Ullage
3	Fuel in Tank
4	Ox in Tank
5	Regulator Downstream
6	Combustion Chamber
7	Nozzle Exit
![image](https://github.com/icl-rocketry/Bipropellant-Propulsion-System-Model/assets/87128082/f13c83ce-71e8-467a-83e3-e66735c71131)

Element Number	Element Description
0	Ox Tank
1	Fuel Tank
2	Pressure Regulator
3	Ox Check Valve
4	Fuel Check Valve
5	Ox Injector
6	Fuel Injector
7	Nozzle
![image](https://github.com/icl-rocketry/Bipropellant-Propulsion-System-Model/assets/87128082/8b7ee27d-b53a-436b-90af-5d62cc926d68)



Some initial examples of pressure and mass flow distributions at a single timestep:
![image](https://github.com/icl-rocketry/Bipropellant-Propulsion-System-Model/assets/87128082/5fb3d786-7690-4d93-90ec-557f88d0cfa4)
![image](https://github.com/icl-rocketry/Bipropellant-Propulsion-System-Model/assets/87128082/1f77837c-3f14-4c27-81da-3326b6472d9f)


