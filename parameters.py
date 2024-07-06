import numpy as np
#Universal Constants
G = 6.674e-11 #Gravitational Constant

#Planet Parameters, by default this set to Earth
Planet_Mass = 5.972e24
Planet_Radius = 6371e3
Planet_mu = G*Planet_Mass #Mass of vessels are negligible


#Target Vessel Orbtial Parameters
Target_rA = 400e3 + Planet_Radius #Apoapsis
Target_rP = 400e3 + Planet_Radius #Periapsis
Target_alpha = (Target_rA + Target_rP)/2 #Semi major axis
Target_ecc = (Target_rA - Target_rP)/(Target_rA + Target_rP) #Eccentricity
Target_i = 100*(np.pi/180) #Inclination
Target_omega = 4*(np.pi/180) #Argument of Periapsis
Target_Omega = 70*(np.pi/180) #Longitude of Ascending Node


#Spacecraft Parameters
Spacecraft_DryMass = 5000 #5 Tonnes without fuel
Spacecraft_Fuel = 30000 #30 Tonnes of fuel on board at start of simulation
Spacecraft_SI = 250 #Specific impulse (s)
Spacecraft_q = 500 #Mass flow rate (kg/s)

#Spacecraft Orbtial Parameters
Spacecraft_rA = 300e3 + Planet_Radius #Apoapsis
Spacecraft_rP = 140e3 + Planet_Radius #Periapsis
Spacecraft_alpha = (Spacecraft_rA + Spacecraft_rP)/2 #Semi major axis
Spacecraft_ecc = (Spacecraft_rA - Spacecraft_rP)/(Spacecraft_rA + Spacecraft_rP) #Eccentricity
Spacecraft_i = 110*(np.pi/180) #Inclination
Spacecraft_omega = 14*(np.pi/180) #Argument of Periapsis
Spacecraft_Omega = 87*(np.pi/180) #Longitude of Ascending Node

