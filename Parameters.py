import numpy as np
from scipy.constants import G
#Simulation Parameters
TimeScale = 10 #This scales relative timescale of the simulation. Default value is 100, meaning one real time second is equivalent to 100 seconds in the simulation. 
TargetOptimizationVariable =  0  #0 to optimize for time (faster intercept), 1 to optimize for fuel (smaller deltaV)

#Universal Constants
G #Gravitational Constant

#Planet Parameters, by default this set to Earth
Planet_Mass = 5.972e24 #(kg)
Planet_Radius = 6371e3 #Assuming a perfectly spherical planet for this simulation, use equatorial rasius (m)
Planet_StandardGravity = 9.807 #Acceleration due to gravity at planet surface (m/s^2)
Planet_mu = G*Planet_Mass #Standard gravitational parameter (m^3s^-2), mass of vessels are negligible 


#Target Vessel Orbtial Parameters
Target_rA = 5000e3 + Planet_Radius #Apoapsis (m)
Target_rP = 500e3 + Planet_Radius  #Periapsis (m)
Target_alpha = (Target_rA + Target_rP)/2 #Semi major axis (m)
Target_ecc = (Target_rA - Target_rP)/(Target_rA + Target_rP) #Eccentricity
Target_i = 45*(np.pi/180) #Inclination (Degree input, converted to radians)
Target_omega = 0*(np.pi/180) #Argument of Periapsis (Degree input, converted to radians)
Target_Omega = 0*(np.pi/180) #Longitude of Ascending Node (Degree input, converted to radians)
Target_T = 0 #Time to periapsis passage (s)
Target_Period = np.sqrt(((4*np.pi**2)*Target_alpha**3)/Planet_mu)

#Spacecraft Parameters
Spacecraft_DryMass = 5000 #Mass without fuel (kg)
Spacecraft_Fuel = 30000 #Fuel on board at start of simulation (kg)
Spacecraft_SI = 250 #Specific impulse (s)
Spacecraft_q = 500 #Mass flow rate (kg/s)
Spacecraft_Thrust = Spacecraft_SI*Planet_StandardGravity*Spacecraft_q #Thrust (N). 

#Spacecraft Orbtial Parameters
Spacecraft_rA = 400e3 + Planet_Radius #Apoapsis (m)
Spacecraft_rP = 3000e3 + Planet_Radius  #Periapsis (m)
Spacecraft_alpha = (Spacecraft_rA + Spacecraft_rP)/2 #Semi major axis (m)
Spacecraft_ecc = (Spacecraft_rA - Spacecraft_rP)/(Spacecraft_rA + Spacecraft_rP) #Eccentricity
Spacecraft_i = 250*(np.pi/180) #Inclination (Degree input, converted to radians) 
Spacecraft_omega = 6*(np.pi/180) #Argument of Periapsis (Degree input, converted to radians)
Spacecraft_Omega = 0*(np.pi/180) #Longitude of Ascending Node (Degree input, converted to radians)
Spacecraft_T = 1000 #Time to periapsis passage (s)
Spacecraft_Period = np.sqrt(((4*np.pi**2)*Spacecraft_alpha**3)/Planet_mu)


def main():
    print(Spacecraft_alpha,Spacecraft_ecc,Spacecraft_i,Spacecraft_omega,Spacecraft_Omega)
    return
