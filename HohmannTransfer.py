from vpython import sphere, box, curve, arrow, vector, rate, color, scene
import numpy as np
import parameters
import sympy as sp

scene.width = 1920
scene.height = 1080

#scale_factor = 1 #1 simulation linear size unit is equivalent to 100km

def KeplerToCartesian(a,e,i,omega,Omega,T,t,mu): #Semi-Major axis, Eccentricity, Inclination, Argument of Periapsis, Longitude of Ascending node, Epoch
    #Mean anomaly
    n = np.sqrt(mu / a**3)
    MA = n * (t - T)  # Mean anomaly
    
    #Define symbolic variables
    EA = sp.Symbol('EA', real=True)
    #Solving for EA symbolically
    EA_solution = sp.nsolve(MA - EA + e * sp.sin(EA), EA, 0)
    #Convert to numerical value
    EA_numeric = float(EA_solution)
    
    #True anomaly (numerical)
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(EA_numeric / 2), np.sqrt(1 - e) * np.cos(EA_numeric / 2))
    
    #Orbital parameters (numerical)
    #Radius
    r = a * (1 - e**2) / (1 + e * np.cos(nu)) 
    #Semi-latus rectum
    p = r * (1 + e * np.cos(nu)) 
    #Angular momentum
    h = np.sqrt(mu * a * (1 - e**2))
    
    #Cartesian coordinates (numerical)
    rx = r * (np.cos(Omega) * np.cos(omega + nu) - np.sin(Omega) * np.sin(omega + nu) * np.cos(i))
    ry = r * (np.sin(Omega) * np.cos(omega + nu) + np.cos(Omega) * np.sin(omega + nu) * np.cos(i))
    rz = r * (np.sin(i) * np.sin(omega + nu))
    
    vx = ((rx * h * e) / (r * p)) * np.sin(nu) - (h / r) * (np.cos(Omega) * np.sin(omega + nu) + np.sin(Omega) * np.cos(omega + nu) * np.cos(i))
    vy = ((ry * h * e) / (r * p)) * np.sin(nu) - (h / r) * (np.sin(Omega) * np.sin(omega + nu) - np.cos(Omega) * np.cos(omega + nu) * np.cos(i))
    vz = ((rz * h * e) / (r * p)) * np.sin(nu) + (h / r) * (np.sin(i) * np.cos(omega + nu))

    return [rx,ry,rz,vx,vy,vz]

def CircularizationBurn():
    #print("Burn",np.linalg.norm([SpacecraftPos[0],SpacecraftPos[1],SpacecraftPos[2]]) - parameters.Planet_Radius)
    return 



#Defining initial conditions for simulation objects
t = 0

arrow_length = 2 * parameters.Planet_Radius 
x_axis = arrow(pos=vector(0, 0, 0), axis=vector(arrow_length, 0, 0), color=color.red, shaftwidth=0.01 * arrow_length)
y_axis = arrow(pos=vector(0, 0, 0), axis=vector(0, arrow_length, 0), color=color.green, shaftwidth=0.01 * arrow_length)
z_axis = arrow(pos=vector(0, 0, 0), axis=vector(0, 0, arrow_length), color=color.blue, shaftwidth=0.01 * arrow_length)

Planet = sphere(pos = vector(0,0,0), radius = parameters.Planet_Radius, color = color.blue)

TargetPosInitial = KeplerToCartesian(parameters.Target_alpha, parameters.Target_ecc, parameters.Target_i, parameters.Target_omega, parameters.Target_Omega,0,t,parameters.Planet_mu)
Target = sphere(pos = vector(TargetPosInitial[0],TargetPosInitial[1],TargetPosInitial[2]), radius = 200000, color = color.red) #Radius is arbitrary

SpacecraftPosInitial = KeplerToCartesian(parameters.Spacecraft_alpha, parameters.Spacecraft_ecc, parameters.Spacecraft_i, parameters.Spacecraft_omega, parameters.Spacecraft_Omega,0,t,parameters.Planet_mu)
Spacecraft = box(pos = vector(SpacecraftPosInitial[0],SpacecraftPosInitial[1],SpacecraftPosInitial[2]), length = 200000, height = 200000, width = 200000, color = color.green) 

TargetPath = curve(colour = color.white)
SpacecraftPath = curve(colour = color.white)



while True:
    rate(60)
    t += 10 #times dont make any sense times dont make any sense times dont make any sense times dont make any sense PLOT IT AND SEE PERIOD

    TargetPos = KeplerToCartesian(parameters.Target_alpha, parameters.Target_ecc, parameters.Target_i, parameters.Target_omega, parameters.Target_Omega,0,t,parameters.Planet_mu)

    SpacecraftPos = KeplerToCartesian(parameters.Spacecraft_alpha, parameters.Spacecraft_ecc, parameters.Spacecraft_i, parameters.Spacecraft_omega, parameters.Spacecraft_Omega,0,t,parameters.Planet_mu)

    if parameters.Spacecraft_ecc > 0.01: #Eccentricity less than 0.01 is arbitrarily set to be the bound for a circular orbit
        if parameters.Spacecraft_rP - 250 <= np.linalg.norm([SpacecraftPos[0],SpacecraftPos[1],SpacecraftPos[2]]) <= parameters.Spacecraft_rP + 250: #250m tolerance on when to burn
            CircularizationBurn()
            continue

    Target.pos.x, Target.pos.y, Target.pos.z = TargetPos[0],TargetPos[1],TargetPos[2]

    Spacecraft.pos.x, Spacecraft.pos.y, Spacecraft.pos.z = SpacecraftPos[0],SpacecraftPos[1],SpacecraftPos[2]
    #print(np.linalg.norm([SpacecraftPos[3],SpacecraftPos[4],SpacecraftPos[5]]))

    TargetPath.append(pos=vector(TargetPos[0],TargetPos[1],TargetPos[2]))

    SpacecraftPath.append(pos=vector(SpacecraftPos[0],SpacecraftPos[1],SpacecraftPos[2]))

    #print(t)

    

    








