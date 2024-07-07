from vpython import sphere, box, curve, arrow, vector, rate, color, scene, button
import numpy as np
import sympy as sp
import parameters
import Graphs

scene.width = 1920
scene.height = 1080

#scale_factor = 1 #1 simulation linear size unit is equivalent to 100km

#Import parameters that will be changed
Target_alpha, Target_ecc, Target_i, Target_omega, Target_Omega, Target_rA, Target_rP = parameters.Target_alpha, parameters.Target_ecc, parameters.Target_i, parameters.Target_omega, parameters.Target_Omega, parameters.Target_rA, parameters.Target_rP
if (Target_rA == Target_rP):
    Target_omega = 0

Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega, Spacecraft_rA, Spacecraft_rP = parameters.Spacecraft_alpha, parameters.Spacecraft_ecc, parameters.Spacecraft_i, parameters.Spacecraft_omega, parameters.Spacecraft_Omega, parameters.Spacecraft_rA, parameters.Spacecraft_rP
if (Spacecraft_rA == Spacecraft_rP):
    Spacecraft_omega = 0

def KeplerToCartesian(a,e,i,omega,Omega,T,t,mu): #Semi-Major axis, Eccentricity, Inclination, Argument of Periapsis, Longitude of Ascending node, Epoch, Standard gravitational parameter
    #Mean anomaly
    n = np.sqrt(mu / a**3)
    MA = n * (t - T)  # Mean anomaly
    
    #Define symbolic variables
    #Eccentric anomaly
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

def CartesianToKepler(r_vec,v_vec):
    #Position and Velocity
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)
    
    #Angular Momentum
    h_vec = np.cross(r_vec,v_vec) 
    h = np.linalg.norm(h_vec) 

    #Eccentricity
    e_vec = np.cross(v_vec,h_vec)/parameters.Planet_mu - r_vec/r
    e = np.linalg.norm(e_vec)

    #Vector pointing towards the ascending node
    n_vec = [-h_vec[1],h_vec[0],0]
    n = np.linalg.norm(n_vec)

    #Inclination
    i = np.arccos(h_vec[2]/h)

    #Lonigitude of Ascending Node
    if (i > 4.8481e-6): #An orbit inclination smaller than 1 arcsecond (1/3600 of a degree = 4.8481e-6 radians) is considered small enough to be planar
        if (n_vec[1] >= 0):
            Omega = np.arccos(n_vec[0]/n)
        else:
            Omega = 2*np.pi - np.arccos(n_vec[0]/n)
    else:
        Omega = 0

    #Argument of Periapsis and True Anomaly 
    if (e > 1e-15): #An orbit with eccentricity less than 1e-15 is considered circular. There are too many computational errors that arise when computing an argument of periapsis for orbits that are extremely close to perfectly circular.
        if (e_vec[2] >= 0):
            omega = np.arccos(np.dot(n_vec,e_vec)/(n*e))
        else:
            omega = 2*np.pi - np.arccos(np.dot(n_vec,e_vec)/(n*e))
        if (np.dot(r_vec,v_vec) >= 0):
            nu = np.arccos(np.dot(e_vec,r_vec)/(e*r))
        else:
            nu = 2*np.pi - np.arccos(np.dot(e_vec,r_vec)/(e*r))
    else:
        omega = 0
        if (i > 4.8481e-6):
            if (r_vec[2] >= 0):
                nu = np.arccos(np.dot(n_vec,r_vec)/(n*r))
            else:
                nu = 2*np.pi - np.arccos(np.dot(n_vec,r_vec)/(n*r))
        else:
            if (v_vec[0 <= 0]):
                nu = np.arccos(r_vec[0]/r)
            else:
                nu = np.pi - np.arccos(r_vec[0]/r)
    
    #Semi-major axis
    alpha = 1/((2/r)-(v**2/parameters.Planet_mu)) 

    #Eccentric anomaly
    EA = np.arctan([((1-e)/(1+e))**(1/2)*np.tan(nu/2)])*2

    #Epoch
    n = np.sqrt(parameters.Planet_mu/alpha**3)
    T = t - 1/n*(EA-e*np.sin(EA))

    return [alpha,e,i,omega,Omega,nu,T]

def CalculateTargetV(r_vec,v_vec,Kep_i):

    return

def CircularizeBurn(r_vec, v_vec, Kep_i):
    #print("Burn",np.linalg.norm([SpacecraftPos[0],SpacecraftPos[1],SpacecraftPos[2]]) - parameters.Planet_Radius)

    return 

running = True
def ExitProgram():
    global running
    running = False
    return

#Defining initial conditions for simulation objects
t = 0

quit_button = button(text = "Quit", bind = ExitProgram) 

arrow_length = 2 * parameters.Planet_Radius 
x_axis = arrow(pos=vector(0, 0, 0), axis=vector(arrow_length, 0, 0), color=color.red, shaftwidth=0.01 * arrow_length)
y_axis = arrow(pos=vector(0, 0, 0), axis=vector(0, arrow_length, 0), color=color.green, shaftwidth=0.01 * arrow_length)
z_axis = arrow(pos=vector(0, 0, 0), axis=vector(0, 0, arrow_length), color=color.blue, shaftwidth=0.01 * arrow_length)

Planet = sphere(pos = vector(0,0,0), radius = parameters.Planet_Radius, color = color.blue)

TargetPosInitial = KeplerToCartesian(Target_alpha, Target_ecc, Target_i, Target_omega, Target_Omega,parameters.Target_T,t,parameters.Planet_mu) #need to fix the 0 at some point
Target = sphere(pos = vector(TargetPosInitial[0],TargetPosInitial[1],TargetPosInitial[2]), radius = 200000, color = color.red) #Size is arbitrary

SpacecraftPosInitial = KeplerToCartesian(Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega,parameters.Spacecraft_T,t,parameters.Planet_mu)
Spacecraft = box(pos = vector(SpacecraftPosInitial[0],SpacecraftPosInitial[1],SpacecraftPosInitial[2]), length = 200000, height = 200000, width = 200000, color = color.green) #Size is arbitrary

TargetPath = curve(colour = color.white)
SpacecraftPath = curve(colour = color.white)

#Simulation runtime loop
while running: 
    rate(60)
    t += 1/60*parameters.TimeScale # PLOT IT AND SEE PERIOD

    TargetPos = KeplerToCartesian(Target_alpha, Target_ecc, Target_i, Target_omega, Target_Omega,parameters.Target_T,t,parameters.Planet_mu)
    SpacecraftPos = KeplerToCartesian(Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega,parameters.Spacecraft_T,t,parameters.Planet_mu)

    if Spacecraft_ecc > 1e-15: #Eccentricity less than 1e-15 is arbitrarily set to be the bound for a circular orbit
        if Spacecraft_rP - 1 <= np.linalg.norm([SpacecraftPos[0],SpacecraftPos[1],SpacecraftPos[2]]) <= Spacecraft_rP + 1: #10m tolerance on when to burn
            #CircularizationBurn()
            #print("Burning")
            continue
        #Update position after burn
        TargetPos = KeplerToCartesian(Target_alpha, Target_ecc, Target_i, Target_omega, Target_Omega,0,t,parameters.Planet_mu)
        SpacecraftPos = KeplerToCartesian(Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega,0,t,parameters.Planet_mu)
    
    Target.pos.x, Target.pos.y, Target.pos.z = TargetPos[0],TargetPos[1],TargetPos[2]
    Spacecraft.pos.x, Spacecraft.pos.y, Spacecraft.pos.z = SpacecraftPos[0],SpacecraftPos[1],SpacecraftPos[2]

    TargetPath.append(pos=vector(TargetPos[0],TargetPos[1],TargetPos[2]))

    SpacecraftPath.append(pos=vector(SpacecraftPos[0],SpacecraftPos[1],SpacecraftPos[2]))

    Graphs.v_x.append(SpacecraftPos[3])
    Graphs.v_y.append(SpacecraftPos[4])
    Graphs.v_z.append(SpacecraftPos[5])
    Graphs.r.append(np.linalg.norm([SpacecraftPos[0],SpacecraftPos[1],SpacecraftPos[2]]))
    Graphs.v.append(np.linalg.norm([SpacecraftPos[3],SpacecraftPos[4],SpacecraftPos[5]]))
    Graphs.t.append(t)

    #print(t)
    # print(Target_alpha, Target_ecc, Target_i, Target_omega, Target_Omega,"K")
    # print(TargetPosCK[0],TargetPosCK[1],TargetPosCK[2],TargetPosCK[3],TargetPosCK[4])

print("Exiting Simulation.")
Graphs.main()




    


    

    








