from vpython import sphere, box, curve, arrow, vector, rate, color, scene, button
import numpy as np
import sympy as sp
from scipy.optimize import minimize
import Parameters
import Graphs

#scale_factor = 1 #1 simulation linear size unit is equivalent to 100km

#Import Parameters that will be changed
Target_alpha, Target_ecc, Target_i, Target_omega, Target_Omega, Target_rA, Target_rP, Target_T = Parameters.Target_alpha, Parameters.Target_ecc, Parameters.Target_i, Parameters.Target_omega, Parameters.Target_Omega, Parameters.Target_rA, Parameters.Target_rP, Parameters.Target_T
if (Target_rA == Target_rP):
    Target_omega = 0

Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega, Spacecraft_rA, Spacecraft_rP, Spacecraft_T = Parameters.Spacecraft_alpha, Parameters.Spacecraft_ecc, Parameters.Spacecraft_i, Parameters.Spacecraft_omega, Parameters.Spacecraft_Omega, Parameters.Spacecraft_rA, Parameters.Spacecraft_rP, Parameters.Spacecraft_T
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
    
    #Orbital Parameters (numerical)
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
    e_vec = np.cross(v_vec,h_vec)/Parameters.Planet_mu - r_vec/r
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
    if (e > 1e-10 and n*e != 0): #An orbit with eccentricity less than 1e-10 is considered circular. There are too many computational errors that arise when computing an argument of periapsis for orbits that are extremely close to perfectly circular.
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
            if (v_vec[0] <= 0):
                nu = np.arccos(r_vec[0]/r)
            else:
                nu = np.pi - np.arccos(r_vec[0]/r)
    
    #Semi-major axis
    alpha = 1/((2/r)-(v**2/Parameters.Planet_mu)) 

    #Eccentric anomaly (IF CIRCULAR, THIS IS DIFFERENT)
    EA = np.arctan([((1-e)/(1+e))**(1/2)*np.tan(nu/2)])*2

    #Epoch
    n = np.sqrt(Parameters.Planet_mu/alpha**3)
    T = t - 1/n*(EA-e*np.sin(EA))

    while (nu > 2*np.pi):
        nu-=2*np.pi

    return [alpha,e,i,omega,Omega,nu,T[0]]

def CalculateTargetV(R, v_initial):
    # Objective function for optimization
    def objective(v):
        V = v_initial + v  # Modify initial velocity by v
        kepler_elements = CartesianToKepler(R, V)
        return kepler_elements[1]  # Return eccentricity
    
    # Initial guess for the optimization (small perturbation)
    initial_guess = np.zeros(3)
    
    # Optimization bounds (could be adjusted based on physical constraints)
    bounds = [(-100e100, 100e100), (-100e100, 100e100), (-100e100, 100e100)]  # Example bounds
    
    # Perform optimization
    result = minimize(objective, initial_guess, bounds=bounds)
    #print(result.x)
    # Extract optimized velocity vector
    v_optimized = result.x
    
    # Return the optimized velocity vector
    return v_initial + v_optimized


def CircularizeBurn(r_vec, v_vec, Kep_i):
    #print("Burn",np.linalg.norm([SpacecraftPos[0],SpacecraftPos[1],SpacecraftPos[2]]) - Parameters.Planet_Radius)

    return 

running = True
def ExitProgram():
    global running
    running = False
    return

#Defining initial conditions for simulation objects
scene.width = 1920
scene.height = 1080

t = 0

quit_button = button(text = "Quit", bind = ExitProgram) 

arrow_length = 2 * Parameters.Planet_Radius 
x_axis = arrow(pos=vector(0, 0, 0), axis=vector(arrow_length, 0, 0), color=color.red, shaftwidth=0.001 * arrow_length)
y_axis = arrow(pos=vector(0, 0, 0), axis=vector(0, arrow_length, 0), color=color.green, shaftwidth=0.001 * arrow_length)
z_axis = arrow(pos=vector(0, 0, 0), axis=vector(0, 0, arrow_length), color=color.blue, shaftwidth=0.001 * arrow_length)

Planet = sphere(pos = vector(0,0,0), radius = Parameters.Planet_Radius, color = color.blue)

TargetPosInitial = KeplerToCartesian(Target_alpha, Target_ecc, Target_i, Target_omega, Target_Omega,Target_T,t,Parameters.Planet_mu) #need to fix the 0 at some point
Target = sphere(pos = vector(TargetPosInitial[0],TargetPosInitial[1],TargetPosInitial[2]), radius = 200000, color = color.red) #Size is arbitrary

SpacecraftPosInitial = KeplerToCartesian(Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega, Spacecraft_T,t,Parameters.Planet_mu)
Spacecraft = box(pos = vector(SpacecraftPosInitial[0],SpacecraftPosInitial[1],SpacecraftPosInitial[2]), length = 200000, height = 200000, width = 200000, color = color.green) #Size is arbitrary

TargetPath = curve(colour = color.white)
SpacecraftPath = curve(colour = color.white)
burn = 0

#Simulation runtime loop
while running: 
    rate(60)
    t += 1/60*Parameters.TimeScale # PLOT IT AND SEE PERIOD

    TargetCart = KeplerToCartesian(Target_alpha, Target_ecc, Target_i, Target_omega, Target_Omega,Target_T,t,Parameters.Planet_mu)
    SpacecraftCart = KeplerToCartesian(Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega,Spacecraft_T,t,Parameters.Planet_mu)

    if Spacecraft_ecc > 1e-10: #Eccentricity less than 1e-10 is arbitrarily set to be the bound for a circular orbit
        if Spacecraft_rA - 100 <= np.linalg.norm([SpacecraftCart[0],SpacecraftCart[1],SpacecraftCart[2]]) <= Spacecraft_rA + 100: #100m tolerance on when to burn
            #Parameters.main()
            #print(SpacecraftCart)
            v_T = CalculateTargetV([SpacecraftCart[0],SpacecraftCart[1],SpacecraftCart[2]],[SpacecraftCart[3],SpacecraftCart[4],SpacecraftCart[5]])
            new_kep = CartesianToKepler([SpacecraftCart[0],SpacecraftCart[1],SpacecraftCart[2]],v_T)
            #print(new_kep[0], new_kep[1], new_kep[2], new_kep[3], new_kep[4], new_kep[6])
            Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega, Spacecraft_T = new_kep[0], new_kep[1], new_kep[2], new_kep[3], new_kep[4], new_kep[6]
            if burn == 0: 
                SpacecraftPath.clear()
                burn = 1
            #continue
        #Update position after burn
        TargetCart = KeplerToCartesian(Target_alpha, Target_ecc, Target_i, Target_omega, Target_Omega, Target_T,t,Parameters.Planet_mu)
        SpacecraftCart = KeplerToCartesian(Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega, Spacecraft_T,t,Parameters.Planet_mu)
    
    Target.pos.x, Target.pos.y, Target.pos.z = TargetCart[0],TargetCart[1],TargetCart[2]
    Spacecraft.pos.x, Spacecraft.pos.y, Spacecraft.pos.z = SpacecraftCart[0],SpacecraftCart[1],SpacecraftCart[2]

    TargetPath.append(pos=vector(TargetCart[0],TargetCart[1],TargetCart[2]))
    SpacecraftPath.append(pos=vector(SpacecraftCart[0],SpacecraftCart[1],SpacecraftCart[2]))

    Graphs.v_x.append(SpacecraftCart[3])
    Graphs.v_y.append(SpacecraftCart[4])
    Graphs.v_z.append(SpacecraftCart[5])
    Graphs.r.append(np.linalg.norm([SpacecraftCart[0],SpacecraftCart[1],SpacecraftCart[2]]))
    Graphs.v.append(np.linalg.norm([SpacecraftCart[3],SpacecraftCart[4],SpacecraftCart[5]]))
    Graphs.t.append(t)

    #print(t)
    #print(SpacecraftCart)
    

print("Exiting Simulation.")
Graphs.main()




    


    

    








