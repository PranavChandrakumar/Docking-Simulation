from vpython import sphere, box, curve, arrow, vector, rate, color, scene, button, label
import numpy as np
import sympy as sp
from scipy.optimize import minimize
import Parameters
import Graphs
import Maneuvers

#Import Parameters that will not be changed
TimeScale = Parameters.TimeScale
Planet_mu = Parameters.Planet_mu
Planet_Radius = Parameters.Planet_Radius

Target_alpha, Target_ecc, Target_i, Target_omega, Target_Omega, Target_rA, Target_rP, Target_T, Target_Period = Parameters.Target_alpha, Parameters.Target_ecc, Parameters.Target_i, Parameters.Target_omega, Parameters.Target_Omega, Parameters.Target_rA, Parameters.Target_rP, Parameters.Target_T, Parameters.Target_Period
if (Target_rA == Target_rP):
    Target_omega = 0
if (Target_T >= Target_Period):
    Target_T = np.mod(Target_T,Target_Period)

Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega, Spacecraft_rA, Spacecraft_rP, Spacecraft_T, Spacecraft_Period = Parameters.Spacecraft_alpha, Parameters.Spacecraft_ecc, Parameters.Spacecraft_i, Parameters.Spacecraft_omega, Parameters.Spacecraft_Omega, Parameters.Spacecraft_rA, Parameters.Spacecraft_rP, Parameters.Spacecraft_T, Parameters.Spacecraft_Period
if (Spacecraft_rA == Spacecraft_rP):
    Spacecraft_omega = 0
if (Spacecraft_T >= Spacecraft_Period):
    Spacecraft_T = np.mod(Spacecraft_T,Spacecraft_Period)

def KeplerToCartesian(a,e,i,omega,Omega,T,t,mu): #Semi-Major axis, Eccentricity, Inclination, Argument of Periapsis, Longitude of Ascending node, Epoch, Standard gravitational parameter
    #Mean anomaly
    n = np.sqrt(mu / a**3)
    MA = n * (t - T)  

    #Eccentric anomaly
    EA = sp.Symbol('EA', real=True)
    #Solving for EA symbolically
    EA_solution = sp.nsolve(MA - EA + e * sp.sin(EA), EA, 0,tol=1e-5)
    #Convert to numerical value
    EA_numeric = float(EA_solution)
    
    #True anomaly (numerical)
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(EA_numeric / 2), np.sqrt(1 - e) * np.cos(EA_numeric / 2))
  
    #Orbital Parameters (numerical)
    #Radius
    if (1 + e * np.cos(nu)) == 0:
        r = a
    else:
        r = a * (1 - e**2) / (1 + e * np.cos(nu)) 
    if r == 0:
        r = a
    
    #Semi-latus rectum
    p = r * (1 + e * np.abs(np.cos(nu)))
    # if p == 0:
    #     p = 1

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
    e_vec = np.cross(v_vec,h_vec)/Planet_mu - r_vec/r
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
    alpha = 1/((2/r)-(v**2/Planet_mu)) 

    #Eccentric anomaly (IF CIRCULAR, THIS IS DIFFERENT)
    
    EA = np.arctan([((1-e)/(1+e))**(1/2)*np.tan(nu/2)])*2

    #Time to periapsis 
    n = np.sqrt(Planet_mu/alpha**3)
    T = t - 1/n*(EA-e*np.sin(EA))

    while (nu > 2*np.pi):
        nu-=2*np.pi

    return [alpha,e,i,omega,Omega,nu,T[0]]

def TimeScaleChange(factor):
    global TimeScale
    if (1 <= TimeScale*factor <= 10000):
        TimeScale = int(TimeScale*factor)
    elif (1 == TimeScale):
        print ("TimeScale is at", TimeScale,"cannot siumate lower than 1 second per second.")
    else:
        print ("TimeScale is at", TimeScale,"cannot siumate higher than 10000 seconds per second.")
    # time_display.text = 'Time Scale: {} seconds per second'.format(TimeScale)
    return

running = True
def ExitProgram():
    global running
    running = False
    return

#Defining initial conditions for simulation objects
scene.width = 1905
scene.height = 960

t = 0

quit_button = button(text = "Quit", bind = ExitProgram)

time_display = label(
    text='Time Scale: {} seconds per second'.format(TimeScale),
    xoffset=-scene.width/2 + 280,  
    yoffset=scene.height/2 - 25,   
    height=15,    
    box=True,    
)

slower_button =  button(text = "Decrease Simulation Speed by 10x", bind = lambda: TimeScaleChange(0.1))
faster_button =  button(text = "Increase Simulation Speed by 10x", bind = lambda: TimeScaleChange(10))


arrow_length = 2 * Planet_Radius 
x_axis = arrow(pos=vector(0, 0, 0), axis=vector(arrow_length, 0, 0), color=color.red, shaftwidth=0.001 * arrow_length)
y_axis = arrow(pos=vector(0, 0, 0), axis=vector(0, arrow_length, 0), color=color.green, shaftwidth=0.001 * arrow_length)
z_axis = arrow(pos=vector(0, 0, 0), axis=vector(0, 0, arrow_length), color=color.blue, shaftwidth=0.001 * arrow_length)

Planet = sphere(pos = vector(0,0,0), radius = Planet_Radius/100, color = color.blue)

TargetPosInitial = KeplerToCartesian(Target_alpha, Target_ecc, Target_i, Target_omega, Target_Omega,Target_T,t,Planet_mu)
Target = sphere(pos = vector(TargetPosInitial[0],TargetPosInitial[1],TargetPosInitial[2]), radius = 200000, color = color.red) #Size is arbitrary

SpacecraftPosInitial = KeplerToCartesian(Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega, Spacecraft_T,t,Planet_mu)
Spacecraft = box(pos = vector(SpacecraftPosInitial[0],SpacecraftPosInitial[1],SpacecraftPosInitial[2]), length = 200000, height = 200000, width = 200000, color = color.green) #Size is arbitrary

TargetPath = curve(colour = color.white)
SpacecraftPath = curve(colour = color.white)

r0 = [[SpacecraftPosInitial[0], SpacecraftPosInitial[1], SpacecraftPosInitial[2]]]
v0 = [[SpacecraftPosInitial[3], SpacecraftPosInitial[4], SpacecraftPosInitial[5]]]
r1 = [[TargetPosInitial[0],TargetPosInitial[1], TargetPosInitial[2]]]
v1 = [[TargetPosInitial[3], TargetPosInitial[4], TargetPosInitial[5]]]
burn = 0

dv = []
dvmag = []
pos = []
C_deltaV = np.array([])
C_pos = np.array([])

#Simulation runtime loop
while running: 
    rate(60)
    t += 1/60*TimeScale 

    TargetCart = KeplerToCartesian(Target_alpha, Target_ecc, Target_i, Target_omega, Target_Omega,Target_T,t,Planet_mu)
    SpacecraftCart = KeplerToCartesian(Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega,Spacecraft_T,t,Planet_mu)

    if t <= Target_Period:
        r1.append([TargetCart[0],TargetCart[1],TargetCart[2]])
        v1.append([TargetCart[3],TargetCart[4],TargetCart[5]])
    if t <= Spacecraft_Period:
        r0.append([SpacecraftCart[0],SpacecraftCart[1],SpacecraftCart[2]])
        v0.append([SpacecraftCart[3],SpacecraftCart[4],SpacecraftCart[5]])
    # if t > max(Target_Period, Spacecraft_Period) and burn == 0:
    #     burn = 1
    #     print("Computing optimal delta V...")
    #     r0_len = len(r0)
    #     r1_len = len(r1)
    #     for i in range(0,min(r0_len,r1_len)):
    #         for j in range (0,max(r0_len,r1_len)):
    #             if r0_len < r1_len:
    #                 deltaV = Maneuvers.CalculateDeltaV(r0[i],r1[j],v0[i],v1[j])
    #             else:
    #                 deltaV = Maneuvers.CalculateDeltaV(r0[j],r1[i],v0[j],v1[i])
    #             vdiff = deltaV[1]-deltaV[0]
    #             dv.append(vdiff)
    #             dvmag.append(np.linalg.norm(vdiff))
    #             pos.append(deltaV[2])
    #             if len(deltaV) == 4:
    #                 burn = 2
    #                 solution = 1
    #                 break
    #         if burn == 2:
    #             break
        
    #     deltaV = min(dvmag)
    #     index = dvmag.index(deltaV)
    #     position = pos[index]
    #     deltaV = dv[index]
    #     C_deltaV = np.append(C_deltaV,deltaV)
    #     C_pos = np.append(C_pos,position)
    #     maneuver = [C_pos,C_deltaV]  
    #     burn = 1
    
    # if burn == 1 and solution == 1 and (maneuver[0][0] - 20000 <= SpacecraftCart[0] <= maneuver[0][0] + 20000) and (maneuver[0][1] - 20000 <= SpacecraftCart[1] <= maneuver[0][1] + 20000) and (maneuver[0][2] - 20000 <= SpacecraftCart[2] <= maneuver[0][2] + 20000):
    #     TimeScale = 1
    #     if (maneuver[0][0] - 1000 <= SpacecraftCart[0] <= maneuver[0][0] + 1000) and (maneuver[0][1] - 1000 <= SpacecraftCart[1] <= maneuver[0][1] + 1000) and (maneuver[0][2] - 1000 <= SpacecraftCart[2] <= maneuver[0][2] + 1000):
    #         new_kep = CartesianToKepler(maneuver[0],maneuver[1])
    #         Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega, Spacecraft_T = new_kep[0], new_kep[1], new_kep[2], new_kep[3], new_kep[4], new_kep[6]            
    #         SpacecraftPath.color = color.red
    #         SpacecraftCart = KeplerToCartesian(Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega,Spacecraft_T,t,Planet_mu)
    #         burn = 2
    # elif burn == 1 and solution == 0:
    #     print("No maneouver calculated.")
    #     solution = 2

########################################################### 
    # if t > max(Target_Period, Spacecraft_Period) and burn == 0:
    #     burn = 1
    #     print("Computing optimal delta V...")
    #     r0_len = len(r0)
    #     r1_len = len(r1)
    #     for j in range (0,max(r0_len,r1_len)):
    #         if r0_len < r1_len:
    #             deltaV = Maneuvers.CalculateDeltaV(r0[0],r1[j],v0[0],v1[j])
    #         else:
    #             deltaV = Maneuvers.CalculateDeltaV(r0[j],r1[0],v0[j],v1[0])
    #         vdiff = deltaV[1]-deltaV[0]
    #         dv.append(vdiff)
    #         dvmag.append(np.linalg.norm(vdiff))
    #         pos.append(deltaV[2])
        
    #     deltaV = min(dvmag)
    #     index = dvmag.index(deltaV)
    #     position = pos[index]
    #     deltaV = dv[index]
    #     C_deltaV = np.append(C_deltaV,deltaV)
    #     C_pos = np.append(C_pos,position)
    #     maneuver = [C_pos,C_deltaV]
          
    # if burn == 1 and (maneuver[0][0] - 20000 <= SpacecraftCart[0] <= maneuver[0][0] + 20000) and (maneuver[0][1] - 20000 <= SpacecraftCart[1] <= maneuver[0][1] + 20000) and (maneuver[0][2] - 20000 <= SpacecraftCart[2] <= maneuver[0][2] + 20000):
    #     TimeScale = 1
    #     if (maneuver[0][0] - 1000 <= SpacecraftCart[0] <= maneuver[0][0] + 1000) and (maneuver[0][1] - 1000 <= SpacecraftCart[1] <= maneuver[0][1] + 1000) and (maneuver[0][2] - 1000 <= SpacecraftCart[2] <= maneuver[0][2] + 1000):
    #         print("fail")
    #         new_kep = CartesianToKepler(maneuver[0],maneuver[1])
    #         Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega, Spacecraft_T = new_kep[0], new_kep[1], new_kep[2], new_kep[3], new_kep[4], new_kep[6]        
    #         print(Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega, Spacecraft_T,t)   
    #         SpacecraftPath.color = color.red
    #         SpacecraftCart = KeplerToCartesian(Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega,Spacecraft_T,t,Planet_mu)
    #         burn = 2
################################################
    if t > max(Target_Period, Spacecraft_Period) and burn == 0:
        burn = 1
        print("Computing optimal delta V...")
        r0_len = len(r0)
        r1_len = len(r1)
        for i in range(0,min(r0_len,r1_len)) :
            for j in range (0,max(r0_len,r1_len)):
                if r0_len < r1_len:
                    deltaV = Maneuvers.CalculateDeltaV(r0[i],r1[j],v0[i],v1[j])
                else:
                    deltaV = Maneuvers.CalculateDeltaV(r0[j],r1[i],v0[j],v1[i])
                if len(deltaV == 3):
                    vdiff = deltaV[1]-deltaV[0]
                    dv.append(vdiff)
                    dvmag.append(np.linalg.norm(vdiff))
                    pos.append(deltaV[2])
                else:
                    print("no")

        deltaV = min(dvmag)
        index = dvmag.index(deltaV)
        position = pos[index]
        deltaV = dv[index]
        C_deltaV = np.append(C_deltaV,deltaV)
        C_pos = np.append(C_pos,position)
        maneuver = [C_pos,C_deltaV]
        print(maneuver)
        
    if burn == 1 and (maneuver[0][0] - 20000 <= SpacecraftCart[0] <= maneuver[0][0] + 20000) and (maneuver[0][1] - 20000 <= SpacecraftCart[1] <= maneuver[0][1] + 20000) and (maneuver[0][2] - 20000 <= SpacecraftCart[2] <= maneuver[0][2] + 20000):
        TimeScale = 1
        if (maneuver[0][0] - 1000 <= SpacecraftCart[0] <= maneuver[0][0] + 1000) and (maneuver[0][1] - 1000 <= SpacecraftCart[1] <= maneuver[0][1] + 1000) and (maneuver[0][2] - 1000 <= SpacecraftCart[2] <= maneuver[0][2] + 1000):
            print("fail")
            new_kep = CartesianToKepler(maneuver[0],maneuver[1])
            Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega, Spacecraft_T = new_kep[0], new_kep[1], new_kep[2], new_kep[3], new_kep[4], new_kep[6]        
            print(Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega, Spacecraft_T,t)   
            SpacecraftPath.color = color.red
            SpacecraftCart = KeplerToCartesian(Spacecraft_alpha, Spacecraft_ecc, Spacecraft_i, Spacecraft_omega, Spacecraft_Omega,Spacecraft_T,t,Planet_mu)
            burn = 2


    Target.pos.x, Target.pos.y, Target.pos.z = TargetCart[0],TargetCart[1],TargetCart[2]
    Spacecraft.pos.x, Spacecraft.pos.y, Spacecraft.pos.z = SpacecraftCart[0],SpacecraftCart[1],SpacecraftCart[2]

    TargetPath.append(pos=vector(TargetCart[0],TargetCart[1],TargetCart[2]))
    SpacecraftPath.append(pos=vector(SpacecraftCart[0],SpacecraftCart[1],SpacecraftCart[2]))

    time_display.text = 'Time Scale: {} seconds per second'.format(TimeScale)

    Graphs.v_x.append(SpacecraftCart[3])
    Graphs.v_y.append(SpacecraftCart[4])
    Graphs.v_z.append(SpacecraftCart[5])
    Graphs.r.append(np.linalg.norm([SpacecraftCart[0],SpacecraftCart[1],SpacecraftCart[2]]))
    Graphs.v.append(np.linalg.norm([SpacecraftCart[3],SpacecraftCart[4],SpacecraftCart[5]]))
    Graphs.t.append(t)

print("Exiting Simulation.")

Graphs.main()





    


    

    








