import numpy as np
import sympy as sp
from scipy.optimize import newton
from scipy.linalg import norm
import Parameters
mu = Parameters.Planet_mu

def CalculateDeltaV(ri,rf,vi,vf):
    ri,rf,vi,vf = np.array(ri), np.array(rf), np.array(vi), np.array(vf)
    c = np.linalg.norm(rf-ri)
    v_mean = (np.linalg.norm(vi)+np.linalg.norm(vf))/2
    alpha = c/2
    #ecc = np.sqrt(1-((c*v_mean)/mu)**2)
    #T = np.sqrt(((4*np.pi**2)*alpha**3)/mu)
    #print(T)
    #deltaV = np.array(LambertsProblem(ri,rf,T,1,mu,-4*np.pi,0.01,4*np.pi**2,200,1e-5))
    #deltaV = np.array(LambertsProblem1(ri,rf,vi,T,mu))

    deltaV = np.array(LambertsProblem2(ri,rf,vi,mu))
    return deltaV

def C2(psi):
    return ((1-np.cos(np.sqrt(psi)))/psi)

def C3(psi):
    return ((np.sqrt(psi)-np.sin(np.sqrt(psi)))/(psi*np.sqrt(psi)))

def LambertsProblem(ri,rf,dt,tm,mu,psiL,psiI,psiU,max,tol): #Arabian Paper
    mu_sqrt = np.sqrt(mu)
    
    #Norms of position vectors
    ri_norm = np.linalg.norm(ri)
    rf_norm = np.linalg.norm(rf)
    
    gamma = np.dot(rf,ri)/(ri_norm*rf_norm)

    Beta = tm*np.sqrt((1-gamma**2))
    
    A = tm*np.sqrt((rf_norm*ri_norm*(1+gamma)))
     #print(A)
    if A == 0:
        return [[0,0,0],[(2**32)-1,(2**32)-1,(2**32)-1],[0,0,0]] 
    
    c2 = 0.5
    c3 = 0.6

    Solution = False

    for i in range(0,max):
        psi = psiI
        
        B = ri_norm + rf_norm + (1/np.sqrt(c2))*(A*(psi*c3-1))

        if (A > 0.0 and B < 0.0):
            psiL += np.pi
            B = -B
        
        chi = np.sqrt(B/c2)

        dt_computed = (1/mu_sqrt)*((chi**3)*c3 + A*np.sqrt(B))
        if i == max - 10:
            dt = dt_computed
        if abs(dt - dt_computed) < tol:
            Solution = True
            break

        if dt_computed <= dt:
            psiL = psi
        else:
            psiU = psi
        
        psi = (psiU + psiL)/2
        c2 = C2(psi)
        c3 = C3(psi)
        
    if Solution == True:
        #print("Solution Found")
        F = 1 - B/ri_norm
        G = A*np.sqrt(B/mu)
        GPrime = 1 - B/rf_norm
    
        vi = (rf - F*ri)/G
        vf = (GPrime*rf-ri)/G
        return [vi,vf,ri,[0,0,0]]
    else:
        #print("No Solution")
        return [[0,0,0],[(2**32)-1,(2**32)-1,(2**32)-1],[0,0,0]] #Maybe difference?

def LambertsProblem1(r1,r2,v1,delta_t, mu): #My own implementation
    r1_norm = norm(r1)
    r2_norm = norm(r2)
    c = r2-r1
    r1r2 = np.dot(r1, r2)

    if c[2] >= 0:
        theta = np.arccos(r1r2 / (r1_norm * r2_norm))
    else:
        theta = 2 * np.pi - np.arccos(r1r2 / (r1_norm * r2_norm))

    A = np.sqrt(r1_norm * r2_norm * (1 + np.cos(theta)))

    def func(z):
        psi = A * (z - np.sin(z)) / np.sqrt(mu)
        c2 = (1 - np.cos(np.sqrt(mu) * delta_t)) / 2
        return psi - np.sqrt(c2)

    z_initial_guess = 0.0
    z_solution = newton(func, z_initial_guess,tol = 1e-2)

    psi = A * (z_solution - np.sin(z_solution)) / np.sqrt(mu)

    f = 1 - psi / r1_norm
    g_dot = 1 - psi / r2_norm
    g = delta_t - psi / np.sqrt(mu)

    v2 = (r2 - f * r1) / g
    v2 += g_dot * v1

    return v1,v2,r1

def LambertsProblem2(r1,r2,v1,mu): #Uni Lecture
    r1mag = np.linalg.norm(r1)
    r2mag = np.linalg.norm(r2)

    c = np.linalg.norm(r2-r1)
    s = (c + r1mag + r2mag)/2
    a = (r1mag + r2mag + c)/4

    alpha = 2*np.arcsin(np.sqrt(s/(2*a)))
    beta = 2*np.arcsin(np.sqrt((s-c)/(2*a)))

    deltaT = np.sqrt((a**3)/mu)*(alpha-beta-(np.sin(alpha)-np.sin(beta)))
    deltaT_min = (np.sqrt(2)/3)*np.sqrt((s**3)/mu)*(1-((s-c)/s)**(3/2))
    #print(deltaT_min,deltaT)
    if deltaT_min > deltaT:
        return [0,0,0,0]
    
    A = np.sqrt((mu/(4*a)))*(1/np.tan(alpha/2))
    B = np.sqrt((mu/(4*a)))*(1/np.tan(beta/2))

    u1 = r1/r1mag
    u2 = r2/r2mag

    uc = (r2 - r1)/c

    v_tar = (B+A)*uc + (B-A)*u1

    return [v1,v_tar,r1]