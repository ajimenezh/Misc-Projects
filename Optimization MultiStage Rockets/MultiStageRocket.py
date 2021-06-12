import math
from numpy import linalg as LA
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
from pyatmos import coesa76

R0 = 6371000
T0 = 264
g_0 = 9.685 / (R0/(T0**2))
mu = 3.986e14 / ((R0**3)/(T0**2))
M0 = 1.814e6 + 1.11e5

pi = math.acos(-1.0)

log = False

from pyatmos import coesa76

def Density(h):
    if (h > 100):
        return 0.0
    rhos_geom, Ts_geom, Ps_geom = coesa76([h])
    return rhos_geom[0]

def Mach(v, h):
    if (h > 100):
        return 10.0
    rhos_geom, Ts_geom, Ps_geom = coesa76([h])
    c = math.sqrt(1.4*287*Ts_geom[0])
    return v/c

def DragCoef(M):
    cd = [1.15, 0.35, 0.9, 0.5, 0.3, 0.2, 0.18, 0.18]
    m = [0.0, 0.4, 1.8, 3, 6, 8, 10, 20]
    return np.interp(M, m, cd)

def ang_to_rad(ang):
    return ang/180*pi

def rad_to_ang(ang):
    return ang/pi*180

class BoosterStage:
    def __init__(self, tau):
        self.tf = tau / T0

        self.T = 3.3362e7 / (M0*R0/(T0**2))
        self.Isp = 264 / T0
        self.beta = self.T / self.Isp / g_0

        self.mp = 1.11e5 / M0
        self.mf = self.tf * self.beta

        # 1.80e6
        tmp = 1.80e6 / M0
        
        self.k = 0.030

        self.m0 = (2.72e6) / M0

        self.r0 = 1
        self.u0 = 0
        self.w0 = math.sqrt(mu / (self.r0**3))
        self.phi = 0.0

        self.t_kick = 15 / T0

        self.A = 79.43
    
    def simulate(self, t, x):
        u = x[0]
        w = x[1]
        r = x[2]
        phi = x[3]
        m = x[4]
        lamb1 = x[5]
        lamb2 = x[6]
        lamb3 = x[7]
        lamb4 = x[8]

        #psi = pi/2
##        if t > self.t_kick:
##           #psi = self.alpha
##           T = 0
##           beta = 0.0
##           psi = math.atan2(lamb1, lamb2)
##           gravity_turn = 1
##           if psi <= self.alpha:
##               # print (t, rad_to_ang(self.alpha), rad_to_ang(psi))
##               self.tf = min(t, self.tf)
##        else:
##            T = self.T
##            beta = self.beta
##            psi = pi / 2
##            gravity_turn = 0

        gravity_turn = 0
        psi = pi/2 + (self.alpha - pi/2) * (t - 0) / (self.tf - 0)
        
        T = self.T
        beta = self.beta

        rho = Density((r*R0 - R0)/1000)
        v = math.sqrt(u*u + (r*w)**2)*(R0/T0)
        M = Mach(v, (r*R0 - R0)/1000)
        Cd = DragCoef(M)
        D = 0.5*rho*(v**2)*self.A*Cd
        
        D /= (M0*R0/(T0**2))

        ang = math.atan2(u, r*w)
        
        w_dot = 1.0/r*(-2*u*w + T/m*math.cos(psi) - D/m*math.cos(ang))

        #print (t, x)
        
        return np.asarray([(-mu/(r**2) + (w**2)*r + T/m * math.sin(psi) -
                            D/m*math.sin(ang)),
                           w_dot,
                           u,
                           w,
                           -beta,
                           gravity_turn * (2*w*lamb2 - lamb3),
                           gravity_turn * (-2*w*lamb1 + (u/r)*lamb2 - lamb4/r),
                           gravity_turn * (-(2*mu/(r**3) + w**2)*lamb1 + lamb2*(w_dot)),
                           0.0])

    def solve(self, delta_t, alpha):
        x = np.asarray([self.u0, self.w0, self.r0, self.phi, self.m0, 1.0, 0.0, 0.0, 0.0])
        self.alpha = alpha
        return self.integrate(x)

    def integrate(self, x):
        solver = scipy.integrate.RK45(self.simulate, 0.0, x, self.tf, rtol=1e-12, atol=1e-12)

        sol = []
        while solver.status == "running":
            solver.step()
            y = solver.y
            psi = pi/2 + (self.alpha - pi/2) * (solver.t - 0) / (self.tf - 0)
            y[5] = math.sin(psi)
            y[6] = math.cos(psi)
            sol.append([solver.t, y])

        if solver.status == "failed":
            print ("Solver RK45 Error")
            return []

        return sol
    
    def is_valid(self, x):
        return x[2] >= 0.8
    

class RocketMultiStage:
    def __init__(self, n_stages, tau, k, T, Isp, m0, y, booster_stage):
        self.tf = tau / T0
        
        self.T_stages = T / (M0*R0/(T0**2))
        self.Isp_stages = Isp / T0
        self.beta_stages = self.T_stages / self.Isp_stages / g_0
        self.mf = self.tf * self.beta_stages

        self.tau = tau
        self.k = k
        self.m0 = m0

        self.booster_stage = booster_stage

        #mp = 2.72e6
        #mf1 = 1.80e6
        #mf2 = 453e3
        #mf3 = 213.7e3
        #mh1 = 111e3
        #mh2 = 31750
        #mh3 = 15870
        #k1 = 0.03
        #k2 = 0.033
        #k3 = 0.12

        #T2 = 6.672e6
        #T3 = 1.12e6

        #Isp2 = 428
        #Isp3 = 850

        # tau = mf / beta
        # beta = T / Isp / g_0

        #beta2 = 1590.69
        #beta3 = 134.45
        
        #tau2 = 284.78
        #tau3 = 1589.43
        
        self.r0 = y[2]
        self.u0 = y[0]
        self.w0 = y[1]
        self.phi = y[3]
        
        self.m0 = y[4]

        self.n_stages = n_stages
    
    def simulate(self, t, x):
        u = x[0]
        w = x[1]
        r = x[2]
        phi = x[3]
        m = x[4]
        lamb1 = x[5]
        lamb2 = x[6]
        lamb3 = x[7]
        lamb4 = x[8]

        #psi = pi / 2.0*0.98
        psi = math.atan2(lamb1, lamb2)
        T = self.T
        beta = self.beta
        
        w_dot = 1.0/r*(-2*u*w + T/m*math.cos(psi))
        
        return np.asarray([-mu/(r**2) + (w**2)*r + T/m * math.sin(psi),
                           w_dot,
                           u,
                           w,
                           -beta,
                           2*w*lamb2 - lamb3,
                           -2*w*lamb1 + (u/r)*lamb2 - lamb4/r,
                           -(2*mu/(r**3) + w**2)*lamb1 + lamb2*(w_dot),
                           0.0])

    def solve2(self, lamb, delta_t):
        x = np.asarray([self.u0, self.w0, self.r0, self.phi, self.m0,
                        lamb[0], lamb[1], 0.0, 0.0])
        return rungeKutta(x, 0.0, self.tf, delta_t, self.simulate, self.is_valid)

    def solve(self, y, t0):
        lamb0 = y[0]
        psi = y[1]
        psi_dot = y[2]
        
        lamb1 = lamb0*math.sin(psi)
        lamb2 = lamb0*math.cos(psi)
        lamb3 = lamb0/math.cos(psi)*(2*self.w0 - psi_dot -
                                           self.u0/self.r0*math.sin(psi)*math.cos(psi))
        #lamb3 = -lamb3

        print ("-->", psi, lamb0, lamb1, lamb2, math.atan2(lamb1, lamb2))
        
        x = np.asarray([self.u0, self.w0, self.r0, self.phi,
                        self.m0 - self.booster_stage.k*self.booster_stage.mf,
                        lamb1, lamb2, lamb3, 0.0])
        res = []
        total_time = t0
        print (t0)
        for i in range(self.n_stages):
            self.T = self.T_stages[i]
            self.beta = self.beta_stages[i]

            if log:
                print ("Initial solution for stage: ", i, " is ", x, ", at time ", total_time*T0)
            solver = scipy.integrate.RK45(self.simulate, 0.0, x, self.tf[i],
                                          rtol=1e-8, atol=1e-8)
            sol = [[total_time, x]]
            while solver.status != "finished":
                solver.step()
                sol.append([solver.t + total_time, solver.y])
            res.append(sol)
            x = sol[-1][1]
            x[4] -= self.mf[i] * self.k[i]
            total_time += self.tf[i]
        return res
    
    def is_valid(self, x):
        return x[2] >= 0.8

    def cost(self, sol):
        x0 = []
        xf = []
        for i in range(self.n_stages + 1):
            x0.append(sol[i][0][1])
            xf.append(sol[i][-1][1])

        Sf = []
        S0 = []
        beta = self.beta_stages
        beta = np.insert(beta, 0, self.booster_stage.beta)
        T = self.T_stages
        T = np.insert(T, 0, self.booster_stage.T)
        k = self.k
        k = np.insert(k, 0, self.booster_stage.k)
        for i in range(self.n_stages + 1):
            u0 = x0[i][0]
            uf = xf[i][0]
            w0 = x0[i][1]
            wf = xf[i][1]
            r0 = x0[i][2]
            rf = xf[i][2]
            m0 = x0[i][4]
            mf = xf[i][4]
            l10 = x0[i][5] if i != 0 else x0[i+1][5]
            l1f = xf[i][5] if i != 0 else xf[i+1][5]
            l20 = x0[i][6] if i != 0 else xf[i+1][6]
            l2f = xf[i][6] if i != 0 else xf[i+1][6]
            l30 = x0[i][7] if i != 0 else xf[i+1][7]
            l3f = xf[i][7] if i != 0 else xf[i+1][7]
            l40 = x0[i][8] if i != 0 else xf[i+1][8]
            l4f = xf[i][8] if i != 0 else xf[i+1][8]

            Sf_i = ((mu/(rf**2) - (wf**2)*rf)*l1f + 2*uf*wf*l2f -
                     uf*l3f - wf*l4f -
                     T[i]/mf*math.sqrt((l1f*l1f)+(l2f*l2f))) / beta[i]
            S0_i = ((mu/(r0**2) - (w0**2)*r0)*l10 + 2*u0*w0*l20 -
                     u0*l30 - w0*l40 -
                     T[i]/m0*math.sqrt((l10*l10)+(l20*l20))) / beta[i]

            Sf.append(Sf_i)
            S0.append(S0_i)

        cost = np.zeros((self.n_stages + 3 + 1))

        cost[0] = sol[-1][-1][1][0]
        print (sol[-1][-1][1][0], sol[-1][-1][1][0]*(R0/T0))

        rd = (R0 + 860e3)/R0
        cost[1] = (sol[-1][-1][1][2] - rd)
        cost[2] = sol[-1][-1][1][1] - math.sqrt(mu/(rd**3))
        
        #cost[1] = sol[-1][-1][1][7]
        #cost[2] = sol[-1][-1][1][6]

        #gamma = -self.k[-1] + ((1 + self.k[-1])*
        #                       math.exp(self.delta_v*self.beta[-1]/self.tau[-1]))
        gamma = 1.0
        for i in range(self.n_stages):
            G_tau_i = 0.0
            for j in range(i+1, self.n_stages + 1, 1):
                G_tau_i += (Sf[j] - S0[j])
                
            G_tau_i = (1 + k[i])*beta[i]*(gamma + G_tau_i)
            G_tau_i += beta[i]*Sf[i]

            cost[3 + i] = G_tau_i
            
        cost[3 + self.n_stages] = (gamma + k[-1])*beta[-1] + beta[-1]*Sf[-1]

        print ("Cost ", cost)

        return cost
                    

def rungeKutta(x0, t0, tf, delta_t, f, is_valid):
    n = int((tf - t0) / delta_t)
    t = t0
    t_real = t0
    x = x0

    res = [[t, x]]
 
    for i in range(n):
        k1 = delta_t*f(t, x)
        k2 = delta_t*f(t + 0.5*delta_t, x + 0.5*k1)
        k3 = delta_t*f(t + 0.5*delta_t, x + 0.5*k2)
        k4 = delta_t*f(t + delta_t, x + k3)

        x = x + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        t = t + delta_t
        t_real = t_real + delta_t
        if not is_valid(x):
            return []

        res.append([t_real, x])
 
    return res

def NewtonRhapson(y0, f, n, m):
    y = y0
    J = np.zeros((n, m))
    delta = 1.0e-1
    for i in range(30):
        f_res = f(y)
        sum = 0
        for x in f_res:
            sum += abs(x)
        print ("Cost ", f_res, sum) 
        for j in range(m):
            y_tmp = y.copy()
            delta = 1e-6
            y_tmp[j] += delta
            f_tmp = f(y_tmp) - f_res
            for k in range(n):
                J[k][j] = f_tmp[k] / delta
        
        #y = y - np.matmul(np.linalg.inv(J), f_res)
        G = np.matmul(np.linalg.inv(np.matmul(np.transpose(J), J)),
                      np.transpose(J))

        try:
            step = 0.1
            while True:
                y_new = y - np.matmul(G, f_res)*step
                print ("New step ", y_new)
                
                if np.sum(abs(f(y_new))) - np.sum(abs(f(y))) > 0:
                    step /= 2
                else:
                    if np.linalg.norm(y_new - y) < 1.0e-6:
                        return y_new

                    #print (i, step, f(y_new))
                    #print (J, np.matmul(np.linalg.inv(J), f_res))
                    y = y_new
                    
                    break
        except np.linalg.LinAlgError:
            print (J)
            return y

    return y

rd = 6.4e6 / R0
wd = math.sqrt(mu / (rd**3))
    
def optimize(y, return_output=False):
    solver = Rocket1Stage()
    s = solver.solve2(y, 1e-4)

    sol = []
    while s.status != "finished":
        s.step()
        sol.append([s.t, s.y])

    if len(sol) == 0:
        return None

    if return_output:
        return sol
    else:
        x = sol[-1][1]
        return np.asarray([(x[0] - 0.0),
                           (x[7] - 0.0),
                           (x[6] - 0.0)])

def booster_stage(y, return_output=False):
    solver = BoosterStage()
    s = solver.solve(y, 1e-1, ang_to_rad(89.5), [sosolvelver.tf])

    sol = []
    while s.status != "finished":
        s.step()
        sol.append([s.t, s.y])

    if len(sol) == 0:
        return None

    if return_output:
        return sol
    else:
        x = sol[-1][1]
        return np.asarray([(x[0] - 0.0),
                           (x[7] - 0.0),
                           (x[6] - 0.0)])

def multi_stage(y, n_stages = 2, return_output=False, lamb0 = 0.4):
    # Initial conditions
    #y = np.insert(y, 0, lamb0)
    lambda0 = y[0]
    alpha = y[2]
    psi = ang_to_rad(y[1])
    tau1 = y[3]
    #psi_dot = y[1]
    tau = []
    for i in range(n_stages):
        tau.append(y[4+i])

    #print (y)

    if log:
        print ("Calc booster stage")
    booster_stage = BoosterStage(tau1)
    
    booster_sol = booster_stage.solve(1e-1, ang_to_rad(alpha))

    if log:
        print ("Calc booster stage derivatives")
    delta_alpha = 1e-1
    booster_stage = BoosterStage(tau1)
    booster_sol_1 = booster_stage.solve(1e-1, ang_to_rad(alpha) + delta_alpha)

    booster_stage = BoosterStage(tau1)
    booster_sol_2 = booster_stage.solve(1e-1, ang_to_rad(alpha) - delta_alpha)

    #print (booster_sol_1[-1], booster_sol_2[-1])

    du_da = (booster_sol_1[-1][1][0] - booster_sol_2[-1][1][0]) / (2*delta_alpha)
    dw_da = (booster_sol_1[-1][1][1] - booster_sol_2[-1][1][1]) / (2*delta_alpha)
    dr_da = (booster_sol_1[-1][1][2] - booster_sol_2[-1][1][2]) / (2*delta_alpha)
    dphi_da = (booster_sol_1[-1][1][3] - booster_sol_2[-1][1][3]) / (2*delta_alpha)

    w0 = booster_sol[-1][1][1]
    u0 = booster_sol[-1][1][0]
    r0 = booster_sol[-1][1][2]

    psi_dot = (2*w0 - u0/r0*math.sin(psi)*math.cos(psi) +
               math.cos(psi)/dr_da*(math.sin(psi)*du_da + r0*math.cos(psi)*dw_da))
    #print (psi_dot)

    if log:
        print ("Initial value for psi dot: ", psi_dot)

    if log:
        print ("Calc multi stage")
    multi_stage_solver = RocketMultiStage(n_stages = 2, tau = np.asarray(tau),
                                          k = np.asarray([0.033, 0.12]),
                                          T = np.asarray([6.672e6, 2.12e6]),
                                          Isp = np.asarray([428, 850]),
                                          m0 = 31.75e3, y = booster_sol[-1][1],
                                          booster_stage = booster_stage)

    #psi_dot = 0.0
    multi_stage_sol = multi_stage_solver.solve([lambda0, psi, psi_dot],
                                               booster_stage.tf)

    if log:
        print (multi_stage_solver.cost(multi_stage_sol))

    if return_output:
        sol = []
        for y in booster_sol:
            sol.append(y)
        for temp in multi_stage_sol:
            for y in temp:
                sol.append(y)
        return sol
    else:
        multi_stage_sol.insert(0, booster_sol)
        #print (multi_stage_solver.cost(multi_stage_sol))
        return multi_stage_solver.cost(multi_stage_sol)

def booster_stage_test(y, lamb0 = 0.4):
    # Initial conditions
    #y = np.insert(y, 0, lamb0)
    lambda0 = y[0]
    alpha = y[1]
    tau1 = y[2]
    tau = []

    #print (y)

    if log:
        print ("Calc booster stage")
    booster_stage = BoosterStage(tau1)

    booster_stage.t_kick = booster_stage.tf
    booster_stage.tf = 1.5*booster_stage.tf
    
    booster_sol = booster_stage.solve(1e-1, ang_to_rad(alpha))

    if log:
        print (multi_stage_solver.cost(multi_stage_sol))

    sol = []
    for y in booster_sol:
        sol.append(y)
    return sol

best_res = 1.e8
best_sol = []

def callback(x, f):
    print (x)
    print (f)

    sum = 0.0
    for t in f:
        sum += t

    global best_res
    global best_sol

    if sum < best_res:
        best_sol = x
        best_res = sum

pi = math.acos(-1.0)

##for i in range(10):
##    x = 181.5 + i*0.1
##    y = np.asarray([0.0077, 92.9,  92, 134.955,
##        182.1,  600])
##    print (x, multi_stage(y))
##sol = optimize(lamb_, True)
##
from scipy.optimize import fsolve
###y = fsolve(optimize, y0, maxfev=150)
###sol = optimize(y, True)
##
from scipy.optimize import newton_krylov
from numpy import cosh, zeros_like, mgrid, zeros
##
###y = newton_krylov(optimize, y0, method='lgmres', verbose=1)    
###sol = optimize(y, True)
##
###y = NewtonRhapson(y0, optimize, 3, 3)
####for i in range(20):
####    tmp = 260 + i*50.0/20
####    print (tmp, multi_stage(
####        np.asarray([pi/2, 89.5, 138, 278, 690]),
####        return_output = False,
####        lamb0 = 0.4))
####sol = multi_stage(np.asarray([pi/2, 89.5, 138, 284, 550]), return_output = True)
##
##x = fsolve(multi_stage, np.asarray([0.2577, 8.8e+01,  138,
##        172,  700]),
##           maxfev=150, full_output=True)
####y = x[0]
##
##y = NewtonRhapson(np.asarray([-6.55e-03,  1.62244414e+00,  8.6661365e+01,
##                              1.40927147e+02, 2.09221506e+02,  7.79295242e+02]), multi_stage, 6, 6)
####-5.24905745e-03  1.43954485e+00  8.66614285e+01  1.41443196e+02
####  2.04166277e+02  7.98903252e+02
#### np.asarray([0.4, pi/2, 89.5, 138, 278, 690])
##
#y = newton_krylov(multi_stage, np.asarray([0.2577, 8.8e+01,  138,
#        172,  700]),
#                  method='lgmres', verbose=1, callback = callback) 
##
y = np.asarray([0.0077, 92.9,  92, 134.955,
        182.1,  600])
#y = np.asarray([-6.86e-2, 8.785e1, 1.4244e2, 2.326e2, 8.398e2])
##sol = best_sol
print (multi_stage(y))
sol = multi_stage(y, return_output = True)
##print (multi_stage(y, return_output = False))
####print (best_sol, best_ans)

#y = np.asarray([0.2577, 8.8e+01,  138,
#        172,  700])
#sol = booster_stage_test(y)

t = []
h = []
v = []
w = []
w2 = []
l1 = []
l2 = []
m = []
phi = []
psi = []

for x in sol:
    t.append(x[0])
    v.append(x[1][0]*(R0/T0))
    h.append((x[1][2]*R0 - R0)/1000)
    w2.append(math.sqrt(mu / (x[1][2]**3)))
    w.append(x[1][1]/T0)
    m.append(x[1][4]*M0)
##    l1.append(x[1][5])
##    l2.append(x[1][6])
    phi.append(x[1][3] / pi * 180)
    psi.append(math.atan2(x[1][5], x[1][6]) / pi * 180)

##plt.plot(t, l1)
##plt.plot(t, l2)
##plt.show()
#plt.plot(t, h)
plt.title('Altitude')
plt.ylabel('h (km)')
plt.xlabel('t (s)')
plt.plot(t, h)
plt.show()
plt.title('Angular velocity')
plt.ylabel('w (rad/s)')
plt.xlabel('t (s)')
plt.plot(t, w)
#plt.plot(t, w2)
plt.show()
plt.title('Velocity')
plt.ylabel('v (m/s)')
plt.xlabel('t (s)')
plt.plot(t, v)
plt.show()
plt.title('Polar angle')
plt.ylabel('phi (rad)')
plt.xlabel('t (s)')
plt.plot(t, phi)
plt.show()
plt.title('Thrust profile')
plt.ylabel('psi (rad)')
plt.xlabel('t (s)')
plt.plot(t, psi)
plt.show()
plt.title('Mass')
plt.ylabel('m (kg)')
plt.xlabel('t (s)')
plt.plot(t, m)
plt.show()
