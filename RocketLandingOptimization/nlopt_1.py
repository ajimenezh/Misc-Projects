import nlopt
from numpy import *
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from sympy.integrals.quadrature import gauss_legendre

PI = math.acos(-1.0)
R0 = 6378.137e3
g0 = 9.8
m0 = 55000
Sref = (2**2) * PI
rho = 1.225

r0 = 3.e3 / R0
V0 = 280 / math.sqrt(g0*R0)
s0 = 0.0
gamma0 = -65 * PI / 180

rf = 1.0
sf = 0.0
Vf = 10 / math.sqrt(g0*R0)
gammaf = -PI/2

Isp = 443 / math.sqrt(R0/g0)
Tmax = 1375.6e3 / (m0*g0)
eps_max = 10 * PI / 180
eps_min = -10 * PI / 180
eps_max_rate = 5 / 180 * PI * math.sqrt(R0/g0)

n = 100
nx = 5
nu = 3
points  = [-1.] + gauss_legendre(n-1, 3)[0] + [1.]

t0 = 0
tf = 25 / math.sqrt(R0/g0)

dl_coef = {}

def alternative_dl(j,t,i):
    if (j,i) in dl_coef:
        return dl_coef[(j, i)]
    y = 0
    t_j = points[j]
    for l in range(n+1):
        if not l==j:
            t_l = points[l]
            k = 1/(t_j-t_l);
            for m in range(n+1):
                if not(m==j) and not(m==l):
                    t_m = points[m]
                    k = k*(t-t_m)/(t_j-t_m)
            y = y + k
    dl_coef[(j, i)] = y
    return y

def objective(x, grad):
    return -x[(n+1)*nx - 1]

it = 0

def constraint(result, x, grad, k, sgn): 
    result[:] = np.zeros(nx)
    t = points[k]
    for idx in range(nx):
        for i in range(n+1):
            result[idx] += np.float64(alternative_dl(i, t, k))*x[i*nx + idx]

    r = x[k*nx]
    V = x[k*nx+2]
    gamma = x[k*nx+3]
    m = x[k*nx+4]
    L0 = 0.5 * rho * V * V * Sref * g0 * R0 / (m0*g0)
    L = x[(n+1)*nx + k*nu]*L0
    D = 0.3*L
    T = x[(n+1)*nx + k*nu + 1]*Tmax
    eps = x[(n+1)*nx + k*nu + 2]
    
    result[0] -= (tf - t0)/2.0*(V*math.sin(gamma))
    result[1] -= (tf - t0)/2.0*(V*math.cos(gamma))
    result[2] -= (tf - t0)/2.0*((-T*math.cos(eps) - D)/m - math.sin(gamma)/(r*r))
    result[3] -= (tf - t0)/2.0*((-T*math.cos(eps) + L)/(m*V) - math.cos(gamma)/(r*r*V))
    result[4] -= (tf - t0)/2.0*(-T/Isp)                       

    return sgn*result

def boundary(result, x, grad, sgn):
    result[:] = np.zeros(2*nx - 1)

    r = x[0]
    s = x[1]
    V = x[2]
    gamma = x[3]
    m = x[4]

    result[0] = r - (1. + r0)
    result[1] = s - s0
    result[2] = V - V0
    result[3] = gamma - gamma0
    result[4] = m - 1.0

    r = x[n*nx]
    s = x[n*nx + 1]
    V = x[n*nx + 2]
    gamma = x[n*nx + 3]

    result[nx] = r - rf
    result[nx+1] = s - sf
    result[nx+2] = V - Vf
    result[nx+3] = gamma - gammaf

    return sgn*result

#opt = nlopt.opt(nlopt.GN_ISRES, 2*n+2)
opt = nlopt.opt(nlopt.LN_COBYLA, (n+1)*(nx+nu))

lower_bounds = []
for i in range(n+1):
    lower_bounds.append(np.float64(1.0))
    lower_bounds.append(np.float64(-2.0))
    lower_bounds.append(np.float64(0.0))
    lower_bounds.append(np.float64(-PI/2))
    lower_bounds.append(np.float64(0.5))
for i in range(n+1):
    lower_bounds.append(np.float64(0.0))
    lower_bounds.append(np.float64(0.2))
    lower_bounds.append(np.float64(eps_min))
    
opt.set_lower_bounds(lower_bounds)

upper_bounds = []
for i in range(n+1):
    upper_bounds.append(np.float64(2.0))
    upper_bounds.append(np.float64(2.0))
    upper_bounds.append(np.float64(2.0))
    upper_bounds.append(np.float64(PI/2))
    upper_bounds.append(np.float64(1.0))
for i in range(n+1):
    upper_bounds.append(np.float64(1.0))
    upper_bounds.append(np.float64(1.0))
    upper_bounds.append(np.float64(eps_max))
    
opt.set_upper_bounds(upper_bounds)

opt.set_min_objective(objective)

tol = [1e-3, 1e-2, 1e-3, 1e-3, 1e-5]

for i in range(n):
    opt.add_inequality_mconstraint(lambda res,x,grad,idx=i: \
                                  constraint(res,x,grad,idx, 1), tol)
    opt.add_inequality_mconstraint(lambda res,x,grad,idx=i: \
                                  constraint(res,x,grad,idx, -1), tol)

tol_2 = [1e-3, 1e-2, 1e-3, 1e-3, 1e-4, 1e-3, 1e-2, 1e-3, 1e-3]

opt.add_inequality_mconstraint(lambda res,x,grad: boundary(res,x,grad,1), tol_2)
opt.add_inequality_mconstraint(lambda res,x,grad: boundary(res,x,grad,-1), tol_2)

opt.set_xtol_rel(1e-4)

in_sol = []
for i in range((n+1)):
    r = 1.0 + r0*(n-i)/n
    s = 0.0
    V = Vf + (V0-Vf)*(n-i)/n
    gamma = gammaf + (gamma0-gammaf)*(n-i)/n
    m = 1.0
    
    in_sol.append(np.float64(r))
    in_sol.append(np.float64(s))
    in_sol.append(np.float64(V))
    in_sol.append(np.float64(gamma))
    in_sol.append(np.float64(m))
for i in range((n+1)):
    in_sol.append(np.float64(0.5))
    in_sol.append(np.float64(0.5))
    in_sol.append(np.float64(0.0))

opt.set_maxtime(60*5)

x = opt.optimize(in_sol)
minf = opt.last_optimum_value()

np.savetxt('data.csv', x, delimiter=',')

print("optimum at ", x)
print("minimum value = ", minf)
print("result code = ", opt.last_optimize_result())
