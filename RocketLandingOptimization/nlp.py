import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from sympy.integrals.quadrature import gauss_legendre

# initial guesses
n = 30
x0 = np.zeros(2*n + 2)
points  = [-1.] + gauss_legendre(n-1, 3)[0] + [1.]

def lagrange(i, t):
    res = 1.0
    t_i = points[i]
    for j in range(n+1):
        if j != i:
            t_j = points[j]
            res *= (t - t_j) / (t_i - t_j)
    return res

def d_lagrange(i, t):
    res = 0.0
    lag = lagrange(i, t)
    t_i = points[i]
    for j in range(n+1):
        if j != i:
            t_j = points[j]
            res += 1.0/(t_j - t_i)
    return res*lag

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

def objective(x):
    return -x[n]

def cons_i(x, k):
    sum_eq = 0.0
    t = points[k]
    for i in range(n+1):
        sum_eq += alternative_dl(i, t, k)*x[i]

##        if k == n-1 and i == n:
##            print (alternative_dl(i, t, k), x[i], sum_eq, 5.0/2.0*(-x[k] + x[k]*x[n+1+k] - x[n+1+k]**2))
##            print (x[i - 1], alternative_dl(i - 1, t, k), alternative_dl(i - 1, t, k)*x[i - 1])
##            print (x[i - 2], alternative_dl(i - 2, t, k), alternative_dl(i - 2, t, k)*x[i - 2])

    sum_eq -= 5.0/2.0*(-x[k] + x[k]*x[n+1+k] - x[n+1+k]**2)

    return sum_eq

def init_cons(x):
    return x[0] - 1.0

def callback(xs):
    print(xs[:n+1])

# show initial objective
print('Initial SSE Objective: ' + str(objective(x0)))

# optimize
b = (0.0,5.0)
bnds = []
for i in range(n+1):
    bnds.append(b)
for i in range(n+1):
    bnds.append((-1.0,1.0))
#con1 = {'type': 'ineq', 'fun': constraint1}
#con2 = {'type': 'eq', 'fun': constraint2}
cons_per_i = [{'type':'eq', 'fun': cons_i, 'args': (i,)} for i in np.arange(n)]
cons_per_i.append({'type': 'eq', 'fun': init_cons})
solution = minimize(objective,x0,method='SLSQP',\
                    bounds=bnds,constraints=cons_per_i,options={'disp': True}, callback=callback)

x = solution.x

# show final objective
print('Final SSE Objective: ' + str(objective(x)))

# print solution
print('Solution')
print(x)

plt.plot(x[:n+1])
plt.show()
