import math
from numpy import linalg as LA
import numpy as np
import matplotlib.pyplot as plt

g_0 = 9.685
mu = 3.986e14

class KeplerianOrbitalElements:
    def __init__(self, r, v):
        h = np.cross(r, v)
        k = np.asarray([0, 0, 1])
        n = np.cross(k, h)
        self.e = ((norm(v)**2 - mu/norm(r))*r - np.dot(r, v)*v) / mu
        self.ee = norm(self.e)
        if self.ee == 1.0:
            self.p = norm(h)**2 / mu
            self.a = p
        else:
            self.E = norm(v)**2/2 - mu/norm(r)
            self.a = -mu / (2*self.E)
            self.p = self.a*(1 - self.ee**2)

        self.i = math.acos(h[2] / norm(h))

def flight_path_angle(r, v):
    return math.asin(np.dot(r, v) / (norm(r) * norm(v)))

def ang_to_rad(x):
    return x / 180 * math.acos(-1.0)

def norm(x):
    return LA.norm(x)

def g(r):
    return - mu/(norm(r)**3)*r

t_m = 337.15 - 155

def T(t):
    if t < t_m:
        return 17437041
    else: 
        return 1306665

def Isp(t):
    if t < t_m:
        return 414.7
    else:
        return 448

def u_T(lag_mul_0, x, t):
    w_0 = math.sqrt(mu/(norm(x[0])**3))

    lag_mul_r = (lag_mul_0[0]*math.cos(w_0*t) +
                 lag_mul_0[1]*w_0*math.sin(w_0*t))
    lag_mul_v = (-(lag_mul_0[0]/w_0)*math.sin(w_0*t) +
                 lag_mul_0[1]*math.cos(w_0*t))

    return lag_mul_v / norm(lag_mul_v)

def LagMul(lag_mul_0, r, t):
    w_0 = math.sqrt(mu/(norm(r)**3))

    lag_mul_r = (lag_mul_0[0]*math.cos(w_0*t) +
                 lag_mul_0[1]*w_0*math.sin(w_0*t))
    lag_mul_v = (-(lag_mul_0[0]/w_0)*math.sin(w_0*t) +
                 lag_mul_0[1]*math.cos(w_0*t))

    return np.concatenate((lag_mul_r, lag_mul_v))

def simulate(x, t, t_real, lag_mul_0):
    return np.asarray([x[1],
                       g(x[0]) + T(t_real)/x[2]*u_T(lag_mul_0, x, t),
                       -T(t_real)/(g_0*Isp(t_real))])

def rungeKutta(x0, t0, tf, delta_t, f, lag_mul_0, lag_mul_f):
    if tf > 1000.0:
        tf = 1000.0
    elif tf < 300.0:
        tf = 300.0
    
    n = int((tf - t0) / delta_t)
    t = t0
    t_real = t0
    x = x0

    res = [[t, x]]

    updated = False
 
    for i in range(n):
        k1 = delta_t*f(x, t, t_real, lag_mul_0)
        k2 = delta_t*f(x + 0.5*k1, t + 0.5*delta_t,
                       t_real + 0.5*delta_t, lag_mul_0)
        k3 = delta_t*f(x + 0.5*k2, t + 0.5*delta_t,
                       t_real + 0.5*delta_t, lag_mul_0)
        k4 = delta_t*f(x + k3, t + delta_t,
                       t_real + delta_t, lag_mul_0)

        x = x + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        t = t + delta_t
        t_real = t_real + delta_t

        if t > t_m and not updated:
            updated = True
            x[2] = 307837
            lag_mul__new = lag_mul_0
            lag_mul__new[1] = LagMul(lag_mul_0, x[0], t)[1]
            lag_mul__new[0] = lag_mul_f
            #lag_mul_0 = lag_mul__new
            #t = 0

        res.append([t_real, x])
 
    return res

def NewtonRhapson(y0, f, delta, n, m):
    y = y0
    J = np.zeros((n, m))
    for i in range(10):
        f_res = f(y)
        for j in range(m):
            y_tmp = y
            y_tmp[j] += delta[j]
            f_tmp = f(y_tmp) - f_res
            for k in range(n):
                J[k][j] = f_tmp[k] / delta[j]
        
        #y = y - np.matmul(np.linalg.inv(J), f_res)
        G = np.matmul(np.linalg.inv(np.matmul(np.transpose(J), J)),
                      np.transpose(J))
        y = y - np.matmul(G, f_res)

        print (y)

        if y[-1] > 1000.0:
            y[-1] = 1000.0
        if y[-1] < 400.0:
            y[-1] = 400.0

        print (f(y))

    return y

def solve_red(y, return_output=False):
    r0 = np.asarray([6.4e6, 1.164e4, 4435])
    v0 = np.asarray([1592, 180, 284])

    m0 = 1227734

    x0 = np.asarray([r0, v0, m0])

    l = v0 / norm(v0)
    lag_mul_0 = np.asarray([[y[0], y[1], y[2]], [l[0], l[1], l[2]]])
    lag_mul_f = np.asarray([y[0], y[1], y[2]])

    h = norm(np.concatenate(lag_mul_0))
    lag_mul_0 = lag_mul_0 / h

    t0 = 0.0

    delta_t = 1.0

    sol = rungeKutta(x0, t0, y[3], delta_t, simulate, lag_mul_0, lag_mul_f)

    if return_output:
        return sol
    else:
        x = sol[-1][1]
        p = KeplerianOrbitalElements(x[0], x[1])
        fpa = flight_path_angle(x[0], x[1])
        #print (y)
        return np.asarray([(norm(x[0]) - 6518889),
                           (norm(x[1]) - 7840.675),
                           p.i - ang_to_rad(28.5),
                           fpa - ang_to_rad(1.0064)])

def solve(y, return_output=False):
    r0 = np.asarray([-5255758, 2044694, 3058330])

    ##gamma = 10.6/180*math.acos(-1.0)
    ##k = np.asarray([0.0, 0.0, 7.3e-5])
    ##beta = 0.0
    ##v2 = ((np.cross(r0, np.cross(k, r0))*math.cos(beta) +
    ##      np.cross(k, r0)*math.sin(beta))*math.cos(gamma) + r0*math.sin(gamma))
    ##print (v2 / norm(v2))

    v0 = np.asarray([413.93, -1839.6, 202.55])

    r0 = np.asarray([6.4e6, 1.164e4, 4435])
    v0 = np.asarray([1592, 180, 284])

    m0 = 1227734

    x0 = np.asarray([r0, v0, m0])

    lag_mul_0 = np.asarray([[y[0], y[1], y[2]], [y[3], y[4], y[5]]])
    lag_mul_f = np.asarray([y[6], y[7], y[8]])

    h = norm(np.concatenate(lag_mul_0))
    lag_mul_0 = lag_mul_0 / h

    t0 = 0.0

    delta_t = 1.0

    sol = rungeKutta(x0, t0, y[9], delta_t, simulate, lag_mul_0, lag_mul_f)

    if return_output:
        return sol
    else:
        x = sol[-1][1]
        p = KeplerianOrbitalElements(x[0], x[1])
        fpa = flight_path_angle(x[0], x[1])
        #print (y)
        return np.asarray([(norm(x[0]) - 6518889) / 6518889,
                           (norm(x[1]) - 7840.675) / 7840.675,
                           p.i - ang_to_rad(28.5),
                           fpa - ang_to_rad(1.0064)])
        return -np.asarray([(norm(x[0]) - 6518889) / 6518889,
                           (norm(x[1]) - 7840.675) / 7840.675,
                           p.i - ang_to_rad(28.5),
                           fpa - ang_to_rad(1.0064),
                           p.ee - 1.1,
                           (p.a - 13237778)/13237778,
                           (norm(LagMul(lag_mul_0, x[0], y[-1])[0]) - 1.0) / norm(LagMul(lag_mul_0, x[0], y[-1])[0]),
                           (norm(LagMul(lag_mul_0, x[0], y[-1])[1]) - 1.0) / norm(LagMul(lag_mul_0, x[0], y[-1])[1]),
                            (norm(LagMul(lag_mul_0, x[0], y[0])[0]) - 1.0) / norm(LagMul(lag_mul_0, x[0], y[0])[0]),
                            (norm(LagMul(lag_mul_0, x[0], y[0])[1]) - 1.0) / norm(LagMul(lag_mul_0, x[0], y[0])[1])])

def solve_2():
    m = 1e15
    result = []
    for i in range(500):
        y = np.zeros(10)
        for j in range(9):
            y[j] = np.random.uniform(0.0,1.0)
        y[9] = 400 + np.random.uniform(0.0,400.0)

        res = solve(y, False)
        err = math.sqrt(res[0]**2 + res[1]**2)
        if err < m:
            m = err
            print ("---->", res)
            print (y, res)
            result = y
    return result

def solve_3():
    r0 = np.asarray([-5255758, 2044694, 3058330])
    v0 = np.asarray([413.93, -1839.6, 202.55])
    
    r0 = np.asarray([6.4e6, 1.164e4, 4435])
    v0 = np.asarray([1592, 180, 284])
    
    l = v0 / norm(v0)
    y0 = np.asarray([1e-4, 1e-4, 1e-4, l[0], l[1], l[2], 1e-4, 1e-4, 1e-4, 600.0])
    delta = np.asarray([1e-5, 1e-5, 1e-5, 1e-0])

    y = y0

    y0 = np.asarray([1e-4, 1.0, 1e-4, 600.0])
    #y = NewtonRhapson(y0, solve_red, delta, 4, 4)

    import scipy.optimize

    def cb(x, f):
        print (x, f)

    # y0 = np.asarray([5.27923474e+01, -7.73258816e+00, -4.35254007e-03, 2.23521737e-01,
    #                -9.69140171e-01, 3.40878432e+03, 4.50308110e+02])

    #y = scipy.optimize.broyden1(solve_red, y0, verbose=True, callback=cb)

    ##y = np.asarray([5.27923474e+01, -7.73258816e+00, -4.35254007e-03, 2.23521737e-01,
    ##                -9.69140171e-01, 3.40878432e+03, 3.50308110e+02])

    from scipy.optimize import newton_krylov
    from numpy import cosh, zeros_like, mgrid, zeros

    y = newton_krylov(solve_red, y0, method='lgmres', verbose=1)

    from scipy.optimize import fsolve

    y0 = np.asarray([1e-4, 1e-4, 1e-4, 600.0])
    #y = fsolve(solve_red, y0, maxfev=150)

    from scipy.optimize import least_squares

    bounds =([0, 0.0, 0.0, 400.0],
             [1, 1.0, 1.0, 900.0])
    
    y0 = np.asarray([1e-4, 1e-4, 1e-4, 600.0])
    #res = least_squares(solve_red, y0, verbose=2, max_nfev=200, bounds=bounds)
    #y = res.x

##    y = [ -9462.60339886,   358.49624092, -2511.6145781,   9994.97059323,
##   158.21711496,   452.79769984,  -148.25226318 , -148.25226318,
##  -299.21317447 ,  300.78672553]

    print (y)
    return y

y = solve_3()
sol = solve_red(y, True)   
                       
t = []
h = []
v = []

for x in sol:
    t.append(x[0])
    h.append(norm(x[1][0]))
    v.append(norm(x[1][1]))

#plt.plot(t, h)
plt.plot(t, h)
plt.show()
plt.plot(t, v)
plt.show()

