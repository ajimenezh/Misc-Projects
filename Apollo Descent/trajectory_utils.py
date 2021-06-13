import numpy as np
import math
import matplotlib.pyplot as plt
from constants import *

N = 5

def derivative(order, n):
    if n == 0:
        return 1
    return order * derivative(order - 1, n - 1)

def f(t, d):
    v = []
    for i in range(d):
        v.append(0)
    for i in range(N - d):
        val = (t**i) * derivative(i+d, d)
        v.append(val)

    return np.asarray(v)

def solve(x, t):
    val = np.zeros((3))
    for i in range(N):
        val = val + (t**i)*x[i]
    return val

def solve2(x, t):
    val = np.zeros((3))
    for i in range(N-1):
        val = val + (t**i)*x[i+1]*derivative(i+1, 1)
    return val

def solve3(x, t):
    val = np.zeros((3))
    for i in range(N-2):
        val = val + (t**i)*x[i+1]*derivative(i+2, 2)
    return val

def calcP63TargetTrajectory():
    ## TODO: Check time sign
    t0 = -11.0*60 - 31.3
    tf = 0

    h0 = 15e3
    z0 = 492e3 - 4416

    # v = sqrt(Gm / r)
    m = M_moon
    vz0 = -math.sqrt(G*m / (h0 + R_moon))

    vx0 = 0

##    hf =R_moon
##    zf = 4416 + 100
##    vxf = 0
##    vzf = 0

    tf = 0
    hf = -541
    zf = 4416

    a = np.array([f(t0, 0), f(t0, 1), f(tf, 2), f(tf, 0), f(tf, 1)])
    b = np.array([h0, vx0, 0, hf, 0])

    px = np.linalg.solve(a, b)

    a = np.array([f(t0, 0), f(t0, 1), f(tf, 2), f(tf, 0), f(tf, 1)])
    b = np.array([z0, vz0, 0, zf, 0])
    
    pz = np.linalg.solve(a, b)

    v = np.zeros((N, 3))

    for i in range(N):
        v[i][0] = px[i]
        v[i][2] = pz[i]

    return v

def calcP64TargetTrajectory():
    ## TODO: Check time sign
    t0 = -177
    tf = 0

    h0 = 2231
    z0 = 7471

    pi = math.acos(-1.0)

    # v = sqrt(Gm / r)
    m = M_moon
    vz0 = -192.3038

    vx0 = -13.0157

    t3 = -31
    h3 = 30
    z3 = 11
    vxf = 0
    vzf = 0

    tf = 0
    hf = 29
    zf = 4.8

    a = np.array([f(t0, 0), f(t0, 1), f(tf, 2), f(tf, 0), f(tf, 1)])
    b = np.array([h0, vx0, 0, hf, 0])

    px = np.linalg.solve(a, b)

    a = np.array([f(t0, 0), f(t0, 1), f(tf, 2), f(tf, 0), f(tf, 1)])
    b = np.array([z0, vz0, 0, zf, 0])
    
    pz = np.linalg.solve(a, b)

    v = np.zeros((N, 3))

    for i in range(N):
        v[i][0] = px[i]
        v[i][2] = pz[i]

    return v

if __name__ == "__main__":
    v = calcP63TargetTrajectory()
    
    v = calcP64TargetTrajectory()
    px = np.zeros((N))
    pz = np.zeros((N))

    print ("r", solve(v, -31))
    print ("v", solve2(v, -31))
    print ("a", solve3(v, -31))

    for i in range(N):
        px[i] = v[i][0]
        pz[i] = v[i][2]

    print (px, pz)

    t0 = -11.0*60 - 31.3
    t0 = -177

    n = 500
    deltat = (0.0 - t0) / n
    deltat = (0.0 - t0) / n

    x = []
    z = []
    for i in range(n + 1):
        xx = solve(v, t0 + deltat*i)[0]
        zz = solve(v, t0 + deltat*i)[2]

        x.append(xx)
        z.append(zz)

    fig, ax = plt.subplots()
    ax.plot(z, x)
    ax.set_xlim(z[0], z[-1] - 50)
    ax.set_xlabel('z (m)')
    ax.set_ylabel('x (m)')
    ax.set_title('Breaking phase target')
    ax.set_title('Approach phase target')

    plt.show()
