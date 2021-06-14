import matplotlib.pyplot as plt
import math

f = open("C:\\Users\\Alex\\CMakeBuilds\\b297fea1-41c7-2032-b5eb-b1e7fa1b40fb\\build\\x64-Release\\out\\output_.txt")

n, nx, nu = [int(x) for x in f.readline()[:-1].split(' ')]

x = []
for i in range(n+1):
    x.append([float(x) for x in f.readline()[:-1].split(' ')[:-1]])

u = []
for i in range(n+1):
    u.append([float(x) for x in f.readline()[:-1].split(' ')[:-1]])

PI = math.acos(-1.0)
R0 = 6378.137e3
g0 = 9.8
m0 = 55000
Sref = (2**2) * PI
rho = 1.225

Isp = 443 / math.sqrt(R0/g0)
Tmax = 1375.6e3 / (m0*g0)
eps_max = 10 * PI / 180
eps_min = -10 * PI / 180
eps_max_rate = 5 / 180 * PI * math.sqrt(R0/g0)
tf = 25

r = [(v[0]-1)*R0 for v in x]
V = [v[2]*math.sqrt(R0*g0) for v in x]
gamma = [v[3]/PI*180 for v in x]
m = [v[4]*55000 for v in x]

T = [v[1]*(m0*g0) for v in u]
eps = [v[2]/PI*180 for v in u]

t = [1.0*tf/n*i for i in range(n+1)]

plt.ylabel('Altitude (km)')
plt.xlabel('Time (s)')
plt.plot(t, r)
plt.show()

plt.ylabel('Velocity (m/s)')
plt.xlabel('Time (s)')
plt.plot(t, V)
plt.show()

plt.ylabel('Flight-path angle (º)')
plt.xlabel('Time (s)')
plt.plot(t, gamma)
plt.show()

plt.ylabel('Mass (Kg)')
plt.xlabel('Time (s)')
plt.plot(t, m)
plt.show()

plt.ylabel('T (KN)')
plt.xlabel('Time (s)')
plt.plot(t, T)   
plt.show()

plt.ylabel('Thrust direction (º)')
plt.xlabel('Time (s)')
plt.plot(t, eps)   
plt.show()

plt.show()

