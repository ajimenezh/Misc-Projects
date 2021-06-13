import utils
from constants import *
import numpy as np
import matplotlib.pyplot as plt

def gravity(r):
    return np.asarray([-G*M_moon/(r**2), 0.0, 0.0])

class P66Solver:
    def __init__(self):
         self.NROD = 0
         self.v_cp_x = 0.0

    def CBP(self):
        return np.asarray([[1.0, 0.0, 0.0],
                           [0.0, 1.0, 0.0],
                           [0.0, 0.0, 1.0]])
         
    ## Input
    ##  - Current velocity: vp
    ##  - Trhust acceleration measurement
    ##  - Moon gravity
    ##  - Thrust correction increment
    ##  - Current mass
    def calc(self, v_p, a_fp, g_p, delta_fa, m):
        CBP = self.CBP()
        
        ## Compute existing vertical acceleration at ROD sample instant
        a_p_x = a_fp[0] + CBP[0][0]*delta_fa / m - g_p[0]

        ## Extrapolate sample-instant measured velocity by effective transport lag
        v_p_x = v_p[0] + a_p_x*0.35

        ## Update Commanded Vertical Velocity incorporating ROD inputs
        self.v_cp_x = self.v_cp_x + self.NROD*0.3
        self.NROD = 0

        ## Compute Thrust-Acceleration Command for Throttle Routing
        a_fcp = -((v_p_x - self.v_cp_x)/1.5 + g_p[0]) / CBP[0][0]

        print (v_p[0]/2, v_p_x, g_p[0], a_fcp)
        
        ## Restrict Trhust acceleration command to produce Trhust within
        ## permitted Thrust region
        a_fcp = max(a_fcp, 0.0933*46706 / m)
        a_fcp = min(a_fcp, 0.6*46706 / m)

        return np.asarray([a_fcp, 0, 0])        


solver = P66Solver()

r_p = np.asarray([59.89, 0.0, -6.667])
v_p = np.asarray([-2.859, 0.0, 0.75644])
a_p = np.asarray([-13.05263, 0.0, 8.6782])
g_p = gravity(r_p[0] + R_moon)

x = []
z = []

delta_t = 2
i = 1

print (i*delta_t)
print (r_p)
print (v_p)
print (a_p)
print ("-----------")

propellant_mass = 8248

while r_p[0] > 0.0:
    g_p = gravity(r_p[0] + R_moon)
    
    a_fcp = solver.calc(v_p, a_p - g_p, g_p, 0.0, M_0 - propellant_mass)

    print (a_fcp)

    a_p = a_fcp + g_p
    v_p += delta_t*a_p
    r_p += delta_t*v_p

    z.append(r_p[2])
    x.append(r_p[0])

    print (i*delta_t)
    print (r_p)
    print (v_p)
    print (a_p)
    print ("-----------")

    i += 1

fig, ax = plt.subplots()
ax.plot(z, x)
#ax.set_xlim(pos_z[0], pos_z[-1] - 10e3)
#ax.set_xlabel('z (m)')
#ax.set_ylabel('x (m)')
#ax.set_title('Breaking phase target')

plt.show()
        
