import utils
from constants import *
from trajectory_utils import *
import numpy as np
import matplotlib.pyplot as plt

T0 = -11.0*60 - 31.3

class P63Solver:
    def __init__(self, l_p, P63_target_trajectory):
        # Landing site
        self.l_p = l_p
        self.t_old = T0
        self.T = T0

        self.P63_target_trajectory = P63_target_trajectory

    def CGP(self, v):
        x_g = utils.toUnit(self.l_p)
        z_g = utils.toUnit(np.asarray([-x_g[2], self.l_p[1], x_g[0]]))
        y_g = utils.toUnit(np.cross(z_g, x_g))
        c_gp = np.asarray([x_g, y_g, z_g])
        return np.matmul(c_gp, v)

    def getTarget(self, t, x):
        v = np.zeros((5, 3))
        for i in range(5):
            for j in range(5-i):
                v[i] = v[i] + x[j+i]*(t**j)*derivative(j+i, i)
        return v

    ## Input
    ##  - Current state: rp, vp
    ##  - Clock time current state: t
    ##  - Current gravity: gp
    def calc(self, r_p, v_p, t, g_p):

        # Update landing site vector for lunar rotation
        self.l_p = utils.norm(self.l_p) * utils.toUnit(
            self.l_p + (t - self.t_old)*np.cross(w_moon_p, self.l_p))

        # Compute current state vector in guidance coordinates
        #print ("pos", r_p, self.l_p, r_p - self.l_p, self.CGP(r_p - self.l_p))
        r_g = self.CGP(r_p - self.l_p)

        v_g = self.CGP(v_p - np.cross(w_moon_p, r_p))

        # Update target referenced time t
        self.T = self.T + (t - self.t_old)
        self.t_old = t

        delta_crit = self.T / 128

        # Target values
        s_tg = self.P63_target_trajectory[4]
        j_tg = self.P63_target_trajectory[3]
        a_tg = self.P63_target_trajectory[2]
        v_tg = self.P63_target_trajectory[1]
        r_tg = self.P63_target_trajectory[0]

        s_tg_z = s_tg[2]*24
        j_tg_z = j_tg[2]*6
        a_tg_z = a_tg[2]*2
        v_tg_z = v_tg[2]
        r_tg_z = r_tg[2]

        v_g_z = v_g[2]
        r_g_z = r_g[2]

##        target = self.getTarget(self.T, self.P63_target_trajectory)
##        s_tg = target[4]
##        j_tg = target[3]*6
##        a_tg = target[2]*2
##        v_tg = target[1]
##        r_tg = target[0]
##
##        s_tg_z = s_tg[2]
##        j_tg_z = j_tg[2]
##        a_tg_z = a_tg[2]
##        v_tg_z = v_tg[2]
##        r_tg_z = r_tg[2]
##
##        print target

        #print ("start", self.T, self.t_old, t)

        # Compute T to satisfy terminal jerk constraint
        while True:
            ## Check T sign, and value
            delta_t = -((j_tg_z * (self.T**3) + 6*a_tg_z * (self.T**2) +
                       (18*v_tg_z + 6*v_g_z)*self.T + 24*(r_tg_z - r_g_z)) /
                       (3*j_tg_z * (self.T**2) + 12*a_tg_z*self.T + 18*v_tg_z + 6*v_g_z))

            ## TODO: Check if this T should be local or global
            self.T = self.T + delta_t

            #print (self.T, delta_t, delta_crit)

            if abs(delta_t) <= abs(delta_crit):
                break

        #print (((j_tg_z * (self.T**3) + 6*a_tg_z * (self.T**2) +
        #               (18*v_tg_z + 6*v_g_z)*self.T + 24*(r_tg_z - r_g_z))))

        # Guidance Equation
        T_leadtime = 0.0
        T_p = self.T + T_leadtime

        T_frac = T_p / self.T

        #print (t, self.T)

        #print ("target", self.P63_target_trajectory, r_g)
        #print (12*(r_tg - r_g), 1/(self.T**2))
        #print (6*v_tg/self.T)
        #print (6*v_g/self.T)
        #print (a_tg)

        ## TODO: Calc target
        self.a_cg = ((3*(T_frac**2) - 2*T_frac)*12*(r_tg - r_g)/(self.T**2) +
                (4*(T_frac**2) - 3*T_frac)*6*v_tg/self.T +
                (2*(T_frac**2) - T_frac)*6*v_g/self.T +
                (6*(T_frac**2) - 6*T_frac + 1)*a_tg)

        #print (12*(r_tg - r_g)/(self.T**2))
        #print (6*v_tg/self.T)
        #print (6*v_g/self.T )
        #print (a_tg)
        #print ("a_cg", self.a_cg)

        # Compute thrust-acceleration command and unit thrust command
        self.a_fcp = self.CGP(self.a_cg - g_p)
        self.u_nfcp = utils.toUnit(self.a_fcp)
                

solver = P63Solver(np.asarray([R_moon, 0, 0], dtype='float64'),
                   calcP63TargetTrajectory())

h = 15e3 + R_moon
z = 492e3
vx = 0
vz = -math.sqrt(G*M_moon / h)
T = T0
delta_t = 5

print (solve(calcP63TargetTrajectory(), T0))
print (solve(calcP63TargetTrajectory(), -1172.76))
print (h, z)
print (vx, vz)

pos_x = []
pos_z = []

for i in range(100):
    T += delta_t
    solver.calc(np.asarray([h, 0, z]), np.asarray([vx, 0, vz]),
                T, np.asarray([0, 0, 0]))

    h += vx*delta_t
    z += vz*delta_t

    vx += solver.a_cg[0]*delta_t
    vz += solver.a_cg[2]*delta_t

    print (h - R_moon, pos_z)
    print (vx, vz)

    pos_x.append(h - R_moon)
    pos_z.append(z)

fig, ax = plt.subplots()
ax.plot(pos_z, pos_x)
ax.set_xlim(pos_z[0], pos_z[-1] - 10e3)
ax.set_xlabel('z (m)')
ax.set_ylabel('x (m)')
ax.set_title('Breaking phase target')

plt.show()
