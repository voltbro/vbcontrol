import vbcontrolpy as vb
import numpy as np
import matplotlib.pyplot as plt


# PMSM Motor Control Test
Ts = 0.0001

# Physical params
R = 1.0
U_max = 24
L = 0.000344
J = 1.04e-2
Kt = 0.66
Kw = 0.969
b=0.005
pole_pairs = 14
gear_ratio = 1

# Cogging torque params
Z = 24
Tk = [0.5/16.0, 0.2/16.0, 0.1/16.0, 0.03/16.0]
k = [1, 2, 3, 4]
alpha = [0, 0.01, 0.017, 0.017]

x = [0.0, 0.0, 0.0, 0.0]
u = [0.0, 5.0]

cur_time = 0.0
cnt = 0

theta = []
d_theta = []
cur = []
time = []

# motor = vb.dynamics.dc_motor(Ts)
# motor.set_physical_params(b, J, K, R, L)

motor = vb.dynamics.pmsm(Ts)
motor.set_physical_params(R, U_max, L, J, Kt, Kw, b, pole_pairs, gear_ratio)
motor.set_cogging_torque_params(Tk, k, alpha, Z)

for i in range(int(1/Ts)):
    x_ = motor.nonlinear_dynamics(x, u)
    
    if cnt % 10/Ts == 0:
        theta.append(x[2])
        d_theta.append(x[3])
        cur.append(x[1])
        time.append(cur_time)

    cur_time += Ts
    x = x_

print(f"Size: {len(time)}")

plt.title("PMSM Motor Nonlinear Dynamics", loc = 'center', size=30)

plt.subplot(3, 1, 1)
plt.plot(time, theta, linewidth = '4')
plt.ylabel("Angle [rad]", size=20)
plt.grid()

plt.subplot(3, 1, 2)
plt.plot(time, d_theta, linewidth = '4')
plt.ylabel("Velocity [rad/s]", size=20)
plt.grid()

plt.subplot(3, 1, 3)
plt.plot(time, cur, linewidth = '4')
plt.xlabel("Current [A]", size=20)
plt.ylabel("time [s]", size=20)
plt.grid()

plt.show()

