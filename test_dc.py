import vbcontrolpy as vb
import numpy as np
import matplotlib.pyplot as plt


print(vb.math.sign(10))

arr = np.array([0.0, 1.0, 9.0, 2.0, -5.0], dtype=np.float64)
print(vb.math.clip(arr, 0.0, 8.0))

# DC Motor Control Test
Ts = 0.01
b=0.001
J=0.00001
K=0.04
R=29.0
L=1.28

x = [0.0, 0.0, 0.0]
u = [1.0]

cur_time = 0.0
cnt = 0

theta = []
d_theta = []
cur = []
time = []

motor = vb.dynamics.dc_motor(Ts)
motor.set_physical_params(b, J, K, R, L)

for i in range(int(1/Ts)):
    x_ = motor.nonlinear_dynamics(x, u)
    
    if cnt % 10/Ts == 0:
        theta.append(x[0])
        d_theta.append(x[1])
        cur.append(x[2])
        time.append(cur_time)

    cur_time += Ts
    x = x_

print(f"Size: {len(time)}")

plt.title("DC Motor Nonlinear Dynamics", loc = 'center', size=30)

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

