VBControl - это C++ библиотека для решения задач теории автоматического управления. Также имеет python-обертку.
Возможности:
 - Linear Kalman Filter
 - Extended Kalman Filter
 - Unscented Kalman Filter
 - Low Pass Filter
 - PID Controller
 - Full State Controller
 - Discrete LQR Solver (Newton Method)
 - Convert from continuous to discrete state space systems
 - Fcition Compensation
 - Math operations: clip, sign
 - DC Motor Dynamical Model
 - Two Wheeled Self-Balancing Robot Dynamics

## Установка
```
cd ~
git clone https://github.com/voltbro/vbcontrol.git
cd vbcontrol
mkdir build
cd build
cmake ..
make
```