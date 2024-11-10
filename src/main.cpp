#include "GyrobroDynamics.hpp"
#include "LinearKalmanFilter.hpp"
#include "ExtendedKalmanFilter.hpp"
#include "UnscentedKalmanFilter.hpp"
#include "ss.hpp"
#include "PID.hpp"
#include "FullStateController.hpp"
#include "DLQR.hpp"
#include "DCMotorDynamics.hpp"
#include "PMSMDynamics.hpp"
#include "PMSMCascadeFOC.hpp"

// #include "LowPassFilter.hpp"

#include "iostream"

#include <matplot/matplot.h>

// #include <stdlib.h>
#include <string>
#include <cmath>
#include <Eigen/Dense>
// #include <matplotlib-cpp/matplotlibcpp.h>
using namespace std;
using namespace Eigen;
using namespace matplot;


void create_lti(MatrixXd &A, MatrixXd &B, MatrixXd &C, double Ts)
{
    // double Ts = 1.0/500.0;
    // int dim_x = 7;
    // int dim_u = 2;
    // int dim_z = 5;

    double g = 9.81;

    // Wheel parameters
    double m = 0.977;						// wheel weight [kg]
    double hw = 0.045;           // wheel width [m]
    double R = 0.0795;						// wheel radius [m]

    // Body parameters
    double M = 4.414;                       // body weight [kg]
    double W = 0.37;						// body width [m]
    double D = 0.09;					    // body depth [m]
    double h = 0.4;						// body height [m]

    double fm = 0.0022;					    // friction coefficient between body & DC motor
    double fw = 0.011;          			    // friction coefficient between wheel & floor

    // Motors parameters
    double M_stator = 0.548;
    double R_stator = 0.0493;                // stator radius [m]
    double hw_stator = 0.02;
    double Jm = 0.00001;						// DC motor inertia moment [kgm**2]
    double Rm = 0.4; //0.8 // motor resistance
    double n = 1.0;							// Gear ratio
    double Kt = 0.42; // Motor torque constant
    double Kb = 0.37;

    MatrixXd AA;
    MatrixXd BB;
    MatrixXd CC;
    // GyrobroDynamics dyn(Ts);
    GyrobroDynamics dyn;
    dyn.set_timestep(Ts);
    dyn.set_physical_params(g, m, hw, R, M, W, D, h, fm, fw, M_stator, R_stator, hw_stator, Rm, Jm, Kt, Kb, n);
    dyn.ss_model(AA, BB, CC);
    A = AA;
    B = BB;
    C = CC;
}

void check_lkf()
{
    double Ts = 0.002;
    int dim_x = 7;
    int dim_u = 2;
    int dim_z = 5;

    VectorXd v1(dim_z);
    VectorXd v2(dim_x);
    VectorXd v3(dim_x);

    MatrixXd A;
    MatrixXd B;
    MatrixXd C;
    MatrixXd H(dim_z, dim_x);
    MatrixXd Q(dim_z, dim_x);
    MatrixXd R(dim_z, dim_x);
    MatrixXd P0(dim_x, dim_x);
    VectorXd x(dim_x);
    VectorXd u(dim_u);
    VectorXd z(dim_z);
    MatrixXd F(dim_x, dim_x);
    MatrixXd G(dim_x, dim_u);

    create_lti(A, B, C, Ts);
    auto [Ad, Bd] = ss::c2d(A, B, Ts);
    cout << "-----------" << endl;
    cout << "Dynamics Matrices" << endl;
    cout << "-----------" << endl;
    cout << "A:" << endl;
    cout << Ad << endl;
    cout << "B:" << endl;
    cout << Bd << endl;

    v1 << 2.0, 2.0, 0.05, 0.1, 0.02;
    R = v1.array().matrix().asDiagonal();

    v2 << pow(Ts,4)/4, 10*pow(Ts,2), 1.0, 4000.0, 100.0, pow(Ts,2), 1.0;
    Q = v2.array().matrix().asDiagonal();
    Q(0,1) = pow(Ts,3)/2;
    Q(0,3) = pow(Ts,2)/2;
    Q(1,0) = pow(Ts,3)/2;
    Q(1,3) = Ts;
    Q(2,4) = Ts;
    Q(3,0) = pow(Ts,2)/2;
    Q(3,1) = Ts;
    Q(4,2) = Ts;
    Q(5,6) = Ts;
    Q(6,5) = Ts;

    cout << "-----------" << endl;
    cout << "Covariance Matrices" << endl;
    cout << "-----------" << endl;
    cout << "R:" << endl;
    cout << R << endl;
    cout << "Q:" << endl;
    cout << Q << endl;

    v3 << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
    P0 = v3.array().matrix().asDiagonal();
    
    x.setZero();

    F = Ad;
    G = Bd;
    H = C;

    LinearKalmanFilter lkf(F, G, H, Q, R);
    lkf.set_initial_estimate(x, P0);

    u << 0.1, -0.1;
    z << 0.1, 0.2, 0.3, 0.4, 0.5;
    cout << "-----------" << endl;
    cout << "LKF Output" << endl;
    cout << "-----------" << endl;
    for (int i=0; i<10; i++)
    {
        x = lkf.step(u, z);

        cout << "i=" << i << endl;
        cout << x << endl;
        cout << "-----------" << endl;
    } 

}

void check_ekf()
{
    double Ts = 0.002;
    int dim_x = 7;
    int dim_u = 2;
    int dim_z = 5;

    VectorXd v1(dim_z);
    VectorXd v2(dim_x);
    VectorXd v3(dim_x);

    MatrixXd HH(dim_z, dim_x);
    MatrixXd QQ(dim_z, dim_x);
    MatrixXd RR(dim_z, dim_x);
    MatrixXd P0(dim_x, dim_x);
    VectorXd x(dim_x);
    VectorXd u(dim_u);
    VectorXd z(dim_z);
    MatrixXd F(dim_x, dim_x);
    MatrixXd G(dim_x, dim_u);

    // EKF params: Q. R. P0, x0
    v1 << 2.0, 2.0, 0.05, 0.1, 0.02;
    RR = v1.array().matrix().asDiagonal();

    v2 << pow(Ts, 4)/4, 10*pow(Ts,2), 1.0, 4000.0, 100.0, pow(Ts,2), 1.0;
    QQ = v2.array().matrix().asDiagonal();
    QQ(0,1) = pow(Ts,3)/2;
    QQ(0,3) = pow(Ts,2)/2;
    QQ(1,0) = pow(Ts,3)/2;
    QQ(1,3) = Ts;
    QQ(2,4) = Ts;
    QQ(3,0) = pow(Ts,2)/2;
    QQ(3,1) = Ts;
    QQ(4,2) = Ts;
    QQ(5,6) = Ts;
    QQ(6,5) = Ts;

    cout << "-----------" << endl;
    cout << "Covariance Matrices" << endl;
    cout << "-----------" << endl;
    cout << "R:" << endl;
    cout << RR << endl;
    cout << "Q:" << endl;
    cout << QQ << endl;

    v3 << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
    P0 = v3.array().matrix().asDiagonal();
    
    x.setZero();

    HH <<      0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    // dynamics model
    double g = 9.81;

    // Wheel parameters
    double m = 0.977;						// wheel weight [kg]
    double hw = 0.045;           // wheel width [m]
    double R = 0.0795;						// wheel radius [m]

    // Body parameters
    double M = 4.414;                       // body weight [kg]
    double W = 0.37;						// body width [m]
    double D = 0.09;					    // body depth [m]
    double h = 0.4;						// body height [m]

    double fm = 0.0022;					    // friction coefficient between body & DC motor
    double fw = 0.011;          			    // friction coefficient between wheel & floor

    // Motors parameters
    double M_stator = 0.548;
    double R_stator = 0.0493;                // stator radius [m]
    double hw_stator = 0.02;
    double Jm = 0.00001;						// DC motor inertia moment [kgm**2]
    double Rm = 0.4; //0.8 // motor resistance
    double n = 1.0;							// Gear ratio
    double Kt = 0.42; // Motor torque constant
    double Kb = 0.37;

    GyrobroDynamics dyn(Ts);
    dyn.set_physical_params(g, m, hw, R, M, W, D, h, fm, fw, M_stator, R_stator, hw_stator, Rm, Jm, Kt, Kb, n);
    
    // // VectorXd(*dyn)(VectorXd, VectorXd);
    // // MatrixXd(*jacobian)(VectorXd);

    ExtendedKalmanFilter ekf;
    ekf.set_params(HH, QQ, RR);
    // ekf.set_dynamics(dyn);
    ekf.set_initial_estimate(x, P0);

    u << 0.1, -0.1;
    z << 0.1, 0.2, 0.3, 0.4, 0.5;

    cout << "-----------" << endl;
    cout << "EKF Output" << endl;
    cout << "-----------" << endl;
    for (int i=0; i<10; i++)
    {
        x = ekf.step(dyn, u, z);

        cout << "i=" << i << endl;
        cout << x << endl;
        cout << "-----------" << endl;
    } 
    cout << "Finish!" << endl;

}

void check_ukf()
{
    double Ts = 0.002;
    int dim_x = 7;
    int dim_u = 2;
    int dim_z = 5;

    VectorXd v1(dim_z);
    VectorXd v2(dim_x);
    VectorXd v3(dim_x);

    // MatrixXd A(dim_x, dim_x);
    // MatrixXd B(dim_x, dim_u);
    // MatrixXd C(dim_z, dim_x);
    MatrixXd HH(dim_z, dim_x);
    MatrixXd QQ(dim_z, dim_x);
    MatrixXd RR(dim_z, dim_x);
    MatrixXd P0(dim_x, dim_x);
    VectorXd x(dim_x);
    VectorXd u(dim_u);
    VectorXd z(dim_z);
    MatrixXd F(dim_x, dim_x);
    MatrixXd G(dim_x, dim_u);

    // EKF params: Q. R. P0, x0
    v1 << 2.0, 0.1, 0.05, 0.1, 0.02;
    RR = v1.array().matrix().asDiagonal();

    v2 << pow(Ts, 4)/4, 10*pow(Ts,2), 1.0, 1000.0, 100.0, 4*pow(Ts,2), 1.0;
    QQ = v2.array().matrix().asDiagonal();
    QQ(0,1) = pow(Ts,3)/2;
    QQ(0,3) = pow(Ts,2)/2;
    QQ(1,0) = pow(Ts,3)/2;
    QQ(1,3) = Ts;
    QQ(2,4) = Ts;
    QQ(3,0) = pow(Ts,2)/2;
    QQ(3,1) = Ts;
    QQ(4,2) = Ts;
    QQ(5,6) = Ts;
    QQ(6,5) = Ts;

    cout << "-----------" << endl;
    cout << "Covariance Matrices" << endl;
    cout << "-----------" << endl;
    cout << "R:" << endl;
    cout << RR << endl;
    cout << "Q:" << endl;
    cout << QQ << endl;

    v3 << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
    P0 = v3.array().matrix().asDiagonal();
    
    x.setZero();

    HH <<      0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    // dynamics model
    double g = 9.81;

    // Wheel parameters
    double m = 0.977;						// wheel weight [kg]
    double hw = 0.045;           // wheel width [m]
    double R = 0.0795;						// wheel radius [m]

    // Body parameters
    double M = 4.414;                       // body weight [kg]
    double W = 0.37;						// body width [m]
    double D = 0.09;					    // body depth [m]
    double h = 0.4;						// body height [m]

    double fm = 0.0022;					    // friction coefficient between body & DC motor
    double fw = 0.011;          			    // friction coefficient between wheel & floor

    // Motors parameters
    double M_stator = 0.548;
    double R_stator = 0.0493;                // stator radius [m]
    double hw_stator = 0.02;
    double Jm = 0.00001;						// DC motor inertia moment [kgm**2]
    double Rm = 0.4; //0.8 // motor resistance
    double n = 1.0;							// Gear ratio
    double Kt = 0.42; // Motor torque constant
    double Kb = 0.37;

    GyrobroDynamics dyn(Ts);
    dyn.set_physical_params(g, m, hw, R, M, W, D, h, fm, fw, M_stator, R_stator, hw_stator, Rm, Jm, Kt, Kb, n);
    
    // // VectorXd(*dyn)(VectorXd, VectorXd);
    // // MatrixXd(*jacobian)(VectorXd);
    float kappa = 0;
    float alpha = 0.1;
    float beta = 2;

    UnscentedKalmanFilter ukf;
    ukf.set_params(HH, QQ, RR, kappa, alpha, beta);
    // ukf.set_dynamics(dyn);
    ukf.set_initial_estimate(x, P0);

    u << 0.1, -0.1;
    z << 0.1, 0.2, 0.3, 0.4, 0.5;

    cout << "-----------" << endl;
    cout << "UKF Output" << endl;
    cout << "-----------" << endl;
    for (int i=0; i<10; i++)
    {
        x = ukf.step(dyn, u, z);

        cout << "i=" << i << endl;
        cout << x << endl;
        cout << "-----------" << endl;
    } 
    cout << "Finish!" << endl;

}

void check_lti()
{
    double Ts = 1.0/500.0;
    // int dim_x = 7;
    // int dim_u = 2;
    // int dim_z = 5;

    MatrixXd A;
    MatrixXd B;
    MatrixXd C;

    create_lti(A, B, C, Ts);
    cout << "----------------------------------" << endl;
    cout << "                A                 " << endl;
    cout << "----------------------------------" << endl;
    cout << A << endl;
    cout << "----------------------------------" << endl;
    cout << endl;
    cout << endl;
    cout << "----------------------------------" << endl;
    cout << "                B                 " << endl;
    cout << "----------------------------------" << endl;
    cout << B << endl;
    cout << "----------------------------------" << endl;
    cout << endl;
    cout << endl;
    cout << "----------------------------------" << endl;
    cout << "                C                 " << endl;
    cout << "----------------------------------" << endl;
    cout << C << endl;
    cout << "----------------------------------" << endl;
    
}

void check_c2d()
{
    double Ts = 1.0/500.0;
    MatrixXd A;
    MatrixXd B;
    MatrixXd C;

    create_lti(A, B, C, Ts);

    auto [Ad, Bd] = ss::c2d(A, B, Ts);

    cout << "----------------------------------" << endl;
    cout << "                Ad                " << endl;
    cout << "----------------------------------" << endl;
    cout << Ad << endl;
    cout << "----------------------------------" << endl;
    cout << endl;
    cout << endl;
    cout << "----------------------------------" << endl;
    cout << "                Bd                " << endl;
    cout << "----------------------------------" << endl;
    cout << Bd << endl;
    cout << "----------------------------------" << endl;
    
}

void check_nonlinear_dyn()
{
    double Ts = 0.002;

    double g = 9.81;

    // Wheel parameters
    double m = 0.977;						// wheel weight [kg]
    double hw = 0.045;           // wheel width [m]
    double R = 0.0795;						// wheel radius [m]

    // Body parameters
    double M = 4.414;                       // body weight [kg]
    double W = 0.37;						// body width [m]
    double D = 0.09;					    // body depth [m]
    double h = 0.4;						// body height [m]

    double fm = 0.0022;					    // friction coefficient between body & DC motor
    double fw = 0.011;          			    // friction coefficient between wheel & floor

    // Motors parameters
    double M_stator = 0.548;
    double R_stator = 0.0493;                // stator radius [m]
    double hw_stator = 0.02;
    double Jm = 0.00001;						// DC motor inertia moment [kgm**2]
    double Rm = 0.4; //0.8 // motor resistance
    double n = 1.0;							// Gear ratio
    double Kt = 0.42; // Motor torque constant
    double Kb = 0.37;

    GyrobroDynamics dyn(Ts);
    dyn.set_physical_params(g, m, hw, R, M, W, D, h, fm, fw, M_stator, R_stator, hw_stator, Rm, Jm, Kt, Kb, n);

    VectorXd x(7), x_(7);
    VectorXd u(2);
    MatrixXd J(7,7);

    cout << "----------------------------------" << endl;
    cout << "Test 1:" << endl;
    cout << "----------------------------------" << endl;
    x << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    u << 0.0, 0.0;
    x_ = dyn.nonlinear_dynamics(x, u);
    J = dyn.jacobian(x);
    cout << "x" << endl;
    cout << x_ << endl;
    cout << "----------------------------------" << endl;
    cout << "J" << endl;
    cout << J << endl;
    cout << "----------------------------------" << endl;
    cout << endl;
    cout << endl;
    cout << "----------------------------------" << endl;
    cout << "Test 2:" << endl;
    cout << "----------------------------------" << endl;
    x << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0;
    u << 8.0, 9.0;
    x_ = dyn.nonlinear_dynamics(x, u);
    J = dyn.jacobian(x_);
    cout << "x" << endl;
    cout << x_ << endl;
    cout << "----------------------------------" << endl;
    cout << "J" << endl;
    cout << J << endl;
    cout << "----------------------------------" << endl;
    cout << endl;
    cout << endl;
}

void check_pid()
{
    double Ts = 1.0;
    double Kp = 10;
    double Kd = 1;
    double Ki = 0.1;
    double u_max = 100;
    double u_min = -50;
    double u;
    double ref, cur;

    PID pid(Kp, Ki, Kd, Ts, u_min, u_max, true);

    ref = 1;
    for(int i=0; i<12; i++)
    {
        cur = (double)i/10;
        u = pid.calculate(ref, cur);
        cout << "ref: " << ref << " | cur: " << cur << " | u: " << u << endl;
    }
}

void check_lqr()
{
    double Ts = 1.0/500.0;
    int dim_x = 7;
    int dim_u = 2;

    MatrixXd A;
    MatrixXd B;
    MatrixXd C;
    MatrixXd Q;
    MatrixXd R;
    MatrixXd N;
    MatrixXd K;

    Q.resize(dim_x, dim_x);
    R.resize(dim_u, dim_u);
    N.resize(dim_x, dim_u);
    Q <<    1,     0,     0,     0,     0,     0,     0,
            0,     1,     0,     0,     0,     0,     0,
            0,     0,     10,    0,     0,     0,     0,
            0,     0,     0,     1,     0,     0,     0,
            0,     0,     0,     0,     1,     0,     0,
            0,     0,     0,     0,     0,     1000,  0,
            0,     0,     0,     0,     0,     0,     100;

    R << 20000, 0,
         0,     20000;

    create_lti(A, B, C, Ts);
    auto [Ad, Bd] = ss::c2d(A, B, Ts);

    DLQR dlqr;
    dlqr.set_params(Ad, Bd, Q, R, N);
    K = dlqr.solveRiccati();

    Eigen::IOFormat fmt(4, 0, ", ", "\n", "[", "]");
    cout << "K:" << endl;
    std::cout << K.format(fmt) << std::endl;
}

void check_dc()
{
    float Ts = 0.01;
    float b=0.001;
    float J=0.00001;
    float K=0.04;
    float R=29.0;
    float L=1.28;

    VectorXd x(3), x_(3);
    VectorXd u(1);

    std::vector<double> theta;
    std::vector<double> d_theta;
    std::vector<double> cur;
    std::vector<double> time;
    double cur_time = 0;
    int cnt = 0;

    DCMotorDynamics motor;
    motor.set_timestep(Ts);
    motor.set_physical_params(b, J, K, R, L);

    x << 0.0, 0.0, 0.0;
    u << 1.0;
    for (int i = 0; i < 1/Ts; i++)
    {
        x_ = motor.nonlinear_dynamics(x, u);
        // cout << "x: " << x_ << endl;
        x = x_;

        if (cnt % 10/Ts == 0)
        {
            theta.push_back(x(0));
            d_theta.push_back(x(1));
            cur.push_back(x(2));
            time.push_back(cur_time);
        }

        cur_time += Ts;
    }

    cout << "Size: " << time.size() << endl;

    auto h = figure(true);
    h->name("DC Motor Nonlinear Dynamics");
    h->size(1500, 1500);
    h->tiledlayout(3, 1);
    h->title("DC Motor Nonlinear Dynamics");
    
    // tiledlayout(3, 1);
    auto ax1 = h->nexttile();
    plot(ax1, time, theta)->line_width(2);
    title(ax1, "Angle");
    ylabel(ax1, "\theta [rad]");
    grid(on);
    
    auto ax2 = h->nexttile();
    plot(ax2, time, d_theta)->line_width(2);
    title(ax2, "Velocity");
    ylabel(ax2, "d_theta [rad/s]");
    grid(on);

    auto ax3 = h->nexttile();
    plot(ax3, time, cur)->line_width(2);
    title(ax3, "Current");
    ylabel(ax3, "I [A]");
    grid(on);

    show();
        
} 

void check_pmsm()
{
    float Ts = 0.0001;

    float R = 1.0;
    float U_max = 24;
    float L = 0.000344;
    float J = 1.04e-2;
    float Kt = 0.66;
    float Kw = 0.969;
    float b=0.005;
    float pole_pairs = 14;
    float gear_ratio = 1;

    float Z = 24;
    VectorXd Tk(4);
    VectorXd k(4);
    VectorXd alpha(4);
    Tk << 0.5, 0.2, 0.1, 0.03;
    Tk = Tk / 16;
    k << 1, 2, 3, 4;
    alpha << 0, 0.01, 0.017, 0.017;
    
    VectorXd x(4), x_(4);
    VectorXd u(2);

    std::vector<double> theta;
    std::vector<double> d_theta;
    std::vector<double> i_d, i_q;
    std::vector<double> time;
    double cur_time = 0;
    int cnt = 0;

    PMSMDynamics motor;
    motor.set_timestep(Ts);
    motor.set_physical_params(R, U_max, L, J, Kt, Kw, b, pole_pairs, gear_ratio);
    motor.set_cogging_torque_params(Tk, k, alpha, Z);

    x << 0.0, 0.0, 0.0, 0.0;
    u << 0.0, 11.39;

    for (int i = 0; i < 0.01/Ts; i++)
    {
        x_ = motor.nonlinear_dynamics(x, u);
        // cout << "x: " << x_ << endl;
        x = x_;

        if (cnt % 100/Ts == 0)
        {
            i_d.push_back(x(0));
            i_q.push_back(x(1));
            theta.push_back(x(2));
            d_theta.push_back(x(3));
            time.push_back(cur_time);
        }

        cur_time += Ts;
    }

    cout << "Size: " << time.size() << endl;

    auto h = figure(true);
    h->name("PMSM Motor Nonlinear Dynamics");
    h->size(1500, 1500);
    h->tiledlayout(4, 1);
    h->title("PMSM Motor Nonlinear Dynamics");
    
    // tiledlayout(3, 1);
    auto ax1 = h->nexttile();
    plot(ax1, time, theta)->line_width(2);
    // title(ax1, "Angle");
    ylabel(ax1, "theta [rad]");
    grid(on);
    
    auto ax2 = h->nexttile();
    plot(ax2, time, d_theta)->line_width(2);
    // title(ax2, "Velocity");
    ylabel(ax2, "d theta [rad/s]");
    grid(on);

    auto ax3 = h->nexttile();
    plot(ax3, time, i_d)->line_width(2);
    // title(ax3, "i_d");
    ylabel(ax3, "i_d [A]");
    grid(on);

    auto ax4 = h->nexttile();
    plot(ax4, time, i_q)->line_width(2);
    // title(ax4, "i_q");
    ylabel(ax4, "i_q [A]");
    grid(on);

    show();
}

void check_pmsm_foc()
{
    double Ts = 0.00001;

    double R = 1.0;
    double U_max = 24;
    double L = 0.000344;
    double J = 1.04e-2;
    double Kt = 0.66;
    double Kw = 0.969;
    double b=0.005;
    double pole_pairs = 14;
    double gear_ratio = 1;

    double Z = 24;
    VectorXd Tk(4);
    VectorXd k(4);
    VectorXd alpha(4);
    Tk << 0.5, 0.2, 0.1, 0.03;
    Tk = Tk / 16;
    k << 1, 2, 3, 4;
    alpha << 0, 0.01, 0.017, 0.017;

    double theta_ref, w_ref, tau_ref, voltage_ref;

    double kp = 8.0;
    double kd = 0.2;
    double kp_id = 1.88;
    double ki_id = 0.0046;
    double kp_iq = 4*1.88;
    double ki_iq = 0.0046;
    
    VectorXd x(4);

    std::vector<double> i_d, i_q;
    std::vector<double> theta;
    std::vector<double> d_theta;
    std::vector<double> time;
    double cur_time = 0;
    int cnt = 0;

    PMSMCascadeFOC foc;
    foc.set_timestep(Ts);
    foc.set_feedback(kp, kd, kp_id, ki_id, kp_iq, ki_iq);
    foc.set_motor_params(R, U_max, L, J, Kt, Kw, b, pole_pairs, gear_ratio, Tk, k, alpha, Z);

    theta_ref = 1.0;
    w_ref = 0;
    tau_ref = 0;
    voltage_ref = 0;

    for (int i = 0; i < 0.5/Ts; i++)
    {
        x = foc.calculate(theta_ref, w_ref, tau_ref, voltage_ref);

        if (cnt % 100/Ts == 0)
        {
            i_d.push_back(x(0));
            i_q.push_back(x(1));
            theta.push_back(x(2));
            d_theta.push_back(x(3));
            time.push_back(cur_time);
        }

        cur_time += Ts;
    } 

    cout << "Size: " << time.size() << endl;

    auto h = figure(true);
    h->name("PMSM Motor FOC Control");
    h->size(1500, 1500);
    h->tiledlayout(4, 1);
    h->title("PMSM Motor FOC Control");
    
    auto ax1 = h->nexttile();
    plot(ax1, time, theta)->line_width(2);
    // title(ax1, "Angle");
    ylabel(ax1, "theta [rad]");
    grid(on);
    
    auto ax2 = h->nexttile();
    plot(ax2, time, d_theta)->line_width(2);
    // title(ax2, "Velocity");
    ylabel(ax2, "d theta [rad/s]");
    grid(on);

    auto ax3 = h->nexttile();
    plot(ax3, time, i_d)->line_width(2);
    // title(ax3, "i_d");
    ylabel(ax3, "i_d [A]");
    grid(on);

    auto ax4 = h->nexttile();
    plot(ax4, time, i_q)->line_width(2);
    // title(ax4, "i_q");
    ylabel(ax4, "i_q [A]");
    grid(on);

    show();
}

int main(int argc, char *argv[])
{
    if (argc > 2)
        throw runtime_error("Too many arguments!");
    if (argc <= 1)
        throw runtime_error("Please use agrument --test=<num> where <num> is a test number from 1 to 11.");

    string arg_str = argv[1];
    arg_str.erase(0, 7);
    int test_num = atoi(arg_str.c_str());

    switch(test_num) {
        case 1:
            // Gyroboy linear dynamics
            cout << "Linear Dynamics checking..." << endl;
            check_lti();
            break;
        case 2:
            // continuous to discrete transform
            cout << "c2d checking..." << endl;
            check_c2d();
            break;
        case 3:
            // Gyroboy nonlinear dynamics
            cout << "Nonlinear Dynamics checking..." << endl;
            check_nonlinear_dyn();
            break;
        case 4:
            // LKF
            cout << "Linear Kalmna Filter checking..." << endl;
            check_lkf();
            break;
        case 5:
            // EKF
            cout << "Extended Kalman Filter checking..." << endl;
            check_ekf();
            break;
        case 6:
            // UKF
            cout << "Unscented Kalman Filter checking..." << endl;
            check_ukf();
            break;
        case 7:
            // PID
            cout << "PID Controller checking..." << endl;
            check_pid();
            break;
        case 8:
            // DLQR
            cout << "LQR Solver checking..." << endl;
            check_lqr();
            break;
        case 9:
            // DC Motor Nonlinear Dynamics
            cout << "DC Motor Nonlinear Dynamics checking..." << endl;
            check_dc();
            break;
        case 10:
            // PMSM Motor Nonlinear Dynamics
            cout << "PMSM Motor Nonlinear Dynamics checking..." << endl;
            check_pmsm();
            break;
        case 11:
            // PMSM Motor FOC Control
            cout << "PMSM Motor FOC Control checking..." << endl;
            check_pmsm_foc();
            break;
        default:
            cout << "Test number must be from 1 to 8" << endl;
    }
    
    

    return 0;
}
