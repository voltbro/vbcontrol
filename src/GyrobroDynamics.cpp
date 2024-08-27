#include "GyrobroDynamics.hpp"

#define DEBUG_MODE (true)

#if DEBUG_MODE
#include <iostream>
using namespace std;
#endif

GyrobroDynamics::GyrobroDynamics() 
    : dim_x(DIM_X_GB), dim_y(DIM_Y_GB), dim_u(DIM_U_GB), A(dim_x, dim_x), B(dim_x, dim_u), 
      C(dim_y, dim_x)
{
    E.resize(2,2);
    F.resize(2,2);
    G.resize(2,2);
    H.resize(2,2);
}

GyrobroDynamics::GyrobroDynamics(double Ts) 
    : dim_x(DIM_X_GB), dim_y(DIM_Y_GB), dim_u(DIM_U_GB), A(dim_x, dim_x), B(dim_x, dim_u), 
      C(dim_y, dim_x)
{
    this->Ts = Ts;

    E.resize(2,2);
    F.resize(2,2);
    G.resize(2,2);
    H.resize(2,2);
    // eye.setZero();
}

//functions
// void GyrobroDynamics::set_timestep(double Ts)
// {
//     this->Ts = Ts;
// }

void GyrobroDynamics::set_physical_params(double g,
                        double m,
                        double hw,
                        double R,
                        double M,
                        double W,
                        double D,
                        double h,
                        double fm,
                        double fw,
                        double M_stator,
                        double R_stator,
                        double hw_stator,
                        double Rm,
                        double Jm,
                        double Kt,
                        double Kb,
                        double n)
{
    Jw = (m * pow(R, 2) / 2.0) - (M_stator * pow(R_stator, 2) / 2.0);
    Jwx = (m * (3 * pow(R, 2) + pow(hw, 2)) / 12.0) - (M_stator * (3 * pow(R_stator, 2) + pow(hw_stator, 2)) / 12.0);
    L = h / 2.0;
    Jpsi = M * pow(L, 2) / 3.0;
    Jphi = M * (pow(W, 2) + pow(D, 2)) / 12.0;
    this->g = g;
    this->m = m;
    this->hw = hw;
    this->R = R;
    this->M = M;
    this->W = W;
    this->D = D;
    this->h = h;
    this->fm = fm;
    this->fw = fw;
    this->M_stator = M_stator;
    this->R_stator = R_stator;
    this->hw_stator = hw_stator;
    this->Rm = Rm;
    this->Jm = Jm;
    this->Kt = Kt;
    this->Kb = Kb;
    this->n = n;

    calc_basic_params();
}

// gives continuous state space model
void GyrobroDynamics::ss_model(MatrixXd &A, MatrixXd &B, MatrixXd &C)
{
    MatrixXd m1(2,2);
    MatrixXd m2(2,2);
    MatrixXd m3(2,2);
    m1 = -invE * G;
    m2 = -invE * F;
    m3 = invE * H;

    this->A << 0.0, 1.0,      0.0,      0.0,      0.0,     0.0, 0.0,
               0.0, 0.0,      0.0,      1.0,      0.0,     0.0, 0.0,
               0.0, 0.0,      0.0,      0.0,      1.0,     0.0, 0.0,
               0.0, m1(0,0),  m1(0,1),  m2(0,0),  m2(0,1), 0.0, 0.0,
               0.0, m1(1,0),  m1(1,1),  m2(1,0),  m2(1,1), 0.0, 0.0,
               0.0, 0.0,      0.0,      0.0,      0.0,     0.0, 1.0,
               0.0, 0.0,      0.0,      0.0,      0.0,     0.0, -J/I; 

    this->B << 0,       0,
               0,       0,
               0,       0,
               m3(0,0), m3(0,1),
               m3(1,0), m3(1,1),
               0,       0,
               -K/I,    K/I;

    this->C << 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    A = this->A;
    B = this->B;
    C = this->C;
}

VectorXd GyrobroDynamics::nonlinear_dynamics(VectorXd &x, VectorXd &u)
{
    if (x.rows() != dim_x)
        throw runtime_error("Size of x must be " + to_string(dim_x) + "!");
    if (u.rows() != dim_u)
        throw runtime_error("Size of u must be " + to_string(dim_u) + "!");
    
    // cout << "GyrobroDynamics::nonlinear_dynamics" << endl;
    VectorXd x_(dim_x);
    double sin_x2 = sin(x(2, 0));
    double cos_x2 = cos(x(2, 0));


    double a = (-F_21*x(3, 0) - F_22*x(4,0) + L2M * pow(x(6,0), 2) * sin_x2*cos(x(2,0)) + LMG*sin_x2 + alpha*(-u(0,0) - u(1,0)))/detE;
    double b = (-F_11*x(3, 0) - F_12*x(4,0) + LMR * pow(x(4,0), 2) * sin_x2 + alpha*(u(0,0) + u(1,0)))/detE;
    
    double eq4 = -E_12*a + E_22*b;
    double eq5 = E_11*a - E_21*b;
    double eq6 = (-J*x(6,0) + K*(-u(0,0) + u(1,0)) - 2*L2M*x(4,0)*x(6,0)*sin_x2*cos_x2)/I;

    x_ << x(0,0) + x(1,0)*Ts,
             x(1,0) + x(3,0)*Ts, 
             x(2,0) + x(4,0)*Ts, 
             x(3,0) + eq4*Ts, 
             x(4,0) + eq5*Ts,
             x(5,0) + x(6,0)*Ts, 
             x(6,0) + eq6*Ts;

    return x_;
}

MatrixXd GyrobroDynamics::jacobian(VectorXd &x)
{
    if (x.rows() != dim_x)
        throw runtime_error("Size of x must be " + to_string(dim_x) + "!");

    MatrixXd jacobian(dim_x, dim_x);
    jacobian.setIdentity();

    jacobian(0,1) = Ts;
    jacobian(1,3) = Ts;
    jacobian(2,4) = Ts;
    jacobian(5,6) = Ts;

    jacobian(3,2) = -Ts*L*M*(E_12*(-2*L*pow(x(6,0),2)*pow(sin(x(2,0)),2) + L*pow(x(6,0),2) + g*cos(x(2,0))) - E_22*R*pow(x(4,0),2)*cos(x(2,0)))/detE;
    jacobian(3,3) = (Ts*(E_12*F_21 - E_22*F_11) + detE)/detE;
    jacobian(3,4) = Ts*(E_12*F_22 - E_22*(F_12 - 2*LMR*x(4,0)*sin(x(2,0))))/detE;
    jacobian(3,6) = -Ts*E_12*L2M*x(6,0)*sin(2*x(2,0))/detE;
    jacobian(4,2) = Ts*L*M*(E_11*(-2*L*pow(x(6,0),2)*pow(sin(x(2,0)),2) + L*pow(x(6,0),2) + g*cos(x(2,0))) - E_21*R*pow(x(4,0),2)*cos(x(2,0)))/detE;
    jacobian(4,3) = Ts*(-E_11*F_21 + E_21*F_11)/detE;
    jacobian(4,4) = (-Ts*(E_11*F_22 - E_21*(F_12 - 2*LMR*x(4,0)*sin(x(2,0)))) + detE)/detE;
    jacobian(4,6) = Ts*E_11*L2M*x(6,0)*sin(2*x(2,0))/detE;

    jacobian(6,2) = -2*Ts*L2M*x(4,0)*x(6,0)*cos(2*x(2,0))/I;
    jacobian(6,4) = -Ts*L2M*x(6,0)*sin(2*x(2,0))/I;
    jacobian(6,6) = (-Ts*(J + L2M*x(4,0)*sin(2*x(2,0))) + I)/I;

    return jacobian;
}

void GyrobroDynamics::calc_basic_params()
{
    LMG = L*M*g;
    L2M = pow(L,2) * M;
    LMR = L*M*R;

    alpha = n * Kt / Rm;
    beta = n * Kt * Kb / Rm + fm;

    F_11 = 2 * (beta + fw);
    F_12 = -2 * beta;
    F_21 = F_12;
    F_22 = 2 * beta;

    E_11 = (2*m + M) * pow(R, 2) + 2*Jw + 2*pow(n, 2)*Jm;
    E_12 = LMR - 2*pow(n,2)*Jm;
    E_21 = E_12;
    E_22 = L2M + Jpsi + 2*pow(n,2)*Jm;
    detE = E_11*E_22 - E_12*E_21;
    

    I = 0.5*m*pow(W,2) + Jphi + (pow(W,2)/(2*pow(R,2)))*(Jw + pow(n,2)*Jm);
    J = (pow(W,2)/(2*pow(R,2)))*(beta + fw);
    K = (W/(2*R))*alpha;

    E << E_11,   E_12,
         E_21,   E_22;

    F << F_11,   F_12,
         F_21,   F_22;

    G << 0,   0,
         0,  -LMG;

    H << alpha,   alpha,
        -alpha,  -alpha;

    invE = E.inverse();
}

