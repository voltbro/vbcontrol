#include "PMSMDynamics.hpp"

PMSMDynamics::PMSMDynamics()
    : dim_x(DIM_X_PMSM), dim_y(DIM_Y_PMSM), dim_u(DIM_U_PMSM),
      A(dim_x, dim_x), B(dim_x, dim_u), C(dim_y, dim_x),
      x_(dim_x)
{

}

PMSMDynamics::PMSMDynamics(double Ts)
    : dim_x(DIM_X_PMSM), dim_y(DIM_Y_PMSM), dim_u(DIM_U_PMSM),
      A(dim_x, dim_x), B(dim_x, dim_u), C(dim_y, dim_x),
      x_(dim_x)
{
    this->Ts = Ts;
}

void PMSMDynamics::set_physical_params( double R, // electric resistance
                                        double U_max, // maximal voltage
                                        double L, // electric inductance
                                        double J, // moment of inertia of the rotor
                                        double Kt, // torque constant
                                        double Kw, // velocity constant
                                        double b, // motor viscous frisction constant
                                        double pole_pairs,
                                        double gear_ratio
                                        ) 
{
    this->b = b;
    this->J = J;
    this->Kt = Kt;
    this->Kw = Kw;
    this->R = R;
    this->L = L;
    this->U_max = U_max;
    this->pole_pairs = pole_pairs;
    this->gear_ratio = gear_ratio;

    //Misc consts
    k_phi = Kw/pole_pairs;
    a11 = -R/L;
    a12 = pole_pairs;

    a21 = -R/L;
    a22 = -pole_pairs;
    a23 = -k_phi*pole_pairs/L;

    a41 = 3*pole_pairs*k_phi/(2*J);
    a42 = -b/J;
    a43 = 1/J;
}

void PMSMDynamics::set_cogging_torque_params (    VectorXd Tk,
                                        VectorXd k,
                                        VectorXd alpha,
                                        double Z
                                        )
{
    this->Tk = Tk;
    this->k = k;
    this->alpha = alpha;
    this->Z = Z;
}

void PMSMDynamics::ss_model(MatrixXd &A, MatrixXd &B, MatrixXd &C)
{
    this->A <<  0, 0, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1,
                0, 0, 0, 0;

    this->B <<  1, 0,
                0, 0,
                0, 0,
                0, 1;

    this->C <<  1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 0;

    A = this->A;
    B = this->B;
    C = this->C;
}

VectorXd PMSMDynamics::cont_nonlinear_dynamics(const VectorXd &x, const VectorXd &u)
{
    x1 = x(0);
    x2 = x(1);
    x3 = x(2);
    x4 = x(3);
    u1 = u(0);
    u2 = u(1);

    if (u(0) > 24)
        u1 = 24;
    else if (u(0) < -24)
        u1 = -24;

    if (u(1) > 24)
        u2 = 24;
    else if (u(1) < -24)
        u2 = -24;
    
    Tcog = 0;
    for (int i = 0; i < 4; i++)
        Tcog += Tk(i) * sin(k(i)*Z*x3 + alpha(i));
    
    d_id = a11*x1 + a12*x2*x4 + u1/L;
    d_iq = a21*x2 + a22*x1*x4 + a23*x4 + u2/L;
    d_theta = x4;
    dd_theta = a41*x2 + a42*x4 + a43*Tcog;

    x_ << d_id, d_iq, d_theta, dd_theta;

    return x_;
}

VectorXd PMSMDynamics::nonlinear_dynamics(VectorXd &x, VectorXd &u)
{
    // RK4 Method
    f1 = Ts * cont_nonlinear_dynamics(x,        u);
    f2 = Ts * cont_nonlinear_dynamics(x + f1/2, u);
    f3 = Ts * cont_nonlinear_dynamics(x + f2/2, u);
    f4 = Ts * cont_nonlinear_dynamics(x + f3,   u);
    
    x_ = x + (f1 + 2*f2 + 2*f3 + f4)/6;
    return x_;
}

MatrixXd PMSMDynamics::jacobian(VectorXd &x)
{
    x1 = x(0);
    x2 = x(1);
    x3 = x(2);
    x4 = x(3);
    MatrixXd A(4,4);

    d_Tcog = 0;
    for (int i = 0; i < 4; i++)
        d_Tcog += Z*k(i)*Tk(i) * cos(k(i)*Z*x3 + alpha(i));

    A << 1 + a11*Ts,  a12*x4*Ts,             0,         a12*x2*Ts,
          a22*x4*Ts, 1 + a21*Ts,             0, (a22*x1 + a23)*Ts,
                  0,          0,             1,                Ts,
                  0,     a41*Ts, a43*d_Tcog*Ts,        1 + a42*Ts;
    return A;
}