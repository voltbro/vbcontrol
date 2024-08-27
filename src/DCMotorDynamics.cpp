#include "DCMotorDynamics.hpp"

DCMotorDynamics::DCMotorDynamics()
    : dim_x(DIM_X_DC), dim_y(DIM_Y_DC), dim_u(DIM_U_DC),
      A(dim_x, dim_x), B(dim_x, dim_u), C(dim_y, dim_x),
      x_(dim_x)
{

}

DCMotorDynamics::DCMotorDynamics(double Ts)
    : dim_x(DIM_X_DC), dim_y(DIM_Y_DC), dim_u(DIM_U_DC),
      A(dim_x, dim_x), B(dim_x, dim_u), C(dim_y, dim_x),
      x_(dim_x)
{
    this->Ts = Ts;
}

void DCMotorDynamics::set_physical_params(double b, // motor viscous frisction constant
                                            double J, // moment of inertia of the rotor
                                            double K, // motor constant
                                            double R, // electric resistance
                                            double L // electric inductance
                                            ) 
{
    this->b = b;
    this->J = J;
    this->K = K;
    this->R = R;
    this->L = L;
}

void DCMotorDynamics::ss_model(MatrixXd &A, MatrixXd &B, MatrixXd &C)
{
    this->A << 0,  1,    0,
         0, -b/J,  K/J,
         0, -K/L, -R/L;

    this->B << 0,
         0,
         1/L;

    this->C << 1, 0, 0;

    A = this->A;
    B = this->B;
    C = this->C;
}

VectorXd DCMotorDynamics::nonlinear_dynamics(VectorXd &x, VectorXd &u)
{
    theta = x(0);
    d_theta = x(1);
    i = x(2);
    V = u(0);

    didt = (V - i*R - K*d_theta)/L;
    dd_theta = (K*i - b*d_theta)/J;

    i += didt * Ts;
    d_theta += dd_theta * Ts;
    theta += d_theta * Ts;

    x_ << theta, d_theta, i;

    return x_;
}

MatrixXd DCMotorDynamics::jacobian(VectorXd &x)
{
    // TODO: create jacobian
    MatrixXd A(2,2);
    A << x(0), 1,
        2,    3;
    return A;
}