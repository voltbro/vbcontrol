#include "PMSMCascadeFOC.hpp"


PMSMCascadeFOC::PMSMCascadeFOC()
    : dim_x(DIM_X_PMSM), dim_u(DIM_U_PMSM), x(dim_x), x_(dim_x), u(dim_u)
{
    i_d_ref = 0;
    x << 0, 0, 0, 0,
    u << 0, 0;

    fb_inited = false;
    motor_inited = false;
    ts_inited = false;
    inited = ts_inited & fb_inited & motor_inited;
}

PMSMCascadeFOC::PMSMCascadeFOC(double Ts)
    : dim_x(DIM_X_PMSM), dim_u(DIM_U_PMSM), x(dim_x), x_(dim_x), u(dim_u)
{
    this->Ts = Ts;

    i_d_ref = 0;
    x << 0, 0, 0, 0,
    u << 0, 0;

    ts_inited = true;
    fb_inited = false;
    motor_inited = false;
    inited = ts_inited & fb_inited & motor_inited;
}

PMSMCascadeFOC::PMSMCascadeFOC(double Ts, double kp, double kd, double kp_id, double ki_id, double kp_iq, double ki_iq)
    : dim_x(DIM_X_PMSM), dim_u(DIM_U_PMSM), x(dim_x), x_(dim_x), u(dim_u)
{
    this->Ts = Ts;

    this->kp = kp;
    this->kd = kd;
    this->kp_id = kp_id;
    this->ki_id = ki_id;
    this->kp_iq = kp_iq;
    this->ki_iq = ki_iq;

    i_d_ref = 0;
    x << 0, 0, 0, 0,
    u << 0, 0;

    ts_inited = true;
    fb_inited = true;
    motor_inited = false;
    inited = ts_inited & fb_inited & motor_inited;
}

void PMSMCascadeFOC::set_timestep(double Ts)
{
    this->Ts = Ts;
    motor.set_timestep(Ts);
    cout << "Ts: " << Ts << endl;

    ts_inited = true;
    inited = ts_inited & fb_inited & motor_inited;
}

void PMSMCascadeFOC::set_feedback(double kp, double kd, double kp_id, double ki_id, double kp_iq, double ki_iq)
{
    this->kp = kp;
    this->kd = kd;
    this->kp_id = kp_id;
    this->ki_id = ki_id;
    this->kp_iq = kp_iq;
    this->ki_iq = ki_iq;

    // pd.init(kp, 0, kd, Ts);
    pi_id.init(kp_id, ki_id, 0, Ts);
    pi_iq.init(kp_iq, ki_iq, 0, Ts);
    cout << "Ts: " << Ts << endl;

    fb_inited = true;
    inited = ts_inited & fb_inited & motor_inited;
}

void PMSMCascadeFOC::set_motor_params(  double R, // electric resistance
                            double U_max, // maximal voltage
                            double L, // electric inductance
                            double J, // moment of inertia of the rotor
                            double Kt, // torque constant
                            double Kw, // velocity constant
                            double b, // motor viscous frisction constant
                            double pole_pairs,
                            double gear_ratio,
                            VectorXd Tk, //coggin torque param
                            VectorXd k, //coggin torque param
                            VectorXd alpha, //coggin torque param
                            double Z //coggin torque param
                            )
{
    
    motor.set_physical_params(R, U_max, L, J, Kt, Kw, b, pole_pairs, gear_ratio);
    motor.set_cogging_torque_params(Tk, k, alpha, Z);
    this->Kt = Kt;
    this->gear_ratio = gear_ratio;

    cout << "Ts: " << Ts << endl;

    motor_inited = true;
    inited = ts_inited & fb_inited & motor_inited;
}

VectorXd PMSMCascadeFOC::calculate(double theta_ref, double w_ref, double tau_ref, double voltage_ref)
{
    if (inited == false)
        throw std::runtime_error("PMSMCascadeFOC is not inited!");
    
    i_d = x(0);
    i_q = x(1);
    theta = x(2);
    w = x(3);

    // Angle P controller
    e_p = theta_ref - theta;
    u_p = kp * e_p;
    // Velocity P controller
    e_w = w_ref - w;
    u_w = kd * e_w;
    // Getting i_q_ref
    i_q_ref = (1/(Kt*gear_ratio)) * (u_p + u_w + tau_ref);
    // PI controllers
    u_d = pi_id.calculate(i_d_ref, i_d);
    u_q = pi_iq.calculate(i_q_ref, i_q) + voltage_ref;
    u << u_d, u_q;

    // Model step
    x_ = motor.nonlinear_dynamics(x, u);
    // cout << "x: " << x_ << endl;
    x = x_;
    
    return x;
}
