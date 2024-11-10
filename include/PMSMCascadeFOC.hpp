#ifndef _PMSMCascadeFOC_hpp_
#define _PMSMCascadeFOC_hpp_

#include <cmath>
#include <Eigen/Dense>
#include <stdexcept>
#include <iostream>

#include "vbmath.hpp"
#include "PID.hpp"
#include "PMSMDynamics.hpp"

using namespace Eigen;
using namespace std;

// PID controller class
class PMSMCascadeFOC
{
public:
	//constructors
    PMSMCascadeFOC();
    PMSMCascadeFOC(double Ts);
	PMSMCascadeFOC(double Ts, double kp, double kd, double kp_id, double ki_id, double kp_iq, double ki_iq);

    void set_timestep(double Ts);
    void set_feedback(double kp, double kd, double kp_id, double ki_id, double kp_iq, double ki_iq);
    void set_motor_params(  double R, // electric resistance
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
                            );

    VectorXd calculate(double theta_ref, double w_ref, double tau_ref, double voltage_ref);

private:
    int dim_x, dim_u;

    VectorXd x, x_;
    VectorXd u;

    bool fb_inited, motor_inited, ts_inited, inited;

    double Ts;

    double kp, kd, kp_id, ki_id, kp_iq, ki_iq;
    double u_min, u_max;
    double Kt, gear_ratio;

    double theta, w, tau, i_d, i_q;
    double theta_ref, w_ref, tau_ref, voltage_ref, i_d_ref, i_q_ref;
    double e_p, e_w, e_d, e_q;
    double u_p, u_w;
    double u_d, u_q;

    PMSMDynamics motor;
    PID pi_id, pi_iq, pd;
};

#endif //_PMSMCascadeFOC_hpp_