
#ifndef _PMSMDynamics_hpp_
#define _PMSMDynamics_hpp_

#define DIM_X_PMSM 4
#define DIM_U_PMSM 2
#define DIM_Y_PMSM 4

// #include <cmath>
#include <Eigen/Dense>
 
#include "SysModel.hpp"
 
using namespace Eigen;
using namespace std;

// Self-Balancing robot dynamics class
class PMSMDynamics : public SysModel 
{
public: 
	//constructors
	PMSMDynamics();
	PMSMDynamics(double Ts); // Ts - timestep [ms]
	//functions
	void set_physical_params(   double R, // electric resistance
                                double U_max, // maximal voltage
                                double L, // electric inductance
                                double J, // moment of inertia of the rotor
                                double Kt, // torque constant
                                double Kw, // velocity constant
                                double b, // motor viscous frisction constant
                                double pole_pairs,
                                double gear_ratio
                                );
    void set_cogging_torque_params (    VectorXd Tk,
                                        VectorXd k,
                                        VectorXd alpha,
                                        double Z
                                        );
	void ss_model(MatrixXd &A, MatrixXd &B, MatrixXd &C) override;
	VectorXd nonlinear_dynamics(VectorXd &x, VectorXd &u) override;
    VectorXd cont_nonlinear_dynamics(const VectorXd &x, const VectorXd &u);
	MatrixXd jacobian(VectorXd &x) override; 
private:
	// void calc_basic_params();

	int dim_x, dim_y, dim_u;
	
	MatrixXd A;
	MatrixXd B;
	MatrixXd C;
    VectorXd x_;

    VectorXd f1, f2, f3, f4;

    // physical params
    double R; // electric resistance
    double U_max; // maximum voltage
    double L; // phase inductance
    double J; // moment of inertia of the rotor
    double Kt; // torque constant
    double Kw; // velocity constant
    double b; // motor viscous frisction constant
    double pole_pairs; // pole pairs number
    double gear_ratio; 
    // cogging torque params
    VectorXd Tk;
    VectorXd k;
    VectorXd alpha;
    double Z;

    // misc constants
    double k_phi;
    double a11, a12, a21, a22, a23, a41, a42, a43;

    double x1, x2, x3, x4, Tcog, d_Tcog, u1, u2;
    double id, iq;
    double d_id, d_iq;
    double theta;
    double d_theta;
    double dd_theta;
    double u_d, u_q;
    
};

#endif //_PMSMDynamics_hpp_