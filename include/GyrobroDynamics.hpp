
#ifndef _GyrobroDynamics_hpp_
#define _GyrobroDynamics_hpp_
 
#define DIM_X_GB 7
#define DIM_U_GB 2
#define DIM_Y_GB 6

#include <cmath>
#include <Eigen/Dense>
#include <string>

#include "SysModel.hpp"

using namespace Eigen;

// Self-Balancing robot dynamics class
class GyrobroDynamics : public SysModel 
{
public: 
	//constructors
	GyrobroDynamics();
	GyrobroDynamics(double Ts);
	std::string get_type() override {return "GyrobroDynamics";}
	//functions
	// void set_timestep(double Ts); // in ms
	void set_physical_params(double g,
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
							double n);
	void ss_model(MatrixXd &A, MatrixXd &B, MatrixXd &C) override;
	VectorXd nonlinear_dynamics(VectorXd &x, VectorXd &u) override;
	MatrixXd jacobian(VectorXd &x) override; 
private:
	void calc_basic_params();

	int dim_x, dim_y, dim_u;
	// double Ts;
	double g;
	double m;
	double hw;
	double R;
	double M;
	double W;
	double D;
	double h;
	double fm;
	double fw;
	double M_stator;
	double R_stator;
	double hw_stator;
	double Rm;
	double Jm;
	double Kt;
	double Kb;
	double n;

	double Jw;
	double Jwx;
	double L;
	double Jpsi;
	double Jphi;

	double alpha, beta;
	double F_11, F_12, F_21, F_22;
	double E_11, E_12, E_21, E_22;
	double detE;
	double I, J, K;
	double LMG, L2M, LMR;

	MatrixXd E, F, G, H, invE;//Matrix <double, 2, 2>

	MatrixXd A;//Matrix <double, 7, 7>
	MatrixXd B;//Matrix <double, 7, 2>
	MatrixXd C;//Matrix <double, 5, 7>
	// MatrixXd eye;//Matrix <double, 7, 7>
};

#endif //_GyrobroDynamics_hpp_