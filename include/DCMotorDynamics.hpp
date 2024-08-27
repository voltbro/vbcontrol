
#ifndef _DCMotorDynamics_hpp_
#define _DCMotorDynamics_hpp_

#define DIM_X_DC 3
#define DIM_U_DC 1
#define DIM_Y_DC 2

#include <cmath>
#include <Eigen/Dense>
 
#include "SysModel.hpp"
 
using namespace Eigen;
using namespace std;

// Self-Balancing robot dynamics class
class DCMotorDynamics : public SysModel 
{
public: 
	//constructors
	DCMotorDynamics();
	DCMotorDynamics(double Ts); // Ts - timestep [ms]
	//functions
	void set_physical_params(double b, // motor viscous frisction constant
                                double J, // moment of inertia of the rotor
                                double K, // motor constant
                                double R, // electric resistance
                                double L // electric inductance
                                );
	void ss_model(MatrixXd &A, MatrixXd &B, MatrixXd &C) override;
	VectorXd nonlinear_dynamics(VectorXd &x, VectorXd &u) override;
	MatrixXd jacobian(VectorXd &x) override; 
private:
	// void calc_basic_params();

	int dim_x, dim_y, dim_u;
	
	MatrixXd A;
	MatrixXd B;
	MatrixXd C;
    VectorXd x_;

    double b; // motor viscous frisction constant
    double J; // moment of inertia of the rotor
    double K; // motor constant
    double R; // electric resistance
    double L; // electric inductance

    double i;
    double didt;
    double theta;
    double d_theta;
    double dd_theta;
    double V;
    

};

#endif //_LowPassFilter_hpp_