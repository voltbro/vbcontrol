
#ifndef _SysModel_hpp_
#define _SysModel_hpp_

#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <string>

using namespace Eigen;
 
// Self-Balancing robot dynamics class
class SysModel{
public:
	//constructors
	SysModel() {}
	SysModel(double Ts): Ts(Ts) {}
	//functions
	void set_timestep(double Ts) // in ms
	{
		this->Ts = Ts;
	}
	virtual std::string get_type() {return "SysModel";}
	// virtual void set_physical_params() = 0
	virtual void ss_model(MatrixXd &A, MatrixXd &B, MatrixXd &C) {
		A.resize(2,2);
		B.resize(2,2);
		C.resize(2,2);

		A << 0, 1,
		     2, 3;

		B << 0, 1,
		     2, 3;

		C << 0, 1,
		     2, 3;
	}
	virtual VectorXd nonlinear_dynamics(VectorXd &x, VectorXd &u) {return x+u;}
	virtual MatrixXd jacobian(VectorXd &x) {
		MatrixXd A(2,2);
		A << x(0), 1,
		     2,    3;
		return A;
		}

	double Ts;
private:
		
};

#endif //_SysModel_hpp_