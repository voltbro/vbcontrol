
#ifndef _ExtendedKalmanFilter_hpp_
#define _ExtendedKalmanFilter_hpp_

#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include "SysModel.hpp"

using namespace Eigen;
using namespace std;

class ExtendedKalmanFilter{
public:
    ExtendedKalmanFilter();
	ExtendedKalmanFilter(
                        // VectorXd(*)(VectorXd, VectorXd),
                        // MatrixXd(*)(VectorXd),
                        const MatrixXd &H,
                        const MatrixXd &Q,
                        const MatrixXd &R);
    
	void set_params(
                    // Eigen::Matrix<double, 7, 1>(*dynamic_model)(Eigen::Matrix<double, 7, 1>, Eigen::Matrix<double, 2, 1>),
                    // MatrixXd(*)(VectorXd),
                    const MatrixXd &H,
                    const MatrixXd &Q,
                    const MatrixXd &R);

    void set_initial_estimate(const VectorXd &x0, const MatrixXd &P0);
    // void set_dynamics(VectorXd(*)(VectorXd, VectorXd));
    // void set_dynamics(SysModel &model);
    VectorXd step(SysModel &model, VectorXd &u, const VectorXd &z);
    MatrixXd get_K();

private:
    void predict(SysModel &model, VectorXd &u);
	VectorXd update(const VectorXd &z);

    // Eigen::Matrix<double, 7, 1>(*dyn)(Eigen::Matrix<double, 7, 1>, Eigen::Matrix<double, 2, 1>);
    // MatrixXd(*jacobian)(VectorXd);
    // SysModel dyn;

	MatrixXd H, Q, R, P, K, P_;
    int dim_x, dim_u, dim_z;
    bool initialized;  
    MatrixXd I;
    MatrixXd dfdx;
    VectorXd x_hat, x_;  
};

#endif //_ExtendedKalmanFilter_hpp_