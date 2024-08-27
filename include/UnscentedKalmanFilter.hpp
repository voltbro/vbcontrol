#ifndef _UnscentedKalmanFilter_hpp_
#define _UnscentedKalmanFilter_hpp_

#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include "SysModel.hpp"

using namespace Eigen;
using namespace std;

class UnscentedKalmanFilter{
public:
    UnscentedKalmanFilter();
	UnscentedKalmanFilter(
                        // VectorXd(*)(VectorXd, VectorXd),
                        // MatrixXd(*)(VectorXd),
                        const MatrixXd &H, // measurement matrix
                        const MatrixXd &Q,
                        const MatrixXd &R,
                        const float kappa,
                        const float alpha,
                        const float beta);
    
	void set_params(
                    // Eigen::Matrix<double, 7, 1>(*dynamic_model)(Eigen::Matrix<double, 7, 1>, Eigen::Matrix<double, 2, 1>),
                    // MatrixXd(*)(VectorXd),
                    const MatrixXd &H,
                    const MatrixXd &Q,
                    const MatrixXd &R,
                    const float kappa,
                    const float alpha,
                    const float beta);

    void set_initial_estimate(const VectorXd &x0, const MatrixXd &P0);
    // void set_dynamics(VectorXd(*)(VectorXd, VectorXd));
    // void set_dynamics(SysModel &model);
    VectorXd step(SysModel &model, VectorXd &u, const VectorXd &z);
    MatrixXd get_K();

private:
    void predict(SysModel &model, VectorXd &u);
	VectorXd update(const VectorXd &z);
    void init_sigma_params();

    // Eigen::Matrix<double, 7, 1>(*dyn)(Eigen::Matrix<double, 7, 1>, Eigen::Matrix<double, 2, 1>);
    // MatrixXd(*jacobian)(VectorXd);
    // SysModel *dyn;

	MatrixXd H, Q, R, P, K, P_;
    int dim_x, dim_u, dim_z;
    size_t dim_chi;
    bool initialized1, initialized2, initialized3;  
    MatrixXd I;
    MatrixXd dfdx;
    VectorXd x_hat, x_;  
    MatrixXd X_;
    float kappa, alpha, beta;
    float lambda;
    VectorXd w_m;
    MatrixXd w_c;
    MatrixXd chi, chi_;
    VectorXd chi_i;
    MatrixXd Z;
    VectorXd mu_z;
    MatrixXd MU_Z;
    MatrixXd Pz, Pxz;
    MatrixXd L_sqrt;
    
    int N;
};

#endif //_UnscentedKalmanFilter_hpp_