
#ifndef _LinearKalmanFilter_hpp_
#define _LinearKalmanFilter_hpp_

#include <cmath>
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

class LinearKalmanFilter{
public:
    LinearKalmanFilter();
	LinearKalmanFilter(const MatrixXd  &F,
                        const MatrixXd &G,
                        const MatrixXd &H,
                        const MatrixXd &Q,
                        const MatrixXd &R);
    
	void set_params(const MatrixXd  &F,
                        const MatrixXd &G,
                        const MatrixXd &H,
                        const MatrixXd &Q,
                        const MatrixXd &R);

    void set_initial_estimate(const VectorXd &x0, const MatrixXd &P0);
    VectorXd step(const VectorXd &u, const VectorXd &z);
    MatrixXd get_K();

private:
    void predict(const VectorXd &u);
	VectorXd update(const VectorXd &z);

	MatrixXd F, G, H, Q, R, P, K, P_;
    int dim_x, dim_u, dim_z;
    bool initialized;  
    MatrixXd I;
    VectorXd x_hat, x_;  
};

#endif //_LinearKalmanFilter_hpp_