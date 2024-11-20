#ifndef _DLQR_hpp_
#define _DLQR_hpp_

#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

#define EPS 0.05

// LQR controller class
class DLQR{
public:
    DLQR();
    DLQR(MatrixXd &A, 
        MatrixXd &B,
        MatrixXd &Q,
        MatrixXd &R,
        MatrixXd &N);

    void set_params(MatrixXd &A, 
        MatrixXd &B,
        MatrixXd &Q,
        MatrixXd &R,
        MatrixXd &N);

    MatrixXd solveRiccati();

private:
    MatrixXd A; 
    MatrixXd B;
    MatrixXd Q;
    MatrixXd R;
    MatrixXd N;
    MatrixXd K;

    int dim_x, dim_u;

    double eps;
};

#endif //_DLQR_hpp_