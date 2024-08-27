#ifndef _FullStateController_hpp_
#define _FullStateController_hpp_

#include <cmath>
#include <Eigen/Dense>
#include <stdexcept>
#include <iostream>

#include "vbmath.hpp"

using namespace Eigen;
using namespace std;

// PID controller class
class FullStateController
{
public:
	//constructors
    FullStateController();
	FullStateController(MatrixXd K, 
        VectorXd u_min,
        VectorXd u_max);

    void init(MatrixXd K, 
        VectorXd u_min,
        VectorXd u_max);

    VectorXd calculate(VectorXd set_val, VectorXd cur_val);

private:
    bool inited;
    MatrixXd K;
    VectorXd u_min, u_max;

    VectorXd e, u;

    int dim_x, dim_u;
};

#endif //_FullStateController_hpp_