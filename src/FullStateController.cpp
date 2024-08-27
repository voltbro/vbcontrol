#include "FullStateController.hpp"


FullStateController::FullStateController()
{
    inited = false;
}

FullStateController::FullStateController(MatrixXd K, 
                                        VectorXd u_min,
                                        VectorXd u_max)
{
    this->K = K;
    this->u_min = u_min;
    this->u_max = u_max;

    dim_x = K.cols();
    dim_u = K.rows();

    if (u_min.rows() != dim_u || u_max.rows() != dim_u)
        throw std::runtime_error("Incorrect dimensions of u_max, u_min or K vectors! ");

    inited = true;
}

void FullStateController::init(MatrixXd K, 
                                VectorXd u_min,
                                VectorXd u_max)
{
    this->K = K;
    this->u_min = u_min;
    this->u_max = u_max;

    if (u_min.rows() != dim_u || u_max.rows() != dim_u)
        throw std::runtime_error("Incorrect dimensions of u_max, u_min or K vectors! ");

    inited = true;
}

VectorXd FullStateController::calculate(VectorXd set_val, VectorXd cur_val)
{
    if (inited == false)
        throw std::runtime_error("FullStateController is not inited! Use FullStateController::init method.");

    if (dim_x != set_val.rows() || dim_x != cur_val.rows())
        throw std::runtime_error("Incorrect dimensions of set_val or cure_val vectors!");

    e = set_val - cur_val;
    u = K * e;
    u = vbmath::clip(u, u_min, u_max);

    return u;
}
