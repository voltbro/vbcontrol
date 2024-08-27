#include "ExtendedKalmanFilter.hpp"

ExtendedKalmanFilter::ExtendedKalmanFilter()
{

}

ExtendedKalmanFilter::ExtendedKalmanFilter(
                                    // VectorXd(*dynamics_model)(VectorXd, VectorXd),
                                    // MatrixXd(*jacobian)(VectorXd),
                                    const MatrixXd &H,
                                    const MatrixXd &Q,
                                    const MatrixXd &R)
  : H(H), Q(Q), R(R),
    dim_x(H.cols()), dim_z(H.rows()), initialized(false),
    I(dim_x, dim_x), x_hat(dim_x), x_(dim_x)
{
    if (Q.cols() != dim_x || Q.rows() != dim_x)
        throw runtime_error("ExtendedKalmanFilter: Q dimensions doesn't equal to x!");
    if (R.cols() != dim_z || R.rows() != dim_z)
        throw runtime_error("ExtendedKalmanFilter: R dimensions doesn't equal to z!");

    I.setIdentity();
    x_.setZero();
    P_.setZero();
}

void ExtendedKalmanFilter::set_params(
                        const MatrixXd &H,
                        const MatrixXd &Q,
                        const MatrixXd &R)
{
    this->H = H;
    this->Q = Q;
    this->R = R;
    dim_x = H.cols();
    dim_z = H.rows();
    initialized = false;
    I.resize(dim_x, dim_x);
    x_hat.resize(dim_x);
    x_.resize(dim_x);

    I.setIdentity();
    x_.setZero();
    P_.setZero();
}

// void ExtendedKalmanFilter::set_dynamics(SysModel &model)
// {
//     // model.get_type();
//     dyn = model;  
// }

void ExtendedKalmanFilter::set_initial_estimate(const VectorXd &x0, const MatrixXd &P0)
{   
    if (x0.rows() != dim_x)
        throw runtime_error("ExtendedKalmanFilter: x0 dimensions doesn't equal to x!");
    if (P0.cols() != dim_x || P0.rows() != dim_x)
        throw runtime_error("ExtendedKalmanFilter: P0 dimensions doesn't equal to x!");

    this->x_hat = x0;
    this->P = P0;
    initialized = true;
}

MatrixXd ExtendedKalmanFilter::get_K()
{
    return K;
}

void ExtendedKalmanFilter::predict(SysModel &model, VectorXd &u)
{
    if(!initialized)
        throw runtime_error("Extended Kalman Filter is not initialized!");
    
    // cout << "It is " << model.get_type() << endl;
    x_ = model.nonlinear_dynamics(x_hat, u);
    // cout << "2" << endl;
    dfdx = model.jacobian(x_hat);
    // cout << "3" << endl;
    P_ = dfdx * P * dfdx.transpose() + Q;
    // cout << "4" << endl;
}

VectorXd ExtendedKalmanFilter::update(const VectorXd &z)
{
    K = P_ * H.transpose() * (H * P_ * H.transpose() + R).inverse();
    x_hat = x_ + K * (z - H * x_);
    P = (I - K*H) * P_ * (I - K*H).transpose() + K * R * K.transpose();
    return x_hat;
}

VectorXd ExtendedKalmanFilter::step(SysModel &model, VectorXd &u, const VectorXd &z)
{
    // cout << "It is" << model.get_type() << endl;
    predict(model, u);
    // cout << "yo" << endl;
    return update(z);
}
