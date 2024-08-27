#include "LinearKalmanFilter.hpp"

LinearKalmanFilter::LinearKalmanFilter()
{
}

LinearKalmanFilter::LinearKalmanFilter(
                                    const MatrixXd &F,
                                    const MatrixXd &G,
                                    const MatrixXd &H,
                                    const MatrixXd &Q,
                                    const MatrixXd &R)
  : F(F), G(G), H(H), Q(Q), R(R),
    dim_x(F.rows()), dim_u(G.rows()), dim_z(H.rows()), initialized(false),
    I(dim_x, dim_x), x_hat(dim_x), x_(dim_x)
{
    
    I.setIdentity();
    x_.setZero();
    P_.setZero();
}

void LinearKalmanFilter::set_params(const MatrixXd  &F,
                        const MatrixXd &G,
                        const MatrixXd &H,
                        const MatrixXd &Q,
                        const MatrixXd &R)
{
    this->F = F;
    this->G = G;
    this->H = H;
    this->Q = Q;
    this->R = R;
    dim_x = F.rows();
    dim_u = G.rows();
    dim_z = H.rows();
    initialized = false;
    I.resize(dim_x, dim_x);
    x_hat.resize(dim_x);
    x_.resize(dim_x);

    I.setIdentity();
    x_.setZero();
    P_.setZero();
}

void LinearKalmanFilter::set_initial_estimate(const VectorXd &x0, const MatrixXd &P0)
{
    this->x_hat = x0;
    this->P = P0;
    initialized = true;
}

MatrixXd LinearKalmanFilter::get_K()
{
    return K;
}

void LinearKalmanFilter::predict(const VectorXd &u)
{
    if(!initialized)
        throw runtime_error("Linear Kalman Filter is not initialized!");

    x_ = F * x_hat + G * u;
    P_ = F * P * F.transpose() + Q;
}

VectorXd LinearKalmanFilter::update(const VectorXd &z)
{
    K = P_ * H.transpose() * (H * P_ * H.transpose() + R).inverse();
    x_hat = x_ + K * (z - H * x_);
    P = (I - K*H) * P_ * (I - K*H).transpose() + K * R * K.transpose();
    return x_hat;
}

VectorXd LinearKalmanFilter::step(const VectorXd &u, const VectorXd &z)
{
    predict(u);
    return update(z);
}
