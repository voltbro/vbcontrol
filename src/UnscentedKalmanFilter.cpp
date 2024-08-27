 #include "UnscentedKalmanFilter.hpp"

UnscentedKalmanFilter::UnscentedKalmanFilter() 
    : initialized1(false), initialized2(false), initialized3(false)
{
    
}

UnscentedKalmanFilter::UnscentedKalmanFilter(
                                    // VectorXd(*dynamics_model)(VectorXd, VectorXd),
                                    // MatrixXd(*jacobian)(VectorXd),
                                    const MatrixXd &H,
                                    const MatrixXd &Q,
                                    const MatrixXd &R,
                                    const float kappa,
                                    const float alpha,
                                    const float beta)
  : H(H), Q(Q), R(R),
    dim_x(H.cols()), dim_z(H.rows()), 
    initialized1(true), initialized2(false), initialized3(false),
    I(dim_x, dim_x), x_hat(dim_x), x_(dim_x),
    kappa(kappa), alpha(alpha), beta(beta)
{
    I.setIdentity();
    x_.setZero();
    P_.setZero();

    init_sigma_params();
}

void UnscentedKalmanFilter::set_params(
                        const MatrixXd &H,
                        const MatrixXd &Q,
                        const MatrixXd &R,
                        const float kappa,
                        const float alpha,
                        const float beta)
{
    this->H = H;
    this->Q = Q;
    this->R = R;
    this->alpha = alpha;
    this->beta = beta;
    this->kappa = kappa;
    dim_x = H.cols();
    dim_z = H.rows();

    
    I.resize(dim_x, dim_x);
    x_hat.resize(dim_x);
    x_.resize(dim_x);

    I.setIdentity();
    x_.setZero();
    P_.setZero();

    init_sigma_params();

    initialized1 = true;
}

// void UnscentedKalmanFilter::set_dynamics(SysModel &model)
// {
//     dyn = &model;  
//     initialized2 = true;
// }

void UnscentedKalmanFilter::set_initial_estimate(const VectorXd &x0, const MatrixXd &P0)
{
    this->x_hat = x0;
    this->P = P0;
    initialized3 = true;
}

void UnscentedKalmanFilter::init_sigma_params()
{
    N = dim_x;
    dim_chi = size_t(2*N+1);

    lambda = pow(alpha, 2) * (N + kappa) - N;
    // cout << lambda << endl;
    
    float w0_m = lambda/(N + lambda);
    // cout << w0_m << endl;

    float w0_c = lambda/(N + lambda) + 1 - pow(alpha, 2) + beta;
    float wi = 1/(2*(N + lambda));
    // cout << w0_c << endl;
    // cout << wi << endl;

    w_m.resize(2*N+1);
    w_m(0) = w0_m;

    VectorXd w_c_diag(2*N+1);
    w_c_diag(0) = w0_c;
    for (size_t i = 1; i < dim_chi; i++)
    {
        w_m(i) = wi;
        w_c_diag(i) = wi;
    }
    
    w_c.resize(2*N+1, 2*N+1);
    w_c = w_c_diag.array().matrix().asDiagonal();

    // cout << w_m << endl;
    // cout << w_c << endl;

    
    chi.resize(dim_x, 2*N+1);
    chi_.resize(dim_x, 2*N+1);
    Z.resize(dim_z, 2*N+1);
    mu_z.resize(dim_z);
    Pz.resize(dim_z, dim_z);
    X_.resize(dim_x, dim_chi);
    MU_Z.resize(dim_z, dim_chi);

}

MatrixXd UnscentedKalmanFilter::get_K()
{
    return K;
}

void UnscentedKalmanFilter::predict(SysModel &model, VectorXd &u)
{
    // sigma transform
    chi.col(0) = x_hat;
    L_sqrt = ((N + lambda) * P).llt().matrixL();
    for (size_t i = 1; i < dim_chi; i++)
    {
        if (i <= size_t(N))
            chi.col(i) = x_hat + L_sqrt.col(i-1);
        else
            chi.col(i) = x_hat - L_sqrt.col(i-N-1);
    }

    // propagation through non-linear function
    for (size_t i = 0; i < dim_chi; i++)
    {
        chi_i = chi.col(i);
        chi_.col(i) = model.nonlinear_dynamics(chi_i, u);
    }

    // mean state and covariance computation
    x_ = chi_ * w_m;
    for (size_t i = 0; i < dim_chi; i++)
        X_.col(i) = x_;
    P_ = (chi_ - X_) * w_c * ((chi_ - X_).transpose()) + Q;
}

VectorXd UnscentedKalmanFilter::update(const VectorXd &z)
{
    // Propagate through observation matrix
    for (size_t i = 0; i < dim_chi; i++)
        Z.col(i) = H * chi_.col(i);

    // Mean and cross-covariance computation
    mu_z = Z * w_m;
    for (size_t i = 0; i < dim_chi; i++)
        MU_Z.col(i) = mu_z;
    Pz = (Z - MU_Z) * w_c * ((Z - MU_Z).transpose()) + R;
    Pxz = (chi_ - X_) * w_c * ((Z - MU_Z).transpose());

    // Kalman Gain computation
    K = Pxz * Pz.inverse();

    // Estimate state and new covariance
    x_hat = x_ + K * (z - mu_z);
    P = P_ - K * Pz * (K.transpose());

    return x_hat;
}

VectorXd UnscentedKalmanFilter::step(SysModel &model, VectorXd &u, const VectorXd &z)
{
    if(!(initialized1 & initialized3))
        throw runtime_error("Unscented Kalman Filter is not initialized!");

    predict(model, u);
    return update(z);
}
