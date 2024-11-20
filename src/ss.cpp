#include "ss.hpp"

std::tuple<MatrixXd, MatrixXd> ss::c2d(MatrixXd A, MatrixXd B, double Ts)
    {
        int dim_x = A.rows();
        int dim_u = B.cols();

        MatrixXd I(dim_x, dim_x);
        I.setIdentity();
        MatrixXd Ad(dim_x, dim_x);
        MatrixXd Bd(dim_x, dim_u);

        Ad = I + A * Ts;
        Bd = B * Ts;

        return std::make_tuple(Ad, Bd);
    }