#ifndef _ss_hpp_
#define _ss_hpp_

#include <Eigen/Dense>
#include <tuple>

using namespace Eigen;

namespace ss
{
    // convert continuous model to discrete
    std::tuple<MatrixXd, MatrixXd> c2d(MatrixXd A, MatrixXd B, double Ts);
    // {
    //     int dim_x = A.rows();
    //     int dim_u = B.cols();

    //     MatrixXd I(dim_x, dim_x);
    //     I.setIdentity();
    //     MatrixXd Ad(dim_x, dim_x);
    //     MatrixXd Bd(dim_x, dim_u);

    //     Ad = I + A * Ts;
    //     Bd = B * Ts;

    //     return std::make_tuple(Ad, Bd);
    // }
}

#endif //_ss_hpp_