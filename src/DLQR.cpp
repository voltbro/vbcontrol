/**
 * @file lqr_solve.cpp
 * @brief LQR solver for discrete time infinite horizon problems
 * @author Jonas Schlagenhauf
 * @version v0.1
 * @date 2016-09-08
 */
#include "DLQR.hpp"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <Eigen/Dense>

using Eigen::Matrix;
using Eigen::Dynamic;
using namespace std;
using namespace Eigen;

DLQR::DLQR() : eps(0.01)
{

}

/**
 * @brief Computes the LQR gain matrix (usually denoted K) for a discrete time
 * infinite horizon problem.
 *
 * @param A State matrix of the underlying system
 * @param B Input matrix of the underlying system
 * @param Q Weight matrix penalizing the state
 * @param R Weight matrix penalizing the controls
 * @param N Weight matrix penalizing state / control pairs
 * @param K Pointer to the generated matrix (has to be a double/dynamic size
 * matrix!)
 * @param eps Delta between iterations that determines when convergence is
 * reached
 */

DLQR::DLQR(MatrixXd &A, 
        MatrixXd &B,
        MatrixXd &Q,
        MatrixXd &R,
        MatrixXd &N) 
  : A(A), B(B), Q(Q), R(R), N(N), 
    dim_x(A.cols()), dim_u(B.cols()),
    eps(EPS)
{
  if (A.rows() != A.cols() || B.rows() != A.rows() || Q.rows() != Q.cols() ||
      Q.rows() != A.rows() || R.rows() != R.cols() || R.rows() != B.cols() ||
      N.rows() != A.rows() || N.cols() != B.cols()) {
    throw runtime_error("One or more matrices have incompatible dimensions. Aborting.");    
  }
}

void DLQR::set_params(MatrixXd &A, 
        MatrixXd &B,
        MatrixXd &Q,
        MatrixXd &R,
        MatrixXd &N)
{
  if (A.rows() != A.cols() || B.rows() != A.rows() || Q.rows() != Q.cols() ||
      Q.rows() != A.rows() || R.rows() != R.cols() || R.rows() != B.cols() ||
      N.rows() != A.rows() || N.cols() != B.cols()) {
    throw runtime_error("One or more matrices have incompatible dimensions. Aborting.");    
  }

  this->A = A; 
  this->B = B;
  this->Q = Q;
  this->R = R;
  this->N = N;
  eps = EPS;

  dim_x = A.cols();
  dim_u = B.cols();
}

MatrixXd DLQR::solveRiccati() { //1e-15) {
  // check if dimensions are compatible
  if (A.rows() != A.cols() || B.rows() != A.rows() || Q.rows() != Q.cols() ||
      Q.rows() != A.rows() || R.rows() != R.cols() || R.rows() != B.cols() ||
      N.rows() != A.rows() || N.cols() != B.cols()) {
    throw runtime_error("One or more matrices have incompatible dimensions. Aborting.");    
  }

  // precompute as much as possible
  MatrixXd B_T = B.transpose();
  MatrixXd Acal = A - B * R.inverse() * N.transpose();
  MatrixXd Acal_T = Acal.transpose();
  MatrixXd Qcal = Q - N * R.inverse() * N.transpose();

  // initialize P with Q
  MatrixXd P = Q;

  // iterate until P converges
  unsigned int numIterations = 0;
  MatrixXd Pold = P;
  while (true) {
    numIterations++;
    
    // compute new P
    P = Acal_T * P * Acal -
        Acal_T * P * B * (R + B_T * P * B).inverse() * B_T * P * Acal + Qcal;

    // update delta
    MatrixXd delta = P - Pold;
    // std::cout << "delta: " << delta.maxCoeff() << std::endl;
    if (fabs(delta.maxCoeff()) < eps) {
      std::cout << "Number of iterations until convergence: " << numIterations
                << std::endl;
      break;
    }
    Pold = P;
  }

  // compute K from P
  K = (R + B_T * P * B).inverse() * (B_T * P * A + N.transpose());

  return K;
}

/**
 * @brief Little piece of code to test the solver
 *
 * @return 1 if something went wrong, 0 otherwise.
 */
// int main() {
//   MatrixXd A;
//   MatrixXd B;
//   MatrixXd Q;
//   MatrixXd R;
//   MatrixXd N;
//   MatrixXd K;

//   int dim_x = 7;
//   int dim_u = 2;

//   A.resize(dim_x, dim_x);
//   B.resize(dim_x, dim_u);
//   Q.resize(dim_x, dim_x);
//   R.resize(dim_u, dim_u);
//   N.resize(dim_x, dim_u);

//   A << 1.00000000e+00,  2.00000000e-03, -2.13027421e-07, 1.99771260e-06,  2.58398739e-10,  0.00000000e+00, 0.00000000e+00,
//        0.00000000e+00,  1.00000000e+00, -2.13027421e-04, 1.99771260e-03,  2.58398739e-07,  0.00000000e+00,  0.00000000e+00,
//        0.00000000e+00,  0.00000000e+00,  1.00013705e+00, 7.19025687e-07,  1.99995919e-03,  0.00000000e+00, 0.00000000e+00,
//        0.00000000e+00,  0.00000000e+00, -2.13027421e-01, 9.97712604e-01,  2.58398739e-04,  0.00000000e+00, 0.00000000e+00,
//        0.00000000e+00,  0.00000000e+00,  1.37052557e-01, 7.19025687e-04,  9.99959188e-01,  0.00000000e+00, 0.00000000e+00,
//        0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00,  0.00000000e+00,  1.00000000e+00, 1.99805084e-03,
//        0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 9.98050837e-01;

//   B << 4.49997698e-08,  4.49997698e-08,
//        4.49997698e-05,  4.49997698e-05,
//        -1.69779492e-05, -1.69779492e-05,
//         4.49997698e-02,  4.49997698e-02,
//        -1.69779492e-02, -1.69779492e-02,
//        -1.33256677e-05,  1.33256677e-05,
//        -1.33256677e-02,  1.33256677e-02;
//   // A << -0.00000000e+00,  1.00000000e+00,  0.00000000e+00, 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00,
//   //       0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 1.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00,
//   //       0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00,  1.00000000e+00,  0.00000000e+00, 0.00000000e+00,
//   //       0.00000000e+00,  0.00000000e+00, -1.06644538e+02, -1.14505385e+00,  2.35994487e-01,  0.00000000e+00, 0.00000000e+00,
//   //       0.00000000e+00,  0.00000000e+00,  6.85660176e+01,  3.59931850e-01, -8.90187746e-02,  0.00000000e+00, 0.00000000e+00,
//   //       0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 1.00000000e+00,
//   //       0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00, -9.75532365e-01;

//   // B << 0.0        ,  0.0        ,
//   //       0.0        ,  0.0        ,
//   //       0.0        ,  0.0        ,
//   //      22.52674653, 22.52674653,
//   //      -8.49724666, -8.49724666,
//   //       0.0        ,  0.0        ,
//   //      -6.66933367,  6.66933367;

//   Q <<  1,     0,     0,     0,     0,     0,     0,
//         0,     1,     0,     0,     0,     0,     0,
//         0,     0,     10,    0,     0,     0,     0,
//         0,     0,     0,     1,     0,     0,     0,
//         0,     0,     0,     0,     1,     0,     0,
//         0,     0,     0,     0,     0,     1000,  0,
//         0,     0,     0,     0,     0,     0,     100;

//   R << 20000, 0,
//        0,     20000;
  
//   DLQR dlqr;
//   dlqr.set_params(A, B, Q, R, N);

//   K = dlqr.solveRiccati();

//   Eigen::IOFormat fmt(4, 0, ", ", "\n", "[", "]");
//   std::cout << "K: " << K.format(fmt) << std::endl;

//   return 0;
// }
