#ifndef _PID_hpp_
#define _PID_hpp_

#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <stdexcept>
#include <iostream>

#include "vbmath.hpp"

using namespace Eigen;
using namespace std;

// PID controller class
class PID{
public:
	//constructors
    PID(); 
	PID(double Kp, 
        double Ki, 
        double Kd,
        double Ts,
        double u_min = -numeric_limits<double>::infinity(),
        double u_max = numeric_limits<double>::infinity(),
        bool use_antiwindup = true);

    void init(double Kp, 
        double Ki, 
        double Kd,
        double Ts,
        double u_min = -numeric_limits<double>::infinity(),
        double u_max = numeric_limits<double>::infinity(),
        bool use_antiwindup = true);

    double calculate(double set_val, double cur_val);
    double calculate(double set_val, double cur_val, double cur_vel);
    // double calculate(double set_val, double cur_val, double ie);
    // double calculate(double set_val, double cur_val, double de, double ie);

private:
    double _calc(double e, double de, double ie);

    bool use_antiwindup;
    bool inited;
    double Kp, Ki, Kd, Ts;
    double u_min, u_max;

    double e, de, ie, prev_e, u;
};

#endif //_PID_hpp_