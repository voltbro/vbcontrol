#include "PID.hpp"


PID::PID()
{
    ie = 0;
    inited = false;
}

PID::PID(double Kp, 
        double Ki, 
        double Kd,
        double Ts,
        double u_min/*  = -numeric_limits<double>::infinity() */,
        double u_max/*  = numeric_limits<double>::infinity() */,
        bool use_antiwindup/*  = true */)
{
    this->Kp = Kp;
    this->Ki = Ki;
    this->Kd = Kd;
    this->Ts = Ts;
    this->u_min = u_min;
    this->u_max = u_max;
    this->use_antiwindup = use_antiwindup;

    ie = 0;

    inited = true;
}

void PID::init(double Kp, 
        double Ki, 
        double Kd,
        double Ts,
        double u_min/*  = -numeric_limits<double>::infinity() */,
        double u_max/*  = numeric_limits<double>::infinity() */,
        bool use_antiwindup/*  = true */)
{
    this->Kp = Kp;
    this->Ki = Ki;
    this->Kd = Kd;
    this->Ts = Ts;
    this->u_min = u_min;
    this->u_max = u_max;
    this->use_antiwindup = use_antiwindup;

    inited = true;
}

double PID::_calc(const double e, const double de, const double ie)
{
    double out = Kp * e + Ki * ie + Kd * de;
    
    out = vbmath::clip(out, u_min, u_max);

    return out;
}

double PID::calculate(double set_val, double cur_val)
{
    if (inited == false)
        throw std::runtime_error("PID Controller is not inited! Use PID::init method.");

    e = set_val - cur_val;
    
    de = (e - prev_e)/Ts;

    if (use_antiwindup)
    {
        if ((e > 0 && u >= u_max) || (e <= 0 && u <= u_min))
        {
            ie += 0;
            // cout << "yo" << endl;
        }
        else    
            ie += (e*Ts);
    }
    
    
    u = _calc(e, de, ie);

    prev_e = e;

    return u;
}

double PID::calculate(double set_val, double cur_val, double cur_vel)
{
    if (inited == false)
        throw std::runtime_error("PID Controller is not inited! Use PID::init method.");

    e = set_val - cur_val;
    
    if (use_antiwindup)
    {
        if ((e > 0 && u >= u_max) || (e <= 0 && u <= u_min))
            ie += 0;
        else    
            ie += (e*Ts);
    }
    
    u = _calc(e, -cur_vel, ie);

    prev_e = e;

    return u;
}