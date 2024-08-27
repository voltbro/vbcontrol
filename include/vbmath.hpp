#ifndef _vbmath_hpp_
#define _vbmath_hpp_

#include <algorithm>
#include <math.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

namespace vbmath
{
    template <typename T>
    inline int sign(T val)
    {
        return (T(0) < val) - (val < T(0));
    }

    template <typename T, typename A>
    inline T clip(T n, A lower, A upper) 
    {
        return std::max(lower, std::min(n, upper));
    }

    template <typename A>
    inline VectorXd clip(VectorXd n,  A lower, A upper) 
    {
        int dim = n.rows();
        // if (dim != lower.rows() || dim != upper.rows())
        //     throw std::runtime_error("Diffrent dimensions of n, lower and upper vectors!");

        VectorXd out(dim);
        for(int i=0; i<dim; i++)
        {
            out(i) = std::max(lower, std::min(A(n(i)), upper));
        }

        return out;
    }

    template <>
    inline VectorXd clip(VectorXd n, VectorXd lower, VectorXd upper) 
    {
        int dim = n.rows();
        if (dim != lower.rows() || dim != upper.rows())
            throw std::runtime_error("Diffrent dimensions of n, lower and upper vectors!");
        // int dim = 5;
        VectorXd out(dim);
        for(int i=0; i<dim; i++)
        {
            out(i) = std::max(lower(i), std::min(n(i), upper(i)));
        }

        return out;
    }

    template <typename T>
    inline T fric_comp(T x, T gain, T offset)
    {
        return sign(x) * (gain * abs(x) + offset);
    }

}

#endif //_vbmath_hpp_