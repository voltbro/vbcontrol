#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
// #include <Eigen/Dense>
namespace py = pybind11;

#include "vbmath.hpp"
#include "DCMotorDynamics.hpp"
#include "SysModel.hpp"
#include "ss.hpp"
#include "PID.hpp"
#include "DLQR.hpp"
#include "GyrobroDynamics.hpp"
#include "ExtendedKalmanFilter.hpp"
#include "LinearKalmanFilter.hpp"
#include "UnscentedKalmanFilter.hpp"
#include "FullStateController.hpp"
#include "LowPassFilter.hpp"

using namespace vbmath;
using namespace pybind11::literals;
using namespace Eigen;

PYBIND11_MODULE(vbcontrolpy, m) 
{
    auto math = m.def_submodule("math");
    math.def("sign", &sign<int>, "This function returns sign of the number");
    math.def("sign", &sign<double>, "This function returns sign of the number");
    math.def("sign", &sign<float>, "This function returns sign of the number");
    math.def("fric_comp", &fric_comp<double>, "Friction compensator");
    math.def("fric_comp", &fric_comp<float>, "Friction compensator");
    math.def("clip", py::overload_cast<float, float, float>(&clip<float, float>), "Clip number");
    math.def("clip", py::overload_cast<double, double, double>(&clip<double, double>), "Clip number");
    math.def("clip", py::overload_cast<int, int, int>(&clip<int, int>), "Clip number");
    math.def("clip", py::overload_cast<VectorXd, float, float>(&clip<float>), "Clip number");
    math.def("clip", py::overload_cast<VectorXd, double, double>(&clip<double>), "Clip number");
    math.def("clip", py::overload_cast<VectorXd, int, int>(&clip<int>), "Clip number");
    math.def("clip", py::overload_cast<VectorXd, VectorXd, VectorXd>(&clip<VectorXd>), "Clip number");


    auto ss = m.def_submodule("ss");
    ss.def("c2d", &ss::c2d);


    py::class_<PID>(m, "pid")
        .def(py::init<>())
        .def(py::init<double, double, double, double, double, double, bool>())
        .def("init", &PID::init)
        .def("calculate", py::overload_cast<double, double>(&PID::calculate))
        .def("calculate", py::overload_cast<double, double, double>(&PID::calculate));

    
    py::class_<DLQR>(m, "dlqr")
        .def(py::init<>())
        .def(py::init<MatrixXd&, MatrixXd&, MatrixXd&, MatrixXd&, MatrixXd&>())
        .def("set_params", &DLQR::set_params)
        .def("solveRiccati", &DLQR::solveRiccati);

    
    py::class_<FullStateController>(m, "full_state_controller")
        .def(py::init<>())
        .def(py::init<MatrixXd, VectorXd, VectorXd>())
        .def("init", &FullStateController::init)
        .def("calculate", &FullStateController::calculate);

    
    py::class_<ExtendedKalmanFilter>(m, "ekf")
        .def(py::init<>())
        .def(py::init<MatrixXd&, MatrixXd&, MatrixXd&>())
        .def("set_params", &ExtendedKalmanFilter::set_params)
        .def("set_initial_estimate", &ExtendedKalmanFilter::set_initial_estimate)
        .def("step", &ExtendedKalmanFilter::step)
        .def("get_K", &ExtendedKalmanFilter::get_K);


    py::class_<LinearKalmanFilter>(m, "lkf")
        .def(py::init<>())
        .def(py::init<MatrixXd&, MatrixXd&, MatrixXd&, MatrixXd&, MatrixXd&>())
        .def("set_params", &LinearKalmanFilter::set_params)
        .def("set_initial_estimate", &LinearKalmanFilter::set_initial_estimate)
        .def("step", &LinearKalmanFilter::step)
        .def("get_K", &LinearKalmanFilter::get_K);

    
    py::class_<UnscentedKalmanFilter>(m, "ukf")
        .def(py::init<>())
        .def(py::init<MatrixXd&, MatrixXd&, MatrixXd&, float, float, float>())
        .def("set_params", &UnscentedKalmanFilter::set_params)
        .def("set_initial_estimate", &UnscentedKalmanFilter::set_initial_estimate)
        .def("step", &UnscentedKalmanFilter::step)
        .def("get_K", &UnscentedKalmanFilter::get_K);


    py::class_<LowPassFilter>(m, "low_pass_filter")
        .def(py::init<>())
        .def(py::init<float, float>())
        .def("getOutput", &LowPassFilter::getOutput)
        .def("reconfigureFilter", &LowPassFilter::reconfigureFilter)
        .def("update", py::overload_cast<float>(&LowPassFilter::update))
        .def("update", py::overload_cast<float, float, float>(&LowPassFilter::update));


    auto dyn = m.def_submodule("dynamics");
    py::class_<SysModel>(dyn, "sys_model");
    py::class_<DCMotorDynamics, SysModel>(dyn, "dc_motor")
        .def(py::init<>())
        .def(py::init<double>())
        .def("set_timestep", &DCMotorDynamics::set_timestep)
        .def("set_physical_params", &DCMotorDynamics::set_physical_params)
        .def("ss_model", &DCMotorDynamics::ss_model)
        .def("nonlinear_dynamics", &DCMotorDynamics::nonlinear_dynamics)
        .def("jacobian", &DCMotorDynamics::jacobian);

    py::class_<GyrobroDynamics, SysModel>(dyn, "gyrobro")
        .def(py::init<>())
        .def(py::init<double>())
        .def("set_timestep", &GyrobroDynamics::set_timestep)
        .def("set_physical_params", &GyrobroDynamics::set_physical_params)
        .def("ss_model", &GyrobroDynamics::ss_model)
        .def("nonlinear_dynamics", &GyrobroDynamics::nonlinear_dynamics)
        .def("jacobian", &GyrobroDynamics::jacobian);

}
