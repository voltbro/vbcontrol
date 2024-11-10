cmake_minimum_required(VERSION 3.14)

# specify the C++ standard

set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall  -O3 -march=native ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -fPIC")
# add_definitions(-DCOMPILEDWITHC11)
message(STATUS "Using flag -std=c++17")

# set the project name
project(VBControl VERSION 1.0)

find_package(Eigen3 REQUIRED)
# find_package(matplotplusplus REQUIRED)
# find_package(pybind11 REQUIRED)

include_directories("include/")

set(SOURCES src/main.cpp 
            src/LinearKalmanFilter.cpp 
            src/ExtendedKalmanFilter.cpp 
            src/UnscentedKalmanFilter.cpp 
            src/LowPassFilter.cpp 
            src/GyrobroDynamics.cpp 
            src/DCMotorDynamics.cpp
            src/PMSMDynamics.cpp
            src/PMSMCascadeFOC.cpp
            src/PID.cpp 
            src/FullStateController.cpp
            src/DLQR.cpp)

# executable
add_subdirectory(matplotplusplus)
add_executable(vbctrl ${SOURCES})
target_link_libraries(vbctrl PUBLIC matplot)

# library
add_library(vbcontrol STATIC ${SOURCES})
target_link_libraries(vbcontrol PUBLIC matplot)

# python wrap
add_subdirectory(pybind11)
pybind11_add_module(vbcontrolpy "./src/pybingings.cpp")
target_link_libraries(vbcontrolpy PUBLIC vbcontrol matplot)
