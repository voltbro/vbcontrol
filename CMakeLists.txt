cmake_minimum_required(VERSION 3.10)

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
# find_package(pybind11 REQUIRED)

include_directories("include/")
set(SOURCES src/main.cpp 
            src/LinearKalmanFilter.cpp 
            src/ExtendedKalmanFilter.cpp 
            src/UnscentedKalmanFilter.cpp 
            src/LowPassFilter.cpp 
            src/GyrobroDynamics.cpp 
            src/DCMotorDynamics.cpp
            src/PID.cpp 
            src/FullStateController.cpp
            src/DLQR.cpp)
# add the executable
# add_executable(vbctrl ${SOURCES})

add_library(vbcontrol STATIC ${SOURCES})

# find_package(pybind11 REQUIRED)
add_subdirectory(pybind11)
pybind11_add_module(vbcontrolpy "./src/pybingings.cpp")
target_link_libraries(vbcontrolpy PUBLIC vbcontrol)
