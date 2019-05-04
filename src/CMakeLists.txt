cmake_minimum_required(VERSION 3.1)
project(diffflame)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 11)

find_package(Threads REQUIRED)

include_directories("/usr/local/include" "/usr/local/opt/lapack/include")
link_directories("/usr/local/lib" "/usr/local/opt/lapack/lib")

add_executable(diffflame main.cpp)
target_link_libraries(diffflame cantera sundials_cvodes sundials_ida sundials_nvecserial fmt lapack blas Threads::Threads)

