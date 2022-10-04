cmake_minimum_required(VERSION 3.10)
project(Proyecto1_Gravedad)
set(CMAKE_CXX_STANDARD 17)

SET(CMAKE_BUILD_TYPE Release)
SET(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -frecord-gcc-switches -march=haswell -pipe -fomit-frame-pointer -mmmx")
find_package(OpenMP REQUIRED)

add_executable(sim-soa-s sim-soa-s.cpp)
add_executable(sim-aos-s sim-aos-s.cpp)
add_executable(sim-soa-p sim-soa-p.cpp)
add_executable(sim-aos-p sim-aos-p.cpp)

target_link_libraries (sim-aos-p PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries (sim-soa-p PUBLIC OpenMP::OpenMP_CXX)