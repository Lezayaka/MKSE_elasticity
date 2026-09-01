cmake_minimum_required(VERSION 3.20)

project(
    mkse_elasticity
    VERSION 1.0.0
    DESCRIPTION "Finite superelement solver for elastic contact problems"
    LANGUAGES CXX
)

add_executable(
    mkse-elasticity
    src/main.cpp
    src/FEM.cpp
    src/FSEM.cpp
    src/LinearAlgebra.cpp
)

target_include_directories(
    mkse-elasticity
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_compile_features(mkse-elasticity PRIVATE cxx_std_20)

if(MSVC)
    target_compile_options(mkse-elasticity PRIVATE /W4 /permissive-)
else()
    target_compile_options(mkse-elasticity PRIVATE -Wall -Wextra -Wpedantic)
endif()
