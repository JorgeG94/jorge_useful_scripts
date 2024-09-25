cmake_minimum_required(VERSION 3.18)
project(hello_world_project LANGUAGES CXX CUDA)

find_package(CUDAToolkit REQUIRED)
find_package(OpenMP REQUIRED)
# Add CPU and CUDA libraries
add_subdirectory(cpu)
add_subdirectory(gpu)
add_subdirectory(main)

