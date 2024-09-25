cmake_minimum_required(VERSION 3.18)
project(hello_world_project LANGUAGES CXX CUDA)

find_package(CUDAToolkit REQUIRED)
find_package(OMP REQUIRED)
# Add CPU and CUDA libraries
add_subdirectory(cpu)
add_subdirectory(cuda)
add_subdirectory(main)

