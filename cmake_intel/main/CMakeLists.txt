# Create the main executable
add_executable(main_exe main.cpp)

# Link the CPU and CUDA libraries
target_link_libraries(main_exe PRIVATE cpu_hello cuda_hello)

