# Create the CPU library
add_library(cpu_hello STATIC cpu.cpp)

# Optionally, specify any specific C++ standards or features
set_target_properties(cpu_hello PROPERTIES CXX_STANDARD 11)

