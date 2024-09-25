# Create the CPU library
add_library(cpu_hello STATIC cpu.cpp)
target_link_libraries(cpu_hello PUBLIC OpenMP::OpenMP_CXX)

# Optionally, specify any specific C++ standards or features
set_target_properties(cpu_hello PROPERTIES CXX_STANDARD 11)

