# Create the CUDA library
add_library(cuda_hello STATIC gpu.cu)
target_link_libraries(cuda_hello PUBLIC CUDA::cudart) 

# Set any CUDA-specific properties or flags
set_target_properties(cuda_hello PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

