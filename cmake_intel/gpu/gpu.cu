#include <stdio.h>
#include <cuda_runtime.h>
#include "gpu.h"

__global__ void hello_from_gpu_kernel() {
    printf("Hello World from the GPU!\n");
}

void hello_from_gpu() {
    hello_from_gpu_kernel<<<1, 1>>>();
    cudaDeviceSynchronize();
}

