#include <cuda_runtime.h>
#include <iostream>
#include <cstdio> // 用于 printf

// 定义 atomicMin 函数
inline __device__
float atomicMin(float *addr, float value) {
    float old = *addr, assumed;
    if (old <= value) return old;
    do {
        assumed = old;
        old = __int_as_float(atomicCAS((unsigned int*)addr, __float_as_int(assumed), __float_as_int(value)));
        value = min(value, old);
    } while (old != assumed);
    return old;
}

// CUDA 内核函数
__global__ void updateMinKernel(float *globalMin, float *localMins, int size) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < size) {
        // 使用 atomicMin 函数更新全局最小值
        atomicMin(globalMin, localMins[idx]);
    }
}

int main() {
    // 分配设备内存
    float *d_globalMin;
    float *d_localMins;
    constexpr size_t size = 1024; // 假设有1024个局部最小值
    cudaMalloc(&d_globalMin, sizeof(float));
    cudaMalloc(&d_localMins, size * sizeof(float));

    // 初始化全局最小值
    float globalMin = FLT_MAX;
    cudaMemcpy(d_globalMin, &globalMin, sizeof(float), cudaMemcpyHostToDevice);

    // 初始化局部最小值数组（例如，使用随机值）
    float localMins[size];
    for (int i = 0; i < size; ++i) {
        // localMins[i] = static_cast<float>(i) / size;
        localMins[i] = static_cast<float>(i) + 0.1;
    }
    cudaMemcpy(d_localMins, localMins, size * sizeof(float), cudaMemcpyHostToDevice);

    // 启动内核
    int threadsPerBlock = 256;
    int blocksPerGrid = (size + threadsPerBlock - 1) / threadsPerBlock;
    updateMinKernel<<<blocksPerGrid, threadsPerBlock>>>(d_globalMin, d_localMins, size);

    // 复制结果回主机
    float h_globalMin;
    cudaMemcpy(&h_globalMin, d_globalMin, sizeof(float), cudaMemcpyDeviceToHost);

    // 打印结果
    printf("Global minimum: %f\n", h_globalMin);

    // 释放设备内存
    cudaFree(d_globalMin);
    cudaFree(d_localMins);

    return 0;
}