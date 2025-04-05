#include "main.h"
#include <iostream>

// test __device__ __host__ void form_product(const double *a, const double *b, int la, int lb, double *d);
// test __host__ __device__ void horizontal_shift(double ae, int l, double *cfs);


__global__ void test_form_product() {
    const int la = 3, lb = 3;
    double a[la] = {1.0, 2.0, 3.0};
    double b[lb] = {4.0, 5.0, 6.0};
    double d[la + lb - 1] = {0.0};

    form_product(a, b, la, lb, d);

    printf("form_product result: ");
    for (int i = 0; i < la + lb - 1; ++i) {
        printf("%f ", d[i]);
    }
    printf("\n");
}

__global__ void test_horizontal_shift() {
    const int l = 5;
    double cfs[l] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double ae = 2.0;

    horizontal_shift(ae, l, cfs);

    printf("horizontal_shift result: ");
    for (int i = 0; i < l; ++i) {
        printf("%f ", cfs[i]);
    }
    printf("\n");
}

int main() {
    std::cout << "Running test_form_product..." << std::endl;
    test_form_product<<<1, 1>>>();
    cudaDeviceSynchronize();

    std::cout << "Running test_horizontal_shift..." << std::endl;
    test_horizontal_shift<<<1, 1>>>();
    cudaDeviceSynchronize();

    return 0;
}