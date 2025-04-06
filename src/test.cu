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

__global__ void test_shift_operator() {
    // Input data
    double vec[3] = {1.0, 2.0, 3.0};
    double s = 1.5;
    double di[3] = {0.5, 1.0, 1.5};
    double qi[6] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    double ds[3] = {0.05, 0.1, 0.15};
    double ddi[3][3] = {{0.1, 0.2, 0.3}, {0.4, 0.5, 0.6}, {0.7, 0.8, 0.9}};
    double dqi[3][6] = {
        {0.01, 0.02, 0.03, 0.04, 0.05, 0.06},
        {0.07, 0.08, 0.09, 0.10, 0.11, 0.12},
        {0.13, 0.14, 0.15, 0.16, 0.17, 0.18}
    };

    // Output data
    double ddj[3][3] = {0};
    double dqj[3][6] = {0};

    // Call the shift_operator function
    shift_operator(vec, s, di, qi, ds, ddi, dqi, ddj, dqj);

    // Print results
    printf("ddj:\n");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            printf("%f ", ddj[i][j]);
        }
        printf("\n");
    }

    printf("dqj:\n");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 6; ++j) {
            printf("%f ", dqj[i][j]);
        }
        printf("\n");
    }
}



__global__ void test_multipole_grad_3d() {
    // Input parameters
    double rpi[3] = {1.0, 2.0, 3.0};
    double rpj[3] = {4.0, 5.0, 6.0};
    double ai = 1.5, aj = 2.5;
    int li[3] = {1, 2, 3};
    int lj[3] = {2, 1, 0};
    double s1d[MAXL2] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2};

    // Output parameters
    double s3d = 0.0;
    double d3d[3] = {0.0};
    double q3d[6] = {0.0};
    double ds3d[3] = {0.0};
    double dd3d[3][3] = {0.0};
    double dq3d[3][6] = {0.0};

    // Call the function
    multipole_grad_3d(rpi, rpj, ai, aj, li, lj, s1d, s3d, d3d, q3d, ds3d, dd3d, dq3d);

    // Print results
    printf("s3d: %f\n", s3d);
    printf("d3d: %f %f %f\n", d3d[0], d3d[1], d3d[2]);
    printf("q3d: %f %f %f %f %f %f\n", q3d[0], q3d[1], q3d[2], q3d[3], q3d[4], q3d[5]);
    printf("ds3d: %f %f %f\n", ds3d[0], ds3d[1], ds3d[2]);

    printf("dd3d:\n");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            printf("%f ", dd3d[i][j]);
        }
        printf("\n");
    }

    printf("dq3d:\n");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 6; ++j) {
            printf("%f ", dq3d[i][j]);
        }
        printf("\n");
    }
}

int main() {
    std::cout << "Running test_form_product..." << std::endl;
    test_form_product<<<1, 1>>>();
    cudaDeviceSynchronize();

    std::cout << "Running test_horizontal_shift..." << std::endl;
    test_horizontal_shift<<<1, 1>>>();
    cudaDeviceSynchronize();

    std::cout << "Running test_shift_operator..." << std::endl;
    test_shift_operator<<<1, 1>>>();
    cudaDeviceSynchronize();


    std::cout << "Running test_multipole_grad_3d..." << std::endl;
    test_multipole_grad_3d<<<1, 1>>>();
    cudaDeviceSynchronize();
    return 0;

    return 0;
}