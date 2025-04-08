#include "main.h"
#include <iostream>
#include <cassert>

// test __device__ __host__ void form_product(const double *a, const double *b, int la, int lb, double *d);
// test __host__ __device__ void horizontal_shift(double ae, int l, double *cfs);
__device__ void assert_close(const double *a, const double *b,
                             const size_t size, const double epsilon)
{
    for (size_t i = 0; i < size; ++i)
    {
        if (fabs(a[i] - b[i]) > epsilon)
        {
            printf("Assertion failed at index %zu: a[%zu] = %f, b[%zu] = %f\n",
                   i, i, a[i], i, b[i]);
            assert(false);
        }
    }
}

__global__ void test_form_product()
{
    const int la = 3, lb = 3;
    double a[la] = {1.0, 2.0, 3.0};
    double b[lb] = {4.0, 5.0, 6.0};
    double d[la + lb - 1] = {0.0};

    form_product(a, b, la, lb, d);

    printf("form_product result: ");
    for (int i = 0; i < la + lb - 1; ++i)
    {
        printf("%f ", d[i]);
    }
    printf("\n");
}

__global__ void test_horizontal_shift()
{
    const int l = 5;
    double cfs[l] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double ae = 2.0;

    horizontal_shift(ae, l, cfs);

    printf("horizontal_shift result: ");
    for (int i = 0; i < l; ++i)
    {
        printf("%f ", cfs[i]);
    }
    printf("\n");
}

__global__ void test_shift_operator()
{
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
        {0.13, 0.14, 0.15, 0.16, 0.17, 0.18}};

    // Output data
    double ddj[3][3] = {0};
    double dqj[3][6] = {0};

    // Call the shift_operator function
    shift_operator(vec, s, di, qi, ds, ddi, dqi, ddj, dqj);

    // Print results
    printf("ddj:\n");
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            printf("%f ", ddj[i][j]);
        }
        printf("\n");
    }

    printf("dqj:\n");
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            printf("%f ", dqj[i][j]);
        }
        printf("\n");
    }
}

__global__ void test_multipole_grad_3d()
{
    // Input parameters
    double rpi[3] = {0.083114717636212157, -0.39033516689727543, -0.25671966537784952};
    double rpj[3] = {-0.804297356879628, 3.7772557251143746, 2.4842645706404505};
    double ai = 5.9405971566712825, aj = 0.61389118221497996;
    int li[3] = {0, 0, 0};
    int lj[3] = {0, 0, 0};
    double s1d[MAXL2] = {1, 0, 0.076283605088381307, 0};

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
    assert(s3d == 1.0);
    printf("s3d: %f\n", s3d);
    double d3d_true[3] = {-0.804297356879628, 3.7772557251143746, 2.4842645706404505};
    assert_close(d3d, d3d_true, 3, 1e-4);
    printf("d3d: %f %f %f\n", d3d[0], d3d[1], d3d[2]);

    double q3d_true[6] = {0.72317784337193691, -3.0380367959679342, 14.343944417997701, -1.9980874279558183, 9.3837025721504457, 6.2478540620277627};
    printf("q3d: %f %f %f %f %f %f\n", q3d[0], q3d[1], q3d[2], q3d[3], q3d[4], q3d[5]);
    assert_close(q3d, q3d_true, 6, 1e-4);

    printf("ds3d: %f %f %f\n", ds3d[0], ds3d[1], ds3d[2]);

    printf("dd3d:\n");
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            printf("%f ", dd3d[i][j]);
        }
        printf("\n");
    }

    printf("dq3d:\n");
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            printf("%f ", dq3d[i][j]);
        }
        printf("\n");
    }
}

__global__ void test_transform0() {
    int li = 0, lj = 2;
    matrix cart;
    double cart_data[6] = {0.6451174322, 0.6451174322, 0.6451174322, 0, 0, 0};
    cart.at = cart_data;
    cart.R = 6, cart.C = 1;

    matrix sphr;
    double sphr_data[5];
    sphr.at = sphr_data;
    sphr.R = 5, sphr.C = 1;

    transform0(lj, li, cart, sphr);

    for (size_t i = 0; i < sphr.R; i++)
    {
        assert(sphr.at[i*sphr.C+0] == 0.0);
    }
    
    // printf("transform0 result: ");
    // for (int i = 0; i < r * c; ++i) {
    //     printf("%f \n", sphr[i]);
    // }
    // printf("\n");
}


__global__ void test_fill_matrix() {
    Matrix2D<float, 128, 128> mat;
    
    for (size_t i = 0; i < mat.rowCount(); i++)
    {
        for (size_t j = 0; j < mat.colCount(); j++)
        {
            mat.at(i, j) = i*4 + j;
        }
    }

    for (size_t i = 0; i < mat.rowCount(); i++)
    {
        for (size_t j = 0; j < mat.colCount(); j++)
        {
            auto val = mat.at(i,j);
            assert(val == i * 4 + j);
        }
    }
}

// #define MSAO 3;


int main()
{
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

    std::cout << "Running test_transform0..." << std::endl;
    test_transform0<<<1, 1>>>();
    cudaDeviceSynchronize();    

    std::cout << "Running test_fill_matrix..." << std::endl;
    test_fill_matrix<<<1, 1>>>();
    cudaDeviceSynchronize();

    std::cout << "Running test_fill_matrix..." << std::endl;
    test_call_multipole_grad_cgto<<<1, 1>>>();
    cudaDeviceSynchronize();
        
    return 0;
}