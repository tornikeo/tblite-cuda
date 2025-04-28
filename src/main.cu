#include "main.h"
#include "xtb.h"
#include <iostream>
#include <assert.h>
#include "lapack/sygvd.h"
#include "blas/level2.h"

int main()
{

    // std::cout << "Running test_form_product..." << std::endl;
    // test_form_product<<<1, 1>>>();
    // cudaDeviceSynchronize();

    // std::cout << "Running test_horizontal_shift..." << std::endl;
    // test_horizontal_shift<<<1, 1>>>();
    // cudaDeviceSynchronize();

    // std::cout << "Running test_shift_operator..." << std::endl;
    // test_shift_operator<<<1, 1>>>();
    // cudaDeviceSynchronize();

    // std::cout << "Running test_multipole_grad_3d..." << std::endl;
    // test_multipole_grad_3d<<<1, 1>>>();
    // cudaDeviceSynchronize();

    // std::cout << "Running test_transform0..." << std::endl;
    // test_transform0<<<1, 1>>>();
    // cudaDeviceSynchronize();    

    // std::cout << "Running test_fill_matrix..." << std::endl;
    // test_fill_matrix<<<1, 1>>>();
    // cudaDeviceSynchronize();

    // std::cout << "Running test_fill_matrix..." << std::endl;
    // test_call_multipole_grad_cgto<<<1, 1>>>();
    // cudaDeviceSynchronize();
    // printf("Hello, world!\n");
    
    test_sygvd();
    test_blas();
    test_xtb();

    return 0;
}