#include <cuda_runtime.h>
#include <stdio.h>
#include "xtb.h"

__device__ 
void xtb_singlepoint(
  const structure_type mol, 
  const xtb_calculator calc, 
  const wavefunction_type wfn,
  float accuracy,
  float energy
) {
  printf("We are here.");
}

__global__
void test_xtb_singlepoint()
{
  structure_type mol = {0};
  xtb_calculator calc = {0};
  wavefunction_type wfn = {0};
  float accuracy = 0.01f;
  float energy = 0.0f;

  xtb_singlepoint(mol, calc, wfn, accuracy, energy);
  printf("Energy: %f\n", energy);
}


void xtb_test() {
  test_xtb_singlepoint<<<1,1>>>();
  cudaDeviceSynchronize();
}