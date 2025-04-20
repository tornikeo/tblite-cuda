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