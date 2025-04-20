#include <cuda_runtime.h>
#include <stdio.h>
#include <cassert>
#include <cstdlib>
#include "xtb.h"



__device__
inline void* calloc(size_t nmem, size_t size) {
  void * ptr = malloc(nmem * size);
  memset(ptr, 0, nmem * size);
  return ptr;
} 

__device__ 
void xtb_singlepoint(
  const structure_type mol, 
  const xtb_calculator calc, 
  const wavefunction_type wfn,
  const float accuracy,
  
  float energy,
  float (&gradient)[][3],
  float (&sigma)[][3],
  const int verbosity
) {
  bool grad, converged, econverged, pconverged;
  int prlevel;
  float econv, pconv, cutoff, elast, dpmom[3], qpmom[6], nel;
  // void *ptr;
  // real(wp), allocatable :: energies(:), edisp(:), erep(:), exbond(:), eint(:), eelec(:)
  // real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), dEdcn(:)
  // real(wp), allocatable :: selfenergy(:), dsedcn(:), lattr(:, :), wdensity(:, :, :)
  // type(integral_type) :: ints
  // real(wp), allocatable :: tmp(:)
  // type(potential_type) :: pot
  // type(container_cache) :: ccache, dcache, icache, hcache, rcache
  // type(broyden_mixer) :: mixer
  // type(timer_type) :: timer
  // type(error_type), allocatable :: 
  
  float *energies = (float*)malloc(mol.nat * sizeof(float));
  assert(energies != NULL);
  memset(energies, 0, mol.nat * sizeof(float));
  printf("Thread %d got pointer: %p\n", threadIdx.x, energies);
  // float *energies = (float*)calloc(mol.nat, sizeof(float));
  // float *edisp = (float*)calloc(mol.nat, sizeof(float));
  // float *erep = (float*)calloc(mol.nat, sizeof(float));
  // float *exbond = (float*)calloc(mol.nat, sizeof(float));
  // float *eint = (float*)calloc(mol.nat, sizeof(float));
  // float *eelec = (float*)calloc(mol.nat, sizeof(float));
  // float *cn = (float*)calloc(mol.nat, sizeof(float));
  // float *cdndr = (float*)calloc(mol.nat * mol.nat * 3, sizeof(float));
  
  // float *dEdcn = (float*)calloc(mol.nat, sizeof(float));
  // float *selfenergy = (float*)calloc(mol.nat, sizeof(float));
  // float *tmp = (float*)calloc(mol.nat, sizeof(float));
  free(energies);
}

__global__
void test_xtb_singlepoint()
{
  structure_type mol = {0};
  xtb_calculator calc = {0};
  wavefunction_type wfn = {0};
  
  mol.nat = 9;

  float accuracy = 0.01f;
  float energy = 0.0f;

  float gradient[9][3] = {0};
  float sigma[3][3] = {0.0};
  int verbosity = 1;

  xtb_singlepoint(mol, calc, wfn, accuracy, energy, gradient, sigma, verbosity);
  printf("Energy: %f\n", energy);
}


void xtb_test() {
  cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
  cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128*1024*1024);
  test_xtb_singlepoint<<<10,1>>>();
  cudaDeviceSynchronize();
}