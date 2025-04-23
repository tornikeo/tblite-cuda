#ifndef POTENTIAL_H
#define POTENTIAL_H
#include "limits.h"
#include "structure.h"
#include "basis/type.h"

// !> Container for density dependent potential-shifts
// type :: potential_type
//    !> Atom-resolved charge-dependent potential shift
//    real(wp), allocatable :: vat(:, :)
//    !> Shell-resolved charge-dependent potential shift
//    real(wp), allocatable :: vsh(:, :)
//    !> Orbital-resolved charge-dependent potential shift
//    real(wp), allocatable :: vao(:, :)

//    !> Atom-resolved dipolar potential
//    real(wp), allocatable :: vdp(:, :, :)
//    !> Atom-resolved quadrupolar potential
//    real(wp), allocatable :: vqp(:, :, :)
// contains
class potential_type {
public:
  // Atom-resolved charge-dependent potential shift
  float vat[MAX_NSPIN][MAX_NAT];
  // Shell-resolved charge-dependent potential shift
  float vsh[MAX_NSPIN][MAX_NSH];
  // Orbital-resolved charge-dependent potential shift
  float vao[MAX_NSPIN][MAX_NAO];

  // Atom-resolved dipolar potential
  float vdp[MAX_NAT][MAX_NSPIN][3];
  // Atom-resolved quadrupolar potential
  float vqp[MAX_NAT][MAX_NSPIN][6];

  // Member function to reset all values to 0.0
  __device__ void reset();
};

__device__
void new_potential(
  potential_type &self, 
  const structure_type &mol, 
  const basis_type &bas,
  const int nspin
);

#endif