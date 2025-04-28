#ifndef WAVEFUNCTION_FERMI_H
#define WAVEFUNCTION_FERMI_H
#include "../limits.h"
#include <math.h>


__device__
void get_electronic_entropy(
  float (&occ)[MAX_NAO],
  const float kt,
  float (&s)
);
__device__
void get_fermi_filling_(
  const int homo,
  const float kt,
  const float (&emo)[MAX_NAO],
  float (&occ)[MAX_NAO],
  float &e_fermi
);

__device__ 
void get_aufbau_filling(
  const float nel,
  int &homo,
  float (&occ)[MAX_NAO]
);

__device__
void get_fermi_filling(
  const float nel,
  const float kt,
  const float (&emo)[MAX_NAO],
  int &homo,
  float (&focc)[MAX_NAO],
  float &e_fermi,
  float &ts
);

#endif