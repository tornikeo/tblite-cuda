#include "fermi.h"
#include "../limits.h"
#include <math.h>


__device__
void get_electronic_entropy(
  float (&occ)[MAX_NAO],
  const float kt,
  float (&s)
)
{
  const float thr = sqrtf(1e-5);
  s = 0.0f;

  for (int iao = 0; iao < MAX_NAO; ++iao) {
    if (occ[iao] > thr && (1.0f - occ[iao]) > thr) {
      s += occ[iao] * logf(occ[iao]) + (1.0f - occ[iao]) * logf(1.0f - occ[iao]);
    }
  }

  s *= kt;
}

__device__
void get_fermi_filling_(
  const int homo,
  const float kt,
  const float (&emo)[MAX_NAO],
  float (&occ)[MAX_NAO],
  float &e_fermi
)
{
    const int max_cycle = 200;
    const float thr = sqrtf(1e-5);

    float occt = static_cast<float>(homo);
    e_fermi = 0.5f * (emo[max(homo - 1, 0)] + emo[min(homo, MAX_NAO - 1)]);

    for (int ncycle = 0; ncycle < max_cycle; ++ncycle) {
        float total_number = 0.0f;
        float total_dfermi = 0.0f;

        for (int iao = 0; iao < MAX_NAO; ++iao) {
            float fermifunct = 0.0f;
            float dfermifunct = 0.0f;

            if ((emo[iao] - e_fermi) / kt < 50.0f) {
                float exp_val = expf((emo[iao] - e_fermi) / kt);
                fermifunct = 1.0f / (exp_val + 1.0f);
                dfermifunct = exp_val / (kt * powf(exp_val + 1.0f, 2));
            }

            occ[iao] = fermifunct;
            total_number += fermifunct;
            total_dfermi += dfermifunct;
        }

        float change_fermi = (occt - total_number) / total_dfermi;
        e_fermi += change_fermi;

        if (fabsf(occt - total_number) <= thr) {
            break;
        }
    }
}


__device__ 
void get_aufbau_filling(
  const float nel,
  int &homo,
  float (&occ)[MAX_NAO]
)
{
    // Initialize occupation array to zero
    for (int i = 0; i < MAX_NAO; ++i) {
        occ[i] = 0.0f;
    }

    // Determine the highest occupied molecular orbital (homo)
    homo = static_cast<int>(floorf(nel));

    // Fill occupation array up to homo
    for (int i = 0; i < min(homo, MAX_NAO); ++i) {
        occ[i] = 1.0f;
    }

    // Handle fractional occupation for the next orbital
    if (homo < MAX_NAO) {
        occ[homo] = fmodf(nel, 1.0f);
    }

    // Adjust homo based on fractional occupation
    if (fmodf(nel, 1.0f) > 0.5f) {
        homo += 1;
    }
}

__device__
void get_fermi_filling(
  const float nel,
  const float kt,
  const float (&emo)[MAX_NAO],
  int &homo,
  float (&focc)[MAX_NAO],
  float &e_fermi,
  float &ts
)
{
    float etmp = 0.0f;
    float stmp = 0.0f;

    ts = 0.0f;
    e_fermi = 0.0f;

    get_aufbau_filling(nel, homo, focc);

    if (homo > 0) {
        get_fermi_filling_(homo, kt, emo, focc, etmp);
        get_electronic_entropy(focc, kt, ts);
        e_fermi = 0.5f * etmp;
    }
}
