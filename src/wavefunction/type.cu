#include "type.h"

__device__ void get_alpha_beta_occupation(
  float nocc, float nuhf, float &nalp, float &nbet)
{
  float ntmp, diff;

  // Ensure we cannot get a negative occupation here
  diff = fminf(nuhf, nocc);
  ntmp = nocc - diff;

  nalp = ntmp / 2.0f + diff;
  nbet = ntmp / 2.0f;
}