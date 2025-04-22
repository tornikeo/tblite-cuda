#include "trafo.h"

__device__
void transform0(
  const int lj, 
  const int li, 
  const float (&cart)[mlao(MAXL)][mlao(MAXL)], 
  float (&sphr)[msao(MAXL)][msao(MAXL)]
) {}

template <size_t D>
__device__
void transform1(
  const int lj, 
  const int li, 
  const float (&cart)[mlao(MAXL)][mlao(MAXL)][D], 
  float (&sphr)[msao(MAXL)][msao(MAXL)][D]
) {}