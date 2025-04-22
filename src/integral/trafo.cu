#include "trafo.h"
/* Moved to multipole. See header as to why */

// template <int A, int B>
// __device__
// void transform0(
//   const int lj, 
//   const int li, 
//   const float (&cart)[A][A], 
//   float (&sphr)[B][B]
// ) {
//   for (size_t i = 0; i < B; i++)
//   {
//    for (size_t j = 0; j < B; j++)
//    {
//       sphr[i][j] = 0;
//    }
//   }
// }

// template <int A, int B, int C>
// __device__
// void transform1(
//   const int lj, 
//   const int li, 
//   const float (&cart)[A][A][C], 
//         float (&sphr)[B][B][C]
// ) {
//   float tmpmat[B][B];

//   for (int i = 0; i < B; i++)
//   {  
//     for (int j = 0; j < B; j++)
//     {
//       tmpmat[i][j] = sphr[i][j][0];
//     }
//   }
//   transform0(lj, li, tmpmat, tmpmat);
// }