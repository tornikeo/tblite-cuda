#ifndef INTEGRAL_OVERLAP_H
#define INTEGRAL_OVERLAP_H
#define SQRTPI 1.77245385
#define SQRTPI3 5.568328

/*
integer, parameter :: maxl = 6
integer, parameter :: maxl2 = maxl*2
integer, parameter :: msao(0:maxl) = [1, 3, 5, 7, 9, 11, 13]
integer, parameter :: mlao(0:maxl) = [1, 3, 6, 10, 15, 21, 28]
integer, parameter :: lmap(0:maxl) = [0, 1, 4, 10, 20, 35, 56]
real(wp), parameter :: sqrtpi = sqrt(pi)
real(wp), parameter :: sqrtpi3 = sqrtpi**3 
*/

__device__
/* {1, 3, 5, 7, 9, 11, 13} */
inline constexpr int msao(int idx) 
{
  constexpr int values[] = {1, 3, 5, 7, 9, 11, 13}; // TODO: Any better way to do this?
  return values[idx];
}

__device__
/* {1, 3, 6, 10, 15, 21, 28}; */
inline constexpr int mlao(int idx)
{
  constexpr int values[] = {1, 3, 6, 10, 15, 21, 28};
  return values[idx];
}

__device__
/* {0, 1, 4, 10, 20, 35, 56} */
inline constexpr int lmap(int idx)
{
  constexpr int values[] = {0, 1, 4, 10, 20, 35, 56};
  return values[idx];
}



// constexpr int maxl = 6;
// constexpr int arraySize = maxl + 1;

// __constant__ int msao[arraySize];
// __constant__ int mlao[arraySize];
// __constant__ int lmap[arraySize];

/*

   ! x (+1), y (-1), z (0) in [-1, 0, 1] sorting
   integer, parameter :: lx(3, 84) = reshape([&
      & 0, &
      & 0,0,1, &
      & 2,0,0,1,1,0, &
      & 3,0,0,2,2,1,0,1,0,1, &
      & 4,0,0,3,3,1,0,1,0,2,2,0,2,1,1, &
      & 5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1, &
      & 6,0,0,3,3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2, &
      & 0, &
      & 1,0,0, &
      & 0,2,0,1,0,1, &
      & 0,3,0,1,0,2,2,0,1,1, &
      & 0,4,0,1,0,3,3,0,1,2,0,2,1,2,1, &
      & 0,5,0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2, &
      & 0,6,0,3,0,3,1,0,0,1,5,5,2,0,0,2,4,4,2,1,3,1,3,2,1,4,1,2, &
      & 0, &
      & 0,1,0, &
      & 0,0,2,0,1,1, &
      & 0,0,3,0,1,0,1,2,2,1, &
      & 0,0,4,0,1,0,1,3,3,0,2,2,1,1,2, &
      & 0,0,5,0,2,0,3,2,3,0,1,0,1,4,4,3,1,1,1,2,2, &
      & 0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,1,1,4,2], &
      & shape(lx), order=[2, 1])

*/

#endif