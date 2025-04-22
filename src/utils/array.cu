#include "array.h"
/* TODO: priority high These don't work yet, due to 'not being defined' somehow */

// 1D array
template <typename T, size_t A>
__device__ 
inline void zero(T (&arr)[A])
{
  for (size_t i = 0; i < A; ++i)
  {
    arr[i] = static_cast<T>(0);
  }
}

// 2D array
template <typename T, size_t A, size_t B>
__device__ 
inline void zero(T (&arr)[A][B])
{
  for (size_t i = 0; i < A; ++i)
  {
    for (size_t j = 0; j < B; ++j)
    {
      arr[i][j] = static_cast<T>(0);
    }
  }
}

// 3D array
template <typename T, size_t A, size_t B, size_t C>
__device__ 
inline void zero(T (&arr)[A][B][C])
{
  for (size_t i = 0; i < A; ++i)
  {
    for (size_t j = 0; j < B; ++j)
    {
      for (size_t k = 0; k < C; ++k)
      {
        arr[i][j][k] = static_cast<T>(0);
      }
    }
  }
}

// 4D array (already implemented)
template <typename T, size_t A, size_t B, size_t C, size_t D>
__device__ 
inline void zero(T (&arr)[A][B][C][D])
{
  for (size_t i = 0; i < A; ++i)
  {
    for (size_t j = 0; j < B; ++j)
    {
      for (size_t k = 0; k < C; ++k)
      {
        for (size_t l = 0; l < D; ++l)
        {
          arr[i][j][k][l] = static_cast<T>(0);
        }
      }
    }
  }
}
