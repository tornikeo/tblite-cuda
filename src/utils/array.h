#ifndef UTILS_ARRAY_H
#define UTILS_ARRAY_H

#define SIZEOF_ARRAY(arr) (sizeof(arr) / sizeof((arr)[0]))
#include <stdio.h>

/* TODO: priority high These don't work yet, 
    due to 'not being defined' somehow 
  in general, templates within in other files don't seem to 
  work. -rdc=true doesn't help this.
*/

template <int A>
__device__ __host__ inline 
void arange(
  float (&arr)[A]
)
{
  int state = 0;
  for (size_t i = 0; i < A; i++)
  {
    arr[i] = state;
    state++;
  }  
}

template <int A, int B>
__device__ __host__ inline 
void arange(
  float (&arr)[A][B]
)
{
  int state = 0;
  for (size_t i = 0; i < A; i++)
  {
    for (size_t j = 0; j < B; j++)
    {
      arr[i][j] = state;
      state++;
    }
  }  
}

template <int A, int B, int C>
__device__ __host__ inline 
void arange(
  float (&arr)[A][B][C]
)
{
  int state = 0;
  for (size_t i = 0; i < A; i++)
  {
    for (size_t j = 0; j < B; j++)
    {
      for (size_t k = 0; k < C; k++)
      {
        arr[i][j][k] = state;
        state++;
      }
    }
  }  
}

template <typename T, int A, int B>
__device__ __host__ inline 
void printr(const T (&arr)[A][B])
{
  printf("\n");
  for (size_t i = 0; i < A; i++) {
    printf("[");
    for (size_t j = 0; j < B; j++) {
      printf("%.2f", static_cast<float>(arr[i][j]));
      if (j < B - 1) {
        printf(", ");
      }
    }
    printf("]");
    if (i < A - 1) {
      printf(",\n");
    }
  }
  printf("\n");
}

template <typename T, int A, int B, int C>
__device__ __host__ inline 
void printr(const T (&arr)[A][B][C])
{
  printf("\n[");
  for (size_t i = 0; i < A; i++) {
    printf("[");
    for (size_t j = 0; j < B; j++) {
      printf("[");
      for (size_t k = 0; k < C; k++) {
        printf("%.2f", static_cast<float>(arr[i][j][k]));
        if (k < C - 1) {
          printf(", ");
        }
      }
      printf("]");
      if (j < B - 1) {
        printf(",\n ");
      }
    }
    printf("]");
    if (i < A - 1) {
      printf(",\n\n ");
    }
  }
  printf("]\n");
}

template <typename T, int A>
__device__ __host__ inline 
void printr(const T (&arr)[A], char fmt[5] = "%.2f\0")
{
  printf("\n");
  printf("[");
  for (size_t i = 0; i < A; i++) {
    printf(fmt, static_cast<float>(arr[i]));
    if (i < A - 1) {
      printf(", ");
    }
  }
  printf("]");
  printf("\n");
}

template <typename T, int A, int B>
__device__ __host__ inline
void assert_isclose(
  const T (&expected_yvec)[A][B], const T (&yvec)[A][B])
{
  for (int i = 0; i < A; i++) {
    for (int j = 0; j < B; j++) {
      if (fabs(static_cast<float>(expected_yvec[i][j]) - static_cast<float>(yvec[i][j])) > 1e-5) {
        printf("Assertion failed at [%i][%i]: expected %.6f, got %.6f\n",
               i, j, static_cast<float>(expected_yvec[i][j]), static_cast<float>(yvec[i][j]));
        return;
      }
    }
  }
  printf("All values are close within tolerance.\n");
}

// 1D array
template <typename T, size_t A>
__device__ 
inline void zero(T (&arr)[A]);

// 2D array
template <typename T, size_t A, size_t B>
__device__ 
inline void zero(T (&arr)[A][B]);

// 3D array
template <typename T, size_t A, size_t B, size_t C>
__device__ 
inline void zero(T (&arr)[A][B][C]);

// 4D array
template <typename T, size_t A, size_t B, size_t C, size_t D>
__device__ 
inline void zero(T (&arr)[A][B][C][D]);

#endif