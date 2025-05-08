#ifndef UTILS_ARRAY_H
#define UTILS_ARRAY_H

#define SIZEOF_ARRAY(arr) (sizeof(arr) / sizeof((arr)[0]))
#include <stdio.h>
#include <cassert>

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

template <int A, int B, int C, int D>
__device__ __host__ inline 
void arange(
  float (&arr)[A][B][C][D]
)
{
  int state = 0;
  for (size_t i = 0; i < A; i++)
  {
    for (size_t j = 0; j < B; j++)
    {
      for (size_t k = 0; k < C; k++)
      {
        for (size_t l = 0; l < D; l++)
        {
          arr[i][j][k][l] = state;
          state++;
        }
      }
    }
  }  
}

template <typename T, int A>
__device__ __host__ inline
void fill(T (&arr)[A], T value) {
  for (int i = 0; i < A; i++) {
    arr[i] = value;
  }
}

template <typename T, int A, int B>
__device__ __host__ inline
void fill(T (&arr)[A][B], T value) {
  for (int i = 0; i < A; i++) {
    for (int j = 0; j < B; j++) {
      arr[i][j] = value;
    }
  }
}

template <typename T, int A, int B, int C>
__device__ __host__ inline
void fill(T (&arr)[A][B][C], T value) {
  for (int i = 0; i < A; i++) {
    for (int j = 0; j < B; j++) {
      for (int k = 0; k < C; k++) {
        arr[i][j][k] = value;
      }
    }
  }
}

template <typename T, int A, int B, int C, int D>
__device__ __host__ inline
void fill(T (&arr)[A][B][C][D], T value) {
  for (int i = 0; i < A; i++) {
    for (int j = 0; j < B; j++) {
      for (int k = 0; k < C; k++) {
        for (int l = 0; l < D; l++) {
          arr[i][j][k][l] = value;
        }
      }
    }
  }
}

__device__ inline
void printr(int n, int m, const float *arr, const char *fmt = "%.2f")
{
  printf("\n");
  for (size_t i = 0; i < n; i++) {
    printf("[");
    for (size_t j = 0; j < m; j++) {
      printf(fmt, static_cast<float>(arr[i * m + j]));
      if (j < m - 1) {
        printf(", ");
      }
    }
    printf("]");
    if (i < n - 1) {
      printf(",\n");
    }
  }
  printf("\n");
}

template <typename T, int A, int B>
__device__ __host__ inline 
void printr(const T (&arr)[A][B], const char *fmt = "%.2f")
{
  printf("\n");
  for (size_t i = 0; i < A; i++) {
    printf("[");
    for (size_t j = 0; j < B; j++) {
      printf(fmt, static_cast<float>(arr[i][j]));
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
void printr(const T (&arr)[A][B][C], const char *fmt = "%.2f")
{
  printf("\n[");
  for (size_t i = 0; i < A; i++) {
    printf("[");
    for (size_t j = 0; j < B; j++) {
      printf("[");
      for (size_t k = 0; k < C; k++) {
        printf(fmt, static_cast<float>(arr[i][j][k]));
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

template <typename T, int A, int B, int C, int D>
__device__ __host__ inline 
void printr(const T (&arr)[A][B][C][D], const char *fmt = "%.2f")
{
  printf("\n[");
  for (size_t i = 0; i < A; i++) {
    printf("[");
    for (size_t j = 0; j < B; j++) {
      printf("[");
      for (size_t k = 0; k < C; k++) {
        printf("[");
        for (size_t l = 0; l < D; l++) {
          printf(fmt, static_cast<float>(arr[i][j][k][l]));
          if (l < D - 1) {
            printf(", ");
          }
        }
        printf("]");
        if (k < C - 1) {
          printf(",\n   ");
        }
      }
      printf("]");
      if (j < B - 1) {
        printf(",\n\n  ");
      }
    }
    printf("]");
    if (i < A - 1) {
      printf(",\n\n\n ");
    }
  }
  printf("]\n");
}

template <typename T, int A>
__device__ __host__ inline 
void printr(const T (&arr)[A], const char *fmt = "%.2f")
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

template <typename T, int A>
__device__ __host__ inline
void assert_isclose(
  const T (&expected_yvec)[A], const T (&yvec)[A])
{
  for (int i = 0; i < A; i++) {
    if (fabs(static_cast<float>(expected_yvec[i]) - static_cast<float>(yvec[i])) > 1e-5) {
      printf("Assertion failed at [%i]: expected %.6f, got %.6f\n",
             i, static_cast<float>(expected_yvec[i]), static_cast<float>(yvec[i]));
      assert(false);
    }
  }
  printf("All values are close within tolerance.\n");
}

template <typename T, int A, int B>
__device__ __host__ inline
void assert_isclose(
  const T (&expected_yvec)[A][B], const T (&yvec)[A][B], const float atol=1e-5)
{
  for (int i = 0; i < A; i++) {
    for (int j = 0; j < B; j++) {
      if (fabs(static_cast<float>(expected_yvec[i][j]) - static_cast<float>(yvec[i][j])) > atol) {
        printf("Assertion failed at [%i][%i]: expected %.6f, got %.6f\n",
               i, j, static_cast<float>(expected_yvec[i][j]), static_cast<float>(yvec[i][j]));
        assert(false);
      }
    }
  }
  printf("All values are close within tolerance.\n");
}

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

template<int N>
__device__ __host__ inline float sum(const float (&arr)[N]) {
  float val = 0;
  #pragma unroll
  for (int i = 0; i < N; i++)
  { 
    val += arr[i];
  }
  return val;
}


template <typename T, int N>
__device__ __host__ inline void set(T (&arr)[N], T (&source)[N])
{
  for (size_t i = 0; i < N; i++)
  {
    arr[i] = source[i];
  }
}


template <typename T, int N, int M>
__device__ __host__ inline void set(T (&arr)[N][M], T (&source)[N][M])
{
  for (size_t i = 0; i < N; i++)
  {
    for (size_t j = 0; j < M; j++)
    {
      arr[i][j] = source[i][j];
    }
  }
}

#endif