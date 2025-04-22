#ifndef UTILS_ARRAY_H
#define UTILS_ARRAY_H

/* TODO: priority high These don't work yet, due to 'not being defined' somehow */

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