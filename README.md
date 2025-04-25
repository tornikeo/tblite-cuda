# tblite-cuda (work-in-progress)

This is a work-in-progress translation of [tblite](https://github.com/tblite/tblite/) from from FORTRAN to [CUDA C++](https://docs.nvidia.com/cuda/cuda-c-programming-guide/).

# Building from source

To build tblite-cuda, you need to have a CUDA compiler (`nvcc`) and a CUDA-capable GPU, with [compute capability](https://developer.nvidia.com/cuda-gpus) of 8.0 or more. Building is done through makefile:

```sh
make
```

This will build an executable into `build/main.bin`, which can be used, like:

```sh
compute-sanitizer ./build/main.bin
```

# Current work

Implementing `xtb_singlepoint` as a CUDA `__device__` function. Once created, this should allow running `xtb_singlepoint` in a batched fashion on a GPU. 

# Worklog

- Implemented next_scf up to `get_density` call, which requires a blas `sygvd` call, to solve a generalized eigenvalue `A*x=lamd*B*x` problem. I would like to use cuSOLVER for this due to perf and simplicity reasons.