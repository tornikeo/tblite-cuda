# tblite-cuda (unfinished sample)

A rather sad attempt at porting tblite singlepoint to CUDA C++. Not useful as is.

This repo contains bottom-up translation efforts from tblite to CUDA C++. Ongoing work is finding a translation for fortran array expressions to C++.

[src/main.cu](./src/main.cu) contains most of the main code. [src/test.cu](./src/test.cu) contains some really simple driver code to test that the code doesn't throw exceptions. Current *unachieved* goal is to build up to `get_hamiltonian_gradient` function. We are not yet there.