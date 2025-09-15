# MPI Communication Patterns

Basic MPI communication patterns.

## SUMMA in 1D

A 1D SUMMA over the $k$-dimension (column of A and B) to compute $\mathbf{C} = \mathbf{A} \mathbf{B}^T$ with a direct-exchange (pairwise) of $B$, allowing each rank to fully compute $C$:

* $A$ is row-distributed (each rank has some subset of rows $m$ and all $k$ columns), which stay local.
* $B$ is column-distributed along $k$ (each rank has $n \times k$).
* $C$ is row-distributed like $A$ (each rank has $m \times n$).
* Implemented in `summa_1d.f90`
* Tested for 2, 3 and 4 processes

```shell
cmake -B cmake_build --fresh  
cmake --build cmake_build
mpirun -np 2 cmake_build/summa 
```
