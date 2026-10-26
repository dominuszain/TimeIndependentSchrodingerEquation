# Fortran solver

`main.f90` builds the finite-difference Hamiltonian for the Coulomb potential
`V(x) = -1/x` on `x ∈ [1, 100]` with 1000 grid points. It calls LAPACK's `dsyev`
through the Intel MKL to diagonalise the symmetric matrix and writes the
eigenvalues and eigenvectors to `eigenout.txt`.

## Build and run

With Intel `ifx` and MKL:

```bash
bash run.sh
```

Or compile directly:

```bash
ifx -O3 -r8 -qmkl main.f90 -o main.exe
./main.exe
```

`run.sh` sets `FORT_FMT_RECL=50000` so the formatted output can write the long
eigenvector rows.
