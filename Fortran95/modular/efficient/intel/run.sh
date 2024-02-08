ifx -O3 -r8 -lblas -llapack main.f90 -o main.exe
FORT_FMT_RECL=10000 ./main.exe

#ifx -O3 -r8 -qmkl main.f90 -o main.exe
#intel math kernel library is a replacement for blas and lapack, it is faster with intel compilers.
#install it using the command: sudo apt install intel-mkl
#see the intel documentation for more information.
