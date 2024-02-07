ifx -O3 -r8 -lblas -llapack main.f90 -o main.exe
FORT_FMT_RECL=10000 ./main.exe
