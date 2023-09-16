FC = gfortran
FFLAGS = -O2  -fdefault-real-8 -fdefault-double-8 
LFLAGS = -lblas -llapack

%.o : %.f95
	$(FC) $(FFLAGS) -c $<

main: main.o
	$(FC) main.o $(LFLAGS) -o main.exe

run: main.exe
	./main.exe

clean:
	rm -f *.o main.exe
