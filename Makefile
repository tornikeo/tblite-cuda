all: 
	mkdir -p bin
	nvfortran -cuda src/type.f90 -o bin/type.o -c
	nvfortran -cuda src/kernel.f90 -o bin/kernel.o -c
	nvfortran -cuda src/h0.f90 -o bin/h0.o -c
	nvfortran -cuda src/main.f90 bin/type.o bin/kernel.o bin/h0.o -o bin/main.bin
	./bin/main.bin

clean: 
	rm -rf bin