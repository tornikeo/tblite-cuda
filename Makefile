bin/main: bin/h0.o
	mkdir -p bin
	nvfortran src/main.f90 bin/type.o bin/h0.o -o bin/main

bin/type.o: src/type.f90
	nvfortran src/type.f90 -o bin/type.o -c

bin/h0.o: bin/type.o
	nvfortran src/h0.f90 -o bin/h0.o -c

clean: 
	rm -rf bin
	rm -rf obj