main: main.cpp helpers.cpp
	g++ -o main main.cpp helpers.cpp

main.o: main.cpp header.hpp constants.hpp AminoAcidMasses.h
	g++ -c main.cpp

helpers.o: helpers.cpp header.hpp constants.hpp
	g++ -c helpers.cpp 

run:
	./main test.fasta

clean:
	rm *.o main
