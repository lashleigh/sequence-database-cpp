main: main.cpp helpers.cpp peptide.cpp protein.cpp
	g++ -o main main.cpp helpers.cpp peptide.cpp protein.cpp

main.o: main.cpp
	g++ -c main.cpp

helpers.o: helpers.cpp
	g++ -c helpers.cpp

peptide.o: peptide.cpp
	g++ -c peptide.cpp

protein.o: protein.cpp
	g++ -c protein.cpp

run:
	./main

clean:
	rm *.o main
