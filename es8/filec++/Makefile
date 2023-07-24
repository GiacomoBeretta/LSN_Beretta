CXXFLAGS = -Wall -O3 -pedantic -g -std=c++11

#CC = g++
#CFLAGS = -Wall -O3 --std=c++11

main: main.o random.o metropolis.o simAnnealing.o energy.o
	g++ -o main main.o random.o metropolis.o simAnnealing.o energy.o

main.o: main.cpp random.hpp metropolis.hpp simAnnealing.hpp energy.hpp
	g++ ${CXXFLAGS} -c main.cpp -o main.o

random.o: random.cpp random.hpp
	g++ ${CXXFLAGS} -c random.cpp -o random.o

metropolis.o: metropolis.cpp metropolis.hpp
	g++ ${CXXFLAGS} -c metropolis.cpp -o metropolis.o

simAnnealing.o: simAnnealing.cpp simAnnealing.hpp
	g++ ${CXXFLAGS} -c simAnnealing.cpp -o simAnnealing.o

energy.o: energy.cpp energy.hpp
	g++ ${CXXFLAGS} -c energy.cpp -o energy.o
	
clean :
	rm *.o main seed.out
