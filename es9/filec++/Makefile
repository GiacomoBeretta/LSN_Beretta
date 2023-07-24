CXXFLAGS = -Wall -O3 -pedantic -g -std=c++11

#CC = g++
#CFLAGS = -Wall -O3 --std=c++11

main: main.o random.o position.o path.o GenAlgTSP.o
	g++ -o main main.o random.o position.o path.o GenAlgTSP.o

main.o: main.cpp random.hpp position.hpp path.hpp GenAlgTSP.hpp
	g++ ${CXXFLAGS} -c main.cpp -o main.o

random.o: random.cpp random.hpp
	g++ ${CXXFLAGS} -c random.cpp -o random.o

position.o: position.cpp position.hpp
	g++ ${CXXFLAGS} -c position.cpp -o position.o

path.o: path.cpp path.hpp
	g++ ${CXXFLAGS} -c path.cpp -o path.o

GenAlgTSP.o: GenAlgTSP.cpp GenAlgTSP.hpp
	g++ ${CXXFLAGS} -c GenAlgTSP.cpp -o GenAlgTSP.o
		
clean :
	rm *.o main seed.out
