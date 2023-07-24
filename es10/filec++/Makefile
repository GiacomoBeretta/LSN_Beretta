CXXFLAGS = -Wall -Wextra -O3 -pedantic -g -std=c++11

CC = mpicxx
#CFLAGS = -Wall -O3 --std=c++11

main: main.o random.o position.o path.o ParallelTempering.o
	${CC} main.o random.o position.o path.o ParallelTempering.o -o main

main.o: main.cpp random.hpp position.hpp path.hpp ParallelTempering.hpp
	${CC} ${CXXFLAGS} -c main.cpp -o main.o

random.o: random.cpp random.hpp
	${CC} ${CXXFLAGS} -c random.cpp -o random.o

position.o: position.cpp position.hpp
	${CC} ${CXXFLAGS} -c position.cpp -o position.o

path.o: path.cpp path.hpp
	${CC} ${CXXFLAGS} -c path.cpp -o path.o

ParallelTempering.o: ParallelTempering.cpp ParallelTempering.hpp
	${CC} ${CXXFLAGS} -c ParallelTempering.cpp -o ParallelTempering.o
		
clean :
	rm *.o main
