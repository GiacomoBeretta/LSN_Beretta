CXXFLAGS = -Wall -O3 -pedantic -g -std=c++11

#CC = g++
#CFLAGS = -Wall -O3 --std=c++11

main: main.o random.o StatsCalculator.o
	g++ -o main main.o random.o  StatsCalculator.o

main.o: main.cpp random.hpp StatsCalculator.hpp
	g++ ${CXXFLAGS} -c main.cpp -o main.o

random.o: random.cpp random.hpp
	g++ ${CXXFLAGS} -c random.cpp -o random.o

StatisticsCalculator.o: StatsCalculator.cpp  StatsCalculator.hpp
	g++ ${CXXFLAGS} -c  StatsCalculator.cpp -o  StatsCalculator.o
	
#main.exe : main.o random.o
#	$(CC) random.o main.o -o main.exe
#main.o : main.cpp
#	$(CC) -c main.cpp -o main.o $(CFLAGS)
#random.o : random.cpp random.h
#	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main seed.out
