CXXFLAGS = -Wall -O3 -pedantic -g -std=c++11

#CC = g++
#CFLAGS = -Wall -O3 --std=c++11

MD_MC : MD_MC.o random.o
	g++ random.o MD_MC.o -o MD_MC
MD_MC.o : MD_MC.cpp
	g++ ${CXXFLAGS} -c MD_MC.cpp -o MD_MC.o

random.o : random.cpp random.hpp
	g++ ${CXXFLAGS} -c random.cpp -o random.o

	
#main.exe : main.o random.o
#	$(CC) random.o main.o -o main.exe
#main.o : main.cpp
#	$(CC) -c main.cpp -o main.o $(CFLAGS)
#random.o : random.cpp random.h
#	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o MD_MC seed.out
