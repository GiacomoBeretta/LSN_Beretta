CC = g++
CFLAGS = -Wall -O3 --std=c++11

ising : ising.o random.o
	$(CC) random.o ising.o -o ising
ising.o : ising.cpp ising.hpp
	$(CC) -c ising.cpp -o ising.o $(CFLAGS)
random.o : random.cpp random.hpp
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o ising
