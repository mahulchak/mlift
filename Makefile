#to compile the executable for mlift, the liftover program

CC = g++
CFLAGS = -g -Wall -std=c++0x 

default: mlift
mlift: mllib.o mlift.o
	$(CC) $(CFLAGS) -o mlift mllib.o mlift.o

mllib.o: mllib.cpp ml.h
	$(CC) $(CFLAGS) -c mllib.cpp

mlift.o: mlift.cpp ml.h
	$(CC) $(CFLAGS) -c mlift.cpp

clean:
	$(RM) *.o 
