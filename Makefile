CC=g++
CFLAGS= -march=native 
LIB= -O3  
SOURCES= PMAC.cpp 
all: 
	$(CC) -o test $(SOURCES) $(LIB) $(CFLAGS) 
clean: 
	rm *.o 