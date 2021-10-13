CC = gcc
MPICC = mpicc

default: all

serial:
	$(CC) -o main main.c helper.c

.PHONY: clean

all: serial

clean:
	rm -f main
