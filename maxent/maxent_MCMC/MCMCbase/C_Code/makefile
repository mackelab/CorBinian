clean:
	rm -f *.o

all:
	make clean
	make learn
	make shuffle
	make fake

compile_all:
	make clean
	make learn
	make shuffle
	make fake

CC = gcc
#CC = gccmemstat

CFLAGS = -Wall -ansi -pedantic
#CFLAGS = 

fake: fake_data.o sampler.o mt19937ar.o
	$(CC) $(CFLAGS) -o fake_data.exe fake_data.o sampler.o mt19937ar.o -lm

learn: learn_stop.o learn_iterative.o sampler.o expec_err.o mt19937ar.o OutputFunctions.o ReadInputs.o
	$(CC) $(CFLAGS) -o learn_stop.exe learn_stop.o learn_iterative.o sampler.o expec_err.o mt19937ar.o OutputFunctions.o ReadInputs.o -lm

learn_iterative.o: learn_iterative.c learn_iterative.h learn_stop.h
	$(CC) $(CFLAGS) -c learn_iterative.c

error: expec_err.o
	$(CC) $(CFLAGS) -o expec_err.exe expec_err.o -lm

fake_data.o: fake_data.c learn_stop.h sampler.h
	$(CC) $(CFLAGS) -c fake_data.c

learn_stop.o: learn_stop.c learn_stop.h learn_iterative.h sampler.h expec_err.h
	$(CC) $(CFLAGS) -c learn_stop.c

expec_err.o: expec_err.c learn_stop.h expec_err.h
	$(CC) $(CFLAGS) -c expec_err.c

mt19937ar.o: mt19937ar.c mt19937ar.h
	$(CC) $(CFLAGS) -c mt19937ar.c

sampler: sampler.o mt19937ar.o
	$(CC) $(CFLAGS) -o sampler.exe sampler.o mt19937ar.o -lm

sampler.o: sampler.c learn_stop.h sampler.h mt19937ar.h
	$(CC) $(CFLAGS) -c sampler.c

shuffle: shuffle.o mt19937ar.h
	$(CC) $(CFLAGS) -o shuffle.exe shuffle.o mt19937ar.o -lm

shuffle.o: shuffle.c mt19937ar.h
	$(CC) $(CFLAGS) -c shuffle.c

OutputFunctions.o: OutputFunctions.c OutputFunctions.h expec_err.h
	$(CC) $(CFLAGS) -c OutputFunctions.c

ReadInputs.o: ReadInputs.c  
	$(CC) $(CFLAGS) -c ReadInputs.c
