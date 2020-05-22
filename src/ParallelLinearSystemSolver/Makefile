CC=mpicc
override CFLAGS +=-O3 -std=c99
LIBS=-lmpi -lm

ParallelLinearSystemSolver: main.c
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

clean: 
	rm -f ParallelLinearSystemSolver
