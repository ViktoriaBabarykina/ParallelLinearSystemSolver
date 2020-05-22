CC = mpicc
override CFLAGS += -O3 -std=c99
LIBS = -lmpi -lm
OBJS = main.o matrix.o

ParallelLinearSystemSolver: $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean: 
	rm -f ParallelLinearSystemSolver $(OBJS)

.PHONY: clean
