CC = mpicc
override CFLAGS += -O3 -std=c99
LIBS = -lmpi -lm
OBJS = main.o matrix.o common_math.o conjugate_gradient.o conjugate_gradient_parallel.o debugging.o

ParallelLinearSystemSolver: $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean: 
	rm -f ParallelLinearSystemSolver $(OBJS)

.PHONY: clean
