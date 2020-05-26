FLAGS = -O3 -Wall -fopenmp
SRCS =   main.cpp IO.cpp conj_grad_solve.cpp
OBJS = main.o IO.o conj_grad_solve.o

ParallelLinearSystemSolver: $(OBJS)  
	mpic++ ${FLAGS} -o $@ $(OBJS) 
	export OMP_NUM_THREADS=2 

%.o: %.cpp 
	mpic++ $(FLAGS) -c $<

clean:
	rm -f $(OBJS) ./ParallelLinearSystemSolver

.PHONY: clean
