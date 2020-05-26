FLAGS = -O3 -Wall -fopenmp
SRCS =   main.cpp IO.cpp conj_grad_solve.cpp
OBJS = main.o IO.o conj_grad_solve.o

ParallelLinearSystemSolver: $(OBJS)  
	mpic++ ${FLAGS} -o $@ $(OBJS)  

%.o: %.cpp 
	mpic++ $(FLAGS) -c $<

clean:
	rm -f $(OBJS) ./ParallelLinearSystemSolver

run:
	./run_default.sh

.PHONY: clean run
