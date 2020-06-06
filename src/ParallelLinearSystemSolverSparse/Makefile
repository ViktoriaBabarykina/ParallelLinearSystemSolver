FLAGS = -O3 -Wall -fopenmp
SRCS =   main.cpp IO.cpp conj_grad_solve.cpp SparseMatrix.cpp Vector.cpp
OBJS = main.o IO.o conj_grad_solve.o SparseMatrix.o Vector.o

ParallelLinearSystemSolverSparse: $(OBJS)  
	mpic++ ${FLAGS} -o $@ $(OBJS)  

%.o: %.cpp 
	mpic++ $(FLAGS) -c $<

clean:
	rm -f $(OBJS) ./ParallelLinearSystemSolverSparse

run:
	./run_default.sh

.PHONY: clean run
