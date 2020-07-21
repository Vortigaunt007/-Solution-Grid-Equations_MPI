all: main
#
#CommandLineArg.o: CommandLineArg.cpp CommandLineArg.h
#		mpicxx CommandLineArg.cpp 

#BasicOperation.o: BasicOperation.cpp BasicOperation.h
#		mpicxx BasicOperation.cpp

#Grid.o: Grid.cpp Grid.h
#		mpicxx Grid.cpp

#MatrixCSR.o: MatrixCSR.cpp MatrixCSR.h
#		mpicxx MatrixCSR.cpp

#Solver.o: Solver.cpp Solver.h
#		mpicxx Solver.cpp

#Vector.o: Vector.cpp Vector.h
#		mpicxx Vector.cpp


#main.o: main.cpp CommandLineArg.o Grid.o MatrixCSR.o Vector.o Solver.o
#		mpicxx main.cpp 


#run: main.o BasicOperation.o CommandLineArg.o Grid.o MatrixCSR.o Solver.o Vector.o
#		mpirun -np 4 main.o BasicOperation.o CommandLineArg.o Grid.o MatrixCSR.o Solver.o Vector.o -o mpitest

main.o: main.cpp CommandLineArg.cpp BasicOperation.cpp Grid.cpp MatrixCSR.cpp Vector.cpp Solver.cpp MPI_Exchange.cpp
		mpicxx main.cpp CommandLineArg.cpp BasicOperation.cpp Grid.cpp MatrixCSR.cpp Vector.cpp Solver.cpp MPI_Exchange.cpp CommandLineArg.h BasicOperation.h Grid.h MatrixCSR.h Vector.h Solver.h MPI_Exchange.h -o main.o 
run: main.o 
		mpirun -np 4 main.o

clean:
		rm -rf *.gch 