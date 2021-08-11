File Description:

graph.cpp : this file generates graph matrix(graph.txt) for the given number of cities (n=500) with maximum distance between two 		cities(1000kms)

graph.txt : contains input map(matrix) for the program

random.cpp : this file generates the number of random numbers given as an input (argument = number of random numbers)
 
random.txt :  contains input random numbers for choosing random city in each iteration
 
helpers.h : This is the header file which contains common functions for serial and parallel program. Some of the functions are finding the next city to visit, calculating best path, calculating the probability etc.

serial.cpp : This file contains the implementation of serial version of ant colony algorithm. The input to execute this file would be graph.txt + random.txt(both mentioned above) + number of ants + number of iteration 

parallel.cpp : This file contains the implementation of parallel version of ant colony algorithm with MPI. The input to execute this file would be graph.txt + random.txt(both mentioned above) + number of ants + number of external iteration + number of iteration over each node 

Slurm Folder : This folder contains the slurm scripts to run parallel.cpp over 1,2,4,8 and 16 nodes.


To compile:
To compile and run this project first need to run instruction 1) and 2) that will generates the input files for the program.
3) instruction is used for the compilation and execution of the serial version of program. the 4) instruction is for the compilation of parallel code and submitting the batch jobs on spartan.


Compile and generate inputs for the program:

1)	g++ graph.cpp -o graph.out
	./graph.out

2)	g++ random.cpp -o random.out
	./random.out 100000

3)	g++ serial.cpp -o serial.out
	./serial.out graph.txt random.txt 32 15000

4)	mpic++ parallel.cpp -o parallel.out
	sbatch node1.slurm
	sbatch node2.slurm
	sbatch node4.slurm
	sbatch node8.slurm
	sbatch node16.slurm
