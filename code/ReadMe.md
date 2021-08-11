## Folder Structure:

```
├── graph.cpp: generates random graph matrix (graph.txt) 
├── graph.txt: contains input map(matrix) for the program
├── helpers.h: This is the header file which contains common functions for serial and parallel program
├── parallel.cpp: This file contains the implementation of parallel version of ant colony algorithm with MPI 
├── random.cpp: generates the number of random numbers
├── random.txt: contains input random numbers for choosing random city in each iteration
├── serial.cpp: This file contains the implementation of serial version of ant colony algorithm
└── slurms: This folder contains the slurm scripts to run parallel.cpp over 1,2,4,8 and 16 nodes
    ├── node1.slurm: runs TSP on a single node
    ├── node2.slurm: runs TSP on a 2 nodes
    ├── node4.slurm: runs TSP on a 4 nodes
    ├── node8.slurm: runs TSP on a 8 nodes
    └── node16.slurm: runs TSP on a 16 nodes
```

## Compile and Run

### Generate input files

To compile and run this project first need to run graph.cpp and random.cpp that will generates the input files for the program.

Generates graph.txt file for the given number of cities (n=500) with maximum distance between two cities(1000kms)
```
g++ graph.cpp -o graph.out
./graph.out
```

Generates random.txt that utilizes for choosing random city in each iteration
```
g++ random.cpp -o random.out
./random.out <random_number>
```

### Run serial version of program

Given number of ants travel on all the city and most used edges are passed to the next iteration. After all iterations optimum path is provided. 
```
g++ serial.cpp -o serial.out
./serial.out <graph_file> <random_number_file> <number_of_ants> <number_of_iteration> 
e.g. ./serial.out graph.txt random.txt 32 15000
```

### Run parallel version of program

In each external iteration every node iterates locally and provides optimum path to the master process. Finally master process decide optimum path after all external iterations.
```
mpic++ parallel.cpp -o parallel.out
./parallel.out <graph_file> <random_number_file> <number_of_ants> <number_of_external_iteration> <number_of_iteration_on_each_node>
e.g. ./parallel.out graph.txt random.txt 32 15000 1
```


### Run slurm files

This code submit the jobs on spartan (High Performance Computing at University of Melbourne)

```
mpic++ parallel.cpp -o parallel.out
sbatch node1.slurm
sbatch node2.slurm
sbatch node4.slurm
sbatch node8.slurm
sbatch node16.slurm
```
