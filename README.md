# Traveling Salesman Problem using Ant colony optimization with MPI

Travelling salesman problem is one of the best-researched algorithms in the stream of computer science. TSP problem is the NP-hard problem and different algorithms are used to solve that problem. Ant colony optimization is one of them, this algorithm uses the action of ants to find the shortest path for food using pheromone.

The idea behind the ant colony algorithm is that each ant must visit all the cities at most once and find the shortest path. The ant chose the next city based on the distance and the pheromone intensity. The probabilistic method is used to choose the next city. After all the ants finished with the finding the path, the most used edges have higher pheromone intensity than other edges. That path would be considered as the best path for the next iteration. The same process is performed until the termination condition is met. The ultimate aim is to find the shortest path between all the cities and the starting and ending city should be the same.

## Parallel Algorithm

Message Passing Interface(MPI) is used to achieve parallelism of the ant-colony algorithm. This provides the facility to perform operation simultaneously on different computational nodes using distributed memory.

The parallel algorithm used for ant-colony optimization is described in figure. The figure only represents about single loop iteration without including the local computation.

![alt text](https://github.com/avpatel26/Parallel_TSP_using_AntColonyOptimization/blob/main/images/parallel_algo.png?raw=true)

The key factor for the algorithm is the distribution of the ants over nodes. Every node in the program calculates the distribution of ants as it is necessary to have information that how many ants every node has. This also helps at the end of the algorithm for the comparison result. Additionally, every node shares not only the best path of iteration but also the best Cost and pheromones matrix on the best path. Later the receiving node compares the best path with its path and updates the cost if the receiving path is better than the local one.

## Comparison of Speedup


 The result obtained by running the MPI program on spartan.

![alt text](https://github.com/avpatel26/Parallel_TSP_using_AntColonyOptimization/blob/main/images/graph.png?raw=true)

The absolute speedup is achieved but that is not so good because of the high cost of the communication between nodes. As the number of nodes increases the communication time between the nodes is also increased because of that the speedup is decreased compared to optimal speedup. The speedup with 16 cores is far worst that the speedup with 8 cores because of the high communication rate so it can be derived than the given parallel algorithm does not scale well.
