#include "helpers.h" 
#include <mpi.h>
int main(int argc, char* argv[]) {
  int rank, size;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (rank == 0) {
    printf("Number of Nodes: %d\n", size);
  }
  int i, j, ant_counter, cities_counter;
  long loop_counter;
  long counter = 0;
  int* graph = NULL;
  double *pheromones;
  double* pheromonesUpdate;
  int* bestPath;
  int* otherBestPath;
  int* currentPath;
  long bestCost = 999999999;
  long otherBestCost;
  double* localPheromonesPath;
  double* otherPheromonesPath;
  long* randomNumbers;
  long nRandomNumbers = 0;
  long random_counter = 0;
  char* graphFile;
  char* randomFile;
  int nAnts;
  int totalNAnts;
  long externalIterations;
  long onNodeIteration;
  double alpha = 1;
  double beta = 1;
  double evaporation = 0.9;
  int nCities = 0;
  int* nAntsPerNode = (int*) malloc(size*sizeof(int));
  
  if (rank == 0) {
    std::ifstream in;
    randomFile = argv[2];
    char out[20];
	in.open(randomFile);
    in >> out;
    nRandomNumbers = atol(out);
    randomNumbers = (long*) malloc(nRandomNumbers*sizeof(long));
    for (i = 0; i < nRandomNumbers; i++) {
      randomNumbers[i] = 0;
    }
    i = 0;
    while(!in.eof()) {
      in >> out; //loading random.txt
      randomNumbers[i] = atol(out);
      i++;
    }
    in.close();
  }
  // Broadcasting number of random numbers
  MPI_Bcast(&nRandomNumbers, 1, MPI_LONG, 0, MPI_COMM_WORLD) ;
  if (rank != 0) {
    randomNumbers = (long*) malloc(nRandomNumbers*sizeof(long));
    for (i = 0; i < nRandomNumbers; i++) {
      randomNumbers[i] = 0;
    }
  }
  //Broadcasting all random numbers
  MPI_Bcast(&randomNumbers[0], nRandomNumbers, MPI_LONG, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    start_time = time();
  }
  
 // argument loading by master process
  if (rank == 0) {
    graphFile = argv[1];
    totalNAnts = atoi(argv[3]);
    externalIterations = atoi(argv[4]);
    onNodeIteration = atoi(argv[5]);
    printf("Total Iterations : %ld\n", externalIterations*onNodeIteration);
    printf("Number of Ants: %d\n", totalNAnts);
    std::ifstream in;
    in.open(graphFile);
    char out[12];
    in >> out;
    nCities = atoi(out);
    printf("Number of Cities: %d\n", nCities);
    graph = (int*) malloc(nCities*nCities*sizeof(int));
    in.close();
    loadGraph(graphFile, graph); // load graph.txt
  }
  
  // broadcasting all the variables
  MPI_Bcast(&totalNAnts, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&onNodeIteration, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&externalIterations, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&evaporation, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nCities, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank != 0) {
    graph = (int*) malloc(nCities*nCities*sizeof(int));
  }
  MPI_Bcast(&graph[0], nCities * nCities, MPI_INT, 0, MPI_COMM_WORLD);
  pheromones = (double*) malloc(nCities*nCities*sizeof(double));
  pheromonesUpdate = (double*) malloc(nCities*nCities*sizeof(double));
  localPheromonesPath = (double*) malloc(nCities*sizeof(double));
  otherBestPath = (int*) malloc(nCities*sizeof(int));
  otherPheromonesPath = (double*) malloc(nCities*sizeof(double));
  for (i = 0; i < nCities; i++) {
    otherBestPath[i] = -1;
  }
  bestPath = (int*) malloc(nCities*sizeof(int));
  currentPath = (int*) malloc(nCities*sizeof(int));
  for (i = 0; i < nCities; i++) {
    currentPath[i] = -1;
    bestPath[i] = -1;
  }
  for (i = 0; i < nCities * nCities; i++) {
    pheromones[i] = 0.1;
  }
	
  // Distributing the ants on different nodes
  int antsPerNode = totalNAnts / size;
  int restAnts = totalNAnts - antsPerNode * size;
  for (i = 0; i < size; i++) {
    nAntsPerNode[i] = antsPerNode;
    if (restAnts > i) {
      nAntsPerNode[i]++;
    }
  }
  nAnts = nAntsPerNode[rank];
  //number of ants in nodes before process
  int nAntsBeforeMe = 0;
  for (i = 0; i < size; i++) {
    if (i < rank) {
      nAntsBeforeMe += nAntsPerNode[i];
    } else {
      i = size;
    }
  }
  random_counter = (random_counter + (onNodeIteration * nAntsBeforeMe * nCities)) % nRandomNumbers;
  long antsBestCost = 999999999;
  
  while (counter < externalIterations) {
    loop_counter = 0;
    while (loop_counter < onNodeIteration) {
      for (ant_counter = 0; ant_counter < nAnts; ant_counter++) {
        for (i = 0; i < nCities; i++) {
          currentPath[i] = -1;
        }
        
        long rand = randomNumbers[random_counter];
        int currentCity = rand % nCities;	// selection of random start city
        random_counter = (random_counter + 1) % nRandomNumbers;
        currentPath[currentCity] = 0;
        for (cities_counter = 1; cities_counter < nCities; cities_counter++) {
          rand = randomNumbers[random_counter];
          random_counter = (random_counter + 1) % nRandomNumbers;
          currentCity = findNextCity(currentCity, currentPath, graph, nCities, pheromones, alpha, beta, rand); // next city to visit
          currentPath[currentCity] = cities_counter;
        }
        long oldCost = bestCost;
        bestCost = computeBestCost(bestCost, bestPath, currentPath, graph, nCities);
        if (oldCost > bestCost) {
          copyIntegerArray(currentPath, bestPath, nCities);
        }
      }
      if (bestCost < antsBestCost) {
        antsBestCost = bestCost;
     }
     for (j = 0; j < nCities * nCities; j++) {
        pheromones[j] *= evaporation;
     }
      updatePheromones(pheromones, bestPath, bestCost, nCities);
      loop_counter++;
    }
    pheromonesValue(localPheromonesPath, bestPath, pheromones, nCities); // find locally computed best path's pheromone value
    for (j = 0; j < nCities*nCities; j++) {
      pheromonesUpdate[j] = 1.0;
    }
    // Define temporary values
    long tempBestCost = bestCost;
    int* tempBestPath = (int*) malloc(nCities * sizeof(int));
    double* tempPheromonesPath= (double*) malloc(nCities * sizeof(double));
    copyIntegerArray(bestPath, tempBestPath, nCities);
    for (i = 0; i < size; i++) {
      if (rank == i) {
        copyIntegerArray(bestPath, otherBestPath, nCities);
        copyDoubleArray(localPheromonesPath, otherPheromonesPath, nCities);
        copyDoubleArray(localPheromonesPath, tempPheromonesPath, nCities);
        otherBestCost = bestCost;
      }
      MPI_Bcast(&otherBestPath[0], nCities, MPI_INT, i, MPI_COMM_WORLD);
      MPI_Bcast(&otherPheromonesPath[0], nCities, MPI_DOUBLE, i, MPI_COMM_WORLD);
      MPI_Bcast(&otherBestCost, 1, MPI_LONG, i, MPI_COMM_WORLD) ;
      // checking value of other nodes is better or not
      if (rank != i) {
        if (otherBestCost < tempBestCost) {
          tempBestCost = otherBestCost;
          copyIntegerArray(otherBestPath, tempBestPath, nCities);
          copyDoubleArray(otherPheromonesPath, tempPheromonesPath, nCities);
        } 
      }
    }
    // Update pheromones from best path received from other nodes
    for (j = 0; j < nCities - 1; j++) {
      pheromonesUpdate[getIndex(tempBestPath[j],tempBestPath[j+1],nCities)] += 1.0;
      pheromonesUpdate[getIndex(tempBestPath[j+1],tempBestPath[j],nCities)] += 1.0;
      pheromones[getIndex(tempBestPath[j],tempBestPath[j+1],nCities)] += tempPheromonesPath[j];
      pheromones[getIndex(tempBestPath[j+1],tempBestPath[j],nCities)] += tempPheromonesPath[j];
    }
    pheromonesUpdate[getIndex(tempBestPath[nCities-1],tempBestPath[0],nCities)] += 1.0;
    pheromonesUpdate[getIndex(tempBestPath[0],tempBestPath[nCities-1],nCities)] += 1.0;
    pheromones[getIndex(tempBestPath[nCities-1],tempBestPath[0],nCities)] += tempPheromonesPath[nCities - 1];
    pheromones[getIndex(tempBestPath[0],tempBestPath[nCities-1],nCities)] += tempPheromonesPath[nCities - 1];
    for (j = 0; j < nCities*nCities; j++) {
      pheromones[j] = pheromones[j] / pheromonesUpdate[j];
    }
    // Set own best cost with new cost
    bestCost = tempBestCost;
    copyIntegerArray(tempBestPath, bestPath, nCities);
    counter++;
    random_counter = (random_counter + (onNodeIteration * (totalNAnts - nAnts) * nCities)) % nRandomNumbers;
  }
  // Merging into the master process
  if (rank == 0) {
    for (i = 1; i < size; i++) {
      MPI_Recv(&otherBestPath[0], nCities, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      long oldCost = bestCost;
      bestCost = computeBestCost(bestCost, bestPath, otherBestPath, graph, nCities);
      if (oldCost > bestCost) {
        copyIntegerArray(otherBestPath, bestPath, nCities);
      }
    }
  } else {
    MPI_Send(&bestPath[0], nCities, MPI_INT, 0, 0, MPI_COMM_WORLD);
  }
  if (rank == 0) {
    printf("Best cost: %ld\n", bestCost);
    end_time = time();
    printf("Execution Time: %f\n", (end_time - start_time));
  }
  free(randomNumbers);
  free(graph);
  free(pheromones);
  free(localPheromonesPath);
  free(otherPheromonesPath);
  free(bestPath);
  free(otherBestPath);
  MPI_Finalize();
  return 0;
}


// compilation:  mpic++ parallel.cpp -o parallel.out
// execution: sbatch [slurm_file] 
//slurm => srun ./parallel.out [map)file] [random_file] [number of ants] [number of external iterations] [number of iterations on each node]
// e.g.     srun  ./parallel.out graph.txt random.txt 32 15000 1
