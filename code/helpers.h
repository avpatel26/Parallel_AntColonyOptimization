#include <cmath> 
#include <cstdlib> 
#include <iostream> 
#include <fstream> 
#include <math.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <sys/time.h> 
#include <limits> 
#include <string.h> 
#include <float.h>
double start_time, end_time; 
int getIndex(int i, int j, int n) {
  return (i * n) + j;
}
void copyIntegerArray(int* in, int* out, int n) {
  int i;
  for (i = 0; i < n; i++) {
    out[i] = (int) in[i];
  }
}
void copyDoubleArray(double* in, double* out, int n) {
  int i;
  for (i = 0; i < n; i++) {
    out[i] = (double) in[i];
  }
}
// returns the pheromones value of path 
void pheromonesValue (double* pheromonesPath, int* bestPath, double* pheromones, int n) {
  int i;
  int previousCity = 0;
  int nextCity = 0;
  for (i = 1; i < n; i++) {
    previousCity = bestPath[i - 1];
    nextCity = bestPath[i];
    pheromonesPath[i-1] = pheromones[getIndex(previousCity, nextCity, n)];
  }
  pheromonesPath[n - 1] = pheromones[getIndex(nextCity, bestPath[0], n)];
}
// Load map from graph.txt 
int loadGraph(char* file, int* map) {
  std::ifstream in;
  int matrixFull = 1;
  int size = 0;
  int x = 0;
  int y = 0;
  int first = 1;
  in.open (file);
  char out[8];
  while (!in.eof() && matrixFull) {
    if (first) {
      first = 0;
      in >> out;
      size = atoi(out);
    } else {
      in >> out;
      int index = getIndex(x,y,size);
      map[index]= atoi(out);
      y = (y + 1) % size;
      if (y == 0) {
        x = (x + 1) % size;
        if (x == 0) {
          matrixFull = 0;
        }
      }
    }
  }
  in.close();
  return 0;
}
// find the probability to visit next city 
void probability(int currentCity, double* probabilities, int* path, int* map, int nCities, double* pheromones, double alpha, double beta) {
  int i;
  double total = 0;
  for (i = 0; i < nCities; i++) {
    if (path[i] != -1 || i == currentCity) {
      probabilities[i] = 0.0;
    } else {
      double p = pow(1.0 / map[getIndex(currentCity,i,nCities)],alpha) * pow(pheromones[getIndex(currentCity,i,nCities)], beta);
      probabilities[i] = p;
      total += p;
    }
  }
  if (total == 0) {
    i = 0;
    for (i = 0; i < nCities; i++) {
      if (path[i] == -1 && i != currentCity) {
        probabilities[i] = 1.0;
        total++;
      }
    }
  }
  for (i = 0; i < nCities; i++) {
    probabilities[i] = probabilities[i] / total;
  }
}
// select next city to visit based on probabilities 
int findNextCity(int currentCity, int* path, int* map, int nCities, double* pheromones, double alpha, double 
beta, long rand) {
  int i = 0;
  double *probabilities;
  probabilities = (double*) malloc(nCities*sizeof(double));
  probability(currentCity, probabilities, path, map, nCities, pheromones, alpha, beta);
  int value = (rand % 100) + 1;
  int sum = 0;
  for (i = 0; i < nCities; i++) {
    sum += ceilf(probabilities[i] * 100);
    if (sum >= value) {
      free(probabilities);
      return i;
    }
  }
  free(probabilities);
  return -1;
}
// compute the best cost 
long computeBestCost(long bestCost, int* bestPath, int* currentPath, int* map, int nCities) {
  int i;
  long currentCost = 0;
  int* orderedCities = (int*) malloc(nCities*sizeof(int));
  for (i = 0; i < nCities; i++) {
    orderedCities[currentPath[i]] = i;
  }
  for (i = 0; i < nCities - 1; i++) {
    currentCost += map[getIndex(orderedCities[i],orderedCities[i + 1], nCities)];
  }
  currentCost += map[getIndex(orderedCities[nCities - 1], orderedCities[0], nCities)];
  if (bestCost > currentCost) {
    return currentCost;
  } else {
    return bestCost;
  }
}
//update the pheromone value in pheromone matrix 
void updatePheromones(double* pheromones, int* path, long cost, int n) {
  int i;
  int* orderedCities = (int*) malloc(n*sizeof(int));
  for (i = 0; i < n; i++) {
    int order = path[i];
    orderedCities[order] = i;
  }
  for (i = 0; i < n - 1; i++) {
    pheromones[getIndex(orderedCities[i],orderedCities[i + 1], n)] += 1.0/cost;
    pheromones[getIndex(orderedCities[i + 1],orderedCities[i], n)] += 1.0/cost;
    if (pheromones[getIndex(orderedCities[i],orderedCities[i + 1], n)] > 1) {
      pheromones[getIndex(orderedCities[i],orderedCities[i + 1], n)] = 1.0;
      pheromones[getIndex(orderedCities[i + 1],orderedCities[i], n)] = 1.0;
    }
  }
  pheromones[getIndex(orderedCities[n - 1],orderedCities[0],n)] += 1.0/cost;
  pheromones[getIndex(orderedCities[0],orderedCities[n - 1], n)] += 1.0/cost;
  if (pheromones[getIndex(orderedCities[n - 1],orderedCities[0],n)] > 1.0) {
    pheromones[getIndex(orderedCities[n - 1],orderedCities[0],n)] = 1.0;
    pheromones[getIndex(orderedCities[0],orderedCities[n - 1], n)] = 1.0;
  }
}
// return current time 
double time() {
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv,&tz);
  return ( (double) tv.tv_sec + (double) tv.tv_usec * 1.e-6 );
}
