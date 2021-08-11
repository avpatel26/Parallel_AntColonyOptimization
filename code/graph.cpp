// Random graph generation with specific number of cities

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>

int main(int argc, char* argv[]) {
  std::ofstream file ("graph.txt"); //input file
  int n = 500; //number of cities
  int maximumDistance = 1000; //maximum distance
  int graph[n][n];
  int i, j;
  file << n <<std::endl;
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      if (i == j) {
        graph[i][j] = 0;
      } else {
        graph[i][j] = (rand() % maximumDistance) + 1;
        graph[j][i] = graph[i][j];
      }
    }
 }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      file << graph[i][j] << " "; 
    }
    file << std::endl;
  }
  file.close();
  return 0;
}

//Compile the program: g++ graph.cpp -o graph.out
//execute the program: ./graph.out
