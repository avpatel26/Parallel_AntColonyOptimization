#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include <iostream> 
#include <fstream> 
#include <limits.h>

int main(int argc, char* argv[]) {
  std::ofstream file ("random.txt");
  int n = atoi(argv[1]);
  int array[n];
  int i;
  for (i = 0; i < n; i++) {
    array[i] = (rand() % INT_MAX);
  }
  file << n << std::endl;
  for (i = 0; i < n; i++) {
    file << array[i] << " ";
  }
  file << std::endl;
  file.close();
  return 0;
}

// compilation of program : g++ random.cpp -o random.out
// execution of program: ./random.out [number of random numbers]
// e.g. : ./random.out 100000