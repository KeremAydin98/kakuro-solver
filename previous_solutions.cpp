#include <iostream>
#include <string>

#include <fstream>
#include <sstream>
#include <vector>

#include <bits/stdc++.h>
#include <array>

#include <random>
#include <list>
#include <omp.h>

void random_with_constraint_solution(int** mat, int** sol_mat, vector<sum> sums, int m, int n){

  // Get the starting timepoint
  double start_time = omp_get_wtime();

  while(true){
    for(int i = 0; i <= m; i++){ 

      vector<int> possible_values = {1, 2, 3, 4, 5, 6, 7, 8, 9};

        if(sums[i].dir == 0){ // down direction

        int column = sums[i].start.second;

        for(int j = 1; j < m; j++){ //rows

          // Remove values which are present in the column or larger than sum
          remove_unusable_values(possible_values, sol_mat, sums, column, m, n, i, j);

          random_device rd;
          mt19937 gen(rd());
          uniform_int_distribution<> dis(0, distance(possible_values.begin(), possible_values.end()) - 1);

          auto it = possible_values.begin();
          advance(it, dis(gen));
          int random_element = *it;

          sol_mat[j][column] = random_element;

          }
        }
        else{ // right direction

          int row = sums[i].start.first;

          for(int j = 1; j < m; j++){ //rows

            // Remove used values which are in the row or larger than sum
            remove_unusable_values(possible_values, sol_mat, sums, row, m, n, i, j);

            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<> dis(0, distance(possible_values.begin(), possible_values.end()) - 1);
            auto it = possible_values.begin();
            advance(it, dis(gen));
            int random_element = *it;

            sol_mat[row][j] = random_element;

            }
        }
        if(possible_values.size() == 0){
          break;
        }
      }

      if(check_solution(sol_mat, sums, m, n)){
        break;
      }
      
    }

  double end_time = omp_get_wtime();
  double elapsed_time = end_time - start_time;

  printf("\nElapsed time: %f\n", elapsed_time);

}

void random_without_constraint_solution(int** mat, int** sol_mat, vector<sum> sums, int m, int n){

  // Get the starting timepoint
  double start_time = omp_get_wtime();

  while(true){
   
    initialize_rand_matrix(sol_mat, m, n);

    if(check_solution(sol_mat, sums, m, n)){
      break;
    }

  }

  double end_time = omp_get_wtime();
  double elapsed_time = end_time - start_time;

  printf("\nElapsed time: %f\n", elapsed_time);

}

void initialize_rand_matrix(int** matrix, int m, int n){

  for(int i = 0; i < m; i++){ //rows
    for(int j = 0; j < n; j++){ //cols
      if(i==0 || j==0){
        matrix[i][j] = -1;
      }
      else{
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<int> dis(1, 9);
        matrix[i][j] = dis(gen);
          }
      }
    }

}