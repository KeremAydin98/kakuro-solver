#include <iostream>
#include <string>

#include <fstream>
#include <sstream>
#include <vector>

#include <bits/stdc++.h>
#include <array>

using namespace std;

enum direction {d_down, d_right, none};

#define COORD pair<int, int>

//#define DEBUG

int iter = 0;

//////////////////////////////////////////////
//Auxiliary functions for preparing problem //
//////////////////////////////////////////////

void display_arr(int* arr, int n){

  cout << "arr: ";

  for(int i = 0; i < n; i++){
    cout << arr[i] << " ";
  }

  cout << endl;
  
}

void print_coords(COORD start, COORD end){

  cout << "Start:" << start.first << "," << start.second << endl;
  cout << "End:" << end.first << "," << end.second << endl;
  
}

int find_length(COORD start, COORD end, direction dir){

  if(dir == d_down)
    return end.first - start.first;
  if(dir == d_right)
    return end.second - start.second;

  return -1;
}

void convert_sol(int** mat, int** &sol_mat, int m, int n){

  sol_mat = new int*[m]; //Rows
  for(int i = 0; i < m; i++){
    sol_mat[i] = new int[n]; //Cols
  }

  for(int i = 0; i < m; i++){
    for(int j = 0; j < m; j++){
      if(mat[i][j] == -2)
	sol_mat[i][j] = -2; //Empty value cell
      else
	sol_mat[i][j] = -1; //Hint or empty cell
    }
  }
}

void print_one_matrix(int** matrix, int m, int n){
  cout << "Matrix: " << endl;
  for(int i = 0; i < m; i++){ //rows
    for(int j = 0; j < n; j++){ //cols
      cout << matrix[i][j] << "\t";
    }
    cout << "\n";
  }
}

///Auxiliary functions

struct sum{
  COORD start;
  COORD end;

  int hint;
  int dir;
  int length;
  int* arr;

  void print_sum(){
    cout << "############################" << endl;
    cout << "Creating sum with: " << endl;
    print_coords(start, end);
    cout << "Hint: " << hint << endl;
    cout << "Direction: " << dir << endl;
    cout << "Length: " << length << endl;
    cout << "############################" << endl;
  }
  
  sum(COORD _start, COORD _end, int _hint, direction _dir):
    start(_start), end(_end), hint(_hint), dir(_dir)
  {
    length = find_length(_start, _end, _dir);
    arr = new int[length];
    #ifdef DEBUG
    cout << "############################" << endl;
    cout << "Creating sum with: " << endl;
    print_coords(start, end);
    cout << "Hint: " << hint << endl;
    cout << "Direction: " << dir << endl;
    cout << "Length: " << length << endl;
    cout << "############################" << endl;
    #endif
  }
  
  //~sum(){
  //delete arr;
  //}
};


COORD find_end(int** matrix, int m, int n, int i, int j, direction dir){ //0 down 1 right

  if(dir == d_right){
    for(int jj = j+1; jj < n; jj++){
      if(matrix[i][jj] != -2 || jj == n - 1){
        if(matrix[i][jj] == -2 && jj == n -1)
            jj++;
          COORD END = COORD(i, jj);
          return END;

      }
    }
  }
  else{// if(dir == d_down)
    for(int ii = i+1; ii < m; ii++){
      if(matrix[ii][j] != -2 || ii == m - 1){
        if(matrix[ii][j] == -2 && ii == m - 1)
                ii++;
              COORD END = COORD(ii, j);
              return END;
        }
      }
  }

  // no matching end found, return a default value
  COORD not_found = COORD(-1, -1);
  return not_found;
  
}


vector<sum> get_sums(int** matrix, int m, int n){

  vector<sum> sums;
  
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      int val = matrix[i][j];
      if(val != -1 && val != -2){
	int hint = val;
	hint = hint / 10;

	if((hint%100) == 0){
	  hint = (int)(hint/100);
	  COORD START = COORD(i, j+1); 
	  COORD END = find_end(matrix, m, n, i, j, d_right);
	  sum _sum = sum(START, END, hint, d_right);
	  sums.push_back(_sum);
	}

	else{
	  int div = (int)(hint/100);
	  int rem = (int)(hint%100);
   
	  if(div == 0 && rem != 0){
	    COORD START = COORD(i+1,j);
	    COORD END = find_end(matrix, m, n, i, j, d_down);
	    sum _sum = sum(START, END, rem, d_down);
	    sums.push_back(_sum);
	  }

	  if(div != 0 && rem != 0){
	    COORD START1 = COORD(i+1,j);
	    COORD START2 = COORD(i,j+1);
	    COORD END1 = find_end(matrix, m, n, i, j, d_down);
	    COORD END2 = find_end(matrix, m, n, i, j, d_right);
	    sum _sum1 = sum(START1, END1, rem, d_down);
	    sum _sum2 = sum(START2, END2, div, d_right);
	    sums.push_back(_sum1);
	    sums.push_back(_sum2);
	  }
	}
      }

      
    }
  }
  return sums;
}
  

void read_matrix(int** &matrix, ifstream &afile, int m, int n){

  matrix = new int*[m]; //rows

  for(int i = 0; i < m; i++){
    matrix[i] = new int[n]; //cols
  }

  int val;
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      afile >> val;
      matrix[i][j] = val;
    }
  }
}

void sol_to_file(int** mat, int** sol_mat, int m, int n){

  string fname = "visualize.kakuro";
  ofstream to_write(fname);

  to_write << m << " " << n << "\n";

  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      if(mat[i][j] != -2)
	to_write << mat[i][j] << " ";
      else
	to_write << sol_mat[i][j] << " ";
    }
    to_write << "\n";
  }

  to_write.close();
}

//////////////////////////////////////////////
//Auxiliary functions for preparing problem //
//////////////////////////////////////////////

///////////////////////////////////////////////////
//Auxiliary functions for preparing CUDA setting //
///////////////////////////////////////////////////

void flatten_sums(vector<sum> sums, int* h_sum_starts_x, int* h_sum_starts_y, int* h_sum_ends_x, int* h_sum_ends_y, int* h_sum_hints, int* h_sum_lengths, int* h_sum_dirs, int no_sums){

  for(int i = 0; i < no_sums; i++){
    
    h_sum_starts_x[i] = sums[i].start.first;
    h_sum_starts_y[i] = sums[i].start.second;

    h_sum_ends_x[i] = sums[i].end.first;
    h_sum_ends_y[i] = sums[i].end.second;
    
    h_sum_hints[i] = sums[i].hint;
    h_sum_lengths[i] = sums[i].length;
    
    h_sum_dirs[i] = sums[i].dir;
  }
  
}

void print_flattened(int* h_sum_starts_x, int* h_sum_starts_y, int* h_sum_ends_x, int* h_sum_ends_y, int* h_sum_hints, int* h_sum_lengths, int* h_sum_dirs, int no_sums){

  cout << "###h_sum_starts_x: " << endl;
  for(int i = 0; i < no_sums; i++){
    cout << h_sum_starts_x[i] << " ";
  }
  cout << endl;

  cout << "###h_sum_starts_y: " << endl;
  for(int i = 0; i < no_sums; i++){
    cout << h_sum_starts_y[i] << " ";
  }
  cout << endl;

  cout << "###h_sum_ends_x: " << endl;
  for(int i = 0; i < no_sums; i++){
    cout << h_sum_ends_x[i] << " ";
  }
  cout << endl;

  cout << "###h_sum_ends_y: " << endl;
  for(int i = 0; i < no_sums; i++){
    cout << h_sum_ends_y[i] << " ";
  }
  cout << endl;

  cout << "###h_sum_hints: " << endl;
  for(int i = 0; i < no_sums; i++){
    cout << h_sum_hints[i] << " ";
  }
  cout << endl;

  cout << "###h_sum_lengths: " << endl;
  for(int i = 0; i < no_sums; i++){
    cout << h_sum_lengths[i] << " ";
  }
  cout << endl;

  cout << "###h_sum_dirs: " << endl;
  for(int i = 0; i < no_sums; i++){
    cout << h_sum_dirs[i] << " ";
  }
  cout << endl;
  
}

void flatten_sol_mat(int** sol_mat, int* d_sol_mat, int m, int n){
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      d_sol_mat[i*n+j] = sol_mat[i][j];
    }
  }  
}

void print_flattened_matrix(int* d_sol_mat, int m, int n){

  cout << "###Flattened matrix: " << endl;
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      cout << d_sol_mat[i*n+j] << " ";
    }
    cout << endl;
  }
  cout << endl;
}

///////////////////////////////////////////////////
//Auxiliary functions for preparing CUDA setting //
///////////////////////////////////////////////////


///////////////////
//CUDA FUNCTIONS //
///////////////////

void init_iteration(int** iteration, int** sol_mat, int m, int n){

  for(int i = 0; i < m; i++){
    iteration[i] = new int[n];
  }
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      if(sol_mat[i][j] == -2){
        iteration[i][j] = 0;
      }
      else{
        iteration[i][j] = -1;
      }
    }
  }
}

void flatten_iteration(int** h_iteration, int* d_iteration, int m, int n){
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      d_iteration[i*n+j] = h_iteration[i][j];
    }
  }  
}

/*

__device__ vector<int> remove_unusable_values(int* d_sol_mat, int i, int j, int k,
                                   int* d_sum_starts_x, int* d_sum_starts_y, int* d_sum_ends_x,
                                   int* d_sum_ends_y, int* d_sum_hints,
                                   int m, int n, vector<sum> sums){

  vector<int> possible_values = {1,2,3,4,5,6,7,8,9};

  // smaller than minimum hint
  // different than the values in row and column
  // minimum value according to length
  std::vector<int> will_remove = {};
  for(int kk=d_sum_starts_x[k]; kk<d_sum_ends_x[k]; kk++){
    will_remove.push_back(d_sol_mat[i * n + j]);
  }
  for(int ll=d_sum_starts_x[k]; ll<d_sum_ends_x[k]; ll++){
    will_remove.push_back(d_sol_mat[i * n + j]);
  }
  for(int mm=0; mm<possible_values.size(); mm++)
  {
    if((possible_values[mm] >= d_sum_hints[0] || possible_values[mm] >= d_sum_hints[1])){
      will_remove.push_back(possible_values[mm]);
    }
  }

  vector<int>::iterator it;
  for (it = will_remove.begin(); it != will_remove.end(); it++) {
      possible_values.erase(remove(possible_values.begin(), possible_values.end(), *it), possible_values.end());
  }

  return possible_values;
}

__device__ void fill_sum(int* d_sol_mat, vector<int> possible_values, int* d_sum_hints, int* d_sum_dirs, int i, int j, int k, int m, int n){

  int summation = 0;
  
  if(d_sum_dirs[k] == 0){

    for(int kk=0; kk<j;kk++){
      if(d_sol_mat[kk * n + j] != -1){ 
        summation += d_sol_mat[kk * n + j];
      }
    }
  }
  else{

    for(int kk=0; kk<j;kk++){
      if(d_sol_mat[i * n + kk] != -1){ 
        summation += d_sol_mat[i * n + kk];
      }
    }
  }

  int last_value = d_sum_hints[k] - summation;

  if(find(possible_values.begin(), possible_values.end(), last_value) != possible_values.end()){

    d_sol_mat[i * n + j] = last_value;

  }
  
}

__device__ bool hasRepetitiveValues(const std::vector<int>& vec) {
    std::set<int> uniqueElements;
    for (int element : vec) {
        if (uniqueElements.count(element) > 0) {
            // Element already exists in the set
            return true;
        }
        uniqueElements.insert(element);
    }
    return false;
}

__device__ bool check_solution(int* d_sol_mat, int no_sums, int m, int n,
                    int* d_sum_starts_x, int* d_sum_starts_y, int* d_sum_ends_x,
                    int* d_sum_ends_y, int* d_sum_hints, int* d_sum_lengths, int* d_sum_dirs){
  
  Confirms the solution to see if it is correct or not
  
  for(int k=0; k<no_sums; k++){
    int summation = 0;
    vector<int> repetitive_or_not;
    if(d_sum_dirs[k] == 0)
    {
      for(int j=d_sum_starts_x[k]; j<d_sum_ends_x[k]; j++)
      {
        summation += d_sol_mat[d_sum_starts_y[k] * n + j];
        repetitive_or_not.push_back(d_sol_mat[d_sum_starts_y[k] * n + j]);
      }
    }
    else{
      for(int j=d_sum_starts_y[k]; j<d_sum_ends_y[k]; j++)
      {
        summation += d_sol_mat[j * n + d_sum_starts_x[k]];
        repetitive_or_not.push_back(d_sol_mat[d_sum_starts_y[k] * n + j]);
      }
    }

    if(d_sum_hints[k] != summation){
        return false;
    }

    if (hasRepetitiveValues(repetitive_or_not)){
      return false;
    }
  }

  return true;

}


__global__ void kakuro_kernel(int* d_sum_starts_x, int* d_sum_starts_y, int* d_sum_ends_x,
                              int* d_sum_ends_y, int* d_sum_hints, int* d_sum_lengths,
                              int* d_sum_dirs, int* d_sol_mat, 
                              int m, int n, int no_sums, volatile bool* solved,
                              int* d_iteration, vector<sum> sums){

  
  int tid = blockDim.x * blockIdx.x + threadIdx.x;

  if(tid < no_sums)
  {
    for(int i=d_sum_starts_x[tid]; i<d_sum_ends_x[tid]; i++)
    {

      for(int j=d_sum_starts_y[tid]; j<d_sum_ends_y[tid]; j++)
      {
 
        vector<int> possible_values = remove_unusable_values(d_sol_mat, i, j, tid,
                                                             d_sum_starts_x, d_sum_starts_y, 
                                                             d_sum_ends_x,
                                                             d_sum_ends_y, d_sum_hints,
                                                             m, n, sums);

        if(possible_values.size() == 0){
          d_iteration[i * n + j] += 1;
          continue;
        }

        if((i == d_sum_ends_x[tid]) && (d_sum_dirs[tid] == 0)){
          int first_value = d_sol_mat[i * n + j];
          fill_sum(d_sol_mat, possible_values, d_sum_hints, d_sum_dirs, i, j, tid, m, n);
          if(d_sol_mat[i * n + j] != first_value){
            continue;
          }
        }
        else if((i == d_sum_ends_y[tid]) && (d_sum_dirs[tid] == 1)){
          int first_value = d_sol_mat[i * n + j];
          fill_sum(d_sol_mat, possible_values, d_sum_hints, d_sum_dirs, i, j, tid, m, n);
          
          if(d_sol_mat[i * n + j] != first_value){
            continue;
          }
        }

        int which_value = d_iteration[i * n + j] % possible_values.size();


        d_sol_mat[i * n + j] = possible_values[which_value];

        d_iteration[i * n + j] += 1;

        if(check_solution(d_sol_mat, no_sums, m, n,
                          d_sum_starts_x, d_sum_starts_y, d_sum_ends_x,
                          d_sum_ends_y, d_sum_hints, d_sum_lengths, d_sum_dirs)){
            *solved = true;
        }
      }
    }
  } 
}
  //About volatile bool* solved:
  //You can get idea from https://stackoverflow.com/questions/12505750/how-can-a-global-function-return-a-value-or-break-out-like-c-c-does%5B/url%5D for how to break out of a CUDA kernel
  //You may or may not use it


///////////////////
//CUDA FUNCTIONS //
///////////////////


__global__ void kakuro_kernel(int* d_sum_starts_x, int* d_sum_starts_y, int* d_sum_ends_x,
                              int* d_sum_ends_y, int* d_sum_hints, int* d_sum_lengths,
                              int* d_sum_dirs, int* d_sol_mat, 
                              int m, int n, int no_sums, volatile bool* solved,
                              int* d_iteration, vector<sum> sums){

  
  int tid = blockDim.x * blockIdx.x + threadIdx.x;

  if(tid < no_sums)
  {
    for(int i=d_sum_starts_x[tid]; i<d_sum_ends_x[tid]; i++)
    {

      for(int j=d_sum_starts_y[tid]; j<d_sum_ends_y[tid]; j++)
      {
 
        // ==============================================================
        vector<int> possible_values = {1,2,3,4,5,6,7,8,9};

        // smaller than minimum hint
        // different than the values in row and column
        // minimum value according to length
        std::vector<int> will_remove = {};
        for(int kk=d_sum_starts_x[tid]; kk<d_sum_ends_x[tid]; kk++){
          will_remove.push_back(d_sol_mat[i * n + j]);
        }
        for(int ll=d_sum_starts_x[tid]; ll<d_sum_ends_x[tid]; ll++){
          will_remove.push_back(d_sol_mat[i * n + j]);
        }
        for(int mm=0; mm<possible_values.size(); mm++)
        {
          if((possible_values[mm] >= d_sum_hints[0] || possible_values[mm] >= d_sum_hints[1])){
            will_remove.push_back(possible_values[mm]);
          }
        }

        for (int value : will_remove) {
            possible_values.erase(remove(possible_values.begin(), possible_values.end(), value), possible_values.end());
        }

        // =============================================================

        if(possible_values.size() == 0){
          d_iteration[i * n + j] += 1;
          continue;
        }

        if((i == d_sum_ends_x[tid]) && (d_sum_dirs[tid] == 0)){
          int first_value = d_sol_mat[i * n + j];

          // =============================================================
          int summation = 0;
          
          if(d_sum_dirs[tid] == 0){

            for(int kk=0; kk<j;kk++){
              if(d_sol_mat[kk * n + j] != -1){ 
                summation += d_sol_mat[kk * n + j];
              }
            }
          }
          else{

            for(int kk=0; kk<j;kk++){
              if(d_sol_mat[i * n + kk] != -1){ 
                summation += d_sol_mat[i * n + kk];
              }
            }
          }

          int last_value = d_sum_hints[tid] - summation;

          if(find(possible_values.begin(), possible_values.end(), last_value) != possible_values.end()){

            d_sol_mat[i * n + j] = last_value;

          }
          
          // =============================================================
          if(d_sol_mat[i * n + j] != first_value){
            continue;
          }
        }
        else if((i == d_sum_ends_y[tid]) && (d_sum_dirs[tid] == 1)){
          int first_value = d_sol_mat[i * n + j];
          // =============================================================
          int summation = 0;
          
          if(d_sum_dirs[tid] == 0){

            for(int kk=0; kk<j;kk++){
              if(d_sol_mat[kk * n + j] != -1){ 
                summation += d_sol_mat[kk * n + j];
              }
            }
          }
          else{

            for(int kk=0; kk<j;kk++){
              if(d_sol_mat[i * n + kk] != -1){ 
                summation += d_sol_mat[i * n + kk];
              }
            }
          }

          int last_value = d_sum_hints[tid] - summation;

          if(find(possible_values.begin(), possible_values.end(), last_value) != possible_values.end()){

            d_sol_mat[i * n + j] = last_value;

          }
          // =============================================================
          if(d_sol_mat[i * n + j] != first_value){
            continue;
          }
        }

        int which_value = d_iteration[i * n + j] % possible_values.size();


        d_sol_mat[i * n + j] = possible_values[which_value];

        d_iteration[i * n + j] += 1;

        // =============================================================
        bool condition = true;
        for(int k=0; k<no_sums; k++){
          int summation = 0;
          vector<int> repetitive_or_not;
          if(d_sum_dirs[k] == 0)
          {
            for(int j=d_sum_starts_x[k]; j<d_sum_ends_x[k]; j++)
            {
              summation += d_sol_mat[d_sum_starts_y[k] * n + j];
              repetitive_or_not.push_back(d_sol_mat[d_sum_starts_y[k] * n + j]);
            }
          }
          else{
            for(int j=d_sum_starts_y[k]; j<d_sum_ends_y[k]; j++)
            {
              summation += d_sol_mat[j * n + d_sum_starts_x[k]];
              repetitive_or_not.push_back(d_sol_mat[d_sum_starts_x[k] * n + j]);
            }
          }

          if(d_sum_hints[k] != summation){
              condition = false;
          }

          bool has_repetitive = false;
          set<int> uniqueElements;
          for (int element : repetitive_or_not) {
              if (uniqueElements.count(element) > 0) {
                  // Element already exists in the set
                  has_repetitive = true;
              }
              uniqueElements.insert(element);
          }

          if (has_repetitive){
            condition = false;
          }
        }

        if(condition){
            *solved = true;
        }
        // =============================================================
      }
    }
  } 
}
*/
__global__ void kakuro_kernel(int* d_sum_starts_x, int* d_sum_starts_y, int* d_sum_ends_x,
                              int* d_sum_ends_y, int* d_sum_hints, int* d_sum_lengths,
                              int* d_sum_dirs, int* d_sol_mat, 
                              int m, int n, int no_sums, volatile bool* solved,
                              int* d_iteration) {
  
  int tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < no_sums) {
    for (int i = d_sum_starts_x[tid]; i < d_sum_ends_x[tid]; i++) {

      for (int j = d_sum_starts_y[tid]; j < d_sum_ends_y[tid]; j++) {
 
        // ==============================================================
        int possible_values[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        int possible_values_size = sizeof(possible_values) / sizeof(int);

        // smaller than minimum hint
        // different than the values in row and column
        // minimum value according to length
        int will_remove[9];
        int will_remove_count = 0;
        for (int kk = d_sum_starts_x[tid]; kk < d_sum_ends_x[tid]; kk++) {
          will_remove[will_remove_count++] = d_sol_mat[i * n + j];
        }
        for (int ll = d_sum_starts_x[tid]; ll < d_sum_ends_x[tid]; ll++) {
          will_remove[will_remove_count++] = d_sol_mat[i * n + j];
        }
        for (int mm = 0; mm < possible_values_size; mm++) {
          if ((possible_values[mm] >= d_sum_hints[0] || possible_values[mm] >= d_sum_hints[1])) {
            will_remove[will_remove_count++] = possible_values[mm];
          }
        }

        for (int k = 0; k < will_remove_count; k++) {
          int value = will_remove[k];
          for (int idx = 0; idx < possible_values_size; idx++) {
            if (possible_values[idx] == value) {
              for (int m = idx + 1; m < possible_values_size; m++) {
                possible_values[m - 1] = possible_values[m];
              }
              possible_values_size--;
              break;
            }
          }
        }
        // =============================================================

        if (possible_values_size == 0) {
          atomicAdd(&d_iteration[i * n + j], 1);
          continue;
        }

        if ((i == d_sum_ends_x[tid]) && (d_sum_dirs[tid] == 0)) {
          int first_value = d_sol_mat[i * n + j];

          // =============================================================
          int summation = 0;
          if (d_sum_dirs[tid] == 0) {
            for (int kk = 0; kk < j; kk++) {
              if (d_sol_mat[kk * n + j] != -1) { 
                summation += d_sol_mat[kk * n + j];
              }
            }
          }
          else {
            for (int kk = 0; kk < j; kk++) {
              if (d_sol_mat[i * n + kk] != -1) { 
                summation += d_sol_mat[i * n + kk];
              }
            }
          }

          int last_value = d_sum_hints[tid] - summation;

          bool found = false;
          for (int idx = 0; idx < possible_values_size; idx++) {
            if (possible_values[idx] == last_value) {
              found = true;
              break;
            }
          }

          if (found) {
            d_sol_mat[i * n + j] = last_value;
          }
          // =============================================================
          if (d_sol_mat[i * n + j] != first_value) {
            continue;
          }
        }
        else if ((i == d_sum_ends_y[tid]) && (d_sum_dirs[tid] == 1)) {
          int first_value = d_sol_mat[i * n + j];
          // =============================================================
          int summation = 0;
          if (d_sum_dirs[tid] == 0) {
            for (int kk = 0; kk < j; kk++) {
              if (d_sol_mat[kk * n + j] != -1) { 
                summation += d_sol_mat[kk * n + j];
              }
            }
          }
          else {
            for (int kk = 0; kk < j; kk++) {
              if (d_sol_mat[i * n + kk] != -1) { 
                summation += d_sol_mat[i * n + kk];
              }
            }
          }

          int last_value = d_sum_hints[tid] - summation;

          bool found = false;
          for (int idx = 0; idx < possible_values_size; idx++) {
            if (possible_values[idx] == last_value) {
              found = true;
              break;
            }
          }

          if (found) {
            d_sol_mat[i * n + j] = last_value;
          }
          // =============================================================
          if (d_sol_mat[i * n + j] != first_value) {
            continue;
          }
        }

        int which_value = d_iteration[i * n + j] % possible_values_size;

        d_sol_mat[i * n + j] = possible_values[which_value];

        atomicAdd(&d_iteration[i * n + j], 1);

        // =============================================================
        bool condition = true;
        for (int k = 0; k < no_sums; k++) {
          int summation = 0;
          int repetitive_or_not[100];  // Assuming a maximum size of repetitive_or_not
          int repetitive_or_not_count = 0;
          if (d_sum_dirs[k] == 0) {
            for (int j = d_sum_starts_x[k]; j < d_sum_ends_x[k]; j++) {
              summation += d_sol_mat[d_sum_starts_y[k] * n + j];
              repetitive_or_not[repetitive_or_not_count++] = d_sol_mat[d_sum_starts_y[k] * n + j];
            }
          }
          else {
            for (int j = d_sum_starts_y[k]; j < d_sum_ends_y[k]; j++) {
              summation += d_sol_mat[j * n + d_sum_starts_x[k]];
              repetitive_or_not[repetitive_or_not_count++] = d_sol_mat[d_sum_starts_x[k] * n + j];
            }
          }

          if (d_sum_hints[k] != summation) {
            condition = false;
          }

          bool has_repetitive = false;
          for (int i = 0; i < repetitive_or_not_count; i++) {
            for (int j = i + 1; j < repetitive_or_not_count; j++) {
              if (repetitive_or_not[i] == repetitive_or_not[j]) {
                has_repetitive = true;
                break;
              }
            }
            if (has_repetitive) {
              break;
            }
          }

          if (has_repetitive) {
            condition = false;
          }
        }

        if (condition) {
          *solved = true;
        }
        // =============================================================
      }
    }
  } 
}


int main(int argc, char** argv){
  
  string filename(argv[1]);
  ifstream file;
  file.open(filename.c_str());

  int m, n;
  file >> m;
  file >> n;

  int** mat;
  read_matrix(mat, file, m, n);
  print_one_matrix(mat, m, n);
  
  int** h_sol_mat;
  convert_sol(mat, h_sol_mat, m, n);
  print_one_matrix(h_sol_mat, m, n);
  
  vector<sum> sums = get_sums(mat, m, n);
  
  //CUDA
  cudaDeviceProp prop; // cudaDeviceProp prop; declares a variable prop of type cudaDeviceProp, which is a structure that holds information about a CUDA device.
  cudaGetDeviceProperties(&prop, 0); // retrieves the properties of the CUDA device with the device ID 0 and stores the information in the prop variable
  printf("==prop== Running on device: %d -- %s \n", 0, prop.name);
  printf("==prop== #of SM -- %d \n", prop.multiProcessorCount);
  printf("==prop== Max Threads Per Block: -- %d \n", prop.maxThreadsPerBlock);

  //To DO 
  // =========================================
  int block_dim = 1024; 
  int grid_dim = (int)ceil(n/block_dim);;  

  // =========================================

  int no_sums = sums.size();

  //Flattening sums and matrix
  int* h_sum_starts_x = new int[no_sums];
  int* h_sum_starts_y = new int[no_sums];
  int* h_sum_ends_x = new int[no_sums];
  int* h_sum_ends_y = new int[no_sums];
  int* h_sum_hints = new int[no_sums];
  int* h_sum_lengths = new int[no_sums];
  int* h_sum_dirs = new int[no_sums];

  // Pair to integers
  flatten_sums(sums, h_sum_starts_x, h_sum_starts_y, h_sum_ends_x, h_sum_ends_y, h_sum_hints, h_sum_lengths, h_sum_dirs, no_sums);

  // Print flattened vector
  print_flattened(h_sum_starts_x, h_sum_starts_y, h_sum_ends_x, h_sum_ends_y, h_sum_hints, h_sum_lengths, h_sum_dirs, no_sums);

  int* d_sol_mat;
  int* d_iteration = new int[m*n];
  d_sol_mat = new int[m*n];
  flatten_sol_mat(h_sol_mat, d_sol_mat, m, n);

  print_flattened_matrix(d_sol_mat, m, n);

  //Declare device pointers and copy data into device
  int *d_sum_starts_x, *d_sum_starts_y, *d_sum_ends_x, *d_sum_ends_y, *d_sum_hints, *d_sum_lengths, *d_sum_dirs, *d_t_mats;

  // ITERATION MATRIX
  int** h_iteration = new int*[m];
  init_iteration(h_iteration, h_sol_mat, m, n);
  flatten_iteration(h_iteration, d_iteration, m, n);
  // ==============================

  cudaMalloc(&d_sum_starts_x, no_sums*sizeof(int));
  cudaMalloc(&d_sum_starts_y, no_sums*sizeof(int));
  cudaMalloc(&d_sum_ends_x, no_sums*sizeof(int));
  cudaMalloc(&d_sum_ends_y, no_sums*sizeof(int));
  cudaMalloc(&d_sum_hints, no_sums*sizeof(int));
  cudaMalloc(&d_sum_lengths, no_sums*sizeof(int));
  cudaMalloc(&d_sum_dirs, no_sums*sizeof(int));
  cudaMalloc(&d_sol_mat, (m*n)*sizeof(int));
  cudaMalloc(&d_iteration, (m*n)*sizeof(int));
  cudaMalloc(&d_t_mats, (m * n * grid_dim * block_dim)*sizeof(int)); //Allocating invidual matrix for each GPU thread
  //You may use this array if you will implement a thread-wise solution

  cudaMemcpy(d_sum_starts_x, h_sum_starts_x, no_sums*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sum_starts_y, h_sum_starts_y, no_sums*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sum_ends_x, h_sum_ends_x, no_sums*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sum_ends_y, h_sum_ends_y, no_sums*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sum_hints, h_sum_hints, no_sums*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sum_lengths, h_sum_lengths, no_sums*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sum_dirs, h_sum_dirs, no_sums*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sol_mat, h_sol_mat, (m*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_iteration, h_iteration, (m*n)*sizeof(int), cudaMemcpyHostToDevice);
  
  bool* solved = new bool;
  *solved = false;
  bool* d_solved = new bool;
  
  cudaMalloc(&d_solved, sizeof(bool));
  cudaMemcpy(d_solved, solved, sizeof(bool), cudaMemcpyHostToDevice);
  

  
  
  // CUDA kernel
  kakuro_kernel<<<grid_dim, block_dim>>>(d_sum_starts_x, d_sum_starts_y, d_sum_ends_x,
                                         d_sum_ends_y, d_sum_hints, d_sum_lengths, d_sum_dirs, 
                                         d_sol_mat, m, n, no_sums, d_solved,
                                         d_iteration);
  // ===============================
  cudaDeviceSynchronize();
  //CUDA
  
  
  print_flattened_matrix(d_sol_mat, m, n);
  //TO DO sol_mat_flattened_to_file(mat, d_sol_mat, m, n)
  //Similiar to sol_mat, use hints from mat and values from d_sol_mat
  
  for(int i = 0; i < n; i++){
    delete mat[i];
    delete h_sol_mat[i];
  }

  delete mat;
  delete h_sol_mat;
  delete d_iteration;

  delete h_sum_starts_x;
  delete h_sum_starts_y;
  delete h_sum_ends_x;
  delete h_sum_ends_y;
  delete h_sum_hints;
  delete h_sum_lengths;
  delete h_sum_dirs;
  delete d_sol_mat;

  cudaFree(d_t_mats);
  cudaFree(d_sum_starts_x);
  cudaFree(d_sum_starts_y);
  cudaFree(d_sum_ends_x);
  cudaFree(d_sum_ends_y);
  cudaFree(d_sum_hints);
  cudaFree(d_sum_lengths);
  cudaFree(d_sum_dirs);
  cudaFree(d_sol_mat);
  cudaFree(d_iteration);
  
  
  return 0;
}
