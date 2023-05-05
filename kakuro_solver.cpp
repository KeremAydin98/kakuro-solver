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

///Auxiliary functions

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

void sol_to_file(int** mat, int** sol_mat, int m, int n, string fname){

  //string fname = "visualize.kakuro";
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

void read_matrix(int** &matrix, ifstream &afile, int m, int n){

  matrix = new int*[m]; //rows

  // Allocates memory for m arrays with size n
  for(int i = 0; i < m; i++){
    matrix[i] = new int[n]; //cols
  }

  // Assigns values from the kakuro files
  int val;
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      afile >> val;
      matrix[i][j] = val;
    }
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
    /*
    #ifdef DEBUG
    cout << "############################" << endl;
    cout << "Creating sum with: " << endl;
    print_coords(start, end);
    cout << "Hint: " << hint << endl;
    cout << "Direction: " << dir << endl;
    cout << "Length: " << length << endl;
    cout << "############################" << endl;
    #endif
    */
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

  // Initialize a vector for summations
  vector<sum> sums;
  
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){

      // value on the i. row and j. column
      int val = matrix[i][j];

      // if the value is not -1 or -2, then it is a hint
      if(val != -1 && val != -2){
      
        int hint = val;
        hint = hint / 10;

        if((hint%100) == 0){

          hint = (int)(hint/100);
          // The coordinate of the first cell in the same row 
          COORD START = COORD(i, j+1); 
          // Returns the coordinate of the last cell(on the right direction) in the same row  as the starting cell that is not equal to -2.
          COORD END = find_end(matrix, m, n, i, j, d_right);

          // append the summations to the sums vector
          sum _sum = sum(START, END, hint, d_right);
          sums.push_back(_sum);
        }

        else{
          int div = (int)(hint/100);
          int rem = (int)(hint%100);
        
          if(div == 0 && rem != 0){
            // The coordinate of the first cell in the same column
            COORD START = COORD(i+1,j);
            // Returns the coordinate of the last cell(on the down direction) in the same column as the starting cell that is not equal to -2.
            COORD END = find_end(matrix, m, n, i, j, d_down);

            // append the summations to the sums vector
            sum _sum = sum(START, END, rem, d_down);
            sums.push_back(_sum);
          }

          if(div != 0 && rem != 0){
            // The coordinate of the first cell in the same column
            COORD START1 = COORD(i+1,j);
            // The coordinate of the first cell in the same row
            COORD START2 = COORD(i,j+1);
            // Returns the coordinate of the last cell(on the down direction) in the same column as the starting cell that is not equal to -2.
            COORD END1 = find_end(matrix, m, n, i, j, d_down);
            // Returns the coordinate of the last cell(on the right direction) in the same row as the starting cell that is not equal to -2.
            COORD END2 = find_end(matrix, m, n, i, j, d_right);

            sum _sum1 = sum(START1, END1, rem, d_down);
            sum _sum2 = sum(START2, END2, div, d_right);
            // append the summations to the sums vector
            sums.push_back(_sum1);
            sums.push_back(_sum2);
          }
        }
      }
      
    }
  }
  return sums;
}

bool check_solution(int** sol_mat, vector<sum> sums, int m, int n){
  /*
  Confirms the solution to see if it is correct or not
  */

  for(int i=0; i<=m; i++){
    
    if(sums[i].dir == 0){ // down direction
      int start = sums[i].start.first;
      int end = sums[i].end.first;
      int column = sums[i].start.second;

      int sum = 0;
      for(int j=start; j<end; j++){
        sum += sol_mat[j][column];
      }

      if(sum != sums[i].hint){
        return false;
      }
    }
    else{ // right direction
      int start = sums[i].start.second;
      int end = sums[i].end.second;
      int row = sums[i].start.first;

      int sum = 0;

      for(int j=start; j<end; j++){
        sum += sol_mat[row][j];
      }

      if(sum != sums[i].hint){
        return false;
      }
      }
    }

    return true;

  }




void solution(int** mat, int** sol_mat, vector<sum> sums, int m, int n){

  //TO DO: Write the solution
  //You can use any algorithm and data type
  //Write your solution to file in main function using sol_to_mat() after solving it

  int rand = 0;

  while(true){
   
    for(int i=0; i<m; i++){

      if(check_solution(sol_mat, sums, m, n)){
        return;
      }
      else{
        break;
      }

    }
  
  }

}

int main(int argc, char** argv){
  
  // Open the file kakuro 
  string filename(argv[1]);
  ifstream file;
  file.open(filename.c_str()); // c_str is used to convert filename variable to C-style string

  // The dimensions of the matrix
  int m, n;
  file >> m;
  file >> n;

  // 2D integer array 
  int** mat;

  // Creates m arrays with size n and assign the values of the matrix from kakuro file
  read_matrix(mat, file, m, n);
  // Prints the puzzle matrix without solution
  print_one_matrix(mat, m, n);
  
  // 2D integer array
  int** sol_mat;

  // Creates m arrays with size n and assign the hints and empty cells 
  convert_sol(mat, sol_mat, m, n);
  // Prints the solution matrix without solution
  print_one_matrix(sol_mat, m, n);
  
  // Get the summation of each row and column
  vector<sum> sums = get_sums(mat, m, n);


  // !!!!!!! SOLUTION !!!!!!!!!
  solution(mat, sol_mat, sums, m, n);


  // Prints the solution matrix with solution
  print_one_matrix(sol_mat, m, n);
  // Write solution matrix to a kakuro file
  sol_to_file(mat, sol_mat, m, n, "solution.kakuro");
  
  // Delete allocated values to save memory
  for (int i = 0; i < n; i++){
    delete mat[i];
    delete sol_mat[i];
  }

  delete mat;
  delete sol_mat;
  
  return 0;
}
