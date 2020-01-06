// This file is a lower triangular sparse matrix solver program

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include "mmio.h"
#include "directedGraph.h"

// A data structure to more easily store a ccs format matrix
struct ccs_matrix {
  int dimension;
  std::vector<int> col_ptr;
  std::vector<int> row_index;
  std::vector<double> mat_vals;
};

// A helper function to determine double equality
bool double_equals(double a, double b, double epsilon = 0.00001)
{
    return abs(a - b) < epsilon;
}

/* A helper function to read in the sparse matrix to compute
 * Adapted from https://math.nist.gov/MatrixMarket/mmio/c/example_read.c
 * @parameter: A pointer to the matrix file to read in
 * @return: A ccs_matrix struct containing the matrix information
 * Note: This function assumes the matrix is a square matrix */
struct ccs_matrix read_matrix(FILE *f)
{
  int ret_code;
  MM_typecode matcode;
  int M, N, nz;

  if (mm_read_banner(f, &matcode) != 0)
  {
      printf("Could not process Matrix Market banner.\n");
      exit(1);
  }

  /*  This is how one can screen matrix types if their application */
  /*  only supports a subset of the Matrix Market data types.      */

  if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
          mm_is_sparse(matcode) )
  {
      printf("Sorry, this application does not support ");
      printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
      exit(1);
  }

  /* find out size of sparse matrix .... */
  if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
      exit(1);

  int cur_index = -1;
  int prev_col_ptr = -1;
  int cur_col_ptr = -1;
  vector<double> val(nz, 0);
  vector<int> row_index;
  vector<int> col_ptr;

  // Reading in values from the matrix
  for (int i = 0; i < nz; i++) {
    fscanf(f, "%d %d %lf\n", &cur_index, &cur_col_ptr, &val[i]);
    cur_index--;
    cur_col_ptr--;

    row_index.push_back(cur_index);
    if (cur_col_ptr != prev_col_ptr) {
      col_ptr.push_back(i);
      prev_col_ptr = cur_col_ptr;
    }
  }

  col_ptr.push_back(nz);

  struct ccs_matrix final_matrix;
  final_matrix.dimension = M;
  final_matrix.col_ptr = col_ptr;
  final_matrix.row_index = row_index;
  final_matrix.mat_vals = val;

  return final_matrix;

}

/* A helper function to read in the target vector value
 * Adapted from https://math.nist.gov/MatrixMarket/mmio/c/example_read.c
 * @parameter: A pointer to the matrix file to read in
 * @return: The target vector value stored in a vector */
vector<double> read_res(FILE *f) {
  int ret_code;
  MM_typecode matcode;
  int M, N, nz;

  if (mm_read_banner(f, &matcode) != 0)
  {
      printf("Could not process Matrix Market banner.\n");
      exit(1);
  }

  /*  This is how one can screen matrix types if their application */
  /*  only supports a subset of the Matrix Market data types.      */

  if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
          mm_is_sparse(matcode) )
  {
      printf("Sorry, this application does not support ");
      printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
      exit(1);
  }

  if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
      exit(1);

  vector<double> val(M, 0);

  int cur_row;
  double cur_val;

  for (int i = 0; i < nz; i++) {
    fscanf(f, "%d %*d %lf\n", &cur_row, &cur_val);
    cur_row--;
    val[cur_row] = cur_val;
  }

  return val;
}

/* A helper function to find the potential non-zero entries in the answer vector
 * The algorithm is based on http://faculty.cse.tamu.edu/davis/publications_files/survey_tech_report.pdf
 * @parameter: n - The matrix dimension (This function assumes square matrix)
 *             Lp - The column pointer vector
 *             Li - The index vector
 *             Lx - The value vector
 *             x - The rhs vector value
 * @return: a dictionary containing topological order of the dependencies graph,
 *          ready for wavefront parallelization */
unordered_map<int, vector<int>> find_nonzero_entries(int n, int *Lp, int *Li,
  double *Lx, double *x) {

  // Creating a DAG representing the entry dependencies
  Graph g;
  for (int i = 0; i < n - 1; i++) {
    for (int j = Lp[i]+1; j < Lp[i + 1]; j++) {
      g.addEdge(i, Li[j]);
      // DEBUG: cout << "Adding edge: " << i << "->" << Li[j] << "\n";
    }
  }

  // The initial non-zeros entries (from the rhs vector, b)
  vector<int> nonzeros;
  for (int i = 0; i < n; i++) {
    if (x[i] != 0) {
      nonzeros.push_back(i);
      // DEBUG: cout << "nonzero: " << i << "\n";
    }
  }

  // Perform a topological sort on the graph to prepare for parallelization
  vector<stack<int>> res = g.topologicalSearch(nonzeros);

  return wavefront(res);
}

// The following function is initially provided by the Paramathic lab
/*
* Lower triangular solver Lx=b
* L is stored in the compressed column storage format
* Inputs are:
* n : the matrix dimension
* Lp : the column pointer of L
* Li : the row index of L
* Lx : the values of L
* In/Out:
* x : the right hand-side b at start and the solution x at the end.
*/
int lsolve (int n, int *Lp, int *Li, double *Lx, double *x) {
  int p;

  if (!Lp || !Li || !x) {
    return (0);
  }

  // Parallel version of the triangular solve, based on wavefront parallelization

  // Step 1: Get the potential non-zero entries in the answer vector, stored in
  // a dictionary ready for wavefront parallelization
  unordered_map<int, vector<int>> wavefront =
    find_nonzero_entries(n, Lp, Li, Lx, x);

  // Step 2: Within each wave there is no co-dependencies, therefore we can safely
  // parallelize the code
  for (int i = 0; i < (int) wavefront.size(); i++) {
    vector<int> cur_group = wavefront[i];
    #pragma omp parallel for
    for (int j = 0; j < (int) cur_group.size(); j++) {
      int cur_val = cur_group[j];
      x [cur_val] /= Lx [Lp [cur_val]];
      for (p = Lp [cur_val]+1 ; p < Lp [cur_val+1] ; p++) {
        x [Li [p]] -= Lx [p] * x [cur_val];
      }
    }
  }

  // The non-parallel version of the triangular solve

  // for (int j = 0 ; j < n ; j++) {
  //   x [j] /= Lx [Lp [j]] ;
  //   for (p = Lp [j]+1 ; p < Lp [j+1] ; p++) {
  //     x [Li [p]] -= Lx [p] * x [j];
  //   }
  // }

  return (1) ;
}

// The following function is provided by the Paramathic lab
/*
* Sparse matrix-vector multiply: y = A*x
* A is stored in the compressed column storage format
* Inputs:
* Ap : the column pointer of A
* Ai : the row index of A
* Ax : the values of A
* x : is a dense vector
* Output:
* y : is a dense vector that stores the result of multiplication
*/
int spmv_ccs (int n, int *Ap, int *Ai, double *Ax, double *x, double *y) {
  int p, j;
  if (!Ap || !x || !y) { return (0); }
  for (j = 0 ; j < n ; j++) {
    for (p = Ap [j] ; p < Ap [j+1] ; p++) {
      y [Ai [p]] += Ax [p] * x [j] ;
    }
  }
  return (1);
}

int main(int argc, char *argv[]) {

  FILE *L_file, *b_file;

  if (argc < 3) {
    fprintf(stderr, "Usage: %s [triangular_matrix] [right_hand_side]\n", argv[0]);
    exit(1);
  }
  else {
    if ((L_file = fopen(argv[1], "r")) == NULL) {
      exit(1);
    }
    if ((b_file = fopen(argv[2], "r")) == NULL) {
      exit(1);
    }
  }

  // DEBUG: Test matrix

  // Test case 1:
  // int dim = 4;
  // vector<int> col_ptr = {0, 3, 5, 7, 8};
  // vector<int> row_index = {0, 1, 3, 1, 2, 2, 3, 3};
  // vector<double> mat_vals = {5, 3, 1, 4, 5, 1, 3, 9};
  // struct ccs_matrix L;
  // L.dimension = dim;
  // L.col_ptr = col_ptr;
  // L.row_index = row_index;
  // L.mat_vals = mat_vals;
  // vector<double> b = {35, 53, 42, 67};
  // vector<double> og_b(b.begin(), b.end());
  // vector<double> test_val(b.size(), 0);

  // Test case 2:
  // int dim = 14;
  // vector<int> col_ptr = {0, 2, 6, 9, 11, 13, 16, 19, 21, 24, 28, 31, 33, 35, 37};
  // vector<int> row_index = {0, 2, 1, 3, 6, 8, 2, 4, 7, 3, 8, 4, 7, 5 ,8, 9, 6,
  // 9, 10, 7, 9, 8, 11, 12, 9, 10, 12, 13, 10, 11, 12, 11, 12, 12, 13, 13};
  // vector<double> mat_vals = {1, 3, 1, 4, 2, 7, 1, 5, 6, 1, 9, 1, 8, 1, 10,
  // 11, 1, 14, 2, 1, 3, 1, 8, 5, 1, 2, 4, 1, 1, 8, 3, 1, 7, 1, 4, 1};
  // struct ccs_matrix L;
  // L.dimension = dim;
  // L.col_ptr = col_ptr;
  // L.row_index = row_index;
  // L.mat_vals = mat_vals;
  // vector<double> b = {0, 0, 0, 4, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0};
  // vector<double> og_b(b.begin(), b.end());
  // vector<double> test_val(b.size(), 0);

  // DEBUG: End of test matrix

  // Setting up the computation
  struct ccs_matrix L = read_matrix(L_file);
  vector<double> b = read_res(b_file);
  vector<double> og_b(b.begin(), b.end());
  vector<double> test_val(b.size(), 0);

  // Run the triangular solve
  lsolve(L.dimension, &L.col_ptr[0], &L.row_index[0], &L.mat_vals[0], &b[0]);

  // Multiply the resulting vector with the original matrix to test for correctness
  spmv_ccs(L.dimension, &L.col_ptr[0], &L.row_index[0], &L.mat_vals[0], &b[0], &test_val[0]);

  // Compare the resulting value with the original rhs vector; if there's an
  // errorenous entry, SCREAM
  for (long unsigned int i = 0; i < b.size(); i++) {
    if (!double_equals(test_val[i], og_b[i])) {
        cout << "Calculation error at the " << i << "th entry: ";
        cout << og_b[i] << " vs " << test_val[i] << '\n';
        return (0);
    }
  }

  return (1);

}
