/*
 * Copyright © 2007-2023 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef _SPARSEMATRIX_HH
#define _SPARSEMATRIX_HH

#include <utility>
#include <vector>
#include <map>
#include <tuple>
#include <stack>
#include <fstream>
#include <string>
#include <ctime>

#include "dynumfpack.h"
#include "dynmex.h"

#include "Mem_Mngr.hh"
#include "Evaluate.hh"

using namespace std;

struct t_save_op_s
{
  short int lag, operat;
  int first, second;
};

struct s_plan
{
  string var, exo;
  int var_num, exo_num;
  vector<pair<int, double>> per_value;
  vector<double> value;
};

struct table_conditional_local_type
{
  bool is_cond;
  int var_exo, var_endo;
  double constrained_value;
};
using vector_table_conditional_local_type = vector<table_conditional_local_type>;
using table_conditional_global_type = map<int, vector_table_conditional_local_type>;

constexpr int IFLD = 0, IFDIV = 1, IFLESS = 2, IFSUB = 3, IFLDZ = 4, IFMUL = 5, IFSTP = 6, IFADD = 7;
constexpr double eps = 1e-15, very_big = 1e24;
constexpr int alt_symbolic_count_max = 1;
constexpr double mem_increasing_factor = 1.1;

class dynSparseMatrix
{
public:
  dynSparseMatrix(Evaluate &evaluator_arg, int y_size_arg, int y_kmin_arg, int y_kmax_arg, bool print_it_arg, bool steady_state_arg, bool block_decomposed_arg, int periods_arg, int minimal_solving_periods_arg, const BasicSymbolTable &symbol_table_arg, bool print_error_arg);
  void Simulate_Newton_Two_Boundaries(int blck, int y_size, int y_kmin, int y_kmax, int Size, int periods, bool cvg, int minimal_solving_periods, int stack_solve_algo, const vector_table_conditional_local_type &vector_table_conditional_local);
  void Simulate_Newton_One_Boundary(bool forward);
  void fixe_u(double **u, int u_count_int, int max_lag_plus_max_lead_plus_1);
  void Read_SparseMatrix(const string &file_name, int Size, int periods, int y_kmin, int y_kmax, bool two_boundaries, int stack_solve_algo, int solve_algo);
  void Close_SaveCode();
  void Singular_display(int block, int Size);
  void End_Solver();
  double g0, gp0, glambda2;
  int try_at_iteration;
  static int find_exo_num(const vector<s_plan> &sconstrained_extended_path, int value);
  static int find_int_date(const vector<pair<int, double>> &per_value, int value);

private:
  void Init_GE(int periods, int y_kmin, int y_kmax, int Size, const map<tuple<int, int, int>, int> &IM);
  void Init_Matlab_Sparse(int periods, int y_kmin, int y_kmax, int Size, const map<tuple<int, int, int>, int> &IM, mxArray *A_m, mxArray *b_m, const mxArray *x0_m) const;
  void Init_UMFPACK_Sparse(int periods, int y_kmin, int y_kmax, int Size, const map<tuple<int, int, int>, int> &IM, SuiteSparse_long **Ap, SuiteSparse_long **Ai, double **Ax, double **b, const mxArray *x0_m, const vector_table_conditional_local_type &vector_table_conditional_local, int block_num) const;
  void Init_Matlab_Sparse_Simple(int Size, const map<tuple<int, int, int>, int> &IM, const mxArray *A_m, const mxArray *b_m, bool &zero_solution, const mxArray *x0_m) const;
  void Init_UMFPACK_Sparse_Simple(int Size, const map<tuple<int, int, int>, int> &IM, SuiteSparse_long **Ap, SuiteSparse_long **Ai, double **Ax, double **b, bool &zero_solution, const mxArray *x0_m) const;
  void Simple_Init(int Size, const map<tuple<int, int, int>, int> &IM, bool &zero_solution);
  void End_GE();
  bool mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc);
  bool golden(double ax, double bx, double cx, double tol, double solve_tolf, double *xmin);
  void Solve_ByteCode_Symbolic_Sparse_GaussianElimination(int Size, bool symbolic, int Block_number);
  bool Solve_ByteCode_Sparse_GaussianElimination(int Size, int blck, int it_);
  void Solve_Matlab_Relaxation(mxArray *A_m, mxArray *b_m, unsigned int Size, double slowc_l);
  static void Print_UMFPack(const SuiteSparse_long *Ap, const SuiteSparse_long *Ai, const double *Ax, int n);
  static void Printfull_UMFPack(const SuiteSparse_long *Ap, const SuiteSparse_long *Ai, const double *Ax, const double *b, int n);
  static void PrintM(int n, const double *Ax, const mwIndex *Ap, const mwIndex *Ai);
  void Solve_LU_UMFPack(SuiteSparse_long *Ap, SuiteSparse_long *Ai, double *Ax, double *b, int n, int Size, double slowc_l, bool is_two_boundaries, int it_, const vector_table_conditional_local_type &vector_table_conditional_local);
  void Solve_LU_UMFPack(SuiteSparse_long *Ap, SuiteSparse_long *Ai, double *Ax, double *b, int n, int Size, double slowc_l, bool is_two_boundaries, int it_);

  void End_Matlab_LU_UMFPack();
  void Solve_Matlab_GMRES(mxArray *A_m, mxArray *b_m, int Size, double slowc, int block, bool is_two_boundaries, int it_, mxArray *x0_m);
  void Solve_Matlab_BiCGStab(mxArray *A_m, mxArray *b_m, int Size, double slowc, int block, bool is_two_boundaries, int it_, mxArray *x0_m, int precond);
  void Check_and_Correct_Previous_Iteration(int y_size, int size);
  bool Simulate_One_Boundary(int blck, int y_size, int size);
  bool solve_linear(int block_num, int y_size, int size, int iter);
  void solve_non_linear(int block_num, int y_size, int size);
  string preconditioner_print_out(string s, int preconditioner, bool ss);
  bool compare(int *save_op, int *save_opa, int *save_opaa, int beg_t, int periods, long nop4, int Size);
  void Insert(int r, int c, int u_index, int lag_index);
  void Delete(int r, int c);
  int At_Row(int r, NonZeroElem **first) const;
  int At_Pos(int r, int c, NonZeroElem **first) const;
  int At_Col(int c, NonZeroElem **first) const;
  int At_Col(int c, int lag, NonZeroElem **first) const;
  int NRow(int r) const;
  int NCol(int c) const;
  int Union_Row(int row1, int row2) const;
  void Print(int Size, const int *b) const;
  int Get_u();
  void Delete_u(int pos);
  void Clear_u();
  void Print_u() const;
  void *Symbolic {nullptr}, *Numeric {nullptr};
  void Check_the_Solution(int periods, int y_kmin, int y_kmax, int Size, double *u, int *pivot, int *b);
  int complete(int beg_t, int Size, int periods, int *b);
  void bksub(int tbreak, int last_period, int Size, double slowc_l);
  void simple_bksub(int it_, int Size, double slowc_l);
  // Computes Aᵀ where A is are sparse. The result is sparse.
  static mxArray *Sparse_transpose(const mxArray *A_m);
  // Computes Aᵀ·B where A and B are sparse. The result is sparse.
  static mxArray *Sparse_mult_SAT_SB(const mxArray *A_m, const mxArray *B_m);
  // Computes Aᵀ·B where A is sparse and B is dense. The result is sparse.
  static mxArray *Sparse_mult_SAT_B(const mxArray *A_m, const mxArray *B_m);
  // Computes Aᵀ·B where A is sparse and B is dense. The result is dense.
  static mxArray *mult_SAT_B(const mxArray *A_m, const mxArray *B_m);
  // Computes A−B where A and B are sparse. The result is sparse.
  static mxArray *Sparse_subtract_SA_SB(const mxArray *A_m, const mxArray *B_m);
  // Computes A−B where A and B are dense. The result is dense.
  static mxArray *subtract_A_B(const mxArray *A_m, const mxArray *B_m);
protected:
  const BasicSymbolTable &symbol_table;
  const bool steady_state; // Whether this is a static or dynamic model

  // Whether to use the block-decomposed version of the bytecode file
  bool block_decomposed;

  Evaluate &evaluator;

  stack<double> Stack;
  int nb_prologue_table_u, nb_first_table_u, nb_middle_table_u, nb_last_table_u;
  int nb_prologue_table_y, nb_first_table_y, nb_middle_table_y, nb_last_table_y;
  int middle_count_loop;
  fstream SaveCode;
  string filename;
  int max_u, min_u;
  clock_t time00;

  Mem_Mngr mem_mngr;
  vector<int> u_liste;
  map<pair<int, int>, NonZeroElem *> Mapped_Array;
  int *NbNZRow, *NbNZCol;
  NonZeroElem **FNZE_R, **FNZE_C;
  int u_count_init;

  int *pivot, *pivotk, *pivot_save;
  double *pivotv, *pivotva;
  int *b;
  bool *line_done;
  bool symbolic, alt_symbolic;
  int alt_symbolic_count;
  int *g_save_op;
  int first_count_loop;
  int g_nop_all;
  double markowitz_c_s;
  double res1a;
  long int nop_all, nop1, nop2;
  map<tuple<int, int, int>, int> IM_i;
  vector<double> residual;
  int u_count_alloc, u_count_alloc_save;
  vector<double *> jac;
  double *jcb;
  double slowc, slowc_save, prev_slowc_save, markowitz_c;
  int y_decal;
  int *index_equa;
  int u_count, tbreak_g;
  int iter;
  int start_compare;
  int restart;
  double g_lambda1, g_lambda2, gp_0;
  double lu_inc_tol;

  SuiteSparse_long *Ap_save, *Ai_save;
  double *Ax_save, *b_save;
  mxArray *A_m_save, *b_m_save;

  int stack_solve_algo, solve_algo;
  int minimal_solving_periods;
  bool print_it;
  int Per_u_, Per_y_;
  int maxit_;
  double *direction;
  double solve_tolf;
  double res1, res2, max_res;
  int max_res_idx;
  int *index_vara;

  double *y, *ya;
  int y_size;
  double *T;
  int nb_row_x;
  int y_kmin, y_kmax, periods;
  double *x, *params;
  double *u;
  double *steady_y;
  double *g1, *r, *res;
  vector<mxArray *> jacobian_block, jacobian_exo_block, jacobian_det_exo_block;
  mxArray *GlobalTemporaryTerms;
  int it_;
  map<int, double> TEF;
  map<pair<int, int>, double> TEFD;
  map<tuple<int, int, int>, double> TEFDD;

  // Information about the current block
  int block_num; // Index of the current block
  int size; // Size of the current block
  BlockSimulationType type;
  bool is_linear;
  int u_count_int;
  vector<Block_contain_type> Block_Contain;

  void compute_block_time(int Per_u_, bool evaluate, bool no_derivatives);
  bool compute_complete(bool no_derivatives, double &res1, double &res2, double &max_res, int &max_res_idx);

  bool compute_complete(double lambda, double *crit);

  bool print_error; // Whether to stop processing on floating point exceptions
};

#endif // _SPARSEMATRIX_HH
