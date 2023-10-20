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

#ifndef _INTERPRETER_HH
#define _INTERPRETER_HH

#include <vector>
#include <string>
#include <cstddef>
#include <utility>
#include <map>
#include <tuple>
#include <stack>
#include <fstream>

#include "dynumfpack.h"
#include "dynmex.h"

#include "ErrorHandling.hh"
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

constexpr int NO_ERROR_ON_EXIT {0}, ERROR_ON_EXIT {1};

class Interpreter
{
private:
  double g0, gp0, glambda2;
  int try_at_iteration;

  void *Symbolic {nullptr}, *Numeric {nullptr};

  const BasicSymbolTable &symbol_table;
  const bool steady_state; // Whether this is a static or dynamic model

  // Whether to use the block-decomposed version of the bytecode file
  bool block_decomposed;

  Evaluate &evaluator;

  fstream SaveCode;

  Mem_Mngr mem_mngr;
  vector<int> u_liste;
  int *NbNZRow, *NbNZCol;
  NonZeroElem **FNZE_R, **FNZE_C;
  int u_count_init;

  int *pivot, *pivotk, *pivot_save;
  double *pivotv, *pivotva;
  int *b;
  bool *line_done;
  bool symbolic, alt_symbolic;
  int alt_symbolic_count;
  double markowitz_c_s;
  double res1a;
  long int nop1;
  map<tuple<int, int, int>, int> IM_i;
  int u_count_alloc, u_count_alloc_save;
  double slowc, slowc_save, prev_slowc_save, markowitz_c;
  int *index_equa; // Actually unused
  int u_count, tbreak_g;
  int iter;
  int start_compare;
  int restart;
  double lu_inc_tol;

  SuiteSparse_long *Ap_save, *Ai_save;
  double *Ax_save, *b_save;

  int stack_solve_algo, solve_algo;
  int minimal_solving_periods;
  int Per_u_, Per_y_;
  int maxit_;
  double *direction;
  double solve_tolf;
  // 1-norm error, square of 2-norm error, ∞-norm error
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

  int verbosity; // Corresponds to options_.verbosity

  vector<int> previous_block_exogenous;
  bool global_temporary_terms;
  bool print; // Whether the “print” command is requested
  int col_x, col_y;
  vector<double> residual;

  void evaluate_over_periods(bool forward);
  void solve_simple_one_periods();
  void solve_simple_over_periods(bool forward);
  void compute_complete_2b();
  void initializeTemporaryTerms();
  void evaluate_a_block(bool initialization, bool single_block, const string &bin_base_name);
  int simulate_a_block(const vector_table_conditional_local_type &vector_table_conditional_local, bool single_block, const string &bin_base_name);
  static string elastic(string str, unsigned int len, bool left);
  void check_for_controlled_exo_validity(const vector<s_plan> &sconstrained_extended_path);
  pair<bool, vector<int>> MainLoop(const string &bin_basename, bool evaluate, int block, bool constrained, const vector<s_plan> &sconstrained_extended_path, const vector_table_conditional_local_type &vector_table_conditional_local);
  void Simulate_Newton_Two_Boundaries(bool cvg, const vector_table_conditional_local_type &vector_table_conditional_local);
  void Simulate_Newton_One_Boundary(bool forward);
  void fixe_u();
  void Read_SparseMatrix(const string &file_name, int periods, int y_kmin, int y_kmax, bool two_boundaries);
  void Singular_display();
  void End_Solver();
  static int find_exo_num(const vector<s_plan> &sconstrained_extended_path, int value);
  static int find_int_date(const vector<pair<int, double>> &per_value, int value);
  void Init_GE();
  void Init_Matlab_Sparse(mxArray *A_m, mxArray *b_m, const mxArray *x0_m) const;
  tuple<SuiteSparse_long *, SuiteSparse_long *, double *, double *> Init_UMFPACK_Sparse(const mxArray *x0_m, const vector_table_conditional_local_type &vector_table_conditional_local) const;
  bool Init_Matlab_Sparse_Simple(const mxArray *A_m, const mxArray *b_m, const mxArray *x0_m) const;
  tuple<bool, SuiteSparse_long *, SuiteSparse_long *, double *, double *> Init_UMFPACK_Sparse_Simple(const mxArray *x0_m) const;
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
  pair<int, NonZeroElem *> At_Row(int r) const;
  NonZeroElem *At_Pos(int r, int c) const;
  pair<int, NonZeroElem *> At_Col(int c) const;
  pair<int, NonZeroElem *> At_Col(int c, int lag) const;
  int NRow(int r) const;
  int NCol(int c) const;
  int Get_u();
  void Delete_u(int pos);
  void Clear_u();

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

  void compute_block_time(int Per_u_, bool evaluate, bool no_derivatives);
  bool compute_complete(bool no_derivatives);

  pair<bool, double> compute_complete(double lambda);

public:
  Interpreter(Evaluate &evaluator_arg, double *params_arg, double *y_arg, double *ya_arg, double *x_arg, double *steady_y_arg,
              double *direction_arg, int y_size_arg,
              int nb_row_x_arg, int periods_arg, int y_kmin_arg, int y_kmax_arg,
              int maxit_arg_, double solve_tolf_arg, double markowitz_c_arg,
              int minimal_solving_periods_arg, int stack_solve_algo_arg, int solve_algo_arg,
              bool global_temporary_terms_arg, bool print_arg, mxArray *GlobalTemporaryTerms_arg,
              bool steady_state_arg, bool block_decomposed_arg, int col_x_arg, int col_y_arg, const BasicSymbolTable &symbol_table_arg, int verbosity_arg);
  pair<bool, vector<int>> extended_path(const string &file_name, bool evaluate, int block, int nb_periods, const vector<s_plan> &sextended_path, const vector<s_plan> &sconstrained_extended_path, const vector<string> &dates, const table_conditional_global_type &table_conditional_global);
  pair<bool, vector<int>> compute_blocks(const string &file_name, bool evaluate, int block);
  void Close_SaveCode();

  inline mxArray *
  get_jacob(int block_num) const
  {
    return jacobian_block[block_num];
  }
  inline mxArray *
  get_jacob_exo(int block_num) const
  {
    return jacobian_exo_block[block_num];
  }
  inline mxArray *
  get_jacob_exo_det(int block_num) const
  {
    return jacobian_det_exo_block[block_num];
  }
  inline vector<double>
  get_residual() const
  {
    return residual;
  }
  inline mxArray *
  get_Temporary_Terms() const
  {
    return GlobalTemporaryTerms;
  }
};

#endif // _INTERPRETER_HH
