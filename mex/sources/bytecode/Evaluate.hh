/*
 * Copyright © 2007-2021 Dynare Team
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

#ifndef EVALUATE_HH_INCLUDED
#define EVALUATE_HH_INCLUDED

#include <vector>
#include <string>

#include "dynmex.h"

#define BYTE_CODE
#include "CodeInterpreter.hh"
#include "ErrorHandling.hh"

class Evaluate : public ErrorMsg
{
private:
  unsigned int EQN_dvar1, EQN_dvar2, EQN_dvar3;
  int EQN_lag1, EQN_lag2, EQN_lag3;
protected:
  mxArray *GlobalTemporaryTerms;
  it_code_type start_code, end_code;
  double pow1(double a, double b);
  double divide(double a, double b);
  double log1(double a);
  double log10_1(double a);
  void evaluate_over_periods(bool forward);
  void solve_simple_one_periods();
  void solve_simple_over_periods(bool forward);
  void compute_block_time(int Per_u_, bool evaluate, bool no_derivatives);
  code_liste_type code_liste;
  it_code_type it_code;
  int Block_Count, Per_u_, Per_y_;
  int it_;
  int maxit_, size_of_direction;
  double *direction;
  double solve_tolf;
  bool GaussSeidel;
  int equation, derivative_equation, derivative_variable;
  string filename;
  int stack_solve_algo, solve_algo;
  bool global_temporary_terms;
  bool print, print_error;
  double res1, res2, max_res;
  int max_res_idx;
  vector<Block_contain_type> Block_Contain;

  int size;
  int *index_vara;

  bool print_it, forward;
  int minimal_solving_periods;
  int type, block_num, symbol_table_endo_nbr, Block_List_Max_Lag, Block_List_Max_Lead, u_count_int, block;
  string file_name, bin_base_name;
  bool Gaussian_Elimination, is_linear;
public:
  bool steady_state;
  double slowc;
  Evaluate();
  Evaluate(int y_size_arg, int y_kmin_arg, int y_kmax_arg, bool print_it_arg, bool steady_state_arg, int periods_arg, int minimal_solving_periods_arg, double slowc);
  void set_block(int size_arg, int type_arg, string file_name_arg, string bin_base_name_arg, int block_num_arg,
                 bool is_linear_arg, int symbol_table_endo_nbr_arg, int Block_List_Max_Lag_arg, int Block_List_Max_Lead_arg, int u_count_int_arg, int block_arg);
  void evaluate_complete(bool no_derivatives);
  bool compute_complete(bool no_derivatives, double &res1, double &res2, double &max_res, int &max_res_idx);
  void compute_complete_2b(bool no_derivatives, double *_res1, double *_res2, double *_max_res, int *_max_res_idx);

  bool compute_complete(double lambda, double *crit);
};

#endif
