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

#ifndef _EVALUATE_HH
#define _EVALUATE_HH

#include <vector>
#include <string>
#include <map>
#include <optional>

#define BYTECODE_MEX
#include "Bytecode.hh"
#include "ErrorHandling.hh"
#include "BasicSymbolTable.hh"

using it_code_type = instructions_list_t::const_iterator;

class Evaluate
{
private:
  ExpressionType EQN_type;
  int EQN_equation, EQN_block, EQN_dvar1;
  int EQN_lag1, EQN_lag2, EQN_lag3;
  map<int, double> TEF;
  map<pair<int, int>, double> TEFD;
  map<tuple<int, int, int>, double> TEFDD;

  string error_location(it_code_type expr_begin, it_code_type faulty_op, bool steady_state, int it_) const;

protected:
  BasicSymbolTable &symbol_table;
  int EQN_block_number;
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
  it_code_type start_code, end_code;
  double pow1(double a, double b);
  double divide(double a, double b);
  double log1(double a);
  double log10_1(double a);
  void evaluate_over_periods(bool forward);
  void solve_simple_one_periods();
  void solve_simple_over_periods(bool forward);
  void compute_block_time(int Per_u_, bool evaluate, bool no_derivatives);
  instructions_list_t code_liste;
  it_code_type it_code;
  int Per_u_, Per_y_;
  int it_;
  int maxit_;
  double *direction;
  double solve_tolf;
  bool print_error;
  double res1, res2, max_res;
  int max_res_idx;
  vector<Block_contain_type> Block_Contain;

  int size;
  int *index_vara;

  BlockSimulationType type;
  int block_num, symbol_table_endo_nbr, u_count_int, block;
  string file_name, bin_base_name;
  bool is_linear;

  bool steady_state;

  /* Prints a bytecode expression in human readable form.
     If faulty_op is not default constructed, it should point to a tag withing
     the expression that created a floating point exception, in which case the
     corresponding mathematical operator will be printed within braces.
     The second output argument points to the tag past the expression. */
  pair<string, it_code_type> print_expression(const it_code_type &expr_begin, const optional<it_code_type> &faulty_op = nullopt) const;

public:
  Evaluate(int y_size_arg, int y_kmin_arg, int y_kmax_arg, bool steady_state_arg, int periods_arg, BasicSymbolTable &symbol_table_arg);
  void set_block(int size_arg, BlockSimulationType type_arg, string file_name_arg, string bin_base_name_arg, int block_num_arg,
                 bool is_linear_arg, int symbol_table_endo_nbr_arg, int u_count_int_arg, int block_arg);
  void evaluate_complete(bool no_derivatives);
  bool compute_complete(bool no_derivatives, double &res1, double &res2, double &max_res, int &max_res_idx);
  void compute_complete_2b(bool no_derivatives, double *_res1, double *_res2, double *_max_res, int *_max_res_idx);

  bool compute_complete(double lambda, double *crit);
};

#endif // _EVALUATE_HH
