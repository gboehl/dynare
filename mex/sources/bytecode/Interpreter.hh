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

#include "dynmex.h"

#include "ErrorHandling.hh"
#include "SparseMatrix.hh"
#define BYTECODE_MEX
#include "Bytecode.hh"

using namespace std;

constexpr int NO_ERROR_ON_EXIT {0}, ERROR_ON_EXIT {1};

class Interpreter : public dynSparseMatrix
{
private:
  vector<int> previous_block_exogenous;
protected:
  void evaluate_a_block(bool initialization);
  int simulate_a_block(const vector_table_conditional_local_type &vector_table_conditional_local);
  void print_a_block();
  string elastic(string str, unsigned int len, bool left);
public:
  Interpreter(double *params_arg, double *y_arg, double *ya_arg, double *x_arg, double *steady_y_arg, double *steady_x_arg,
              double *direction_arg, size_t y_size_arg,
              size_t nb_row_x_arg, int periods_arg, int y_kmin_arg, int y_kmax_arg,
              int maxit_arg_, double solve_tolf_arg, size_t size_of_direction_arg, int y_decal_arg, double markowitz_c_arg,
              string &filename_arg, int minimal_solving_periods_arg, int stack_solve_algo_arg, int solve_algo_arg,
              bool global_temporary_terms_arg, bool print_arg, bool print_error_arg, mxArray *GlobalTemporaryTerms_arg,
              bool steady_state_arg, bool block_decomposed_arg, bool print_it_arg, int col_x_arg, int col_y_arg, BasicSymbolTable &symbol_table_arg);
  bool extended_path(const string &file_name, const string &bin_basename, bool evaluate, int block, int &nb_blocks, int nb_periods, const vector<s_plan> &sextended_path, const vector<s_plan> &sconstrained_extended_path, const vector<string> &dates, const table_conditional_global_type &table_conditional_global);
  bool compute_blocks(const string &file_name, const string &bin_basename, bool evaluate, int block, int &nb_blocks);
  void check_for_controlled_exo_validity(FBEGINBLOCK_ *fb, const vector<s_plan> &sconstrained_extended_path);
  bool MainLoop(const string &bin_basename, const CodeLoad &code, bool evaluate, int block, bool last_call, bool constrained, const vector<s_plan> &sconstrained_extended_path, const vector_table_conditional_local_type &vector_table_conditional_local);
  void ReadCodeFile(const string &file_name, CodeLoad &code);

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
