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

#include <sstream>
#include <algorithm>
#include <cstring>
#include <filesystem>

#include "Interpreter.hh"

constexpr double BIG = 1.0e+8, SMALL = 1.0e-5;

Interpreter::Interpreter(double *params_arg, double *y_arg, double *ya_arg, double *x_arg, double *steady_y_arg,
                         double *direction_arg, size_t y_size_arg,
                         size_t nb_row_x_arg, int periods_arg, int y_kmin_arg, int y_kmax_arg,
                         int maxit_arg_, double solve_tolf_arg, size_t size_of_direction_arg, int y_decal_arg, double markowitz_c_arg,
                         string &filename_arg, int minimal_solving_periods_arg, int stack_solve_algo_arg, int solve_algo_arg,
                         bool global_temporary_terms_arg, bool print_arg, bool print_error_arg, mxArray *GlobalTemporaryTerms_arg,
                         bool steady_state_arg, bool block_decomposed_arg, bool print_it_arg, int col_x_arg, int col_y_arg, BasicSymbolTable &symbol_table_arg)
: dynSparseMatrix {y_size_arg, y_kmin_arg, y_kmax_arg, print_it_arg, steady_state_arg, block_decomposed_arg, periods_arg, minimal_solving_periods_arg, symbol_table_arg}
{
  params = params_arg;
  y = y_arg;
  ya = ya_arg;
  x = x_arg;
  steady_y = steady_y_arg;
  direction = direction_arg;
  nb_row_x = nb_row_x_arg;
  periods = periods_arg;
  maxit_ = maxit_arg_;
  solve_tolf = solve_tolf_arg;
  size_of_direction = size_of_direction_arg;
  slowc = 1;
  slowc_save = 1;
  y_decal = y_decal_arg;
  markowitz_c = markowitz_c_arg;
  filename = filename_arg;
  T = nullptr;
  minimal_solving_periods = minimal_solving_periods_arg;
  stack_solve_algo = stack_solve_algo_arg;
  solve_algo = solve_algo_arg;
  global_temporary_terms = global_temporary_terms_arg;
  print = print_arg;
  col_x = col_x_arg;
  col_y = col_y_arg;
  GlobalTemporaryTerms = GlobalTemporaryTerms_arg;
  print_error = print_error_arg;
  print_it = print_it_arg;
}

void
Interpreter::evaluate_a_block(bool initialization)
{
  it_code_type begining;

  switch (type)
    {
    case BlockSimulationType::evaluateForward:
      if (steady_state)
        {
          compute_block_time(0, true, false);
          if (block >= 0)
            for (int j = 0; j < size; j++)
              residual[j] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
          else
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
        }
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, false);
              if (block >= 0)
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+Block_Contain[j].Equation] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
            }
        }
      break;
    case BlockSimulationType::solveForwardSimple:
      g1 = static_cast<double *>(mxMalloc(size*size*sizeof(double)));
      test_mxMalloc(g1, __LINE__, __FILE__, __func__, size*size*sizeof(double));
      r = static_cast<double *>(mxMalloc(size*sizeof(double)));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, false);
          if (block < 0)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, false);
              if (block < 0)
                for (int j = 0; j < size; j++)
                  residual[Per_y_+Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = r[j];
            }
        }
      mxFree(g1);
      mxFree(r);
      break;
    case BlockSimulationType::solveForwardComplete:
      if (initialization)
        {
          fixe_u(&u, u_count_int, u_count_int);
          Read_SparseMatrix(bin_base_name, size, 1, 0, 0, false, stack_solve_algo, solve_algo);
        }
#ifdef DEBUG
      mexPrintf("in SOLVE FORWARD COMPLETE r = mxMalloc(%d*sizeof(double))\n", size);
#endif
      r = static_cast<double *>(mxMalloc(size*sizeof(double)));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, false);
          if (block < 0)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, false);
              if (block < 0)
                for (int j = 0; j < size; j++)
                  residual[it_*y_size+Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = r[j];
            }
        }
      mxFree(r);
      break;
    case BlockSimulationType::evaluateBackward:
      if (steady_state)
        {
          compute_block_time(0, true, false);
          if (block >= 0)
            for (int j = 0; j < size; j++)
              residual[j] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
          else
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
        }
      else
        {
          begining = it_code;
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, false);
              if (block >= 0)
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+Block_Contain[j].Equation] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
            }
        }
      break;
    case BlockSimulationType::solveBackwardSimple:
      g1 = static_cast<double *>(mxMalloc(size*size*sizeof(double)));
      test_mxMalloc(g1, __LINE__, __FILE__, __func__, size*size*sizeof(double));
      r = static_cast<double *>(mxMalloc(size*sizeof(double)));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, false);
          if (block < 0)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          begining = it_code;
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, false);
              if (block < 0)
                for (int j = 0; j < size; j++)
                  residual[Per_y_+Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = r[j];
            }
        }
      mxFree(g1);
      mxFree(r);
      break;
    case BlockSimulationType::solveBackwardComplete:
      if (initialization)
        {
          fixe_u(&u, u_count_int, u_count_int);
          Read_SparseMatrix(bin_base_name, size, 1, 0, 0, false, stack_solve_algo, solve_algo);
        }
      r = static_cast<double *>(mxMalloc(size*sizeof(double)));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, false);
          if (block < 0)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          begining = it_code;
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, false);
              if (block < 0)
                for (int j = 0; j < size; j++)
                  residual[Per_y_+Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = r[j];
            }
        }
      mxFree(r);
      break;
    case BlockSimulationType::solveTwoBoundariesSimple:
    case BlockSimulationType::solveTwoBoundariesComplete:
      if (initialization)
        {
          fixe_u(&u, u_count_int, u_count_int);
          Read_SparseMatrix(bin_base_name, size, periods, y_kmin, y_kmax, true, stack_solve_algo, solve_algo);
        }
      u_count = u_count_int*(periods+y_kmax+y_kmin);
      r = static_cast<double *>(mxMalloc(size*sizeof(double)));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size*sizeof(double));
      begining = it_code;
      for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
        {
          Per_u_ = (it_-y_kmin)*u_count_int;
          Per_y_ = it_*y_size;
          it_code = begining;
          compute_block_time(Per_u_, true, false);
          if (block < 0)
            for (int j = 0; j < size; j++)
              residual[it_*y_size+Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[it_*size+j] = r[j];
        }
      mxFree(r);
      break;
    case BlockSimulationType::unknown:
      throw FatalException{"UNKNOWN block simulation type: impossible case"};
    }
}

int
Interpreter::simulate_a_block(const vector_table_conditional_local_type &vector_table_conditional_local)
{
  it_code_type begining;
  max_res = 0;
  max_res_idx = 0;
  bool cvg;
  double *y_save;
#ifdef DEBUG
  mexPrintf("simulate_a_block type = %d, periods=%d, y_kmin=%d, y_kmax=%d\n", type, periods, y_kmin, y_kmax);
  mexEvalString("drawnow;");
#endif
  switch (type)
    {
    case BlockSimulationType::evaluateForward:
#ifdef DEBUG
      mexPrintf("EVALUATE FORWARD\n");
      mexEvalString("drawnow;");
#endif
      evaluate_over_periods(true);
      break;
    case BlockSimulationType::evaluateBackward:
#ifdef DEBUG
      mexPrintf("EVALUATE BACKWARD\n");
      mexEvalString("drawnow;");
#endif
      evaluate_over_periods(false);
      break;
    case BlockSimulationType::solveForwardSimple:
#ifdef DEBUG
      mexPrintf("SOLVE FORWARD SIMPLE size=%d\n", size);
      mexEvalString("drawnow;");
#endif
      solve_simple_over_periods(true);
      break;
    case BlockSimulationType::solveBackwardSimple:
#ifdef DEBUG
      mexPrintf("SOLVE BACKWARD SIMPLE\n");
      mexEvalString("drawnow;");
#endif
      solve_simple_over_periods(false);
      break;
    case BlockSimulationType::solveForwardComplete:
#ifdef DEBUG
      mexPrintf("SOLVE FORWARD COMPLETE\n");
      mexEvalString("drawnow;");
#endif
      if (vector_table_conditional_local.size())
        {
          auto curr_it_code = it_code;
          evaluate_a_block(true);
          it_code = curr_it_code;
        }
      else
        {
          fixe_u(&u, u_count_int, u_count_int);
          Read_SparseMatrix(bin_base_name, size, 1, 0, 0, false, stack_solve_algo, solve_algo);
        }
      start_code = it_code;
      Per_u_ = 0;

      Simulate_Newton_One_Boundary(true);

      mxFree(u);
      mxFree(index_equa);
      mxFree(index_vara);
      memset(direction, 0, size_of_direction);
      End_Solver();
      break;
    case BlockSimulationType::solveBackwardComplete:
#ifdef DEBUG
      mexPrintf("SOLVE BACKWARD COMPLETE\n");
      mexEvalString("drawnow;");
#endif
      if (vector_table_conditional_local.size())
        {
          auto curr_it_code = it_code;
          evaluate_a_block(true);
          it_code = curr_it_code;
        }
      else
        {
          fixe_u(&u, u_count_int, u_count_int);
          Read_SparseMatrix(bin_base_name, size, 1, 0, 0, false, stack_solve_algo, solve_algo);
        }
      start_code = it_code;
      Per_u_ = 0;

      Simulate_Newton_One_Boundary(false);

      mxFree(index_equa);
      mxFree(index_vara);
      memset(direction, 0, size_of_direction);
      mxFree(u);
      End_Solver();
      break;
    case BlockSimulationType::solveTwoBoundariesSimple:
    case BlockSimulationType::solveTwoBoundariesComplete:
#ifdef DEBUG
      mexPrintf("SOLVE TWO BOUNDARIES\n");
      mexEvalString("drawnow;");
#endif
      if (steady_state)
        {
          mexPrintf("SOLVE TWO BOUNDARIES in a steady state model: impossible case\n");
          return ERROR_ON_EXIT;
        }
      if (vector_table_conditional_local.size())
        {
          auto curr_it_code = it_code;
          evaluate_a_block(true);
          it_code = curr_it_code;
        }
      else
        {
          fixe_u(&u, u_count_int, u_count_int);
          Read_SparseMatrix(bin_base_name, size, periods, y_kmin, y_kmax, true, stack_solve_algo, solve_algo);
        }
      u_count = u_count_int*(periods+y_kmax+y_kmin);
      r = static_cast<double *>(mxMalloc(size*sizeof(double)));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size*sizeof(double));
      res = static_cast<double *>(mxMalloc(size*periods*sizeof(double)));
      test_mxMalloc(res, __LINE__, __FILE__, __func__, size*periods*sizeof(double));
      y_save = static_cast<double *>(mxMalloc(y_size*sizeof(double)*(periods+y_kmax+y_kmin)));
      test_mxMalloc(y_save, __LINE__, __FILE__, __func__, y_size*sizeof(double)*(periods+y_kmax+y_kmin));
      start_code = it_code;
      iter = 0;
      if (!is_linear
          || stack_solve_algo == 4) // On linear blocks, stack_solve_algo=4 may
                                    // need more than one iteration to find the
                                    // optimal (unitary!) path length
        {
          cvg = false;
          glambda2 = g0 = very_big;
          try_at_iteration = 0;
          int u_count_saved = u_count;
          while (!(cvg || (iter > maxit_)))
            {
              res2 = 0;
              res1 = 0;
              max_res = 0;
              max_res_idx = 0;
              copy_n(y, y_size*(periods+y_kmax+y_kmin), y_save);
              if (vector_table_conditional_local.size())
                for (auto & it1 : vector_table_conditional_local)
                  if (it1.is_cond)
                    y[it1.var_endo + y_kmin * size] = it1.constrained_value;
              compute_complete_2b(false, &res1, &res2, &max_res, &max_res_idx);
              end_code = it_code;
              if (!(isnan(res1) || isinf(res1)))
                cvg = (max_res < solve_tolf);
              if (isnan(res1) || isinf(res1) || (stack_solve_algo == 4 && iter > 0))
                copy_n(y_save, y_size*(periods+y_kmax+y_kmin), y);
              u_count = u_count_saved;
              int prev_iter = iter;
              Simulate_Newton_Two_Boundaries(block_num, symbol_table_endo_nbr, y_kmin, y_kmax, size, periods, cvg, minimal_solving_periods, stack_solve_algo, vector_table_conditional_local);
              iter++;
              if (iter > prev_iter)
                {
                  g0 = res2;
                  gp0 = -res2;
                  try_at_iteration = 0;
                  slowc_save = slowc;
                }
            }
          if (!cvg)
            throw FatalException{"In Solve two boundaries, convergence not achieved in block "
                                 + to_string(block_num+1) + ", after "
                                 + to_string(iter) + " iterations"};
        }
      else
        {
          res1 = 0;
          res2 = 0;
          max_res = 0; max_res_idx = 0;

          compute_complete_2b(false, &res1, &res2, &max_res, &max_res_idx);
          end_code = it_code;

          cvg = false;
          Simulate_Newton_Two_Boundaries(block_num, symbol_table_endo_nbr, y_kmin, y_kmax, size, periods, cvg, minimal_solving_periods, stack_solve_algo, vector_table_conditional_local);
          max_res = 0; max_res_idx = 0;
        }
      slowc = 1; // slowc is modified when stack_solve_algo=4, so restore it
      it_code = end_code;
      if (r)
        mxFree(r);
      if (y_save)
        mxFree(y_save);
      if (u)
        mxFree(u);
      if (index_vara)
        mxFree(index_vara);
      if (index_equa)
        mxFree(index_equa);
      if (res)
        mxFree(res);
      memset(direction, 0, size_of_direction);
      End_Solver();
      break;
    default:
      throw FatalException{"In simulate_a_block, Unknown type = " + to_string(static_cast<int>(type))};
      return ERROR_ON_EXIT;
    }
  return NO_ERROR_ON_EXIT;
}

void
Interpreter::print_a_block()
{
  it_code_type begining;
  if (block < 0)
    mexPrintf("\nBlock %d\n", block_num+1);
  else
    mexPrintf("\nBlock %d\n", block+1);
  mexPrintf("----------\n");
  if (steady_state)
    residual = vector<double>(size);
  else
    residual = vector<double>(size*(periods+y_kmin));
  bool go_on = true;
  bool space = false;
  while (go_on)
    {
      if ((*it_code)->op_code == Tags::FENDBLOCK)
        {
          go_on = false;
          it_code++;
        }
      else
        {
          string s;
          tie(s, it_code) = print_expression(it_code);
          if (s == "if (evaluate)" || s == "else")
            space = false;
          if (s.length() > 0)
            {
              if (space)
                mexPrintf("  %s\n", s.c_str());
              else
                mexPrintf("%s\n", s.c_str());
              mexEvalString("drawnow;");
            }
          if (s == "if (evaluate)" || s == "else")
            space = true;
        }
    }
}

void
Interpreter::ReadCodeFile(const string &file_name, CodeLoad &code)
{
  filesystem::path codfile {file_name + "/model/bytecode/" + (block_decomposed ? "block/" : "")
    + (steady_state ? "static" : "dynamic") + ".cod"};

  //First read and store in memory the code
  code_liste = code.get_op_code(codfile);
  EQN_block_number = code.get_block_number();
  if (!code_liste.size())
    throw FatalException{"In compute_blocks, " + codfile.string() + " cannot be opened"};
  if (block >= code.get_block_number())
    throw FatalException{"In compute_blocks, input argument block = " + to_string(block+1)
                         + " is greater than the number of blocks in the model ("
                         + to_string(code.get_block_number()) + " see M_.block_structure_stat.block)"};
}

void
Interpreter::check_for_controlled_exo_validity(FBEGINBLOCK_ *fb, const vector<s_plan> &sconstrained_extended_path)
{
  vector<int> exogenous = fb->get_exogenous();
  vector<int> endogenous = fb->get_endogenous();
  for (auto & it : sconstrained_extended_path)
    {
      if (find(endogenous.begin(), endogenous.end(), it.exo_num) != endogenous.end()
          && find(exogenous.begin(), exogenous.end(), it.var_num) == exogenous.end())
        throw FatalException{"\nThe conditional forecast involving as constrained variable "
                             + symbol_table.getName(SymbolType::endogenous, it.exo_num)
                             + " and as endogenized exogenous " + symbol_table.getName(SymbolType::exogenous, it.var_num)
                             + " that do not appear in block=" + to_string(Block_Count+1)
                             + ")\nYou should not use block in model options"};
      else if (find(endogenous.begin(), endogenous.end(), it.exo_num) != endogenous.end()
               && find(exogenous.begin(), exogenous.end(), it.var_num) != exogenous.end()
               && (fb->get_type() == BlockSimulationType::evaluateForward
                   || fb->get_type() == BlockSimulationType::evaluateBackward))
        throw FatalException{"\nThe conditional forecast cannot be implemented for the block="
                             + to_string(Block_Count+1) + ") that has to be evaluated instead to be solved\nYou should not use block in model options"};
      else if (find(previous_block_exogenous.begin(), previous_block_exogenous.end(), it.var_num)
               != previous_block_exogenous.end())
        throw FatalException{"\nThe conditional forecast involves in the block "
                             + to_string(Block_Count+1) + " the endogenized exogenous "
                             + symbol_table.getName(SymbolType::exogenous, it.var_num)
                             + " that appear also in a previous block\nYou should not use block in model options"};
    }
  for (auto it : exogenous)
    previous_block_exogenous.push_back(it);
}

bool
Interpreter::MainLoop(const string &bin_basename, const CodeLoad &code, bool evaluate, int block, bool last_call, bool constrained, const vector<s_plan> &sconstrained_extended_path, const vector_table_conditional_local_type &vector_table_conditional_local)
{
  int var;
  Block_Count = -1;
  bool go_on = true;
  double max_res_local = 0;
  int max_res_idx_local = 0;

  if (block < 0)
    {
      if (steady_state)
        residual = vector<double>(y_size);
      else
        residual = vector<double>(y_size*(periods+y_kmin));
    }

  while (go_on)
    {
      switch ((*it_code)->op_code)
        {
        case Tags::FBEGINBLOCK:
          Block_Count++;
#ifdef DEBUG
          mexPrintf("---------------------------------------------------------\n");
          if (block < 0)
            mexPrintf("FBEGINBLOCK Block_Count=%d\n", Block_Count+1);
          else
            mexPrintf("FBEGINBLOCK block=%d\n", block+1);
#endif
          //it's a new block
          {
            auto *fb = static_cast<FBEGINBLOCK_ *>(*it_code);
            Block_Contain = fb->get_Block_Contain();
            it_code++;
            if (constrained)
              check_for_controlled_exo_validity(fb, sconstrained_extended_path);
            set_block(fb->get_size(), fb->get_type(), file_name, bin_basename, Block_Count, fb->get_is_linear(), fb->get_endo_nbr(), fb->get_u_count_int(), block);
            if (print)
              print_a_block();
            else if (evaluate)
              {
#ifdef DEBUG
                mexPrintf("jacobian_block=mxCreateDoubleMatrix(%d, %d, mxREAL)\n", fb->get_size(), fb->get_nb_col_jacob());
#endif
                jacobian_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_nb_col_jacob(), mxREAL));
                if (!steady_state)
                  {
#ifdef DEBUG
                    mexPrintf("allocates jacobian_exo_block( %d, %d, mxREAL)\n", fb->get_size(), fb->get_exo_size());
                    mexPrintf("(0) Allocating Jacobian\n");
#endif

                    jacobian_exo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_exo_size(), mxREAL));
                    jacobian_det_exo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_det_exo_size(), mxREAL));
                  }
                if (block >= 0)
                  {
                    if (steady_state)
                      residual = vector<double>(fb->get_size());
                    else
                      residual = vector<double>(fb->get_size()*(periods+y_kmin));
                  }
                evaluate_a_block(true);
              }
            else
              {
#ifdef DEBUG
                mexPrintf("endo in Block_Count=%d, size=%d, type=%d, steady_state=%d, print_it=%d, fb->get_is_linear()=%d, fb->get_endo_nbr()=%d, fb->get_Max_Lag()=%d, fb->get_Max_Lead()=%d, fb->get_u_count_int()=%d\n",
                          Block_Count+1, fb->get_size(), fb->get_type(), steady_state, print_it, fb->get_is_linear(), fb->get_endo_nbr(), fb->get_Max_Lag(), fb->get_Max_Lead(), fb->get_u_count_int());
#endif
                bool result;
                if (sconstrained_extended_path.size())
                  {
                    //mexPrintf("(1) Allocating Jacobian fb->get_size()=%d fb->get_nb_col_jacob()=%d\n", fb->get_size(), fb->get_nb_col_jacob());
                    jacobian_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_nb_col_jacob(), mxREAL));
                    //mexPrintf("mxGetPr(jacobian_block[block_num])=%x\n",mxGetPr(jacobian_block[0]));
                    jacobian_exo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_exo_size(), mxREAL));
                    jacobian_det_exo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_det_exo_size(), mxREAL));
                    residual = vector<double>(fb->get_size()*(periods+y_kmin));
                    result = simulate_a_block(vector_table_conditional_local);

                    mxDestroyArray(jacobian_block.back());
                    jacobian_block.pop_back();
                    mxDestroyArray(jacobian_exo_block.back());
                    jacobian_exo_block.pop_back();
                    mxDestroyArray(jacobian_det_exo_block.back());
                    jacobian_det_exo_block.pop_back();
                  }
                else
                  result = simulate_a_block(vector_table_conditional_local);
                //mexPrintf("OKe\n");
                if (max_res > max_res_local)
                  {
                    max_res_local = max_res;
                    max_res_idx_local = max_res_idx;
                  }
                if (result == ERROR_ON_EXIT)
                  return ERROR_ON_EXIT;
              }
            if (last_call)
              delete fb;
          }
          if (block >= 0)
            go_on = false;

          break;
        case Tags::FEND:
#ifdef DEBUG
          mexPrintf("FEND\n");
#endif
          go_on = false;
          it_code++;
          break;
        case Tags::FDIMT:
#ifdef DEBUG
          mexPrintf("FDIMT size=%d\n", static_cast<FDIMT_ *>(*it_code)->get_size());
#endif
          var = static_cast<FDIMT_ *>(*it_code)->get_size();
          if (T)
            mxFree(T);
          T = static_cast<double *>(mxMalloc(var*(periods+y_kmin+y_kmax)*sizeof(double)));
          test_mxMalloc(T, __LINE__, __FILE__, __func__, var*(periods+y_kmin+y_kmax)*sizeof(double));
          if (block >= 0)
            it_code = code_liste.begin() + code.get_begin_block(block);
          else
            it_code++;
          break;
        case Tags::FDIMST:
#ifdef DEBUG
          mexPrintf("FDIMST size=%d\n", static_cast<FDIMST_ *>(*it_code)->get_size());
#endif
          var = static_cast<FDIMST_ *>(*it_code)->get_size();
          if (T)
            mxFree(T);
          if (global_temporary_terms)
            {
              if (!GlobalTemporaryTerms)
                {
                  mexPrintf("GlobalTemporaryTerms is nullptr\n");
                  mexEvalString("drawnow;");
                }
              if (var != static_cast<int>(mxGetNumberOfElements(GlobalTemporaryTerms)))
                GlobalTemporaryTerms = mxCreateDoubleMatrix(var, 1, mxREAL);
              T = mxGetPr(GlobalTemporaryTerms);
            }
          else
            {
              T = static_cast<double *>(mxMalloc(var*sizeof(double)));
              test_mxMalloc(T, __LINE__, __FILE__, __func__, var*sizeof(double));
            }

          if (block >= 0)
            it_code = code_liste.begin() + code.get_begin_block(block);
          else
            it_code++;
          break;
        default:
          throw FatalException{"In compute_blocks, unknown command "
                               + to_string(static_cast<int>((*it_code)->op_code)) + " (block="
                               + to_string(Block_Count) + ")"};
        }
    }
  max_res = max_res_local;
  max_res_idx = max_res_idx_local;
  Close_SaveCode();
  return true;
}

string
Interpreter::elastic(string str, unsigned int len, bool left)
{
  if (str.length() > len)
    return str;
  else
    {
      int diff = len - str.length();
      if (diff % 2 == 0)
        {
          if (left)
            {
              //mexPrintf("(1) diff=%d\n",diff);
              str.insert(str.end(), diff-1, ' ');
              str.insert(str.begin(), 1, ' ');
            }
          else
            {
              str.insert(str.end(), diff/2, ' ');
              str.insert(str.begin(), diff/2, ' ');
            }
        }
      else
        {
          if (left)
            {
              //mexPrintf("(2) diff=%d\n",diff);
              str.insert(str.end(), diff-1, ' ');
              str.insert(str.begin(), 1, ' ');
            }
          else
            {
              str.insert(str.end(), ceil(diff/2), ' ');
              str.insert(str.begin(), ceil(diff/2+1), ' ');
            }
        }
      return str;
    }
}

bool
Interpreter::extended_path(const string &file_name, const string &bin_basename, bool evaluate, int block, int &nb_blocks, int nb_periods, const vector<s_plan> &sextended_path, const vector<s_plan> &sconstrained_extended_path, const vector<string> &dates, const table_conditional_global_type &table_conditional_global)
{
  CodeLoad code;

  ReadCodeFile(file_name, code);
  it_code = code_liste.begin();
  it_code_type Init_Code = code_liste.begin();
  size_t size_of_direction = y_size*col_y*sizeof(double);
  auto *y_save = static_cast<double *>(mxMalloc(size_of_direction));
  test_mxMalloc(y_save, __LINE__, __FILE__, __func__, size_of_direction);
  auto *x_save = static_cast<double *>(mxMalloc(nb_row_x * col_x *sizeof(double)));
  test_mxMalloc(x_save, __LINE__, __FILE__, __func__, nb_row_x * col_x *sizeof(double));

  vector_table_conditional_local_type vector_table_conditional_local;
  vector_table_conditional_local.clear();

  int endo_name_length_l = static_cast<int>(symbol_table.maxEndoNameLength());
  for (int j = 0; j < col_x* nb_row_x; j++)
    {
      x_save[j] = x[j];
      x[j] = 0;
    }
  for (int j = 0; j < col_x; j++)
    x[y_kmin + j * nb_row_x] = x_save[y_kmin + j * nb_row_x];
  for (int i = 0; i < y_size * col_y; i++)
    y_save[i] = y[i];
  if (endo_name_length_l < 8)
    endo_name_length_l = 8;
  bool old_print_it = print_it;
  print_it = false;
  ostringstream res1;
  res1 << std::scientific << 2.54656875434865131;
  int real_max_length = res1.str().length();
  int date_length = dates[0].length();
  int table_length = 2 + date_length + 3 + endo_name_length_l + 3 + real_max_length + 3 + 3 + 2 + 6 + 2;
  string line;
  line.insert(line.begin(), table_length, '-');
  line.insert(line.length(), "\n");
  if (old_print_it)
    {
      mexPrintf("\nExtended Path simulation:\n");
      mexPrintf("-------------------------\n");
      mexPrintf(line.c_str());
      string title = "|" + elastic("date", date_length+2, false) + "|" + elastic("variable", endo_name_length_l+2, false) + "|" + elastic("max. value", real_max_length+2, false) + "| iter. |" + elastic("cvg", 5, false) + "|\n";
      mexPrintf(title.c_str());
      mexPrintf(line.c_str());
    }
  for (int t = 0; t < nb_periods; t++)
    {
      nb_blocks = 0;
      previous_block_exogenous.clear();
      if (old_print_it)
        {
          mexPrintf("|%s|", elastic(dates[t], date_length+2, false).c_str());
          mexEvalString("drawnow;");
        }
      for (const auto & it : sextended_path)
        x[y_kmin + (it.exo_num - 1) * nb_row_x] = it.value[t];

      it_code = Init_Code;
      vector_table_conditional_local.clear();
      if (auto it = table_conditional_global.find(t); it != table_conditional_global.end())
        vector_table_conditional_local = it->second;
      if (t < nb_periods)
        MainLoop(bin_basename, code, evaluate, block, false, true, sconstrained_extended_path, vector_table_conditional_local);
      else
        MainLoop(bin_basename, code, evaluate, block, true, true, sconstrained_extended_path, vector_table_conditional_local);
      for (int j = 0; j < y_size; j++)
        {
          y_save[j + (t + y_kmin) * y_size] = y[j + y_kmin * y_size];
          if (y_kmin > 0)
            y[j] = y[j + y_kmin * y_size];
        }
      for (int j = 0; j < col_x; j++)
        {
          x_save[t + y_kmin + j * nb_row_x] = x[y_kmin + j * nb_row_x];
          if (t < nb_periods)
            x[y_kmin + j * nb_row_x] = x_save[t + 1 + y_kmin + j * nb_row_x];
        }

      if (old_print_it)
        {
          ostringstream res1;
          res1 << std::scientific << max_res;
          mexPrintf("%s|%s| %4d  |  x  |\n", elastic(symbol_table.getName(SymbolType::endogenous, max_res_idx), endo_name_length_l+2, true).c_str(), elastic(res1.str(), real_max_length+2, false).c_str(), iter);
          mexPrintf(line.c_str());
          mexEvalString("drawnow;");
        }
    }
  print_it = old_print_it;
  for (int i = 0; i < y_size * col_y; i++)
    y[i] = y_save[i];
  for (int j = 0; j < col_x * nb_row_x; j++)
    x[j] = x_save[j];
  if (*Init_Code)
    mxFree(*Init_Code);
  if (y_save)
    mxFree(y_save);
  if (x_save)
    mxFree(x_save);
  nb_blocks = Block_Count+1;
  if (T && !global_temporary_terms)
    mxFree(T);
  return true;
}

bool
Interpreter::compute_blocks(const string &file_name, const string &bin_basename, bool evaluate, int block, int &nb_blocks)
{
  CodeLoad code;
  ReadCodeFile(file_name, code);

  //The big loop on intructions
  it_code = code_liste.begin();
  auto Init_Code = it_code;
  vector<s_plan> s_plan_junk;
  vector_table_conditional_local_type vector_table_conditional_local_junk;

  MainLoop(bin_basename, code, evaluate, block, true, false, s_plan_junk, vector_table_conditional_local_junk);

  mxFree(*Init_Code);
  nb_blocks = Block_Count+1;
  if (T && !global_temporary_terms)
    mxFree(T);
  return true;
}
