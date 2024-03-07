/*
 * Copyright © 2007-2024 Dynare Team
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

#include <algorithm>
#include <array>
#include <cassert>
#include <cfenv>
#include <chrono>
#include <filesystem>
#include <limits>
#include <numeric>
#include <sstream>
#include <type_traits>

#include "Interpreter.hh"

Interpreter::Interpreter(Evaluate& evaluator_arg, double* params_arg, double* y_arg, double* ya_arg,
                         double* x_arg, double* steady_y_arg, double* direction_arg, int y_size_arg,
                         int nb_row_x_arg, int periods_arg, int y_kmin_arg, int y_kmax_arg,
                         int maxit_arg_, double solve_tolf_arg, double markowitz_c_arg,
                         int minimal_solving_periods_arg, int stack_solve_algo_arg,
                         int solve_algo_arg, bool print_arg,
                         const mxArray* GlobalTemporaryTerms_arg, bool steady_state_arg,
                         bool block_decomposed_arg, int col_x_arg, int col_y_arg,
                         const BasicSymbolTable& symbol_table_arg, int verbosity_arg) :
    symbol_table {symbol_table_arg},
    steady_state {steady_state_arg},
    block_decomposed {block_decomposed_arg},
    evaluator {evaluator_arg},
    minimal_solving_periods {minimal_solving_periods_arg},
    y_size {y_size_arg},
    y_kmin {y_kmin_arg},
    y_kmax {y_kmax_arg},
    periods {periods_arg},
    verbosity {verbosity_arg}
{
  pivotva = nullptr;
  mem_mngr.init_Mem();
  symbolic = true;
  alt_symbolic = false;
  alt_symbolic_count = 0;
  res1a = 9.0e60;
  tbreak_g = 0;
  start_compare = 0;
  restart = 0;
  IM_i.clear();
  lu_inc_tol = 1e-10;
  Symbolic = nullptr;
  Numeric = nullptr;
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
  slowc = 1;
  slowc_save = 1;
  markowitz_c = markowitz_c_arg;
  minimal_solving_periods = minimal_solving_periods_arg;
  stack_solve_algo = stack_solve_algo_arg;
  solve_algo = solve_algo_arg;
  print = print_arg;
  col_x = col_x_arg;
  col_y = col_y_arg;

  int ntt {evaluator.getNumberOfTemporaryTerms()};
  if (GlobalTemporaryTerms_arg)
    {
      if (steady_state)
        assert(ntt == static_cast<int>(mxGetNumberOfElements(GlobalTemporaryTerms_arg)));
      else
        assert(periods * ntt == static_cast<int>(mxGetNumberOfElements(GlobalTemporaryTerms_arg)));
      GlobalTemporaryTerms = mxDuplicateArray(GlobalTemporaryTerms_arg);
    }
  else
    {
      if (steady_state)
        GlobalTemporaryTerms = mxCreateDoubleMatrix(ntt, 1, mxREAL);
      else
        GlobalTemporaryTerms = mxCreateDoubleMatrix(periods, ntt, mxREAL);
    }
  T = mxGetPr(GlobalTemporaryTerms);
}

void
Interpreter::evaluate_over_periods(bool forward)
{
  if (steady_state)
    compute_block_time(0, false, false);
  else
    {
      if (forward)
        {
          for (it_ = y_kmin; it_ < periods + y_kmin; it_++)
            compute_block_time(0, false, false);
          it_ = periods + y_kmin - 1; // Do not leave it_ in inconsistent state
        }
      else
        {
          for (it_ = periods + y_kmin - 1; it_ >= y_kmin; it_--)
            compute_block_time(0, false, false);
          it_ = y_kmin; // Do not leave it_ in inconsistent state (see #1727)
        }
    }
}

void
Interpreter::solve_simple_one_periods()
{
  bool cvg = false;
  int iter = 0;
  double ya;
  double slowc = 1;
  res1 = 0;
  while (!(cvg || iter >= maxit_))
    {
      Per_y_ = it_ * y_size;
      ya = y[Block_Contain[0].Variable + Per_y_];
      compute_block_time(0, false, false);
      if (!isfinite(res1))
        {
          res1 = std::numeric_limits<double>::quiet_NaN();
          while ((isinf(res1) || isnan(res1)) && (slowc > 1e-9))
            {
              compute_block_time(0, false, false);
              if (!isfinite(res1))
                {
                  slowc /= 1.5;
                  if (verbosity >= 2)
                    mexPrintf("Reducing the path length in Newton step slowc=%f\n", slowc);
                  feclearexcept(FE_ALL_EXCEPT);
                  y[Block_Contain[0].Variable + Per_y_] = ya - slowc * (r[0] / g1);
                  if (fetestexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW))
                    {
                      res1 = numeric_limits<double>::quiet_NaN();
                      if (verbosity >= 1)
                        mexPrintf("      Singularity in block %d", block_num + 1);
                    }
                }
            }
        }
      double rr;
      rr = r[0];
      cvg = (fabs(rr) < solve_tolf);
      if (cvg)
        continue;
      feclearexcept(FE_ALL_EXCEPT);
      y[Block_Contain[0].Variable + Per_y_] += -slowc * (rr / g1);
      if (fetestexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW))
        {
          res1 = numeric_limits<double>::quiet_NaN();
          if (verbosity >= 1)
            mexPrintf("      Singularity in block %d", block_num + 1);
        }
      iter++;
    }
  if (!cvg)
    throw FatalException {"In Solve Forward simple, convergence not achieved in block "
                          + to_string(block_num + 1) + ", after " + to_string(iter)
                          + " iterations"};
}

void
Interpreter::solve_simple_over_periods(bool forward)
{
  r = static_cast<double*>(mxMalloc(sizeof(double)));
  test_mxMalloc(r, __LINE__, __FILE__, __func__, sizeof(double));
  if (steady_state)
    {
      it_ = 0;
      solve_simple_one_periods();
    }
  else
    {
      if (forward)
        {
          for (it_ = y_kmin; it_ < periods + y_kmin; it_++)
            solve_simple_one_periods();
          it_ = periods + y_kmin - 1; // Do not leave it_ in inconsistent state
        }
      else
        {
          for (it_ = periods + y_kmin - 1; it_ >= y_kmin; it_--)
            solve_simple_one_periods();
          it_ = y_kmin; // Do not leave it_ in inconsistent state (see #1727)
        }
    }
  mxFree(r);
}

void
Interpreter::compute_complete_2b()
{
  res1 = 0;
  res2 = 0;
  max_res = 0;
  for (it_ = y_kmin; it_ < periods + y_kmin; it_++)
    {
      Per_u_ = (it_ - y_kmin) * u_count_int;
      Per_y_ = it_ * y_size;
      int shift = (it_ - y_kmin) * size;
      compute_block_time(Per_u_, false, false);
      if (!(isnan(res1) || isinf(res1)))
        for (int i = 0; i < size; i++)
          {
            double rr;
            rr = r[i];
            res[i + shift] = rr;
            if (max_res < fabs(rr))
              {
                max_res = fabs(rr);
                max_res_idx = i;
              }
            res2 += rr * rr;
            res1 += fabs(rr);
          }
      else
        return;
    }
  it_ = periods + y_kmin - 1; // Do not leave it_ in inconsistent state
}

void
Interpreter::evaluate_a_block(bool initialization, bool single_block, const string& bin_base_name)
{
  switch (type)
    {
    case BlockSimulationType::evaluateForward:
      if (steady_state)
        {
          compute_block_time(0, true, false);
          if (single_block)
            for (int j = 0; j < size; j++)
              residual[j] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
          else
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation]
                  = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
        }
      else
        {
          for (it_ = y_kmin; it_ < periods + y_kmin; it_++)
            {
              Per_y_ = it_ * y_size;
              compute_block_time(0, true, false);
              if (single_block)
                for (int j = 0; j < size; j++)
                  residual[(it_ - y_kmin) * size + j]
                      = y[it_ * y_size + Block_Contain[j].Variable]
                        - ya[it_ * y_size + Block_Contain[j].Variable];
              else
                for (int j = 0; j < size; j++)
                  residual[(it_ - y_kmin) * y_size + Block_Contain[j].Equation]
                      = y[it_ * y_size + Block_Contain[j].Variable]
                        - ya[it_ * y_size + Block_Contain[j].Variable];
            }
        }
      break;
    case BlockSimulationType::solveForwardSimple:
      r = static_cast<double*>(mxMalloc(size * sizeof(double)));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size * sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, false);
          if (!single_block)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          for (it_ = y_kmin; it_ < periods + y_kmin; it_++)
            {
              Per_y_ = it_ * y_size;
              compute_block_time(0, true, false);
              if (!single_block)
                for (int j = 0; j < size; j++)
                  residual[(it_ - y_kmin) * y_size + Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[(it_ - y_kmin) * size + j] = r[j];
            }
        }
      mxFree(r);
      break;
    case BlockSimulationType::solveForwardComplete:
      if (initialization)
        {
          fixe_u();
          Read_SparseMatrix(bin_base_name, false);
        }
#ifdef DEBUG
      mexPrintf("in SOLVE FORWARD COMPLETE r = mxMalloc(%d*sizeof(double))\n", size);
#endif
      r = static_cast<double*>(mxMalloc(size * sizeof(double)));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size * sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, false);
          if (!single_block)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          for (it_ = y_kmin; it_ < periods + y_kmin; it_++)
            {
              Per_y_ = it_ * y_size;
              compute_block_time(0, true, false);
              if (!single_block)
                for (int j = 0; j < size; j++)
                  residual[(it_ - y_kmin) * y_size + Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[(it_ - y_kmin) * size + j] = r[j];
            }
        }
      mxFree(r);
      break;
    case BlockSimulationType::evaluateBackward:
      if (steady_state)
        {
          compute_block_time(0, true, false);
          if (single_block)
            for (int j = 0; j < size; j++)
              residual[j] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
          else
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation]
                  = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
        }
      else
        {
          for (it_ = periods + y_kmin - 1; it_ >= y_kmin; it_--)
            {
              Per_y_ = it_ * y_size;
              compute_block_time(0, true, false);
              if (single_block)
                for (int j = 0; j < size; j++)
                  residual[(it_ - y_kmin) * size + j]
                      = y[it_ * y_size + Block_Contain[j].Variable]
                        - ya[it_ * y_size + Block_Contain[j].Variable];
              else
                for (int j = 0; j < size; j++)
                  residual[(it_ - y_kmin) * y_size + Block_Contain[j].Equation]
                      = y[it_ * y_size + Block_Contain[j].Variable]
                        - ya[it_ * y_size + Block_Contain[j].Variable];
            }
        }
      break;
    case BlockSimulationType::solveBackwardSimple:
      r = static_cast<double*>(mxMalloc(size * sizeof(double)));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size * sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, false);
          if (!single_block)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          for (it_ = periods + y_kmin - 1; it_ >= y_kmin; it_--)
            {
              Per_y_ = it_ * y_size;
              compute_block_time(0, true, false);
              if (!single_block)
                for (int j = 0; j < size; j++)
                  residual[(it_ - y_kmin) * y_size + Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[(it_ - y_kmin) * size + j] = r[j];
            }
        }
      mxFree(r);
      break;
    case BlockSimulationType::solveBackwardComplete:
      if (initialization)
        {
          fixe_u();
          Read_SparseMatrix(bin_base_name, false);
        }
      r = static_cast<double*>(mxMalloc(size * sizeof(double)));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size * sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, false);
          if (!single_block)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          for (it_ = periods + y_kmin - 1; it_ >= y_kmin; it_--)
            {
              Per_y_ = it_ * y_size;
              compute_block_time(0, true, false);
              if (!single_block)
                for (int j = 0; j < size; j++)
                  residual[(it_ - y_kmin) * y_size + Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[(it_ - y_kmin) * size + j] = r[j];
            }
        }
      mxFree(r);
      break;
    case BlockSimulationType::solveTwoBoundariesSimple:
    case BlockSimulationType::solveTwoBoundariesComplete:
      if (initialization)
        {
          fixe_u();
          Read_SparseMatrix(bin_base_name, true);
        }
      u_count = u_count_int * (periods + y_kmax + y_kmin);
      r = static_cast<double*>(mxMalloc(size * sizeof(double)));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size * sizeof(double));
      for (it_ = y_kmin; it_ < periods + y_kmin; it_++)
        {
          Per_u_ = (it_ - y_kmin) * u_count_int;
          Per_y_ = it_ * y_size;
          compute_block_time(Per_u_, true, false);
          if (!single_block)
            for (int j = 0; j < size; j++)
              residual[(it_ - y_kmin) * y_size + Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[(it_ - y_kmin) * size + j] = r[j];
        }
      mxFree(r);
      break;
    }
}

int
Interpreter::simulate_a_block(
    const vector_table_conditional_local_type& vector_table_conditional_local, bool single_block,
    const string& bin_base_name)
{
  max_res = 0;
  max_res_idx = 0;
  bool cvg;
  double* y_save;
#ifdef DEBUG
  mexPrintf("simulate_a_block type = %d, periods=%d, y_kmin=%d, y_kmax=%d\n", type, periods, y_kmin,
            y_kmax);
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
        evaluate_a_block(true, single_block, bin_base_name);
      else
        {
          fixe_u();
          Read_SparseMatrix(bin_base_name, false);
        }
      Per_u_ = 0;

      Simulate_Newton_One_Boundary(true);

      mxFree(u);
      mxFree(index_equa);
      mxFree(index_vara);
      fill_n(direction, y_size * col_y, 0);
      End_Solver();
      break;
    case BlockSimulationType::solveBackwardComplete:
#ifdef DEBUG
      mexPrintf("SOLVE BACKWARD COMPLETE\n");
      mexEvalString("drawnow;");
#endif
      if (vector_table_conditional_local.size())
        evaluate_a_block(true, single_block, bin_base_name);
      else
        {
          fixe_u();
          Read_SparseMatrix(bin_base_name, false);
        }
      Per_u_ = 0;

      Simulate_Newton_One_Boundary(false);

      mxFree(index_equa);
      mxFree(index_vara);
      fill_n(direction, y_size * col_y, 0);
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
          if (verbosity >= 1)
            mexPrintf("SOLVE TWO BOUNDARIES in a steady state model: impossible case\n");
          return ERROR_ON_EXIT;
        }
      if (vector_table_conditional_local.size())
        evaluate_a_block(true, single_block, bin_base_name);
      else
        {
          fixe_u();
          Read_SparseMatrix(bin_base_name, true);
        }
      u_count = u_count_int * (periods + y_kmax + y_kmin);
      r = static_cast<double*>(mxMalloc(size * sizeof(double)));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size * sizeof(double));
      res = static_cast<double*>(mxMalloc(size * periods * sizeof(double)));
      test_mxMalloc(res, __LINE__, __FILE__, __func__, size * periods * sizeof(double));
      y_save
          = static_cast<double*>(mxMalloc(y_size * sizeof(double) * (periods + y_kmax + y_kmin)));
      test_mxMalloc(y_save, __LINE__, __FILE__, __func__,
                    y_size * sizeof(double) * (periods + y_kmax + y_kmin));
      iter = 0;
      if (!is_linear || stack_solve_algo == 4) // On linear blocks, stack_solve_algo=4 may
                                               // need more than one iteration to find the
                                               // optimal (unitary!) path length
        {
          cvg = false;
          glambda2 = g0 = very_big;
          try_at_iteration = 0;
          int u_count_saved = u_count;
          while (!(cvg || (iter >= maxit_)))
            {
              res2 = 0;
              res1 = 0;
              max_res = 0;
              max_res_idx = 0;
              copy_n(y, y_size * (periods + y_kmax + y_kmin), y_save);
              if (vector_table_conditional_local.size())
                for (auto& it1 : vector_table_conditional_local)
                  if (it1.is_cond)
                    y[it1.var_endo + y_kmin * size] = it1.constrained_value;
              compute_complete_2b();
              if (!(isnan(res1) || isinf(res1)))
                cvg = (max_res < solve_tolf);
              if (isnan(res1) || isinf(res1) || (stack_solve_algo == 4 && iter > 0))
                copy_n(y_save, y_size * (periods + y_kmax + y_kmin), y);
              u_count = u_count_saved;
              int prev_iter = iter;
              Simulate_Newton_Two_Boundaries(cvg, vector_table_conditional_local);
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
            throw FatalException {"In Solve two boundaries, convergence not achieved in block "
                                  + to_string(block_num + 1) + ", after " + to_string(iter)
                                  + " iterations"};
        }
      else
        {
          res1 = 0;
          res2 = 0;
          max_res = 0;
          max_res_idx = 0;

          compute_complete_2b();

          cvg = false;
          Simulate_Newton_Two_Boundaries(cvg, vector_table_conditional_local);
          max_res = 0;
          max_res_idx = 0;
        }
      slowc = 1; // slowc is modified when stack_solve_algo=4, so restore it
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
      fill_n(direction, y_size * col_y, 0);
      End_Solver();
      break;
    default:
      throw FatalException {"In simulate_a_block, Unknown type = "
                            + to_string(static_cast<int>(type))};
      return ERROR_ON_EXIT;
    }
  return NO_ERROR_ON_EXIT;
}

void
Interpreter::check_for_controlled_exo_validity(const vector<s_plan>& sconstrained_extended_path)
{
  vector<int> exogenous {evaluator.getCurrentBlockExogenous()};
  vector<int> endogenous {evaluator.getCurrentBlockEndogenous()};
  for (auto& it : sconstrained_extended_path)
    {
      if (find(endogenous.begin(), endogenous.end(), it.exo_num) != endogenous.end()
          && find(exogenous.begin(), exogenous.end(), it.var_num) == exogenous.end())
        throw FatalException {"\nThe conditional forecast involving as constrained variable "
                              + symbol_table.getName(SymbolType::endogenous, it.exo_num)
                              + " and as endogenized exogenous "
                              + symbol_table.getName(SymbolType::exogenous, it.var_num)
                              + " that do not appear in block=" + to_string(block_num + 1)
                              + ")\nYou should not use block in model options"};
      else if (find(endogenous.begin(), endogenous.end(), it.exo_num) != endogenous.end()
               && find(exogenous.begin(), exogenous.end(), it.var_num) != exogenous.end()
               && (type == BlockSimulationType::evaluateForward
                   || type == BlockSimulationType::evaluateBackward))
        throw FatalException {"\nThe conditional forecast cannot be implemented for the block="
                              + to_string(block_num + 1)
                              + ") that has to be evaluated instead to be solved\nYou should not "
                                "use block in model options"};
      else if (find(previous_block_exogenous.begin(), previous_block_exogenous.end(), it.var_num)
               != previous_block_exogenous.end())
        throw FatalException {
            "\nThe conditional forecast involves in the block " + to_string(block_num + 1)
            + " the endogenized exogenous "
            + symbol_table.getName(SymbolType::exogenous, it.var_num)
            + " that appear also in a previous block\nYou should not use block in model options"};
    }
  for (auto it : exogenous)
    previous_block_exogenous.push_back(it);
}

pair<bool, vector<int>>
Interpreter::MainLoop(const string& bin_basename, bool evaluate, int block, bool constrained,
                      const vector<s_plan>& sconstrained_extended_path,
                      const vector_table_conditional_local_type& vector_table_conditional_local)
{
  int nb_blocks {evaluator.getTotalBlockNumber()};

  if (block >= nb_blocks)
    throw FatalException {"Interpreter::MainLoop: Input argument block = " + to_string(block + 1)
                          + " is greater than the number of blocks in the model ("
                          + to_string(nb_blocks) + " see M_.block_structure"
                          + (steady_state ? "_stat" : "") + ".block)"};

  vector<int> blocks;
  if (block < 0)
    {
      blocks.resize(nb_blocks);
      iota(blocks.begin(), blocks.end(), 0);
    }
  else
    blocks.push_back(block);

  jacobian_block.resize(nb_blocks);
  jacobian_exo_block.resize(nb_blocks);
  jacobian_det_exo_block.resize(nb_blocks);

  double max_res_local = 0;
  int max_res_idx_local = 0;

  if (block < 0)
    {
      if (steady_state)
        residual = vector<double>(y_size);
      else
        residual = vector<double>(y_size * periods);
    }

  for (int current_block : blocks)
    {
      evaluator.gotoBlock(current_block);
      block_num = current_block;
      size = evaluator.getCurrentBlockSize();
      type = evaluator.getCurrentBlockType();
      is_linear = evaluator.isCurrentBlockLinear();
      Block_Contain = evaluator.getCurrentBlockEquationsAndVariables();
      u_count_int = evaluator.getCurrentBlockUCount();

      if (constrained)
        check_for_controlled_exo_validity(sconstrained_extended_path);
      if (print)
        {
          if (steady_state)
            residual = vector<double>(size);
          else
            residual = vector<double>(size * periods);
          evaluator.printCurrentBlock();
        }
      else if (evaluate)
        {
#ifdef DEBUG
          mexPrintf("jacobian_block=mxCreateDoubleMatrix(%d, %d, mxREAL)\n", size,
                    getCurrentBlockNbColJacob());
#endif
          jacobian_block[current_block]
              = mxCreateDoubleMatrix(size, evaluator.getCurrentBlockNbColJacob(), mxREAL);
          if (!steady_state)
            {
#ifdef DEBUG
              mexPrintf("allocates jacobian_exo_block( %d, %d, mxREAL)\n", size,
                        evaluator.getCurrentBlockExoSize());
              mexPrintf("(0) Allocating Jacobian\n");
#endif

              jacobian_exo_block[current_block]
                  = mxCreateDoubleMatrix(size, evaluator.getCurrentBlockExoSize(), mxREAL);
              jacobian_det_exo_block[current_block]
                  = mxCreateDoubleMatrix(size, evaluator.getCurrentBlockExoDetSize(), mxREAL);
            }
          if (block >= 0)
            {
              if (steady_state)
                residual = vector<double>(size);
              else
                residual = vector<double>(size * periods);
            }
          evaluate_a_block(true, block >= 0, bin_basename);
        }
      else
        {
#ifdef DEBUG
          mexPrintf("endo in block %d, size=%d, type=%d, steady_state=%d, is_linear=%d, "
                    "endo_nbr=%d, u_count_int=%d\n",
                    current_block + 1, size, type, steady_state, is_linear, symbol_table_endo_nbr,
                    u_count_int);
#endif
          bool result;
          if (sconstrained_extended_path.size())
            {
              jacobian_block[current_block]
                  = mxCreateDoubleMatrix(size, evaluator.getCurrentBlockNbColJacob(), mxREAL);
              jacobian_exo_block[current_block]
                  = mxCreateDoubleMatrix(size, evaluator.getCurrentBlockExoSize(), mxREAL);
              jacobian_det_exo_block[current_block]
                  = mxCreateDoubleMatrix(size, evaluator.getCurrentBlockExoDetSize(), mxREAL);
              residual = vector<double>(size * periods);
              result = simulate_a_block(vector_table_conditional_local, block >= 0, bin_basename);
            }
          else
            result = simulate_a_block(vector_table_conditional_local, block >= 0, bin_basename);
          if (max_res > max_res_local)
            {
              max_res_local = max_res;
              max_res_idx_local = max_res_idx;
            }
          if (result == ERROR_ON_EXIT)
            return {ERROR_ON_EXIT, {}};
        }
    }

  max_res = max_res_local;
  max_res_idx = max_res_idx_local;
  Close_SaveCode();
  return {true, blocks};
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
              // mexPrintf("(1) diff=%d\n",diff);
              str.insert(str.end(), diff - 1, ' ');
              str.insert(str.begin(), 1, ' ');
            }
          else
            {
              str.insert(str.end(), diff / 2, ' ');
              str.insert(str.begin(), diff / 2, ' ');
            }
        }
      else
        {
          if (left)
            {
              // mexPrintf("(2) diff=%d\n",diff);
              str.insert(str.end(), diff - 1, ' ');
              str.insert(str.begin(), 1, ' ');
            }
          else
            {
              str.insert(str.end(), ceil(diff / 2), ' ');
              str.insert(str.begin(), ceil(diff / 2 + 1), ' ');
            }
        }
      return str;
    }
}

pair<bool, vector<int>>
Interpreter::extended_path(const string& file_name, bool evaluate, int block, int nb_periods,
                           const vector<s_plan>& sextended_path,
                           const vector<s_plan>& sconstrained_extended_path,
                           const vector<string>& dates,
                           const table_conditional_global_type& table_conditional_global)
{
  size_t size_of_direction = y_size * col_y * sizeof(double);
  auto* y_save = static_cast<double*>(mxMalloc(size_of_direction));
  test_mxMalloc(y_save, __LINE__, __FILE__, __func__, size_of_direction);
  auto* x_save = static_cast<double*>(mxMalloc(nb_row_x * col_x * sizeof(double)));
  test_mxMalloc(x_save, __LINE__, __FILE__, __func__, nb_row_x * col_x * sizeof(double));

  vector_table_conditional_local_type vector_table_conditional_local;
  vector_table_conditional_local.clear();

  int endo_name_length_l = static_cast<int>(symbol_table.maxEndoNameLength());
  for (int j = 0; j < col_x * nb_row_x; j++)
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
  int old_verbosity {verbosity};
  verbosity = 0;
  ostringstream res1;
  res1 << std::scientific << 2.54656875434865131;
  int real_max_length = res1.str().length();
  int date_length = dates[0].length();
  int table_length
      = 2 + date_length + 3 + endo_name_length_l + 3 + real_max_length + 3 + 3 + 2 + 6 + 2;
  string line;
  line.insert(line.begin(), table_length, '-');
  line.insert(line.length(), "\n");
  if (old_verbosity >= 1)
    {
      mexPrintf("\nExtended Path simulation:\n");
      mexPrintf("-------------------------\n");
      mexPrintf(line.c_str());
      string title = "|" + elastic("date", date_length + 2, false) + "|"
                     + elastic("variable", endo_name_length_l + 2, false) + "|"
                     + elastic("max. value", real_max_length + 2, false) + "| iter. |"
                     + elastic("cvg", 5, false) + "|\n";
      mexPrintf(title.c_str());
      mexPrintf(line.c_str());
    }
  bool r;
  vector<int> blocks;
  for (int t = 0; t < nb_periods; t++)
    {
      previous_block_exogenous.clear();
      if (old_verbosity >= 1)
        {
          mexPrintf("|%s|", elastic(dates[t], date_length + 2, false).c_str());
          mexEvalString("drawnow;");
        }
      for (const auto& it : sextended_path)
        x[y_kmin + (it.exo_num - 1) * nb_row_x] = it.value[t];

      vector_table_conditional_local.clear();
      if (auto it = table_conditional_global.find(t); it != table_conditional_global.end())
        vector_table_conditional_local = it->second;
      tie(r, blocks) = MainLoop(file_name, evaluate, block, true, sconstrained_extended_path,
                                vector_table_conditional_local);
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

      if (old_verbosity >= 1)
        {
          ostringstream res1;
          res1 << std::scientific << max_res;
          mexPrintf("%s|%s| %4d  |  x  |\n",
                    elastic(symbol_table.getName(SymbolType::endogenous, max_res_idx),
                            endo_name_length_l + 2, true)
                        .c_str(),
                    elastic(res1.str(), real_max_length + 2, false).c_str(), iter);
          mexPrintf(line.c_str());
          mexEvalString("drawnow;");
        }
    }
  verbosity = old_verbosity;
  for (int i = 0; i < y_size * col_y; i++)
    y[i] = y_save[i];
  for (int j = 0; j < col_x * nb_row_x; j++)
    x[j] = x_save[j];
  if (y_save)
    mxFree(y_save);
  if (x_save)
    mxFree(x_save);

  return {true, blocks};
}

pair<bool, vector<int>>
Interpreter::compute_blocks(const string& file_name, bool evaluate, int block)
{
  // The big loop on intructions
  vector<s_plan> s_plan_junk;
  vector_table_conditional_local_type vector_table_conditional_local_junk;

  auto [r, blocks] = MainLoop(file_name, evaluate, block, false, s_plan_junk,
                              vector_table_conditional_local_junk);

  return {true, blocks};
}

int
Interpreter::NRow(int r) const
{
  return NbNZRow[r];
}

int
Interpreter::NCol(int c) const
{
  return NbNZCol[c];
}

pair<int, NonZeroElem*>
Interpreter::At_Row(int r) const
{
  return {NbNZRow[r], FNZE_R[r]};
}

NonZeroElem*
Interpreter::At_Pos(int r, int c) const
{
  NonZeroElem* first {FNZE_R[r]};
  while (first->c_index != c)
    first = first->NZE_R_N;
  return first;
}

pair<int, NonZeroElem*>
Interpreter::At_Col(int c) const
{
  return {NbNZCol[c], FNZE_C[c]};
}

pair<int, NonZeroElem*>
Interpreter::At_Col(int c, int lag) const
{
  NonZeroElem* first {FNZE_C[c]};
  int i = 0;
  while (first->lag_index != lag && first)
    first = first->NZE_C_N;
  if (first)
    {
      NonZeroElem* firsta {first};
      if (!firsta->NZE_C_N)
        i++;
      else
        {
          while (firsta->lag_index == lag && firsta->NZE_C_N)
            {
              firsta = firsta->NZE_C_N;
              i++;
            }
          if (firsta->lag_index == lag)
            i++;
        }
    }
  return {i, first};
}

void
Interpreter::Delete(int r, int c)
{
  NonZeroElem *first = FNZE_R[r], *firsta = nullptr;

  while (first->c_index != c)
    {
      firsta = first;
      first = first->NZE_R_N;
    }
  if (firsta)
    firsta->NZE_R_N = first->NZE_R_N;
  if (first == FNZE_R[r])
    FNZE_R[r] = first->NZE_R_N;
  NbNZRow[r]--;

  first = FNZE_C[c];
  firsta = nullptr;
  while (first->r_index != r)
    {
      firsta = first;
      first = first->NZE_C_N;
    }

  if (firsta)
    firsta->NZE_C_N = first->NZE_C_N;
  if (first == FNZE_C[c])
    FNZE_C[c] = first->NZE_C_N;

  u_liste.push_back(first->u_index);
  mem_mngr.mxFree_NZE(first);
  NbNZCol[c]--;
}

void
Interpreter::Insert(int r, int c, int u_index, int lag_index)
{
  NonZeroElem *firstn, *first, *firsta, *a;
  firstn = mem_mngr.mxMalloc_NZE();
  first = FNZE_R[r];
  firsta = nullptr;
  while (first->c_index < c && (a = first->NZE_R_N))
    {
      firsta = first;
      first = a;
    }
  firstn->u_index = u_index;
  firstn->r_index = r;
  firstn->c_index = c;
  firstn->lag_index = lag_index;
  if (first->c_index > c)
    {
      if (first == FNZE_R[r])
        FNZE_R[r] = firstn;
      if (firsta)
        firsta->NZE_R_N = firstn;
      firstn->NZE_R_N = first;
    }
  else
    {
      first->NZE_R_N = firstn;
      firstn->NZE_R_N = nullptr;
    }
  NbNZRow[r]++;
  first = FNZE_C[c];
  firsta = nullptr;
  while (first->r_index < r && (a = first->NZE_C_N))
    {
      firsta = first;
      first = a;
    }
  if (first->r_index > r)
    {
      if (first == FNZE_C[c])
        FNZE_C[c] = firstn;
      if (firsta)
        firsta->NZE_C_N = firstn;
      firstn->NZE_C_N = first;
    }
  else
    {
      first->NZE_C_N = firstn;
      firstn->NZE_C_N = nullptr;
    }

  NbNZCol[c]++;
}

void
Interpreter::Close_SaveCode()
{
  SaveCode.close();
}

void
Interpreter::Read_SparseMatrix(const string& file_name, bool two_boundaries)
{
  unsigned int eq, var;
  int lag;
  mem_mngr.fixe_file_name(file_name);
  if (!SaveCode.is_open())
    {
      filesystem::path binfile {file_name + "/model/bytecode/" + (block_decomposed ? "block/" : "")
                                + (steady_state ? "static" : "dynamic") + ".bin"};
      SaveCode.open(binfile, ios::in | ios::binary);
      if (!SaveCode.is_open())
        throw FatalException {"In Read_SparseMatrix, " + binfile.string() + " cannot be opened"};
    }
  IM_i.clear();
  if (two_boundaries)
    {
      if (stack_solve_algo == 5)
        {
          for (int i = 0; i < u_count_init - size; i++)
            {
              int val;
              SaveCode.read(reinterpret_cast<char*>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char*>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char*>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char*>(&val), sizeof(val));
              IM_i[{eq, var, lag}] = val;
            }
          for (int j = 0; j < size; j++)
            IM_i[{j, size * (periods + y_kmax), 0}] = j;
        }
      else if ((stack_solve_algo >= 0 && stack_solve_algo <= 4) || stack_solve_algo == 6)
        {
          for (int i = 0; i < u_count_init - size; i++)
            {
              int val;
              SaveCode.read(reinterpret_cast<char*>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char*>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char*>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char*>(&val), sizeof(val));
              IM_i[{var - lag * size, -lag, eq}] = val;
            }
          for (int j = 0; j < size; j++)
            IM_i[{size * (periods + y_kmax), 0, j}] = j;
        }
      else
        throw FatalException {"Invalid value for solve_algo or stack_solve_algo"};
    }
  else
    {
      if ((stack_solve_algo == 5 && !steady_state) || (solve_algo == 5 && steady_state))
        {
          for (int i = 0; i < u_count_init; i++)
            {
              int val;
              SaveCode.read(reinterpret_cast<char*>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char*>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char*>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char*>(&val), sizeof(val));
              IM_i[{eq, var, lag}] = val;
            }
        }
      else if (steady_state
               || ((stack_solve_algo >= 0 && stack_solve_algo <= 4) || stack_solve_algo == 6))
        {
          for (int i = 0; i < u_count_init; i++)
            {
              int val;
              SaveCode.read(reinterpret_cast<char*>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char*>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char*>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char*>(&val), sizeof(val));
              IM_i[{var - lag * size, -lag, eq}] = val;
            }
        }
      else
        throw FatalException {"Invalid value for solve_algo or stack_solve_algo"};
    }

  int index_vara_size {size * (two_boundaries ? periods + y_kmin + y_kmax : 1)};
  index_vara = static_cast<int*>(mxMalloc(index_vara_size * sizeof(int)));
  test_mxMalloc(index_vara, __LINE__, __FILE__, __func__, index_vara_size * sizeof(int));
  for (int j = 0; j < size; j++)
    SaveCode.read(reinterpret_cast<char*>(&index_vara[j]), sizeof(*index_vara));
  if (two_boundaries)
    for (int i = 1; i < periods + y_kmin + y_kmax; i++)
      for (int j = 0; j < size; j++)
        index_vara[j + size * i] = index_vara[j + size * (i - 1)] + y_size;
  index_equa = static_cast<int*>(mxMalloc(size * sizeof(int)));
  test_mxMalloc(index_equa, __LINE__, __FILE__, __func__, size * sizeof(int));
  for (int j = 0; j < size; j++)
    SaveCode.read(reinterpret_cast<char*>(&index_equa[j]), sizeof(*index_equa));
}

bool
Interpreter::Simple_Init()
{
  pivot = static_cast<int*>(mxMalloc(size * sizeof(int)));
  test_mxMalloc(pivot, __LINE__, __FILE__, __func__, size * sizeof(int));
  pivot_save = static_cast<int*>(mxMalloc(size * sizeof(int)));
  test_mxMalloc(pivot_save, __LINE__, __FILE__, __func__, size * sizeof(int));
  pivotk = static_cast<int*>(mxMalloc(size * sizeof(int)));
  test_mxMalloc(pivotk, __LINE__, __FILE__, __func__, size * sizeof(int));
  pivotv = static_cast<double*>(mxMalloc(size * sizeof(double)));
  test_mxMalloc(pivotv, __LINE__, __FILE__, __func__, size * sizeof(double));
  pivotva = static_cast<double*>(mxMalloc(size * sizeof(double)));
  test_mxMalloc(pivotva, __LINE__, __FILE__, __func__, size * sizeof(double));
  b = static_cast<int*>(mxMalloc(size * sizeof(int)));
  test_mxMalloc(b, __LINE__, __FILE__, __func__, size * sizeof(int));
  line_done = static_cast<bool*>(mxMalloc(size * sizeof(bool)));
  test_mxMalloc(line_done, __LINE__, __FILE__, __func__, size * sizeof(bool));

  mem_mngr.init_CHUNK_BLCK_SIZE(u_count);
  int i = size * sizeof(NonZeroElem*);
  FNZE_R = static_cast<NonZeroElem**>(mxMalloc(i));
  test_mxMalloc(FNZE_R, __LINE__, __FILE__, __func__, i);
  FNZE_C = static_cast<NonZeroElem**>(mxMalloc(i));
  test_mxMalloc(FNZE_C, __LINE__, __FILE__, __func__, i);
  auto** temp_NZE_R = static_cast<NonZeroElem**>(mxMalloc(i));
  test_mxMalloc(temp_NZE_R, __LINE__, __FILE__, __func__, i);
  auto** temp_NZE_C = static_cast<NonZeroElem**>(mxMalloc(i));
  test_mxMalloc(temp_NZE_C, __LINE__, __FILE__, __func__, i);
  i = size * sizeof(int);
  NbNZRow = static_cast<int*>(mxMalloc(i));
  test_mxMalloc(NbNZRow, __LINE__, __FILE__, __func__, i);
  NbNZCol = static_cast<int*>(mxMalloc(i));
  test_mxMalloc(NbNZCol, __LINE__, __FILE__, __func__, i);
  for (i = 0; i < size; i++)
    {
      line_done[i] = false;
      FNZE_C[i] = nullptr;
      FNZE_R[i] = nullptr;
      temp_NZE_C[i] = nullptr;
      temp_NZE_R[i] = nullptr;
      NbNZRow[i] = 0;
      NbNZCol[i] = 0;
    }
  int u_count1 = size;
  for (auto& [key, value] : IM_i)
    {
      auto& [eq, var, lag] = key;
      if (lag == 0) /*Build the index for sparse matrix containing the jacobian : u*/
        {
          NbNZRow[eq]++;
          NbNZCol[var]++;
          NonZeroElem* first = mem_mngr.mxMalloc_NZE();
          first->NZE_C_N = nullptr;
          first->NZE_R_N = nullptr;
          first->u_index = u_count1;
          first->r_index = eq;
          first->c_index = var;
          first->lag_index = lag;
          if (!FNZE_R[eq])
            FNZE_R[eq] = first;
          if (!FNZE_C[var])
            FNZE_C[var] = first;
          if (temp_NZE_R[eq])
            temp_NZE_R[eq]->NZE_R_N = first;
          if (temp_NZE_C[var])
            temp_NZE_C[var]->NZE_C_N = first;
          temp_NZE_R[eq] = first;
          temp_NZE_C[var] = first;
          u_count1++;
        }
    }
  double cum_abs_sum = 0;
  for (int i = 0; i < size; i++)
    {
      b[i] = i;
      cum_abs_sum += fabs(u[i]);
    }
  bool zero_solution {cum_abs_sum < 1e-20};

  mxFree(temp_NZE_R);
  mxFree(temp_NZE_C);
  u_count = u_count1;

  return zero_solution;
}

bool
Interpreter::Init_Matlab_Sparse_One_Boundary(const mxArray* A_m, const mxArray* b_m,
                                             const mxArray* x0_m) const
{
  double* b = mxGetPr(b_m);
  if (!b)
    throw FatalException {"In Init_Matlab_Sparse_One_Boundary, can't retrieve b vector"};
  double* x0 = mxGetPr(x0_m);
  if (!x0)
    throw FatalException {"In Init_Matlab_Sparse_One_Boundary, can't retrieve x0 vector"};
  mwIndex* Ai = mxGetIr(A_m);
  if (!Ai)
    throw FatalException {"In Init_Matlab_Sparse_One_Boundary, can't allocate Ai index vector"};
  mwIndex* Aj = mxGetJc(A_m);
  if (!Aj)
    throw FatalException {"In Init_Matlab_Sparse_One_Boundary, can't allocate Aj index vector"};
  double* A = mxGetPr(A_m);
  if (!A)
    throw FatalException {"In Init_Matlab_Sparse_One_Boundary, can't retrieve A matrix"};

  for (int i = 0; i < y_size * (periods + y_kmin); i++)
    ya[i] = y[i];
#ifdef DEBUG
  unsigned int max_nze = mxGetNzmax(A_m);
#endif
  unsigned int NZE = 0;
  int last_var = 0;
  double cum_abs_sum = 0;
  for (int i = 0; i < size; i++)
    {
      b[i] = u[i];
      cum_abs_sum += fabs(b[i]);
      x0[i] = y[i];
    }
  bool zero_solution {cum_abs_sum < 1e-20};

  Aj[0] = 0;
  last_var = 0;
  for (auto& [key, index] : IM_i)
    {
      auto& [var, ignore, eq] = key;
      if (var != last_var)
        {
          Aj[1 + last_var] = NZE;
          last_var = var;
        }
#ifdef DEBUG
      if (index < 0 || index >= u_count_alloc || index > size + size * size)
        throw FatalException {"In Init_Matlab_Sparse_One_Boundary, index (" + to_string(index)
                              + ") out of range for u vector max = " + to_string(size + size * size)
                              + " allocated = " + to_string(u_count_alloc)};
      if (NZE >= max_nze)
        throw FatalException {
            "In Init_Matlab_Sparse_One_Boundary, exceeds the capacity of A_m sparse matrix"};
#endif
      A[NZE] = u[index];
      Ai[NZE] = eq;
      NZE++;
#ifdef DEBUG
      if (eq < 0 || eq >= size)
        throw FatalException {"In Init_Matlab_Sparse_One_Boundary, index (" + to_string(eq)
                              + ") out of range for b vector"};
      if (var < 0 || var >= size)
        throw FatalException {"In Init_Matlab_Sparse_One_Boundary, index (" + to_string(var)
                              + ") out of range for index_vara vector"};
      if (index_vara[var] < 0 || index_vara[var] >= y_size)
        throw FatalException {"In Init_Matlab_Sparse_One_Boundary, index ("
                              + to_string(index_vara[var])
                              + ") out of range for y vector max=" + to_string(y_size) + " (0)"};
#endif
    }
  Aj[size] = NZE;

  return zero_solution;
}

tuple<bool, SuiteSparse_long*, SuiteSparse_long*, double*, double*>
Interpreter::Init_UMFPACK_Sparse_One_Boundary(const mxArray* x0_m) const
{
  auto* b = static_cast<double*>(mxMalloc(size * sizeof(double)));
  test_mxMalloc(b, __LINE__, __FILE__, __func__, size * sizeof(double));
  if (!b)
    throw FatalException {"In Init_UMFPACK_Sparse_One_Boundary, can't retrieve b vector"};
  double* x0 = mxGetPr(x0_m);
  if (!x0)
    throw FatalException {"In Init_UMFPACK_Sparse_One_Boundary, can't retrieve x0 vector"};
  SuiteSparse_long* Ap
      = static_cast<SuiteSparse_long*>(mxMalloc((size + 1) * sizeof(SuiteSparse_long)));
  test_mxMalloc(Ap, __LINE__, __FILE__, __func__, (size + 1) * sizeof(SuiteSparse_long));
  if (!Ap)
    throw FatalException {"In Init_UMFPACK_Sparse_One_Boundary, can't allocate Ap index vector"};
  size_t prior_nz = IM_i.size();
  SuiteSparse_long* Ai
      = static_cast<SuiteSparse_long*>(mxMalloc(prior_nz * sizeof(SuiteSparse_long)));
  test_mxMalloc(Ai, __LINE__, __FILE__, __func__, prior_nz * sizeof(SuiteSparse_long));
  if (!Ai)
    throw FatalException {"In Init_UMFPACK_Sparse_One_Boundary, can't allocate Ai index vector"};
  auto* Ax = static_cast<double*>(mxMalloc(prior_nz * sizeof(double)));
  test_mxMalloc(Ax, __LINE__, __FILE__, __func__, prior_nz * sizeof(double));
  if (!Ax)
    throw FatalException {"In Init_UMFPACK_Sparse_One_Boundary, can't retrieve Ax matrix"};
  for (int i = 0; i < size; i++)
    {
      int eq = index_vara[i];
      ya[eq + it_ * y_size] = y[eq + it_ * y_size];
    }
#ifdef DEBUG
  unsigned int max_nze = prior_nz; // mxGetNzmax(A_m);
#endif
  unsigned int NZE = 0;
  int last_var = 0;
  double cum_abs_sum = 0;
  for (int i = 0; i < size; i++)
    {
      b[i] = u[i];
      cum_abs_sum += fabs(b[i]);
      x0[i] = y[i];
    }
  bool zero_solution {cum_abs_sum < 1e-20};

  Ap[0] = 0;
  last_var = 0;
  for (auto& [key, index] : IM_i)
    {
      auto& [var, ignore, eq] = key;
      if (var != last_var)
        {
          Ap[1 + last_var] = NZE;
          last_var = var;
        }
#ifdef DEBUG
      if (index < 0 || index >= u_count_alloc || index > size + size * size)
        throw FatalException {"In Init_UMFPACK_Sparse_One_Boundary, index (" + to_string(index)
                              + ") out of range for u vector max = " + to_string(size + size * size)
                              + " allocated = " + to_string(u_count_alloc)};
      if (NZE >= max_nze)
        throw FatalException {
            "In Init_UMFPACK_Sparse_One_Boundary, exceeds the capacity of A_m sparse matrix"};
#endif
      Ax[NZE] = u[index];
      Ai[NZE] = eq;
      NZE++;
#ifdef DEBUG
      if (eq < 0 || eq >= size)
        throw FatalException {"In Init_UMFPACK_Sparse_One_Boundary, index (" + to_string(eq)
                              + ") out of range for b vector"};
      if (var < 0 || var >= size)
        throw FatalException {"In Init_UMFPACK_Sparse_One_Boundary, index (" + to_string(var)
                              + ") out of range for index_vara vector"};
      if (index_vara[var] < 0 || index_vara[var] >= y_size)
        throw FatalException {"In Init_UMFPACK_Sparse_One_Boundary, index ("
                              + to_string(index_vara[var])
                              + ") out of range for y vector max=" + to_string(y_size) + " (0)"};
#endif
    }
  Ap[size] = NZE;

  return {zero_solution, Ap, Ai, Ax, b};
}

int
Interpreter::find_exo_num(const vector<s_plan>& sconstrained_extended_path, int value)
{
  auto it = find_if(sconstrained_extended_path.begin(), sconstrained_extended_path.end(),
                    [=](auto v) { return v.exo_num == value; });
  if (it != sconstrained_extended_path.end())
    return it - sconstrained_extended_path.begin();
  else
    return -1;
}

int
Interpreter::find_int_date(const vector<pair<int, double>>& per_value, int value)
{
  auto it = find_if(per_value.begin(), per_value.end(), [=](auto v) { return v.first == value; });
  if (it != per_value.end())
    return it - per_value.begin();
  else
    return -1;
}

tuple<SuiteSparse_long*, SuiteSparse_long*, double*, double*>
Interpreter::Init_UMFPACK_Sparse_Two_Boundaries(
    const mxArray* x0_m,
    const vector_table_conditional_local_type& vector_table_conditional_local) const
{
  int n = periods * size;
  auto* b = static_cast<double*>(mxMalloc(n * sizeof(double)));
  if (!b)
    throw FatalException {"In Init_UMFPACK_Sparse_Two_Boundaries, can't retrieve b vector"};
  double* x0 = mxGetPr(x0_m);
  if (!x0)
    throw FatalException {"In Init_UMFPACK_Sparse_Two_Boundaries, can't retrieve x0 vector"};
  SuiteSparse_long* Ap
      = static_cast<SuiteSparse_long*>(mxMalloc((n + 1) * sizeof(SuiteSparse_long)));
  test_mxMalloc(Ap, __LINE__, __FILE__, __func__, (n + 1) * sizeof(SuiteSparse_long));
  if (!Ap)
    throw FatalException {"In Init_UMFPACK_Sparse_Two_Boundaries, can't allocate Ap index vector"};
  size_t prior_nz = IM_i.size() * periods;
  SuiteSparse_long* Ai
      = static_cast<SuiteSparse_long*>(mxMalloc(prior_nz * sizeof(SuiteSparse_long)));
  test_mxMalloc(Ai, __LINE__, __FILE__, __func__, prior_nz * sizeof(SuiteSparse_long));
  if (!Ai)
    throw FatalException {"In Init_UMFPACK_Sparse_Two_Boundaries, can't allocate Ai index vector"};
  auto* Ax = static_cast<double*>(mxMalloc(prior_nz * sizeof(double)));
  test_mxMalloc(Ax, __LINE__, __FILE__, __func__, prior_nz * sizeof(double));
  if (!Ax)
    throw FatalException {"In Init_UMFPACK_Sparse_Two_Boundaries, can't retrieve Ax matrix"};
  for (int i = 0; i < y_size * (periods + y_kmin); i++)
    ya[i] = y[i];
  unsigned int NZE = 0;
  int last_var = 0;
  for (int i = 0; i < periods * size; i++)
    {
      b[i] = 0;
      x0[i] = y[index_vara[size * y_kmin + i]];
    }
  double* jacob_exo;
  int row_x = 0;
#ifdef DEBUG
  int col_x;
#endif
  if (vector_table_conditional_local.size())
    {
      jacob_exo = mxGetPr(jacobian_exo_block[block_num]);
      row_x = mxGetM(jacobian_exo_block[block_num]);
#ifdef DEBUG
      col_x = mxGetN(jacobian_exo_block[block_num]);
#endif
    }
  else
    jacob_exo = nullptr;
#ifdef DEBUG
  int local_index;
#endif

  bool fliped = false;
  bool fliped_exogenous_derivatives_updated = false;
  int flip_exo;
  Ap[0] = 0;
  for (int t = 0; t < periods; t++)
    {
      last_var = -1;
      int var = 0;
      for (auto& [key, value] : IM_i)
        {
          var = get<0>(key);
          int eq = get<2>(key) + size * t;
          int lag = -get<1>(key);
          int index = value + (t - lag) * u_count_init;
          if (var != last_var)
            {
              Ap[1 + last_var + t * size] = NZE;
              last_var = var;
              if (var < size * (periods + y_kmax))
                {
                  if (t == 0 && vector_table_conditional_local.size())
                    {
                      fliped = vector_table_conditional_local[var].is_cond;
                      fliped_exogenous_derivatives_updated = false;
                    }
                  else
                    fliped = false;
                }
              else
                fliped = false;
            }
          if (fliped)
            {
              if (t == 0 && var < (periods + y_kmax) * size && lag == 0
                  && vector_table_conditional_local.size())
                {
                  flip_exo = vector_table_conditional_local[var].var_exo;
#ifdef DEBUG
                  local_index = eq;
#endif
                  if (!fliped_exogenous_derivatives_updated)
                    {
                      fliped_exogenous_derivatives_updated = true;
                      for (int k = 0; k < row_x; k++)
                        {
                          if (jacob_exo[k + row_x * flip_exo] != 0)
                            {
                              Ax[NZE] = jacob_exo[k + row_x * flip_exo];
                              Ai[NZE] = k;
                              NZE++;

#ifdef DEBUG
                              if (local_index < 0 || local_index >= size * periods)
                                throw FatalException {
                                    "In Init_UMFPACK_Sparse_Two_Boundaries, index ("
                                    + to_string(local_index) + ") out of range for b vector"};
                              if (k + row_x * flip_exo < 0 || k + row_x * flip_exo >= row_x * col_x)
                                throw FatalException {
                                    "In Init_UMFPACK_Sparse_Two_Boundaries, index ("
                                    + to_string(var + size * (y_kmin + t + lag))
                                    + ") out of range for jacob_exo vector"};
                              if (t + y_kmin + flip_exo * nb_row_x < 0
                                  || t + y_kmin + flip_exo * nb_row_x >= nb_row_x * this->col_x)
                                throw FatalException {
                                    "In Init_UMFPACK_Sparse_Two_Boundaries, index ("
                                    + to_string(index_vara[var + size * (y_kmin + t + lag)])
                                    + ") out of range for x vector max="
                                    + to_string(nb_row_x * this->col_x)};
#endif
                              u[k] -= jacob_exo[k + row_x * flip_exo]
                                      * x[t + y_kmin + flip_exo * nb_row_x];
                            }
                        }
                    }
                }
            }

          if (var < (periods + y_kmax) * size)
            {
              int ti_y_kmin = -min(t, y_kmin);
              int ti_y_kmax = min(periods - (t + 1), y_kmax);
              int ti_new_y_kmax = min(t, y_kmax);
              int ti_new_y_kmin = -min(periods - (t + 1), y_kmin);
              if (lag <= ti_new_y_kmax && lag >= ti_new_y_kmin) /*Build the index for sparse matrix
                                                                   containing the jacobian : u*/
                {
#ifdef DEBUG
                  if (NZE >= prior_nz)
                    throw FatalException {"In Init_UMFPACK_Sparse_Two_Boundaries, exceeds the "
                                          "capacity of allocated sparse matrix"};
#endif
                  if (!fliped)
                    {
                      Ax[NZE] = u[index];
                      Ai[NZE] = eq - lag * size;
                      NZE++;
                    }
                  else /*if (fliped)*/
                    {
#ifdef DEBUG
                      if (eq - lag * size < 0 || eq - lag * size >= size * periods)
                        throw FatalException {"In Init_UMFPACK_Sparse_Two_Boundaries, index ("
                                              + to_string(eq - lag * size)
                                              + ") out of range for b vector"};
                      if (var + size * (y_kmin + t) < 0
                          || var + size * (y_kmin + t) >= size * (periods + y_kmin + y_kmax))
                        throw FatalException {"In Init_UMFPACK_Sparse_Two_Boundaries, index ("
                                              + to_string(var + size * (y_kmin + t))
                                              + ") out of range for index_vara vector"};
                      if (index_vara[var + size * (y_kmin + t)] < 0
                          || index_vara[var + size * (y_kmin + t)]
                                 >= y_size * (periods + y_kmin + y_kmax))
                        throw FatalException {"In Init_UMFPACK_Sparse_Two_Boundaries, index ("
                                              + to_string(index_vara[var + size * (y_kmin + t)])
                                              + ") out of range for y vector max="
                                              + to_string(y_size * (periods + y_kmin + y_kmax))};
#endif
                      b[eq - lag * size] += u[index] * y[index_vara[var + size * (y_kmin + t)]];
                    }
                }
              if (lag > ti_y_kmax || lag < ti_y_kmin)
                {
#ifdef DEBUG
                  if (eq < 0 || eq >= size * periods)
                    throw FatalException {"In Init_UMFPACK_Sparse_Two_Boundaries, index ("
                                          + to_string(eq) + ") out of range for b vector"};
                  if (var + size * (y_kmin + t + lag) < 0
                      || var + size * (y_kmin + t + lag) >= size * (periods + y_kmin + y_kmax))
                    throw FatalException {"In Init_UMFPACK_Sparse_Two_Boundaries, index ("
                                          + to_string(var + size * (y_kmin + t + lag))
                                          + ") out of range for index_vara vector"};
                  if (index_vara[var + size * (y_kmin + t + lag)] < 0
                      || index_vara[var + size * (y_kmin + t + lag)]
                             >= y_size * (periods + y_kmin + y_kmax))
                    throw FatalException {"In Init_UMFPACK_Sparse_Two_Boundaries, index ("
                                          + to_string(index_vara[var + size * (y_kmin + t + lag)])
                                          + ") out of range for y vector max="
                                          + to_string(y_size * (periods + y_kmin + y_kmax))};
#endif
                  b[eq] += u[index + lag * u_count_init]
                           * y[index_vara[var + size * (y_kmin + t + lag)]];
                }
            }
          else /* ...and store it in the u vector*/
            {
#ifdef DEBUG
              if (index < 0 || index >= u_count_alloc)
                throw FatalException {"In Init_UMFPACK_Sparse_Two_Boundaries, index ("
                                      + to_string(index) + ") out of range for u vector"};
              if (eq < 0 || eq >= (size * periods))
                throw FatalException {"In Init_UMFPACK_Sparse_Two_Boundaries, index ("
                                      + to_string(eq) + ") out of range for b vector"};
#endif
              b[eq] += u[index];
            }
        }
    }
  Ap[size * periods] = NZE;
#ifdef DEBUG
  mexPrintf("Ax = [");
  for (int i = 0; i < static_cast<int>(NZE); i++)
    mexPrintf("%f ", Ax[i]);
  mexPrintf("]\n");

  mexPrintf("Ap = [");
  for (int i = 0; i < n + 1; i++)
    mexPrintf("%d ", Ap[i]);
  mexPrintf("]\n");

  mexPrintf("Ai = [");
  for (int i = 0; i < static_cast<int>(NZE); i++)
    mexPrintf("%d ", Ai[i]);
  mexPrintf("]\n");
#endif

  return {Ap, Ai, Ax, b};
}

void
Interpreter::Init_Matlab_Sparse_Two_Boundaries(const mxArray* A_m, const mxArray* b_m,
                                               const mxArray* x0_m) const
{
  double* b = mxGetPr(b_m);

  if (!b)
    throw FatalException {"In Init_Matlab_Sparse_Two_Boundaries, can't retrieve b vector"};
  double* x0 = mxGetPr(x0_m);
  if (!x0)
    throw FatalException {"In Init_Matlab_Sparse_Two_Boundaries, can't retrieve x0 vector"};
  mwIndex* Aj = mxGetJc(A_m);
  if (!Aj)
    throw FatalException {"In Init_Matlab_Sparse_Two_Boundaries, can't allocate Aj index vector"};
  mwIndex* Ai = mxGetIr(A_m);
  if (!Ai)
    throw FatalException {"In Init_Matlab_Sparse_Two_Boundaries, can't allocate Ai index vector"};
  double* A = mxGetPr(A_m);
  if (!A)
    throw FatalException {"In Init_Matlab_Sparse_Two_Boundaries, can't retrieve A matrix"};

  for (int i = 0; i < y_size * (periods + y_kmin); i++)
    ya[i] = y[i];
  unsigned int NZE = 0;
  int last_var = 0;
  for (int i = 0; i < periods * size; i++)
    {
      b[i] = 0;
      x0[i] = y[index_vara[size * y_kmin + i]];
    }
  Aj[0] = 0;
  for (int t = 0; t < periods; t++)
    {
      last_var = 0;
      for (auto& [key, value] : IM_i)
        {
          int var = get<0>(key);
          if (var != last_var)
            {
              Aj[1 + last_var + t * size] = NZE;
              last_var = var;
            }
          int eq = get<2>(key) + size * t;
          int lag = -get<1>(key);
          int index = value + (t - lag) * u_count_init;
          if (var < (periods + y_kmax) * size)
            {
              int ti_y_kmin = -min(t, y_kmin);
              int ti_y_kmax = min(periods - (t + 1), y_kmax);
              int ti_new_y_kmax = min(t, y_kmax);
              int ti_new_y_kmin = -min(periods - (t + 1), y_kmin);
              if (lag <= ti_new_y_kmax && lag >= ti_new_y_kmin) /*Build the index for sparse matrix
                                                                   containing the jacobian : u*/
                {
#ifdef DEBUG
                  if (index < 0 || index >= u_count_alloc || index > size + size * size)
                    throw FatalException {"In Init_Matlab_Sparse_Two_Boundaries, index ("
                                          + to_string(index) + ") out of range for u vector max = "
                                          + to_string(size + size * size)
                                          + " allocated = " + to_string(u_count_alloc)};
                  if (NZE >= prior_nz)
                    throw FatalException {"In Init_Matlab_Sparse_Two_Boundaries, exceeds the "
                                          "capacity of allocated sparse matrix"};
#endif
                  A[NZE] = u[index];
                  Ai[NZE] = eq - lag * size;
                  NZE++;
                }
              if (lag > ti_y_kmax || lag < ti_y_kmin)
                {
#ifdef DEBUG
                  if (eq < 0 || eq >= size * periods)
                    throw FatalException {"In Init_Matlab_Sparse_Two_Boundaries, index ("
                                          + to_string(eq) + ") out of range for b vector"};
                  if (var + size * (y_kmin + t + lag) < 0
                      || var + size * (y_kmin + t + lag) >= size * (periods + y_kmin + y_kmax))
                    throw FatalException {"In Init_Matlab_Sparse_Two_Boundaries, index ("
                                          + to_string(var + size * (y_kmin + t + lag))
                                          + ") out of range for index_vara vector"};
                  if (index_vara[var + size * (y_kmin + t + lag)] < 0
                      || index_vara[var + size * (y_kmin + t + lag)]
                             >= y_size * (periods + y_kmin + y_kmax))
                    throw FatalException {"In Init_Matlab_Sparse_Two_Boundaries, index ("
                                          + to_string(index_vara[var + size * (y_kmin + t + lag)])
                                          + ") out of range for y vector max="
                                          + to_string(y_size * (periods + y_kmin + y_kmax))};
#endif
                  b[eq] += u[index + lag * u_count_init]
                           * y[index_vara[var + size * (y_kmin + t + lag)]];
                }
            }
          else /* ...and store it in the u vector*/
            {
#ifdef DEBUG
              if (index < 0 || index >= u_count_alloc)
                throw FatalException {"In Init_Matlab_Sparse_Two_Boundaries, index ("
                                      + to_string(index) + ") out of range for u vector"};
              if (eq < 0 || eq >= (size * periods))
                throw FatalException {"In Init_Matlab_Sparse_Two_Boundaries, index ("
                                      + to_string(eq) + ") out of range for b vector"};
#endif
              b[eq] += u[index];
            }
        }
    }
  Aj[size * periods] = NZE;
}

void
Interpreter::Init_Gaussian_Elimination()
{
  double tmp_b = 0.0;
  pivot = static_cast<int*>(mxMalloc(size * periods * sizeof(int)));
  test_mxMalloc(pivot, __LINE__, __FILE__, __func__, size * periods * sizeof(int));
  pivot_save = static_cast<int*>(mxMalloc(size * periods * sizeof(int)));
  test_mxMalloc(pivot_save, __LINE__, __FILE__, __func__, size * periods * sizeof(int));
  pivotk = static_cast<int*>(mxMalloc(size * periods * sizeof(int)));
  test_mxMalloc(pivotk, __LINE__, __FILE__, __func__, size * periods * sizeof(int));
  pivotv = static_cast<double*>(mxMalloc(size * periods * sizeof(double)));
  test_mxMalloc(pivotv, __LINE__, __FILE__, __func__, size * periods * sizeof(double));
  pivotva = static_cast<double*>(mxMalloc(size * periods * sizeof(double)));
  test_mxMalloc(pivotva, __LINE__, __FILE__, __func__, size * periods * sizeof(double));
  b = static_cast<int*>(mxMalloc(size * periods * sizeof(int)));
  test_mxMalloc(b, __LINE__, __FILE__, __func__, size * periods * sizeof(int));
  line_done = static_cast<bool*>(mxMalloc(size * periods * sizeof(bool)));
  test_mxMalloc(line_done, __LINE__, __FILE__, __func__, size * periods * sizeof(bool));
  mem_mngr.init_CHUNK_BLCK_SIZE(u_count);
  int i = (periods + y_kmax + 1) * size * sizeof(NonZeroElem*);
  FNZE_R = static_cast<NonZeroElem**>(mxMalloc(i));
  test_mxMalloc(FNZE_R, __LINE__, __FILE__, __func__, i);
  FNZE_C = static_cast<NonZeroElem**>(mxMalloc(i));
  test_mxMalloc(FNZE_C, __LINE__, __FILE__, __func__, i);
  auto** temp_NZE_R = static_cast<NonZeroElem**>(mxMalloc(i));
  test_mxMalloc(temp_NZE_R, __LINE__, __FILE__, __func__, i);
  auto** temp_NZE_C = static_cast<NonZeroElem**>(mxMalloc(i));
  test_mxMalloc(temp_NZE_C, __LINE__, __FILE__, __func__, i);
  i = (periods + y_kmax + 1) * size * sizeof(int);
  NbNZRow = static_cast<int*>(mxMalloc(i));
  test_mxMalloc(NbNZRow, __LINE__, __FILE__, __func__, i);
  NbNZCol = static_cast<int*>(mxMalloc(i));
  test_mxMalloc(NbNZCol, __LINE__, __FILE__, __func__, i);

  for (int i = 0; i < periods * size; i++)
    {
      b[i] = 0;
      line_done[i] = false;
    }
  for (int i = 0; i < (periods + y_kmax + 1) * size; i++)
    {
      FNZE_C[i] = nullptr;
      FNZE_R[i] = nullptr;
      temp_NZE_C[i] = nullptr;
      temp_NZE_R[i] = nullptr;
      NbNZRow[i] = 0;
      NbNZCol[i] = 0;
    }
  // pragma omp parallel for ordered private(it4, ti_y_kmin, ti_y_kmax, eq, var, lag)
  // schedule(dynamic)
  for (int t = 0; t < periods; t++)
    {
      int ti_y_kmin = -min(t, y_kmin);
      int ti_y_kmax = min(periods - (t + 1), y_kmax);
      int eq = -1;
      // pragma omp ordered
      for (auto& [key, value] : IM_i)
        {
          int var = get<1>(key);
          if (eq != get<0>(key) + size * t)
            tmp_b = 0;
          eq = get<0>(key) + size * t;
          if (var < (periods + y_kmax) * size)
            {
              int lag {get<2>(key)};
              if (lag <= ti_y_kmax && lag >= ti_y_kmin) /*Build the index for sparse matrix
                                                           containing the jacobian : u*/
                {
                  var += size * t;
                  NbNZRow[eq]++;
                  NbNZCol[var]++;
                  NonZeroElem* first = mem_mngr.mxMalloc_NZE();
                  first->NZE_C_N = nullptr;
                  first->NZE_R_N = nullptr;
                  first->u_index = value + u_count_init * t;
                  first->r_index = eq;
                  first->c_index = var;
                  first->lag_index = lag;
                  if (FNZE_R[eq] == nullptr)
                    FNZE_R[eq] = first;
                  if (FNZE_C[var] == nullptr)
                    FNZE_C[var] = first;
                  if (temp_NZE_R[eq] != nullptr)
                    temp_NZE_R[eq]->NZE_R_N = first;
                  if (temp_NZE_C[var] != nullptr)
                    temp_NZE_C[var]->NZE_C_N = first;
                  temp_NZE_R[eq] = first;
                  temp_NZE_C[var] = first;
                }
              else /*Build the additive terms ooutside the simulation periods related to the first
                      lags and the last leads...*/
                {
                  if (lag < ti_y_kmin)
                    tmp_b += u[value + u_count_init * t] * y[index_vara[var + size * (y_kmin + t)]];
                  else
                    tmp_b += u[value + u_count_init * t] * y[index_vara[var + size * (y_kmin + t)]];
                }
            }
          else /* ...and store it in the u vector*/
            {
              b[eq] = value + u_count_init * t;
              u[b[eq]] += tmp_b;
              tmp_b = 0;
            }
        }
    }
  mxFree(temp_NZE_R);
  mxFree(temp_NZE_C);
}

int
Interpreter::Get_u()
{
  if (!u_liste.empty())
    {
      int i = u_liste.back();
      u_liste.pop_back();
      return i;
    }
  else
    {
      if (u_count < u_count_alloc)
        {
          int i = u_count;
          u_count++;
          return i;
        }
      else
        {
          u_count_alloc += 5 * u_count_alloc_save;
          u = static_cast<double*>(mxRealloc(u, u_count_alloc * sizeof(double)));
          if (!u)
            throw FatalException {"In Get_u, memory exhausted (realloc("
                                  + to_string(u_count_alloc * sizeof(double)) + "))"};
          int i = u_count;
          u_count++;
          return i;
        }
    }
}

void
Interpreter::Delete_u(int pos)
{
  u_liste.push_back(pos);
}

void
Interpreter::Clear_u()
{
  u_liste.clear();
}

void
Interpreter::End_Gaussian_Elimination()
{
  mem_mngr.Free_All();
  mxFree(FNZE_R);
  mxFree(FNZE_C);
  mxFree(NbNZRow);
  mxFree(NbNZCol);
  mxFree(b);
  mxFree(line_done);
  mxFree(pivot);
  mxFree(pivot_save);
  mxFree(pivotk);
  mxFree(pivotv);
  mxFree(pivotva);
}

bool
Interpreter::compare(int* save_op, int* save_opa, int* save_opaa, int beg_t, long nop4)
{
  long nop = nop4 / 2;
  double r = 0.0;
  bool OK = true;
  int* diff1 = static_cast<int*>(mxMalloc(nop * sizeof(int)));
  test_mxMalloc(diff1, __LINE__, __FILE__, __func__, nop * sizeof(int));
  int* diff2 = static_cast<int*>(mxMalloc(nop * sizeof(int)));
  test_mxMalloc(diff2, __LINE__, __FILE__, __func__, nop * sizeof(int));
  int max_save_ops_first = -1;
  long j = 0, i = 0;
  while (i < nop4 && OK)
    {
      auto* save_op_s = reinterpret_cast<t_save_op_s*>(&save_op[i]);
      auto* save_opa_s = reinterpret_cast<t_save_op_s*>(&save_opa[i]);
      auto* save_opaa_s = reinterpret_cast<t_save_op_s*>(&save_opaa[i]);
      diff1[j] = save_op_s->first - save_opa_s->first;
      max_save_ops_first = max(max_save_ops_first, save_op_s->first + diff1[j] * (periods - beg_t));
      switch (save_op_s->operat)
        {
        case IFLD:
        case IFDIV:
          OK = (save_op_s->operat == save_opa_s->operat && save_opa_s->operat == save_opaa_s->operat
                && diff1[j] == (save_opa_s->first - save_opaa_s->first));
          i += 2;
          break;
        case IFLESS:
        case IFSUB:
          diff2[j] = save_op_s->second - save_opa_s->second;
          OK = (save_op_s->operat == save_opa_s->operat && save_opa_s->operat == save_opaa_s->operat
                && diff1[j] == (save_opa_s->first - save_opaa_s->first)
                && diff2[j] == (save_opa_s->second - save_opaa_s->second));
          i += 3;
          break;
        default:
          throw FatalException {"In compare, unknown operator = " + to_string(save_op_s->operat)};
        }
      j++;
    }
  // the same pivot for all remaining periods
  if (OK)
    {
      for (int i = beg_t; i < periods; i++)
        for (int j = 0; j < size; j++)
          pivot[i * size + j] = pivot[(i - 1) * size + j] + size;
      if (max_save_ops_first >= u_count_alloc)
        {
          u_count_alloc += max_save_ops_first;
          u = static_cast<double*>(mxRealloc(u, u_count_alloc * sizeof(double)));
          if (!u)
            throw FatalException {"In compare, memory exhausted (realloc("
                                  + to_string(u_count_alloc * sizeof(double)) + "))"};
        }
      for (int t = 1; t < periods - beg_t - y_kmax; t++)
        {
          int i = j = 0;
          while (i < nop4)
            {
              auto* save_op_s = reinterpret_cast<t_save_op_s*>(&save_op[i]);
              double* up = &u[save_op_s->first + t * diff1[j]];
              switch (save_op_s->operat)
                {
                case IFLD:
                  r = *up;
                  i += 2;
                  break;
                case IFDIV:
                  *up /= r;
                  i += 2;
                  break;
                case IFSUB:
                  *up -= u[save_op_s->second + t * diff2[j]] * r;
                  ;
                  i += 3;
                  break;
                case IFLESS:
                  *up = -u[save_op_s->second + t * diff2[j]] * r;
                  i += 3;
                  break;
                }
              j++;
            }
        }
      int t1 = max(1, periods - beg_t - y_kmax);
      int periods_beg_t = periods - beg_t;
      for (int t = t1; t < periods_beg_t; t++)
        {
          int i = j = 0;
          int gap = periods_beg_t - t;
          while (i < nop4)
            {
              if (auto* save_op_s = reinterpret_cast<t_save_op_s*>(&save_op[i]);
                  save_op_s->lag < gap)
                {
                  double* up = &u[save_op_s->first + t * diff1[j]];
                  switch (save_op_s->operat)
                    {
                    case IFLD:
                      r = *up;
                      i += 2;
                      break;
                    case IFDIV:
                      *up /= r;
                      i += 2;
                      break;
                    case IFSUB:
                      *up -= u[save_op_s->second + t * diff2[j]] * r;
                      i += 3;
                      break;
                    case IFLESS:
                      *up = -u[save_op_s->second + t * diff2[j]] * r;
                      i += 3;
                      break;
                    }
                }
              else
                switch (save_op_s->operat)
                  {
                  case IFLD:
                  case IFDIV:
                    i += 2;
                    break;
                  case IFSUB:
                  case IFLESS:
                    i += 3;
                    break;
                  }
              j++;
            }
        }
    }
  mxFree(diff1);
  mxFree(diff2);
  return OK;
}

int
Interpreter::complete(int beg_t)
{
  double yy = 0.0;

  int size_of_save_code = (1 + y_kmax) * size * (size + 1 + 4) / 2 * 4;
  int* save_code = static_cast<int*>(mxMalloc(size_of_save_code * sizeof(int)));
  test_mxMalloc(save_code, __LINE__, __FILE__, __func__, size_of_save_code * sizeof(int));
  int size_of_diff = (1 + y_kmax) * size * (size + 1 + 4);
  int* diff = static_cast<int*>(mxMalloc(size_of_diff * sizeof(int)));
  test_mxMalloc(diff, __LINE__, __FILE__, __func__, size_of_diff * sizeof(int));
  long cal_y = y_size * y_kmin;

  long i = (beg_t + 1) * size - 1;
  long nop = 0;
  for (long j = i; j > i - size; j--)
    {
      long pos = pivot[j];
      NonZeroElem* first;
      long nb_var;
      tie(nb_var, first) = At_Row(pos);
      first = first->NZE_R_N;
      nb_var--;
      save_code[nop] = IFLDZ;
      save_code[nop + 1] = 0;
      save_code[nop + 2] = 0;
      save_code[nop + 3] = 0;
#ifdef DEBUG
      if ((nop + 3) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop + 2, size_of_save_code);
#endif
      nop += 4;
      for (long k = 0; k < nb_var; k++)
        {
          save_code[nop] = IFMUL;
          save_code[nop + 1] = index_vara[first->c_index] + cal_y;
          save_code[nop + 2] = first->u_index;
          save_code[nop + 3] = first->lag_index;
#ifdef DEBUG
          if ((nop + 3) >= size_of_save_code)
            mexPrintf("out of save_code[%d] (bound=%d)\n", nop + 2, size_of_save_code);
#endif
          nop += 4;
          first = first->NZE_R_N;
        }
      save_code[nop] = IFADD;
      save_code[nop + 1] = b[pos];
      save_code[nop + 2] = 0;
      save_code[nop + 3] = 0;
#ifdef DEBUG
      if ((nop + 3) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop + 2, size_of_save_code);
#endif
      nop += 4;
      save_code[nop] = IFSTP;
      save_code[nop + 1] = index_vara[j] + y_size * y_kmin;
      save_code[nop + 2] = 0;
      save_code[nop + 3] = 0;
#ifdef DEBUG
      if ((nop + 2) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop + 2, size_of_save_code);
#endif
      nop += 4;
    }
  i = beg_t * size - 1;
  long nop1 = 0, nopa = 0;
  for (long j = i; j > i - size; j--)
    {
      long pos = pivot[j];
      NonZeroElem* first;
      long nb_var;
      tie(nb_var, first) = At_Row(pos);
      first = first->NZE_R_N;
      nb_var--;
      diff[nopa] = 0;
      diff[nopa + 1] = 0;
      nopa += 2;
      nop1 += 4;
      for (long k = 0; k < nb_var; k++)
        {
          diff[nopa] = save_code[nop1 + 1] - (index_vara[first->c_index] + cal_y);
          diff[nopa + 1] = save_code[nop1 + 2] - (first->u_index);
#ifdef DEBUG
          if ((nop1 + 2) >= size_of_save_code)
            mexPrintf("out of save_code[%d] (bound=%d)\n", nop1 + 2, size_of_save_code);
          if ((nopa + 1) >= size_of_diff)
            mexPrintf("out of diff[%d] (bound=%d)\n", nopa + 2, size_of_diff);
#endif
          nopa += 2;
          nop1 += 4;
          first = first->NZE_R_N;
        }
      diff[nopa] = save_code[nop1 + 1] - (b[pos]);
      diff[nopa + 1] = 0;
#ifdef DEBUG
      if ((nop1 + 3) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop1 + 2, size_of_save_code);
      if ((nopa + 1) >= size_of_diff)
        mexPrintf("out of diff[%d] (bound=%d)\n", nopa + 2, size_of_diff);
#endif
      nopa += 2;
      nop1 += 4;
      diff[nopa] = save_code[nop1 + 1] - (index_vara[j] + y_size * y_kmin);
      diff[nopa + 1] = 0;
#ifdef DEBUG
      if ((nop1 + 4) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop1 + 2, size_of_save_code);
      if ((nopa + 1) >= size_of_diff)
        mexPrintf("out of diff[%d] (bound=%d)\n", nopa + 2, size_of_diff);
#endif
      nopa += 2;
      nop1 += 4;
    }
  long max_var = (periods + y_kmin) * y_size;
  long min_var = y_kmin * y_size;
  for (int t = periods + y_kmin - 1; t >= beg_t + y_kmin; t--)
    {
      int j = 0, k;
      int ti = t - y_kmin - beg_t;
      for (int i = 0; i < nop; i += 4)
        {
          switch (save_code[i])
            {
            case IFLDZ:
              yy = 0;
              break;
            case IFMUL:
              k = save_code[i + 1] + ti * diff[j];
              if (k < max_var && k > min_var)
                yy += y[k] * u[save_code[i + 2] + ti * diff[j + 1]];
              break;
            case IFADD:
              yy = -(yy + u[save_code[i + 1] + ti * diff[j]]);
              break;
            case IFSTP:
              k = save_code[i + 1] + ti * diff[j];
              double err = yy - y[k];
              y[k] += slowc * (err);
              break;
            }
          j += 2;
        }
    }
  mxFree(save_code);
  mxFree(diff);
  return (beg_t);
}

void
Interpreter::bksub(int tbreak, int last_period)
{
  for (int i = 0; i < y_size * (periods + y_kmin); i++)
    y[i] = ya[i];
  if (symbolic && tbreak)
    last_period = complete(tbreak);
  else
    last_period = periods;
  for (int t = last_period + y_kmin - 1; t >= y_kmin; t--)
    {
      int ti = (t - y_kmin) * size;
      int cal = y_kmin * size;
      int cal_y = y_size * y_kmin;
      for (int i = ti - 1; i >= ti - size; i--)
        {
          int j = i + cal;
          int pos = pivot[i + size];
          auto [nb_var, first] = At_Row(pos);
          first = first->NZE_R_N;
          nb_var--;
          int eq = index_vara[j] + y_size;
          double yy = 0;
          for (int k = 0; k < nb_var; k++)
            {
              yy += y[index_vara[first->c_index] + cal_y] * u[first->u_index];
              first = first->NZE_R_N;
            }
          yy = -(yy + y[eq] + u[b[pos]]);
          direction[eq] = yy;
          y[eq] += slowc * yy;
        }
    }
}

void
Interpreter::simple_bksub()
{
  for (int i = 0; i < y_size; i++)
    y[i + it_ * y_size] = ya[i + it_ * y_size];
  for (int i = size - 1; i >= 0; i--)
    {
      int pos = pivot[i];
      auto [nb_var, first] = At_Row(pos);
      first = first->NZE_R_N;
      nb_var--;
      int eq = index_vara[i];
      double yy = 0;
      for (int k = 0; k < nb_var; k++)
        {
          yy += y[index_vara[first->c_index] + it_ * y_size] * u[first->u_index];
          first = first->NZE_R_N;
        }
      yy = -(yy + y[eq + it_ * y_size] + u[b[pos]]);
      direction[eq + it_ * y_size] = yy;
      y[eq + it_ * y_size] += slowc * yy;
    }
}

mxArray*
Interpreter::subtract_A_B(const mxArray* A_m, const mxArray* B_m)
{
  size_t n_A = mxGetN(A_m);
  size_t m_A = mxGetM(A_m);
  double* A_d = mxGetPr(A_m);
  size_t n_B = mxGetN(B_m);
  double* B_d = mxGetPr(B_m);
  mxArray* C_m = mxCreateDoubleMatrix(m_A, n_B, mxREAL);
  double* C_d = mxGetPr(C_m);
  for (int j = 0; j < static_cast<int>(n_A); j++)
    for (unsigned int i = 0; i < m_A; i++)
      {
        size_t index = j * m_A + i;
        C_d[index] = A_d[index] - B_d[index];
      }
  return C_m;
}

mxArray*
Interpreter::Sparse_subtract_SA_SB(const mxArray* A_m, const mxArray* B_m)
{
  size_t n_A = mxGetN(A_m);
  size_t m_A = mxGetM(A_m);
  mwIndex* A_i = mxGetIr(A_m);
  mwIndex* A_j = mxGetJc(A_m);
  size_t total_nze_A = A_j[n_A];
  double* A_d = mxGetPr(A_m);
  size_t n_B = mxGetN(B_m);
  mwIndex* B_i = mxGetIr(B_m);
  mwIndex* B_j = mxGetJc(B_m);
  size_t total_nze_B = B_j[n_B];
  double* B_d = mxGetPr(B_m);
  mxArray* C_m = mxCreateSparse(m_A, n_B, m_A * n_B, mxREAL);
  mwIndex* C_i = mxGetIr(C_m);
  mwIndex* C_j = mxGetJc(C_m);
  double* C_d = mxGetPr(C_m);
  unsigned int nze_B = 0, nze_C = 0, nze_A = 0;
  unsigned int A_col = 0, B_col = 0, C_col = 0;
  C_j[C_col] = 0;
  while (nze_A < total_nze_A || nze_B < total_nze_B)
    {
      while (nze_A >= static_cast<unsigned int>(A_j[A_col + 1]) && (nze_A < total_nze_A))
        A_col++;
      size_t A_row = A_i[nze_A];
      while (nze_B >= static_cast<unsigned int>(B_j[B_col + 1]) && (nze_B < total_nze_B))
        B_col++;
      size_t B_row = B_i[nze_B];
      if (A_col == B_col)
        {
          if (A_row == B_row && (nze_B < total_nze_B && nze_A < total_nze_A))
            {
              C_d[nze_C] = A_d[nze_A++] - B_d[nze_B++];
              C_i[nze_C] = A_row;
              while (C_col < A_col)
                C_j[++C_col] = nze_C;
              C_j[A_col + 1] = nze_C++;
              C_col = A_col;
            }
          else if ((A_row < B_row && nze_A < total_nze_A) || nze_B == total_nze_B)
            {
              C_d[nze_C] = A_d[nze_A++];
              C_i[nze_C] = A_row;
              while (C_col < A_col)
                C_j[++C_col] = nze_C;
              C_j[A_col + 1] = nze_C++;
              C_col = A_col;
            }
          else
            {
              C_d[nze_C] = -B_d[nze_B++];
              C_i[nze_C] = B_row;
              while (C_col < B_col)
                C_j[++C_col] = nze_C;
              C_j[B_col + 1] = nze_C++;
              C_col = B_col;
            }
        }
      else if ((A_col < B_col && nze_A < total_nze_A) || nze_B == total_nze_B)
        {
          C_d[nze_C] = A_d[nze_A++];
          C_i[nze_C] = A_row;
          while (C_col < A_col)
            C_j[++C_col] = nze_C;
          C_j[A_col + 1] = nze_C++;
          C_col = A_col;
        }
      else
        {
          C_d[nze_C] = -B_d[nze_B++];
          C_i[nze_C] = B_row;
          while (C_col < B_col)
            C_j[++C_col] = nze_C;
          C_j[B_col + 1] = nze_C++;
          C_col = B_col;
        }
    }
  while (C_col < n_B)
    C_j[++C_col] = nze_C;
  mxSetNzmax(C_m, nze_C);
  return C_m;
}

mxArray*
Interpreter::mult_SAT_B(const mxArray* A_m, const mxArray* B_m)
{
  size_t n_A = mxGetN(A_m);
  mwIndex* A_i = mxGetIr(A_m);
  mwIndex* A_j = mxGetJc(A_m);
  double* A_d = mxGetPr(A_m);
  size_t n_B = mxGetN(B_m);
  double* B_d = mxGetPr(B_m);
  mxArray* C_m = mxCreateDoubleMatrix(n_A, n_B, mxREAL);
  double* C_d = mxGetPr(C_m);
  for (int j = 0; j < static_cast<int>(n_B); j++)
    for (unsigned int i = 0; i < n_A; i++)
      {
        double sum = 0;
        size_t nze_A = A_j[i];
        while (nze_A < static_cast<unsigned int>(A_j[i + 1]))
          {
            size_t i_A = A_i[nze_A];
            sum += A_d[nze_A++] * B_d[i_A];
          }
        C_d[j * n_A + i] = sum;
      }
  return C_m;
}

mxArray*
Interpreter::Sparse_mult_SAT_B(const mxArray* A_m, const mxArray* B_m)
{
  size_t n_A = mxGetN(A_m);
  mwIndex* A_i = mxGetIr(A_m);
  mwIndex* A_j = mxGetJc(A_m);
  double* A_d = mxGetPr(A_m);
  size_t n_B = mxGetN(B_m);
  size_t m_B = mxGetM(B_m);
  double* B_d = mxGetPr(B_m);
  mxArray* C_m = mxCreateSparse(n_A, n_B, n_A * n_B, mxREAL);
  mwIndex* C_i = mxGetIr(C_m);
  mwIndex* C_j = mxGetJc(C_m);
  double* C_d = mxGetPr(C_m);
  unsigned int nze_C = 0;
  // unsigned int nze_A = 0;
  unsigned int C_col = 0;
  C_j[C_col] = 0;
  // #pragma omp parallel for
  for (unsigned int j = 0; j < n_B; j++)
    for (unsigned int i = 0; i < n_A; i++)
      {
        double sum = 0;
        size_t nze_A = A_j[i];
        while (nze_A < static_cast<unsigned int>(A_j[i + 1]))
          {
            size_t i_A = A_i[nze_A];
            sum += A_d[nze_A++] * B_d[i_A];
          }
        if (fabs(sum) > 1e-10)
          {
            C_d[nze_C] = sum;
            C_i[nze_C] = i;
            while (C_col < j)
              C_j[++C_col] = nze_C;
            nze_C++;
          }
      }
  while (C_col < m_B)
    C_j[++C_col] = nze_C;
  mxSetNzmax(C_m, nze_C);
  return C_m;
}

mxArray*
Interpreter::Sparse_mult_SAT_SB(const mxArray* A_m, const mxArray* B_m)
{
  size_t n_A = mxGetN(A_m);
  mwIndex* A_i = mxGetIr(A_m);
  mwIndex* A_j = mxGetJc(A_m);
  double* A_d = mxGetPr(A_m);
  size_t n_B = mxGetN(B_m);
  mwIndex* B_i = mxGetIr(B_m);
  mwIndex* B_j = mxGetJc(B_m);
  double* B_d = mxGetPr(B_m);
  mxArray* C_m = mxCreateSparse(n_A, n_B, n_A * n_B, mxREAL);
  mwIndex* C_i = mxGetIr(C_m);
  mwIndex* C_j = mxGetJc(C_m);
  double* C_d = mxGetPr(C_m);
  size_t nze_B = 0, nze_C = 0, nze_A = 0;
  unsigned int C_col = 0;
  C_j[C_col] = 0;
  for (unsigned int j = 0; j < n_B; j++)
    for (unsigned int i = 0; i < n_A; i++)
      {
        double sum = 0;
        nze_B = B_j[j];
        nze_A = A_j[i];
        while (nze_A < static_cast<unsigned int>(A_j[i + 1])
               && nze_B < static_cast<unsigned int>(B_j[j + 1]))
          {
            size_t i_A = A_i[nze_A];
            size_t i_B = B_i[nze_B];
            if (i_A == i_B)
              sum += A_d[nze_A++] * B_d[nze_B++];
            else if (i_A < i_B)
              nze_A++;
            else
              nze_B++;
          }
        if (fabs(sum) > 1e-10)
          {
            C_d[nze_C] = sum;
            C_i[nze_C] = i;
            while (C_col < j)
              C_j[++C_col] = nze_C;
            nze_C++;
          }
      }
  while (C_col < n_B)
    C_j[++C_col] = nze_C;
  mxSetNzmax(C_m, nze_C);
  return C_m;
}

mxArray*
Interpreter::Sparse_transpose(const mxArray* A_m)
{
  size_t n_A = mxGetN(A_m);
  size_t m_A = mxGetM(A_m);
  mwIndex* A_i = mxGetIr(A_m);
  mwIndex* A_j = mxGetJc(A_m);
  size_t total_nze_A = A_j[n_A];
  double* A_d = mxGetPr(A_m);
  mxArray* C_m = mxCreateSparse(n_A, m_A, total_nze_A, mxREAL);
  mwIndex* C_i = mxGetIr(C_m);
  mwIndex* C_j = mxGetJc(C_m);
  double* C_d = mxGetPr(C_m);
  unsigned int nze_C = 0, nze_A = 0;
  fill_n(C_j, m_A + 1, 0);
  map<pair<mwIndex, unsigned int>, double> B2;
  for (unsigned int i = 0; i < n_A; i++)
    while (nze_A < static_cast<unsigned int>(A_j[i + 1]))
      {
        C_j[A_i[nze_A] + 1]++;
        B2[{A_i[nze_A], i}] = A_d[nze_A];
        nze_A++;
      }
  for (unsigned int i = 0; i < m_A; i++)
    C_j[i + 1] += C_j[i];
  for (auto& [key, val] : B2)
    {
      C_d[nze_C] = val;
      C_i[nze_C++] = key.second;
    }
  return C_m;
}

void
Interpreter::compute_block_time(int my_Per_u_, bool evaluate, bool no_derivatives)
{
#ifdef DEBUG
  mexPrintf("compute_block_time\n");
#endif
  double *jacob {nullptr}, *jacob_exo {nullptr}, *jacob_exo_det {nullptr};
  if (evaluate)
    {
      jacob = mxGetPr(jacobian_block[block_num]);
      if (!steady_state)
        {
          jacob_exo = mxGetPr(jacobian_exo_block[block_num]);
          jacob_exo_det = mxGetPr(jacobian_det_exo_block[block_num]);
        }
    }

  try
    {
      evaluator.evaluateBlock(it_, y_kmin, y, y_size, x, nb_row_x, params, steady_y, g1, u,
                              my_Per_u_, T, periods, TEF, TEFD, TEFDD, r, jacob, jacob_exo,
                              jacob_exo_det, evaluate, no_derivatives);
    }
  catch (FloatingPointException& e)
    {
      res1 = numeric_limits<double>::quiet_NaN();
      if (verbosity >= 2)
        mexPrintf("%s\n      %s\n", e.message.c_str(), e.location.c_str());
    }
}

bool
Interpreter::compute_complete(bool no_derivatives)
{
  bool result;
  res1 = 0;
  compute_block_time(0, false, no_derivatives);
  if (!(isnan(res1) || isinf(res1)))
    {
      res1 = 0;
      res2 = 0;
      max_res = 0;
      for (int i = 0; i < size; i++)
        {
          double rr;
          rr = r[i];
          if (max_res < fabs(rr))
            {
              max_res = fabs(rr);
              max_res_idx = i;
            }
          res2 += rr * rr;
          res1 += fabs(rr);
        }
      result = true;
    }
  else
    result = false;
  return result;
}

pair<bool, double>
Interpreter::compute_complete(double lambda)
{
  double res2_ = 0, max_res_ = 0;
  int max_res_idx_ = 0;
  if (steady_state)
    {
      it_ = 0;
      for (int i = 0; i < size; i++)
        {
          int eq = index_vara[i];
          y[eq] = ya[eq] + lambda * direction[eq];
        }
      Per_u_ = 0;
      Per_y_ = 0;
      if (compute_complete(true))
        res2_ = res2;
      else
        return {false, numeric_limits<double>::quiet_NaN()};
    }
  else
    {
      for (int it = y_kmin; it < periods + y_kmin; it++)
        for (int i = 0; i < size; i++)
          {
            int eq = index_vara[i];
            y[eq + it * y_size] = ya[eq + it * y_size] + lambda * direction[eq + it * y_size];
          }
      for (it_ = y_kmin; it_ < periods + y_kmin; it_++)
        {
          Per_u_ = (it_ - y_kmin) * u_count_int;
          Per_y_ = it_ * y_size;
          if (compute_complete(true))
            {
              res2_ += res2;
              if (max_res > max_res_)
                {
                  max_res = max_res_;
                  max_res_idx = max_res_idx_;
                }
            }
          else
            return {false, numeric_limits<double>::quiet_NaN()};
        }
      it_ = periods + y_kmin - 1; // Do not leave it_ in inconsistent state
    }
  if (verbosity >= 2)
    mexPrintf("  lambda=%e, res2=%e\n", lambda, res2_);
  double crit {res2_ / 2};
  return {true, crit};
}

tuple<bool, double, double, double, double>
Interpreter::mnbrak(double& ax, double& bx)
{
  constexpr double GOLD = 1.618034;
  constexpr double GLIMIT = 100.0;
  constexpr double TINY = 1.0e-20;

  constexpr tuple failval
      = {false, numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN(),
         numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};

  auto sign = [](double a, double b) { return b >= 0.0 ? fabs(a) : -fabs(a); };

  if (verbosity >= 2)
    mexPrintf("bracketing ax=%f, bx=%f\n", ax, bx);

  auto [success, fa] = compute_complete(ax);
  if (!success)
    return failval;

  auto [success2, fb] = compute_complete(bx);
  if (!success2)
    return failval;

  if (fb > fa)
    {
      swap(ax, bx);
      swap(fa, fb);
    }

  double cx = bx + GOLD * (bx - ax);
  auto [success3, fc] = compute_complete(cx);
  if (!success3)
    return failval;

  while (fb > fc)
    {
      double r = (bx - ax) * (fb - fc);
      double q = (bx - cx) * (fb - fa);
      double u
          = bx - ((bx - cx) * q - (bx - ax) * r) / (2.0 * sign(fmax(fabs(q - r), TINY), q - r));
      double ulim = bx + GLIMIT * (cx - bx);
      double fu;
      if ((bx - u) * (u - cx) > 0.0)
        {
          tie(success, fu) = compute_complete(u);
          if (!success)
            return failval;
          if (fu < fc)
            {
              ax = bx;
              bx = u;
              fa = fb;
              fb = fu;
              goto success;
            }
          else if (fu > fb)
            {
              cx = u;
              fc = fu;
              goto success;
            }
          u = cx + GOLD * (cx - bx);
          tie(success, fu) = compute_complete(u);
          if (!success)
            return failval;
        }
      else if ((cx - u) * (u - ulim) > 0.0)
        {
          tie(success, fu) = compute_complete(u);
          if (!success)
            return failval;
          if (fu < fc)
            {
              bx = cx;
              cx = u;
              u = cx + GOLD * (cx - bx);
              fb = fc;
              fc = fu;
              tie(success, fu) = compute_complete(u);
              if (!success)
                return failval;
            }
        }
      else if ((u - ulim) * (ulim - cx) >= 0.0)
        {
          u = ulim;
          tie(success, fu) = compute_complete(u);
          if (!success)
            return failval;
        }
      else
        {
          u = cx + GOLD * (cx - bx);
          tie(success, fu) = compute_complete(u);
          if (!success)
            return failval;
        }
      ax = bx;
      bx = cx;
      cx = u;
      fa = fb;
      fb = fc;
      fc = fu;
    }

success:
  return {true, cx, fa, fb, fc};
}

pair<bool, double>
Interpreter::golden(double ax, double bx, double cx, double tol)
{
  constexpr pair failval = {false, numeric_limits<double>::quiet_NaN()};
  const double R = 0.61803399;
  const double C = (1.0 - R);
  if (verbosity >= 2)
    mexPrintf("golden\n");
  int iter = 0, max_iter = 100;
  double x1, x2;
  double x0 = ax;
  double x3 = cx;
  if (fabs(cx - bx) > fabs(bx - ax))
    {
      x1 = bx;
      x2 = bx + C * (cx - bx);
    }
  else
    {
      x2 = bx;
      x1 = bx - C * (bx - ax);
    }
  auto [success, f1] = compute_complete(x1);
  if (!success)
    return failval;
  auto [success2, f2] = compute_complete(x2);
  if (!success2)
    return failval;
  while (fabs(x3 - x0) > tol * (fabs(x1) + fabs(x2)) && f1 > solve_tolf && f2 > solve_tolf
         && iter < max_iter && abs(x1 - x2) > 1e-4)
    {
      if (f2 < f1)
        {
          x0 = x1;
          x1 = x2;
          x2 = R * x1 + C * x3;
          f1 = f2;
          tie(success, f2) = compute_complete(x2);
          if (!success)
            return failval;
        }
      else
        {
          x3 = x2;
          x2 = x1;
          x1 = R * x2 + C * x0;
          f2 = f1;
          tie(success, f1) = compute_complete(x1);
          if (!success)
            return failval;
        }
      iter++;
    }
  double xmin {f1 < f2 ? x1 : x2};
  return {true, xmin};
}

void
Interpreter::End_Solver()
{
  if (((stack_solve_algo == 0 || stack_solve_algo == 4) && !steady_state)
      || (solve_algo == 6 && steady_state))
    {
      if (Symbolic)
        {
          umfpack_dl_free_symbolic(&Symbolic);
          Symbolic = nullptr;
        }
      if (Numeric)
        {
          umfpack_dl_free_numeric(&Numeric);
          Numeric = nullptr;
        }
    }
}

void
Interpreter::Solve_LU_UMFPack_Two_Boundaries(
    SuiteSparse_long* Ap, SuiteSparse_long* Ai, double* Ax, double* b,
    const vector_table_conditional_local_type& vector_table_conditional_local)
{
  int n {size * periods};
  SuiteSparse_long sys = 0;
  std::array<double, UMFPACK_CONTROL> Control;
  std::array<double, UMFPACK_INFO> Info;
  std::vector<double> res(n);

  umfpack_dl_defaults(Control.data());
  SuiteSparse_long status = 0;
  if (iter == 0)
    {
      if (Symbolic)
        umfpack_dl_free_symbolic(&Symbolic);
      status = umfpack_dl_symbolic(n, n, Ap, Ai, Ax, &Symbolic, Control.data(), Info.data());
      if (status != UMFPACK_OK)
        {
          umfpack_dl_report_info(Control.data(), Info.data());
          umfpack_dl_report_status(Control.data(), status);
          throw FatalException {"umfpack_dl_symbolic failed"};
        }
    }
  if (Numeric)
    umfpack_dl_free_numeric(&Numeric);
  status = umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control.data(), Info.data());
  if (status != UMFPACK_OK)
    {
      umfpack_dl_report_info(Control.data(), Info.data());
      umfpack_dl_report_status(Control.data(), status);
      throw FatalException {"umfpack_dl_numeric failed"};
    }
  status = umfpack_dl_solve(sys, Ap, Ai, Ax, res.data(), b, Numeric, Control.data(), Info.data());
  if (status != UMFPACK_OK)
    {
      umfpack_dl_report_info(Control.data(), Info.data());
      umfpack_dl_report_status(Control.data(), status);
      throw FatalException {"umfpack_dl_solve failed"};
    }

  if (vector_table_conditional_local.size())
    {
      for (int t = 0; t < periods; t++)
        if (t == 0)
          {
            for (int i = 0; i < size; i++)
              {
                bool fliped = vector_table_conditional_local[i].is_cond;
                if (fliped)
                  {
                    int eq = index_vara[i + size * (y_kmin)];
                    int flip_exo = vector_table_conditional_local[i].var_exo;
                    double yy = -(res[i] + x[y_kmin + flip_exo * nb_row_x]);
                    direction[eq] = 0;
                    x[flip_exo * nb_row_x + y_kmin] += slowc * yy;
                  }
                else
                  {
                    int eq = index_vara[i + size * (y_kmin)];
                    double yy = -(res[i] + y[eq]);
                    direction[eq] = yy;
                    y[eq] += slowc * yy;
                  }
              }
          }
        else
          {
            for (int i = 0; i < size; i++)
              {
                int eq = index_vara[i + size * (t + y_kmin)];
                double yy = -(res[i + size * t] + y[eq]);
                direction[eq] = yy;
                y[eq] += slowc * yy;
              }
          }
    }
  else
    {
      for (int i = 0; i < n; i++)
        {
          int eq = index_vara[i + size * y_kmin];
          double yy = -(res[i] + y[eq]);
          direction[eq] = yy;
          y[eq] += slowc * yy;
        }
    }

  mxFree(Ap);
  mxFree(Ai);
  mxFree(Ax);
  mxFree(b);
}

void
Interpreter::Solve_LU_UMFPack_One_Boundary(SuiteSparse_long* Ap, SuiteSparse_long* Ai, double* Ax,
                                           double* b)
{
  SuiteSparse_long sys = 0;
  std::array<double, UMFPACK_CONTROL> Control;
  std::array<double, UMFPACK_INFO> Info;
  std::vector<double> res(size);

  umfpack_dl_defaults(Control.data());
  SuiteSparse_long status = 0;
  if (iter == 0)
    {
      if (Symbolic)
        umfpack_dl_free_symbolic(&Symbolic);
      status = umfpack_dl_symbolic(size, size, Ap, Ai, Ax, &Symbolic, Control.data(), Info.data());
      if (status != UMFPACK_OK)
        {
          umfpack_dl_report_info(Control.data(), Info.data());
          umfpack_dl_report_status(Control.data(), status);
          throw FatalException {"umfpack_dl_symbolic failed"};
        }
    }
  if (Numeric)
    umfpack_dl_free_numeric(&Numeric);
  status = umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control.data(), Info.data());
  if (status != UMFPACK_OK)
    {
      umfpack_dl_report_info(Control.data(), Info.data());
      umfpack_dl_report_status(Control.data(), status);
      throw FatalException {"umfpack_dl_numeric failed"};
    }
  status = umfpack_dl_solve(sys, Ap, Ai, Ax, res.data(), b, Numeric, Control.data(), Info.data());
  if (status != UMFPACK_OK)
    {
      umfpack_dl_report_info(Control.data(), Info.data());
      umfpack_dl_report_status(Control.data(), status);
      throw FatalException {"umfpack_dl_solve failed"};
    }

  for (int i = 0; i < size; i++)
    {
      int eq = index_vara[i];
      double yy = -(res[i] + y[eq + it_ * y_size]);
      direction[eq] = yy;
      y[eq + it_ * y_size] += slowc * yy;
    }
  mxFree(Ap);
  mxFree(Ai);
  mxFree(Ax);
  mxFree(b);
}

void
Interpreter::Solve_Matlab_GMRES(mxArray* A_m, mxArray* b_m, bool is_two_boundaries, mxArray* x0_m)
{
  size_t n = mxGetM(A_m);
  std::array field_names {"droptol", "type"};
  std::array dims {static_cast<mwSize>(1)};
  mxArray* Setup
      = mxCreateStructArray(dims.size(), dims.data(), field_names.size(), field_names.data());
  mxSetFieldByNumber(Setup, 0, 0, mxCreateDoubleScalar(lu_inc_tol));
  mxSetFieldByNumber(Setup, 0, 1, mxCreateString("ilutp"));
  std::array<mxArray*, 2> lhs0;
  std::array rhs0 {A_m, Setup};
  if (mexCallMATLAB(lhs0.size(), lhs0.data(), rhs0.size(), rhs0.data(), "ilu"))
    throw FatalException("In GMRES, the incomplete LU decomposition (ilu) has failed");
  mxArray* L1 = lhs0[0];
  mxArray* U1 = lhs0[1];
  /*[za,flag1] = gmres(g1a,b,Blck_size,1e-6,Blck_size*periods,L1,U1);*/
  std::array rhs {A_m,
                  b_m,
                  mxCreateDoubleScalar(size),
                  mxCreateDoubleScalar(1e-6),
                  mxCreateDoubleScalar(static_cast<double>(n)),
                  L1,
                  U1,
                  x0_m};
  std::array<mxArray*, 2> lhs;
  mexCallMATLAB(lhs.size(), lhs.data(), rhs.size(), rhs.data(), "gmres");
  mxArray* z = lhs[0];
  mxArray* flag = lhs[1];
  double flag1 {mxGetScalar(flag)};
  mxDestroyArray(rhs0[1]);
  mxDestroyArray(rhs[2]);
  mxDestroyArray(rhs[3]);
  mxDestroyArray(rhs[4]);
  mxDestroyArray(rhs[5]);
  mxDestroyArray(rhs[6]);
  if (flag1 > 0)
    {
      if (flag1 == 1)
        mexWarnMsgTxt(
            ("Error in bytecode: No convergence inside GMRES, in block " + to_string(block_num + 1))
                .c_str());
      else if (flag1 == 2)
        mexWarnMsgTxt(("Error in bytecode: Preconditioner is ill-conditioned, in block "
                       + to_string(block_num + 1))
                          .c_str());
      else if (flag1 == 3)
        mexWarnMsgTxt(("Error in bytecode: GMRES stagnated (Two consecutive iterates were the "
                       "same.), in block "
                       + to_string(block_num + 1))
                          .c_str());
      lu_inc_tol /= 10;
    }
  else
    {
      double* res = mxGetPr(z);
      if (is_two_boundaries)
        for (int i = 0; i < static_cast<int>(n); i++)
          {
            int eq = index_vara[i + size * y_kmin];
            double yy = -(res[i] + y[eq]);
            direction[eq] = yy;
            y[eq] += slowc * yy;
          }
      else
        for (int i = 0; i < static_cast<int>(n); i++)
          {
            int eq = index_vara[i];
            double yy = -(res[i] + y[eq + it_ * y_size]);
            direction[eq] = yy;
            y[eq + it_ * y_size] += slowc * yy;
          }
    }
  mxDestroyArray(A_m);
  mxDestroyArray(b_m);
  mxDestroyArray(x0_m);
  mxDestroyArray(z);
  mxDestroyArray(flag);
}

void
Interpreter::Solve_Matlab_BiCGStab(mxArray* A_m, mxArray* b_m, bool is_two_boundaries,
                                   mxArray* x0_m, int preconditioner)
{
  /* precond = 0  => Jacobi
     precond = 1  => Incomplet LU decomposition*/
  size_t n = mxGetM(A_m);
  mxArray *L1 = nullptr, *U1 = nullptr, *Diag = nullptr;

  if (preconditioner == 0)
    {
      std::array<mxArray*, 1> lhs0;
      std::array rhs0 {A_m, mxCreateDoubleScalar(0)};
      mexCallMATLAB(lhs0.size(), lhs0.data(), rhs0.size(), rhs0.data(), "spdiags");
      mxArray* tmp = lhs0[0];
      double* tmp_val = mxGetPr(tmp);
      Diag = mxCreateSparse(n, n, n, mxREAL);
      mwIndex* Diag_i = mxGetIr(Diag);
      mwIndex* Diag_j = mxGetJc(Diag);
      double* Diag_val = mxGetPr(Diag);
      for (size_t i = 0; i < n; i++)
        {
          Diag_val[i] = tmp_val[i];
          Diag_j[i] = i;
          Diag_i[i] = i;
        }
      Diag_j[n] = n;
    }
  else if (preconditioner == 1)
    {
      /*[L1, U1] = ilu(g1a=;*/
      std::array field_names {"type", "droptol", "milu", "udiag", "thresh"};
      const int type = 0, droptol = 1, milu = 2, udiag = 3, thresh = 4;
      std::array dims {static_cast<mwSize>(1)};
      mxArray* Setup
          = mxCreateStructArray(dims.size(), dims.data(), field_names.size(), field_names.data());
      mxSetFieldByNumber(Setup, 0, type, mxCreateString("ilutp"));
      mxSetFieldByNumber(Setup, 0, droptol, mxCreateDoubleScalar(lu_inc_tol));
      mxSetFieldByNumber(Setup, 0, milu, mxCreateString("off"));
      mxSetFieldByNumber(Setup, 0, udiag, mxCreateDoubleScalar(0));
      mxSetFieldByNumber(Setup, 0, thresh, mxCreateDoubleScalar(1));
      std::array<mxArray*, 2> lhs0;
      std::array rhs0 {A_m, Setup};
      if (mexCallMATLAB(lhs0.size(), lhs0.data(), rhs0.size(), rhs0.data(), "ilu"))
        throw FatalException {"In BiCGStab, the incomplete LU decomposition (ilu) has failed"};
      L1 = lhs0[0];
      U1 = lhs0[1];
      mxDestroyArray(Setup);
    }
  double flags = 2;
  mxArray* z = nullptr;
  if (steady_state) /*Octave BicStab algorihtm involves a 0 division in case of a preconditionner
                       equal to the LU decomposition of A matrix*/
    {
      mxArray* res = mult_SAT_B(Sparse_transpose(A_m), x0_m);
      double* resid = mxGetPr(res);
      double* b = mxGetPr(b_m);
      for (int i = 0; i < static_cast<int>(n); i++)
        resid[i] = b[i] - resid[i];
      std::array<mxArray*, 1> lhs;
      std::array rhs {L1, res};
      mexCallMATLAB(lhs.size(), lhs.data(), rhs.size(), rhs.data(), "mldivide");
      std::array rhs2 {U1, lhs[0]};
      mexCallMATLAB(lhs.size(), lhs.data(), rhs2.size(), rhs2.data(), "mldivide");
      z = lhs[0];
      double* phat = mxGetPr(z);
      double* x0 = mxGetPr(x0_m);
      for (int i = 0; i < static_cast<int>(n); i++)
        phat[i] = x0[i] + phat[i];

      /*Check the solution*/
      res = mult_SAT_B(Sparse_transpose(A_m), z);
      resid = mxGetPr(res);
      double cum_abs = 0;
      for (int i = 0; i < static_cast<int>(n); i++)
        {
          resid[i] = b[i] - resid[i];
          cum_abs += fabs(resid[i]);
        }
      if (cum_abs > 1e-7)
        flags = 2;
      else
        flags = 0;
      mxDestroyArray(res);
    }

  if (flags == 2)
    {
      if (preconditioner == 0)
        {
          /*[za,flag1] = bicgstab(g1a,b,1e-6,Blck_size*periods,L1,U1);*/
          std::array rhs {A_m, b_m, mxCreateDoubleScalar(1e-6),
                          mxCreateDoubleScalar(static_cast<double>(n)), Diag};
          std::array<mxArray*, 2> lhs;
          mexCallMATLAB(lhs.size(), lhs.data(), rhs.size(), rhs.data(), "bicgstab");
          z = lhs[0];
          mxArray* flag = lhs[1];
          flags = mxGetScalar(flag);
          mxDestroyArray(flag);
          mxDestroyArray(rhs[2]);
          mxDestroyArray(rhs[3]);
          mxDestroyArray(rhs[4]);
        }
      else if (preconditioner == 1)
        {
          /*[za,flag1] = bicgstab(g1a,b,1e-6,Blck_size*periods,L1,U1);*/
          std::array rhs {A_m,
                          b_m,
                          mxCreateDoubleScalar(1e-6),
                          mxCreateDoubleScalar(static_cast<double>(n)),
                          L1,
                          U1,
                          x0_m};
          std::array<mxArray*, 2> lhs;
          mexCallMATLAB(lhs.size(), lhs.data(), rhs.size(), rhs.data(), "bicgstab");
          z = lhs[0];
          mxArray* flag = lhs[1];
          flags = mxGetScalar(flag);
          mxDestroyArray(flag);
          mxDestroyArray(rhs[2]);
          mxDestroyArray(rhs[3]);
          mxDestroyArray(rhs[4]);
          mxDestroyArray(rhs[5]);
        }
    }

  if (flags > 0)
    {
      if (flags == 1)
        mexWarnMsgTxt(("Error in bytecode: No convergence inside BiCGStab, in block "
                       + to_string(block_num + 1))
                          .c_str());
      else if (flags == 2)
        mexWarnMsgTxt(("Error in bytecode: Preconditioner is ill-conditioned, in block "
                       + to_string(block_num + 1))
                          .c_str());
      else if (flags == 3)
        mexWarnMsgTxt(("Error in bytecode: BiCGStab stagnated (Two consecutive iterates were the "
                       "same.), in block "
                       + to_string(block_num + 1))
                          .c_str());
      lu_inc_tol /= 10;
    }
  else
    {
      double* res = mxGetPr(z);
      if (is_two_boundaries)
        for (int i = 0; i < static_cast<int>(n); i++)
          {
            int eq = index_vara[i + size * y_kmin];
            double yy = -(res[i] + y[eq]);
            direction[eq] = yy;
            y[eq] += slowc * yy;
          }
      else
        for (int i = 0; i < static_cast<int>(n); i++)
          {
            int eq = index_vara[i];
            double yy = -(res[i] + y[eq + it_ * y_size]);
            direction[eq] = yy;
            y[eq + it_ * y_size] += slowc * yy;
          }
    }
  mxDestroyArray(A_m);
  mxDestroyArray(b_m);
  mxDestroyArray(x0_m);
  mxDestroyArray(z);
}

void
Interpreter::Singular_display()
{
  Simple_Init();
  std::array rhs {mxCreateDoubleMatrix(size, size, mxREAL)};
  double* pind = mxGetPr(rhs[0]);
  for (int j = 0; j < size * size; j++)
    pind[j] = 0.0;
  for (int ii = 0; ii < size; ii++)
    {
      auto [nb_eq, first] = At_Col(ii);
      for (int j = 0; j < nb_eq; j++)
        {
          int k = first->u_index;
          int jj = first->r_index;
          pind[ii * size + jj] = u[k];
          first = first->NZE_C_N;
        }
    }
  std::array<mxArray*, 3> lhs;
  mexCallMATLAB(lhs.size(), lhs.data(), rhs.size(), rhs.data(), "svd");
  mxArray* SVD_u = lhs[0];
  mxArray* SVD_s = lhs[1];
  double* SVD_ps = mxGetPr(SVD_s);
  double* SVD_pu = mxGetPr(SVD_u);
  for (int i = 0; i < size; i++)
    if (abs(SVD_ps[i * (1 + size)]) < 1e-12)
      {
        mexPrintf(" The following equations form a linear combination:\n    ");
        double max_u = 0;
        for (int j = 0; j < size; j++)
          if (abs(SVD_pu[j + i * size]) > abs(max_u))
            max_u = SVD_pu[j + i * size];
        vector<int> equ_list;
        for (int j = 0; j < size; j++)
          {
            double rr = SVD_pu[j + i * size] / max_u;
            if (rr < -1e-10)
              {
                equ_list.push_back(j);
                if (rr != -1)
                  mexPrintf(" - %3.2f*Dequ_%d_dy", abs(rr), j + 1);
                else
                  mexPrintf(" - Dequ_%d_dy", j + 1);
              }
            else if (rr > 1e-10)
              {
                equ_list.push_back(j);
                if (j > 0)
                  if (rr != 1)
                    mexPrintf(" + %3.2f*Dequ_%d_dy", rr, j + 1);
                  else
                    mexPrintf(" + Dequ_%d_dy", j + 1);
                else if (rr != 1)
                  mexPrintf(" %3.2f*Dequ_%d_dy", rr, j + 1);
                else
                  mexPrintf(" Dequ_%d_dy", j + 1);
              }
          }
        mexPrintf(" = 0\n");
      }
  mxDestroyArray(lhs[0]);
  mxDestroyArray(lhs[1]);
  mxDestroyArray(lhs[2]);
  if (block_num > 1)
    throw FatalException {"In Solve_ByteCode_Sparse_GaussianElimination, singular system in block "
                          + to_string(block_num + 1)};
  else
    throw FatalException {"In Solve_ByteCode_Sparse_GaussianElimination, singular system"};
}

bool
Interpreter::Solve_ByteCode_Sparse_GaussianElimination()
{
  int pivj = 0, pivk = 0;
  auto** bc = static_cast<NonZeroElem**>(mxMalloc(size * sizeof(NonZeroElem*)));
  test_mxMalloc(bc, __LINE__, __FILE__, __func__, size * sizeof(*bc));
  auto* piv_v = static_cast<double*>(mxMalloc(size * sizeof(double)));
  test_mxMalloc(piv_v, __LINE__, __FILE__, __func__, size * sizeof(double));
  int* pivj_v = static_cast<int*>(mxMalloc(size * sizeof(int)));
  test_mxMalloc(pivj_v, __LINE__, __FILE__, __func__, size * sizeof(int));
  int* pivk_v = static_cast<int*>(mxMalloc(size * sizeof(int)));
  test_mxMalloc(pivk_v, __LINE__, __FILE__, __func__, size * sizeof(int));
  int* NR = static_cast<int*>(mxMalloc(size * sizeof(int)));
  test_mxMalloc(NR, __LINE__, __FILE__, __func__, size * sizeof(int));

  for (int i = 0; i < size; i++)
    {
      /*finding the max-pivot*/
      double piv = 0, piv_abs = 0;
      auto [nb_eq, first] = At_Col(i);
      int l = 0;
      int N_max = 0;
      bool one = false;
      for (int j = 0; j < nb_eq; j++)
        {
          if (!line_done[first->r_index])
            {
              int k = first->u_index;
              int jj = first->r_index;
              int NRow_jj = NRow(jj);

              piv_v[l] = u[k];
              double piv_fabs = fabs(u[k]);
              pivj_v[l] = jj;
              pivk_v[l] = k;
              NR[l] = NRow_jj;
              if (NRow_jj == 1 && !one)
                {
                  one = true;
                  piv_abs = piv_fabs;
                  N_max = NRow_jj;
                }
              if (!one)
                {
                  if (piv_fabs > piv_abs)
                    piv_abs = piv_fabs;
                  if (NRow_jj > N_max)
                    N_max = NRow_jj;
                }
              else if (NRow_jj == 1)
                {
                  if (piv_fabs > piv_abs)
                    piv_abs = piv_fabs;
                  if (NRow_jj > N_max)
                    N_max = NRow_jj;
                }
              l++;
            }
          first = first->NZE_C_N;
        }
      if (piv_abs < eps)
        {
          mxFree(piv_v);
          mxFree(pivj_v);
          mxFree(pivk_v);
          mxFree(NR);
          mxFree(bc);
          if (steady_state)
            {
              if (verbosity >= 1)
                {
                  if (block_num > 1)
                    mexPrintf("Error: singular system in Simulate_NG in block %d\n", block_num + 1);
                  else
                    mexPrintf("Error: singular system in Simulate_NG");
                }
              return true;
            }
          else
            {
              if (block_num > 1)
                throw FatalException {
                    "In Solve_ByteCode_Sparse_GaussianElimination, singular system in block "
                    + to_string(block_num + 1)};
              else
                throw FatalException {
                    "In Solve_ByteCode_Sparse_GaussianElimination, singular system"};
            }
        }
      if (!one)
        {
          double markovitz = 0, markovitz_max = -9e70;
          for (int j = 0; j < l; j++)
            {
              if (N_max > 0 && NR[j] > 0)
                {
                  if (fabs(piv_v[j]) > 0)
                    {
                      if (markowitz_c > 0)
                        markovitz = exp(
                            log(fabs(piv_v[j]) / piv_abs)
                            - markowitz_c
                                  * log(static_cast<double>(NR[j]) / static_cast<double>(N_max)));
                      else
                        markovitz = fabs(piv_v[j]) / piv_abs;
                    }
                  else
                    markovitz = 0;
                }
              else
                markovitz = fabs(piv_v[j]) / piv_abs;
              if (markovitz > markovitz_max)
                {
                  piv = piv_v[j];
                  pivj = pivj_v[j]; // Line number
                  pivk = pivk_v[j]; // positi
                  markovitz_max = markovitz;
                }
            }
        }
      else
        for (int j = 0; j < l; j++)
          {
            if (NR[j] == 1)
              {
                piv = piv_v[j];
                pivj = pivj_v[j]; // Line number
                pivk = pivk_v[j]; // positi
              }
          }
      pivot[i] = pivj;
      pivotk[i] = pivk;
      pivotv[i] = piv;
      line_done[pivj] = true;

      /*divide all the non zeros elements of the line pivj by the max_pivot*/
      int nb_var;
      tie(nb_var, first) = At_Row(pivj);
      for (int j = 0; j < nb_var; j++)
        {
          u[first->u_index] /= piv;
          first = first->NZE_R_N;
        }
      u[b[pivj]] /= piv;
      /*subtract the elements on the non treated lines*/
      tie(nb_eq, first) = At_Col(i);
      auto [nb_var_piva, first_piva] = At_Row(pivj);
      int nb_eq_todo = 0;
      for (int j = 0; j < nb_eq && first; j++)
        {
          if (!line_done[first->r_index])
            bc[nb_eq_todo++] = first;
          first = first->NZE_C_N;
        }
      for (int j = 0; j < nb_eq_todo; j++)
        {
          first = bc[j];
          int row = first->r_index;
          double first_elem = u[first->u_index];

          int nb_var_piv = nb_var_piva;
          NonZeroElem* first_piv = first_piva;
          auto [nb_var_sub, first_sub] = At_Row(row);
          int l_sub = 0, l_piv = 0;
          int sub_c_index = first_sub->c_index, piv_c_index = first_piv->c_index;
          while (l_sub < nb_var_sub || l_piv < nb_var_piv)
            if (l_sub < nb_var_sub && (sub_c_index < piv_c_index || l_piv >= nb_var_piv))
              {
                first_sub = first_sub->NZE_R_N;
                if (first_sub)
                  sub_c_index = first_sub->c_index;
                else
                  sub_c_index = size;
                l_sub++;
              }
            else if (sub_c_index > piv_c_index || l_sub >= nb_var_sub)
              {
                int tmp_u_count = Get_u();
                Insert(row, first_piv->c_index, tmp_u_count, 0);
                u[tmp_u_count] = -u[first_piv->u_index] * first_elem;
                first_piv = first_piv->NZE_R_N;
                if (first_piv)
                  piv_c_index = first_piv->c_index;
                else
                  piv_c_index = size;
                l_piv++;
              }
            else
              {
                if (i == sub_c_index)
                  {
                    NonZeroElem* firsta = first;
                    NonZeroElem* first_suba = first_sub->NZE_R_N;
                    Delete(first_sub->r_index, first_sub->c_index);
                    first = firsta->NZE_C_N;
                    first_sub = first_suba;
                    if (first_sub)
                      sub_c_index = first_sub->c_index;
                    else
                      sub_c_index = size;
                    l_sub++;
                    first_piv = first_piv->NZE_R_N;
                    if (first_piv)
                      piv_c_index = first_piv->c_index;
                    else
                      piv_c_index = size;
                    l_piv++;
                  }
                else
                  {
                    u[first_sub->u_index] -= u[first_piv->u_index] * first_elem;
                    first_sub = first_sub->NZE_R_N;
                    if (first_sub)
                      sub_c_index = first_sub->c_index;
                    else
                      sub_c_index = size;
                    l_sub++;
                    first_piv = first_piv->NZE_R_N;
                    if (first_piv)
                      piv_c_index = first_piv->c_index;
                    else
                      piv_c_index = size;
                    l_piv++;
                  }
              }
          u[b[row]] -= u[b[pivj]] * first_elem;
        }
    }
  for (int i = 0; i < y_size; i++)
    ya[i + it_ * y_size] = y[i + it_ * y_size];

  slowc_save = slowc;
  simple_bksub();
  End_Gaussian_Elimination();
  mxFree(piv_v);
  mxFree(pivj_v);
  mxFree(pivk_v);
  mxFree(NR);
  mxFree(bc);
  return false;
}

void
Interpreter::Solve_ByteCode_Symbolic_Sparse_GaussianElimination(bool symbolic)
{
  /*Triangularisation at each period of a block using a simple gaussian Elimination*/
  int *save_op = nullptr, *save_opa = nullptr, *save_opaa = nullptr;
  long int nop = 0, nopa = 0;
  bool record = false;

  int pivj = 0, pivk = 0;
  int tbreak = 0, last_period = periods;

  auto* piv_v = static_cast<double*>(mxMalloc(size * sizeof(double)));
  test_mxMalloc(piv_v, __LINE__, __FILE__, __func__, size * sizeof(double));
  int* pivj_v = static_cast<int*>(mxMalloc(size * sizeof(int)));
  test_mxMalloc(pivj_v, __LINE__, __FILE__, __func__, size * sizeof(int));
  int* pivk_v = static_cast<int*>(mxMalloc(size * sizeof(int)));
  test_mxMalloc(pivk_v, __LINE__, __FILE__, __func__, size * sizeof(int));
  int* NR = static_cast<int*>(mxMalloc(size * sizeof(int)));
  test_mxMalloc(NR, __LINE__, __FILE__, __func__, size * sizeof(int));
  auto** bc = static_cast<NonZeroElem**>(mxMalloc(size * sizeof(NonZeroElem*)));
  test_mxMalloc(bc, __LINE__, __FILE__, __func__, size * sizeof(NonZeroElem*));

  for (int t = 0; t < periods; t++)
    {
#ifdef MATLAB_MEX_FILE
      if (utIsInterruptPending())
        throw UserException {};
#endif

      if (record && symbolic)
        {
          save_op = static_cast<int*>(mxMalloc(nop * sizeof(int)));
          test_mxMalloc(save_op, __LINE__, __FILE__, __func__, nop * sizeof(int));
          nopa = nop;
        }
      nop = 0;
      Clear_u();
      int ti = t * size;
      for (int i = ti; i < size + ti; i++)
        {
          /*finding the max-pivot*/
          double piv = 0, piv_abs = 0;
          auto [nb_eq, first] = At_Col(i, 0);
          if ((symbolic && t <= start_compare) || !symbolic)
            {
              int l = 0, N_max = 0;
              bool one = false;
              piv_abs = 0;
              for (int j = 0; j < nb_eq; j++)
                {
                  if (!line_done[first->r_index])
                    {
                      int k = first->u_index;
                      int jj = first->r_index;
                      int NRow_jj = NRow(jj);
                      piv_v[l] = u[k];
                      double piv_fabs = fabs(u[k]);
                      pivj_v[l] = jj;
                      pivk_v[l] = k;
                      NR[l] = NRow_jj;
                      if (NRow_jj == 1 && !one)
                        {
                          one = true;
                          piv_abs = piv_fabs;
                          N_max = NRow_jj;
                        }
                      if (!one)
                        {
                          if (piv_fabs > piv_abs)
                            piv_abs = piv_fabs;
                          if (NRow_jj > N_max)
                            N_max = NRow_jj;
                        }
                      else if (NRow_jj == 1)
                        {
                          if (piv_fabs > piv_abs)
                            piv_abs = piv_fabs;
                          if (NRow_jj > N_max)
                            N_max = NRow_jj;
                        }
                      l++;
                    }
                  first = first->NZE_C_N;
                }
              double markovitz = 0, markovitz_max = -9e70;
              int NR_max = 0;
              if (!one)
                for (int j = 0; j < l; j++)
                  {
                    if (N_max > 0 && NR[j] > 0)
                      {
                        if (fabs(piv_v[j]) > 0)
                          {
                            if (markowitz_c > 0)
                              markovitz = exp(log(fabs(piv_v[j]) / piv_abs)
                                              - markowitz_c
                                                    * log(static_cast<double>(NR[j])
                                                          / static_cast<double>(N_max)));
                            else
                              markovitz = fabs(piv_v[j]) / piv_abs;
                          }
                        else
                          markovitz = 0;
                      }
                    else
                      markovitz = fabs(piv_v[j]) / piv_abs;
                    if (markovitz > markovitz_max)
                      {
                        piv = piv_v[j];
                        pivj = pivj_v[j]; // Line number
                        pivk = pivk_v[j]; // positi
                        markovitz_max = markovitz;
                        NR_max = NR[j];
                      }
                  }
              else
                for (int j = 0; j < l; j++)
                  {
                    if (N_max > 0 && NR[j] > 0)
                      {
                        if (fabs(piv_v[j]) > 0)
                          {
                            if (markowitz_c > 0)
                              markovitz = exp(log(fabs(piv_v[j]) / piv_abs)
                                              - markowitz_c
                                                    * log(static_cast<double>(NR[j])
                                                          / static_cast<double>(N_max)));
                            else
                              markovitz = fabs(piv_v[j]) / piv_abs;
                          }
                        else
                          markovitz = 0;
                      }
                    else
                      markovitz = fabs(piv_v[j]) / piv_abs;
                    if (NR[j] == 1)
                      {
                        piv = piv_v[j];
                        pivj = pivj_v[j]; // Line number
                        pivk = pivk_v[j]; // positi
                        markovitz_max = markovitz;
                        NR_max = NR[j];
                      }
                  }
              if (fabs(piv) < eps && verbosity >= 1)
                mexPrintf(
                    "==> Error NR_max=%d, N_max=%d and piv=%f, piv_abs=%f, markovitz_max=%f\n",
                    NR_max, N_max, piv, piv_abs, markovitz_max);
              if (NR_max == 0 && verbosity >= 1)
                mexPrintf("==> Error NR_max=0 and piv=%f, markovitz_max=%f\n", piv, markovitz_max);
              pivot[i] = pivj;
              pivot_save[i] = pivj;
              pivotk[i] = pivk;
              pivotv[i] = piv;
            }
          else
            {
              pivj = pivot[i - size] + size;
              pivot[i] = pivj;
              first = At_Pos(pivj, i);
              pivk = first->u_index;
              piv = u[pivk];
              piv_abs = fabs(piv);
            }
          line_done[pivj] = true;

          if (record && symbolic)
            {
              if (nop + 1 >= nopa)
                {
                  nopa = static_cast<long>(mem_increasing_factor * static_cast<double>(nopa));
                  save_op = static_cast<int*>(mxRealloc(save_op, nopa * sizeof(int)));
                }
              auto* save_op_s = reinterpret_cast<t_save_op_s*>(&save_op[nop]);
              save_op_s->operat = IFLD;
              save_op_s->first = pivk;
              save_op_s->lag = 0;
              nop += 2;
              if (piv_abs < eps)
                {
                  if (block_num > 1)
                    throw FatalException {"In Solve_ByteCode_Symbolic_Sparse_GaussianElimination, "
                                          "singular system in block "
                                          + to_string(block_num + 1)};
                  else
                    throw FatalException {
                        "In Solve_ByteCode_Symbolic_Sparse_GaussianElimination, singular system"};
                }
              /*divide all the non zeros elements of the line pivj by the max_pivot*/
              int nb_var;
              tie(nb_var, first) = At_Row(pivj);
              for (int j = 0; j < nb_var; j++)
                {
                  u[first->u_index] /= piv;
                  if (nop + j * 2 + 1 >= nopa)
                    {
                      nopa = static_cast<long>(mem_increasing_factor * static_cast<double>(nopa));
                      save_op = static_cast<int*>(mxRealloc(save_op, nopa * sizeof(int)));
                    }
                  save_op_s = reinterpret_cast<t_save_op_s*>(&save_op[nop + j * 2]);
                  save_op_s->operat = IFDIV;
                  save_op_s->first = first->u_index;
                  save_op_s->lag = first->lag_index;
                  first = first->NZE_R_N;
                }
              nop += nb_var * 2;
              u[b[pivj]] /= piv;
              if (nop + 1 >= nopa)
                {
                  nopa = static_cast<long>(mem_increasing_factor * static_cast<double>(nopa));
                  save_op = static_cast<int*>(mxRealloc(save_op, nopa * sizeof(int)));
                }
              save_op_s = reinterpret_cast<t_save_op_s*>(&save_op[nop]);
              save_op_s->operat = IFDIV;
              save_op_s->first = b[pivj];
              save_op_s->lag = 0;
              nop += 2;
              // Subtract the elements on the non treated lines
              tie(nb_eq, first) = At_Col(i);
              auto [nb_var_piva, first_piva] = At_Row(pivj);

              int nb_eq_todo = 0;
              for (int j = 0; j < nb_eq && first; j++)
                {
                  if (!line_done[first->r_index])
                    bc[nb_eq_todo++] = first;
                  first = first->NZE_C_N;
                }
              for (int j = 0; j < nb_eq_todo; j++)
                {
                  t_save_op_s* save_op_s_l;
                  NonZeroElem* first = bc[j];
                  int row = first->r_index;
                  double first_elem = u[first->u_index];
                  if (nop + 1 >= nopa)
                    {
                      nopa = static_cast<long>(mem_increasing_factor * static_cast<double>(nopa));
                      save_op = static_cast<int*>(mxRealloc(save_op, nopa * sizeof(int)));
                    }
                  save_op_s_l = reinterpret_cast<t_save_op_s*>(&save_op[nop]);
                  save_op_s_l->operat = IFLD;
                  save_op_s_l->first = first->u_index;
                  save_op_s_l->lag = abs(first->lag_index);
                  nop += 2;

                  int nb_var_piv = nb_var_piva;
                  NonZeroElem* first_piv = first_piva;
                  auto [nb_var_sub, first_sub] = At_Row(row);
                  int l_sub = 0;
                  int l_piv = 0;
                  int sub_c_index = first_sub->c_index;
                  int piv_c_index = first_piv->c_index;
                  int tmp_lag = first_sub->lag_index;
                  while (l_sub < nb_var_sub /*=NRow(row)*/ || l_piv < nb_var_piv)
                    {
                      if (l_sub < nb_var_sub && (sub_c_index < piv_c_index || l_piv >= nb_var_piv))
                        {
                          /* There is no nonzero element at row pivot for this
                             column ⇒ Nothing to do for the current element got to
                             next column */
                          first_sub = first_sub->NZE_R_N;
                          if (first_sub)
                            sub_c_index = first_sub->c_index;
                          else
                            sub_c_index = size * periods;
                          l_sub++;
                        }
                      else if (sub_c_index > piv_c_index || l_sub >= nb_var_sub)
                        {
                          // There is an nonzero element at row pivot but not at the current row=>
                          // insert a negative element in the current row
                          int tmp_u_count = Get_u();
                          int lag = first_piv->c_index / size - row / size;
                          Insert(row, first_piv->c_index, tmp_u_count, lag);
                          u[tmp_u_count] = -u[first_piv->u_index] * first_elem;
                          if (nop + 2 >= nopa)
                            {
                              nopa = static_cast<long>(mem_increasing_factor
                                                       * static_cast<double>(nopa));
                              save_op = static_cast<int*>(mxRealloc(save_op, nopa * sizeof(int)));
                            }
                          save_op_s_l = reinterpret_cast<t_save_op_s*>(&save_op[nop]);
                          save_op_s_l->operat = IFLESS;
                          save_op_s_l->first = tmp_u_count;
                          save_op_s_l->second = first_piv->u_index;
                          save_op_s_l->lag = max(first_piv->lag_index, abs(tmp_lag));
                          nop += 3;
                          first_piv = first_piv->NZE_R_N;
                          if (first_piv)
                            piv_c_index = first_piv->c_index;
                          else
                            piv_c_index = size * periods;
                          l_piv++;
                        }
                      else /*first_sub->c_index==first_piv->c_index*/
                        {
                          if (i == sub_c_index)
                            {
                              NonZeroElem* firsta = first;
                              NonZeroElem* first_suba = first_sub->NZE_R_N;
                              Delete(first_sub->r_index, first_sub->c_index);
                              first = firsta->NZE_C_N;
                              first_sub = first_suba;
                              if (first_sub)
                                sub_c_index = first_sub->c_index;
                              else
                                sub_c_index = size * periods;
                              l_sub++;
                              first_piv = first_piv->NZE_R_N;
                              if (first_piv)
                                piv_c_index = first_piv->c_index;
                              else
                                piv_c_index = size * periods;
                              l_piv++;
                            }
                          else
                            {
                              u[first_sub->u_index] -= u[first_piv->u_index] * first_elem;
                              if (nop + 3 >= nopa)
                                {
                                  nopa = static_cast<long>(mem_increasing_factor
                                                           * static_cast<double>(nopa));
                                  save_op
                                      = static_cast<int*>(mxRealloc(save_op, nopa * sizeof(int)));
                                }
                              save_op_s_l = reinterpret_cast<t_save_op_s*>(&save_op[nop]);
                              save_op_s_l->operat = IFSUB;
                              save_op_s_l->first = first_sub->u_index;
                              save_op_s_l->second = first_piv->u_index;
                              save_op_s_l->lag = max(abs(tmp_lag), first_piv->lag_index);
                              nop += 3;
                              first_sub = first_sub->NZE_R_N;
                              if (first_sub)
                                sub_c_index = first_sub->c_index;
                              else
                                sub_c_index = size * periods;
                              l_sub++;
                              first_piv = first_piv->NZE_R_N;
                              if (first_piv)
                                piv_c_index = first_piv->c_index;
                              else
                                piv_c_index = size * periods;
                              l_piv++;
                            }
                        }
                    }
                  u[b[row]] -= u[b[pivj]] * first_elem;

                  if (nop + 3 >= nopa)
                    {
                      nopa = static_cast<long>(mem_increasing_factor * static_cast<double>(nopa));
                      save_op = static_cast<int*>(mxRealloc(save_op, nopa * sizeof(int)));
                    }
                  save_op_s_l = reinterpret_cast<t_save_op_s*>(&save_op[nop]);
                  save_op_s_l->operat = IFSUB;
                  save_op_s_l->first = b[row];
                  save_op_s_l->second = b[pivj];
                  save_op_s_l->lag = abs(tmp_lag);
                  nop += 3;
                }
            }
          else if (symbolic)
            {
              nop += 2;
              if (piv_abs < eps)
                {
                  if (block_num > 1)
                    throw FatalException {"In Solve_ByteCode_Symbolic_Sparse_GaussianElimination, "
                                          "singular system in block "
                                          + to_string(block_num + 1)};
                  else
                    throw FatalException {
                        "In Solve_ByteCode_Symbolic_Sparse_GaussianElimination, singular system"};
                }
              // Divide all the non zeros elements of the line pivj by the max_pivot
              int nb_var;
              tie(nb_var, first) = At_Row(pivj);
              for (int j = 0; j < nb_var; j++)
                {
                  u[first->u_index] /= piv;
                  first = first->NZE_R_N;
                }
              nop += nb_var * 2;
              u[b[pivj]] /= piv;
              nop += 2;
              // Subtract the elements on the non treated lines
              tie(nb_eq, first) = At_Col(i);
              auto [nb_var_piva, first_piva] = At_Row(pivj);

              int nb_eq_todo = 0;
              for (int j = 0; j < nb_eq && first; j++)
                {
                  if (!line_done[first->r_index])
                    bc[nb_eq_todo++] = first;
                  first = first->NZE_C_N;
                }
              for (int j = 0; j < nb_eq_todo; j++)
                {
                  NonZeroElem* first = bc[j];
                  int row = first->r_index;
                  double first_elem = u[first->u_index];
                  nop += 2;
                  int nb_var_piv = nb_var_piva;
                  NonZeroElem* first_piv = first_piva;
                  auto [nb_var_sub, first_sub] = At_Row(row);
                  int l_sub = 0;
                  int l_piv = 0;
                  int sub_c_index = first_sub->c_index;
                  int piv_c_index = first_piv->c_index;
                  while (l_sub < nb_var_sub /*= NRow(row)*/ || l_piv < nb_var_piv)
                    {
                      if (l_sub < nb_var_sub && (sub_c_index < piv_c_index || l_piv >= nb_var_piv))
                        {
                          /* There is no nonzero element at row pivot for this
                             column ⇒ Nothing to do for the current element got to
                             next column */
                          first_sub = first_sub->NZE_R_N;
                          if (first_sub)
                            sub_c_index = first_sub->c_index;
                          else
                            sub_c_index = size * periods;
                          l_sub++;
                        }
                      else if (sub_c_index > piv_c_index || l_sub >= nb_var_sub)
                        {
                          /* There is an nonzero element at row pivot but not
                             at the current row ⇒ insert a negative element in the
                             current row */
                          int tmp_u_count = Get_u();
                          int lag = first_piv->c_index / size - row / size;
                          Insert(row, first_piv->c_index, tmp_u_count, lag);
                          u[tmp_u_count] = -u[first_piv->u_index] * first_elem;
                          nop += 3;
                          first_piv = first_piv->NZE_R_N;
                          if (first_piv)
                            piv_c_index = first_piv->c_index;
                          else
                            piv_c_index = size * periods;
                          l_piv++;
                        }
                      else /*first_sub->c_index==first_piv->c_index*/
                        {
                          if (i == sub_c_index)
                            {
                              NonZeroElem* firsta = first;
                              NonZeroElem* first_suba = first_sub->NZE_R_N;
                              Delete(first_sub->r_index, first_sub->c_index);
                              first = firsta->NZE_C_N;
                              first_sub = first_suba;
                              if (first_sub)
                                sub_c_index = first_sub->c_index;
                              else
                                sub_c_index = size * periods;
                              l_sub++;
                              first_piv = first_piv->NZE_R_N;
                              if (first_piv)
                                piv_c_index = first_piv->c_index;
                              else
                                piv_c_index = size * periods;
                              l_piv++;
                            }
                          else
                            {
                              u[first_sub->u_index] -= u[first_piv->u_index] * first_elem;
                              nop += 3;
                              first_sub = first_sub->NZE_R_N;
                              if (first_sub)
                                sub_c_index = first_sub->c_index;
                              else
                                sub_c_index = size * periods;
                              l_sub++;
                              first_piv = first_piv->NZE_R_N;
                              if (first_piv)
                                piv_c_index = first_piv->c_index;
                              else
                                piv_c_index = size * periods;
                              l_piv++;
                            }
                        }
                    }
                  u[b[row]] -= u[b[pivj]] * first_elem;
                  nop += 3;
                }
            }
        }
      if (symbolic)
        {
          if (t > static_cast<int>(periods * 0.35))
            {
              symbolic = false;
              mxFree(save_opaa);
              mxFree(save_opa);
              mxFree(save_op);
            }
          else if (record && nop == nop1)
            {
              if (t > static_cast<int>(periods * 0.35))
                {
                  symbolic = false;
                  if (save_opaa)
                    {
                      mxFree(save_opaa);
                      save_opaa = nullptr;
                    }
                  if (save_opa)
                    {
                      mxFree(save_opa);
                      save_opa = nullptr;
                    }
                  if (save_op)
                    {
                      mxFree(save_op);
                      save_op = nullptr;
                    }
                }
              else if (save_opa && save_opaa)
                {
                  if (compare(save_op, save_opa, save_opaa, t, nop))
                    {
                      tbreak = t;
                      tbreak_g = tbreak;
                      break;
                    }
                }
              if (save_opa)
                {
                  if (save_opaa)
                    {
                      mxFree(save_opaa);
                      save_opaa = nullptr;
                    }
                  save_opaa = save_opa;
                }
              save_opa = save_op;
            }
          else
            {
              if (nop == nop1)
                record = true;
              else
                {
                  record = false;
                  if (save_opa)
                    {
                      mxFree(save_opa);
                      save_opa = nullptr;
                    }
                  if (save_opaa)
                    {
                      mxFree(save_opaa);
                      save_opaa = nullptr;
                    }
                }
            }
          nop1 = nop;
        }
    }
  mxFree(bc);
  mxFree(piv_v);
  mxFree(pivj_v);
  mxFree(pivk_v);
  mxFree(NR);
  if (symbolic)
    {
      if (save_op)
        mxFree(save_op);
      if (save_opa)
        mxFree(save_opa);
      if (save_opaa)
        mxFree(save_opaa);
    }

  // The backward substitution
  for (int i = 0; i < y_size * (periods + y_kmin); i++)
    ya[i] = y[i];
  slowc_save = slowc;
  bksub(tbreak, last_period);
  End_Gaussian_Elimination();
}

void
Interpreter::Check_and_Correct_Previous_Iteration()
{
  if (isnan(res1) || isinf(res1) || (res2 > g0 && iter > 0))
    {
      while (isnan(res1) || isinf(res1))
        {
          prev_slowc_save = slowc_save;
          slowc_save /= 1.1;
          for (int i = 0; i < size; i++)
            {
              int eq = index_vara[i];
              y[eq + it_ * y_size]
                  = ya[eq + it_ * y_size] + slowc_save * direction[eq + it_ * y_size];
            }
          compute_complete(true);
        }

      while (res2 > g0 && slowc_save > 1e-1)
        {
          prev_slowc_save = slowc_save;
          slowc_save /= 1.5;
          for (int i = 0; i < size; i++)
            {
              int eq = index_vara[i];
              y[eq + it_ * y_size]
                  = ya[eq + it_ * y_size] + slowc_save * direction[eq + it_ * y_size];
            }
          compute_complete(true);
        }
      if (verbosity >= 2)
        mexPrintf("Error: Simulation diverging, trying to correct it using slowc=%f\n", slowc_save);
      for (int i = 0; i < size; i++)
        {
          int eq = index_vara[i];
          y[eq + it_ * y_size] = ya[eq + it_ * y_size] + slowc_save * direction[eq + it_ * y_size];
        }
      compute_complete(false);
    }
  else
    for (int i = 0; i < size; i++)
      {
        int eq = index_vara[i];
        y[eq + it_ * y_size] = ya[eq + it_ * y_size] + slowc_save * direction[eq + it_ * y_size];
      }
  slowc_save = slowc;
}

bool
Interpreter::Simulate_One_Boundary()
{
  mxArray *b_m = nullptr, *A_m = nullptr, *x0_m = nullptr;
  SuiteSparse_long *Ap = nullptr, *Ai = nullptr;
  double *Ax = nullptr, *b = nullptr;
  int preconditioner = 1;

  try_at_iteration = 0;
  Clear_u();
  bool singular_system = false;
  u_count_alloc_save = u_count_alloc;

  if (isnan(res1) || isinf(res1))
    {
#ifdef DEBUG
      for (int j = 0; j < y_size; j++)
        {
          bool select = false;
          for (int i = 0; i < size; i++)
            if (j == index_vara[i])
              {
                select = true;
                break;
              }
          if (select)
            mexPrintf("-> variable %s (%d) at time %d = %f direction = %f\n",
                      get_variable(SymbolType::endogenous, j).c_str(), j + 1, it_,
                      y[j + it_ * y_size], direction[j + it_ * y_size]);
          else
            mexPrintf("   variable %s (%d) at time %d = %f direction = %f\n",
                      get_variable(SymbolType::endogenous, j).c_str(), j + 1, it_,
                      y[j + it_ * y_size], direction[j + it_ * y_size]);
        }
#endif
      if (steady_state)
        {
          if (verbosity >= 1)
            {
              if (iter == 0)
                mexPrintf(" the initial values of endogenous variables are too far from the "
                          "solution.\nChange them!\n");
              else
                mexPrintf(
                    " dynare cannot improve the simulation in block %d at time %d (variable %d)\n",
                    block_num + 1, it_ + 1, index_vara[max_res_idx] + 1);
              mexEvalString("drawnow;");
            }
        }
      else
        {
          if (iter == 0)
            throw FatalException {"In Simulate_One_Boundary, The initial values of endogenous "
                                  "variables are too far from the solution. Change them!"};
          else
            throw FatalException {
                "In Simulate_One_Boundary, Dynare cannot improve the simulation in block "
                + to_string(block_num + 1) + " at time " + to_string(it_ + 1) + " (variable "
                + to_string(index_vara[max_res_idx] + 1)};
        }
    }

  if (verbosity >= 1)
    {
      if (steady_state)
        {
          switch (solve_algo)
            {
            case 5:
              mexPrintf("MODEL STEADY STATE: (method=Sparse Gaussian Elimination)\n");
              break;
            case 6:
              mexPrintf("MODEL STEADY STATE: (method=Sparse LU)\n");
              break;
            case 7:
              mexPrintf(preconditioner_print_out("MODEL STEADY STATE: (method=GMRES)\n",
                                                 preconditioner, true)
                            .c_str());
              break;
            case 8:
              mexPrintf(preconditioner_print_out("MODEL STEADY STATE: (method=BiCGStab)\n",
                                                 preconditioner, true)
                            .c_str());
              break;
            }
        }

      mexPrintf("------------------------------------\n");
      mexPrintf("      Iteration no. %d\n", iter + 1);
      mexPrintf("      Inf-norm error = %.3e\n", static_cast<double>(max_res));
      mexPrintf("      2-norm error   = %.3e\n", static_cast<double>(sqrt(res2)));
      mexPrintf("      1-norm error   = %.3e\n", static_cast<double>(res1));
      mexPrintf("------------------------------------\n");
    }
  bool zero_solution;

  if ((solve_algo == 5 && steady_state) || (stack_solve_algo == 5 && !steady_state))
    zero_solution = Simple_Init();
  else
    {
      x0_m = mxCreateDoubleMatrix(size, 1, mxREAL);
      if (!x0_m)
        throw FatalException {"In Simulate_One_Boundary, can't allocate x0_m vector"};
      if (!((solve_algo == 6 && steady_state)
            || ((stack_solve_algo == 0 || stack_solve_algo == 1 || stack_solve_algo == 4
                 || stack_solve_algo == 6)
                && !steady_state)))
        {
          b_m = mxCreateDoubleMatrix(size, 1, mxREAL);
          if (!b_m)
            throw FatalException {"In Simulate_One_Boundary, can't allocate b_m vector"};
          A_m = mxCreateSparse(size, size, min(static_cast<int>(IM_i.size() * 2), size * size),
                               mxREAL);
          if (!A_m)
            throw FatalException {"In Simulate_One_Boundary, can't allocate A_m matrix"};
          zero_solution = Init_Matlab_Sparse_One_Boundary(A_m, b_m, x0_m);
        }
      else
        {
          tie(zero_solution, Ap, Ai, Ax, b) = Init_UMFPACK_Sparse_One_Boundary(x0_m);
          if (Ap_save[size] != Ap[size])
            {
              mxFree(Ai_save);
              mxFree(Ax_save);
              Ai_save
                  = static_cast<SuiteSparse_long*>(mxMalloc(Ap[size] * sizeof(SuiteSparse_long)));
              test_mxMalloc(Ai_save, __LINE__, __FILE__, __func__,
                            Ap[size] * sizeof(SuiteSparse_long));
              Ax_save = static_cast<double*>(mxMalloc(Ap[size] * sizeof(double)));
              test_mxMalloc(Ax_save, __LINE__, __FILE__, __func__, Ap[size] * sizeof(double));
            }
          copy_n(Ap, size + 1, Ap_save);
          copy_n(Ai, Ap[size], Ai_save);
          copy_n(Ax, Ap[size], Ax_save);
          copy_n(b, size, b_save);
        }
    }
  if (zero_solution)
    for (int i = 0; i < size; i++)
      {
        int eq = index_vara[i];
        double yy = -(y[eq + it_ * y_size]);
        direction[eq] = yy;
        y[eq + it_ * y_size] += slowc * yy;
      }
  else
    {
      if ((solve_algo == 5 && steady_state) || (stack_solve_algo == 5 && !steady_state))
        singular_system = Solve_ByteCode_Sparse_GaussianElimination();
      else if ((solve_algo == 7 && steady_state) || (stack_solve_algo == 2 && !steady_state))
        Solve_Matlab_GMRES(A_m, b_m, false, x0_m);
      else if ((solve_algo == 8 && steady_state) || (stack_solve_algo == 3 && !steady_state))
        Solve_Matlab_BiCGStab(A_m, b_m, false, x0_m, preconditioner);
      else if ((solve_algo == 6 && steady_state)
               || ((stack_solve_algo == 0 || stack_solve_algo == 1 || stack_solve_algo == 4
                    || stack_solve_algo == 6)
                   && !steady_state))
        {
          Solve_LU_UMFPack_One_Boundary(Ap, Ai, Ax, b);
          mxDestroyArray(x0_m);
        }
    }
  return singular_system;
}

bool
Interpreter::solve_linear(bool do_check_and_correct)
{
  bool cvg = false;
  compute_complete(false);
  cvg = (max_res < solve_tolf);
  if (!cvg || isnan(res1) || isinf(res1))
    {
      if (do_check_and_correct)
        Check_and_Correct_Previous_Iteration();
      bool singular_system = Simulate_One_Boundary();
      if (singular_system && verbosity >= 1)
        Singular_display();
    }
  return cvg;
}

void
Interpreter::solve_non_linear()
{
  max_res_idx = 0;
  bool cvg = false;
  iter = 0;
  glambda2 = g0 = very_big;
  while (!cvg && iter < maxit_)
    {
      cvg = solve_linear(iter > 0);
      g0 = res2;
      iter++;
    }
  if (!cvg)
    {
      if (steady_state)
        throw FatalException {
            "In Solve Forward/Backward Complete, convergence not achieved in block "
            + to_string(block_num + 1) + ", after " + to_string(iter) + " iterations"};
      else
        throw FatalException {
            "In Solve Forward/Backward Complete, convergence not achieved in block "
            + to_string(block_num + 1) + ", at time " + to_string(it_) + ", after "
            + to_string(iter) + " iterations"};
    }
}

void
Interpreter::Simulate_Newton_One_Boundary(bool forward)
{
  r = static_cast<double*>(mxMalloc(size * sizeof(double)));
  test_mxMalloc(r, __LINE__, __FILE__, __func__, size * sizeof(double));
  iter = 0;
  if ((solve_algo == 6 && steady_state)
      || ((stack_solve_algo == 0 || stack_solve_algo == 1 || stack_solve_algo == 4
           || stack_solve_algo == 6)
          && !steady_state))
    {
      Ap_save = static_cast<SuiteSparse_long*>(mxMalloc((size + 1) * sizeof(SuiteSparse_long)));
      test_mxMalloc(Ap_save, __LINE__, __FILE__, __func__, (size + 1) * sizeof(SuiteSparse_long));
      Ap_save[size] = 0;
      Ai_save = static_cast<SuiteSparse_long*>(mxMalloc(1 * sizeof(SuiteSparse_long)));
      test_mxMalloc(Ai_save, __LINE__, __FILE__, __func__, 1 * sizeof(SuiteSparse_long));
      Ax_save = static_cast<double*>(mxMalloc(1 * sizeof(double)));
      test_mxMalloc(Ax_save, __LINE__, __FILE__, __func__, 1 * sizeof(double));
      b_save = static_cast<double*>(mxMalloc((size) * sizeof(SuiteSparse_long)));
      test_mxMalloc(b_save, __LINE__, __FILE__, __func__, (size) * sizeof(SuiteSparse_long));
    }
  if (steady_state)
    {
      it_ = 0;
      if (!is_linear)
        solve_non_linear();
      else
        solve_linear(false);
    }
  else if (forward)
    {
      if (!is_linear)
        for (it_ = y_kmin; it_ < periods + y_kmin; it_++)
          solve_non_linear();
      else
        for (it_ = y_kmin; it_ < periods + y_kmin; it_++)
          solve_linear(false);
    }
  else
    {
      if (!is_linear)
        for (it_ = periods + y_kmin - 1; it_ >= y_kmin; it_--)
          solve_non_linear();
      else
        for (it_ = periods + y_kmin - 1; it_ >= y_kmin; it_--)
          solve_linear(false);
    }
  if ((solve_algo == 6 && steady_state)
      || ((stack_solve_algo == 0 || stack_solve_algo == 1 || stack_solve_algo == 4
           || stack_solve_algo == 6)
          && !steady_state))
    {
      mxFree(Ap_save);
      mxFree(Ai_save);
      mxFree(Ax_save);
      mxFree(b_save);
    }
  mxFree(r);
}

string
Interpreter::preconditioner_print_out(string s, int preconditioner, bool ss)
{
  int n = s.length();
  string tmp = ", preconditioner=";
  switch (preconditioner)
    {
    case 0:
      if (ss)
        tmp.append("Jacobi on static jacobian");
      else
        tmp.append("Jacobi on dynamic jacobian");
      break;
    case 1:
      if (ss)
        tmp.append("incomplete lutp on static jacobian");
      else
        tmp.append("incomplete lu0 on dynamic jacobian");
      break;
    case 2:
      tmp.append("incomplete lutp on dynamic jacobian");
      break;
    case 3:
      tmp.append("lu on static jacobian");
      break;
    }
  s.insert(n - 2, tmp);
  return s;
}

void
Interpreter::Simulate_Newton_Two_Boundaries(
    bool cvg, const vector_table_conditional_local_type& vector_table_conditional_local)
{
  double top = 0.5;
  double bottom = 0.1;
  int preconditioner = 2;
  if (start_compare == 0)
    start_compare = y_kmin;
  u_count_alloc_save = u_count_alloc;
  auto t1 {chrono::high_resolution_clock::now()};
  nop1 = 0;
  mxArray *b_m = nullptr, *A_m = nullptr, *x0_m = nullptr;
  double *Ax {nullptr}, *b {nullptr};
  SuiteSparse_long *Ap = nullptr, *Ai = nullptr;

  assert(stack_solve_algo == 0 || stack_solve_algo == 2 || stack_solve_algo == 3
         || stack_solve_algo == 4 || stack_solve_algo == 5);

  if (isnan(res1) || isinf(res1) || (res2 > 12 * g0 && iter > 0))
    {
      if (iter == 0 || fabs(slowc_save) < 1e-8)
        {
          if (verbosity >= 2)
            mexPrintf("res1 = %f, res2 = %f g0 = %f iter = %d\n", res1, res2, g0, iter);
          for (int j = 0; j < y_size; j++)
            {
              bool select = false;
              for (int i = 0; i < size; i++)
                if (j == index_vara[i])
                  {
                    select = true;
                    break;
                  }
              if (verbosity >= 2)
                {
                  if (select)
                    mexPrintf("-> variable %s (%d) at time %d = %f direction = %f\n",
                              symbol_table.getName(SymbolType::endogenous, j).c_str(), j + 1, it_,
                              y[j + it_ * y_size], direction[j + it_ * y_size]);
                  else
                    mexPrintf("   variable %s (%d) at time %d = %f direction = %f\n",
                              symbol_table.getName(SymbolType::endogenous, j).c_str(), j + 1, it_,
                              y[j + it_ * y_size], direction[j + it_ * y_size]);
                }
            }
          if (iter == 0)
            throw FatalException {
                "In Simulate_Newton_Two_Boundaries, the initial values of endogenous variables are "
                "too far from the solution. Change them!"};
          else
            throw FatalException {
                "In Simulate_Newton_Two_Boundaries, dynare cannot improve the simulation in block "
                + to_string(block_num + 1) + " at time " + to_string(it_ + 1) + " (variable "
                + to_string(index_vara[max_res_idx] + 1) + " = " + to_string(max_res) + ")"};
        }
      if (!(isnan(res1) || isinf(res1)) && !(isnan(g0) || isinf(g0))
          && (stack_solve_algo == 4 || stack_solve_algo == 5))
        {
          if (try_at_iteration == 0)
            {
              prev_slowc_save = slowc_save;
              slowc_save = max(-gp0 / (2 * (res2 - g0 - gp0)), bottom);
            }
          else
            {
              double t1 = res2 - gp0 * slowc_save - g0;
              double t2 = glambda2 - gp0 * prev_slowc_save - g0;
              double a = (1 / (slowc_save * slowc_save) * t1
                          - 1 / (prev_slowc_save * prev_slowc_save) * t2)
                         / (slowc_save - prev_slowc_save);
              double b = (-prev_slowc_save / (slowc_save * slowc_save) * t1
                          + slowc_save / (prev_slowc_save * prev_slowc_save) * t2)
                         / (slowc_save - prev_slowc_save);
              prev_slowc_save = slowc_save;
              slowc_save = max(min(-b + sqrt(b * b - 3 * a * gp0) / (3 * a), top * slowc_save),
                               bottom * slowc_save);
            }
          glambda2 = res2;
          try_at_iteration++;
          if (slowc_save <= bottom)
            {
              for (int i = 0; i < y_size * (periods + y_kmin); i++)
                y[i] = ya[i] + direction[i];
              g0 = res2;
              gp0 = -res2;
              try_at_iteration = 0;
              iter--;
              return;
            }
        }
      else
        {
          prev_slowc_save = slowc_save;
          slowc_save /= 1.05;
        }
      if (verbosity >= 2)
        {
          if (isnan(res1) || isinf(res1))
            mexPrintf("The model cannot be evaluated, trying to correct it using slowc=%f\n",
                      slowc_save);
          else
            mexPrintf("Simulation diverging, trying to correct it using slowc=%f\n", slowc_save);
        }
      for (int i = 0; i < y_size * (periods + y_kmin); i++)
        y[i] = ya[i] + slowc_save * direction[i];
      iter--;
      return;
    }
  u_count += u_count_init;
  if (stack_solve_algo == 5)
    {
      if (alt_symbolic && alt_symbolic_count < alt_symbolic_count_max)
        {
          if (verbosity >= 2)
            mexPrintf("Pivoting method will be applied only to the first periods.\n");
          alt_symbolic = false;
          symbolic = true;
          markowitz_c = markowitz_c_s;
          alt_symbolic_count++;
        }
      if (res1 / res1a - 1 > -0.3 && symbolic && iter > 0)
        {
          if (restart > 2)
            {
              if (verbosity >= 2)
                mexPrintf("Divergence or slowdown occurred during simulation.\nIn the next "
                          "iteration, pivoting method will be applied to all periods.\n");
              symbolic = false;
              alt_symbolic = true;
              markowitz_c_s = markowitz_c;
              markowitz_c = 0;
            }
          else
            {
              if (verbosity >= 2)
                mexPrintf("Divergence or slowdown occurred during simulation.\nIn the next "
                          "iteration, pivoting method will be applied for a longer period.\n");
              start_compare = min(tbreak_g, periods);
              restart++;
            }
        }
      else
        {
          start_compare = max(y_kmin, minimal_solving_periods);
          restart = 0;
        }
    }
  res1a = res1;
  if (verbosity >= 1)
    {
      if (iter == 0)
        {
          switch (stack_solve_algo)
            {
            case 0:
              mexPrintf("MODEL SIMULATION: (method=Sparse LU solver on stacked system)\n");
              break;
            case 2:
              mexPrintf(
                  preconditioner_print_out("MODEL SIMULATION: (method=GMRES on stacked system)\n",
                                           preconditioner, false)
                      .c_str());
              break;
            case 3:
              mexPrintf(preconditioner_print_out(
                            "MODEL SIMULATION: (method=BiCGStab on stacked system)\n",
                            preconditioner, false)
                            .c_str());
              break;
            case 4:
              mexPrintf("MODEL SIMULATION: (method=Sparse LU solver with optimal path length on "
                        "stacked system)\n");
              break;
            case 5:
              mexPrintf("MODEL SIMULATION: (method=LBJ with Sparse Gaussian Elimination)\n");
              break;
            }
        }
      mexPrintf("------------------------------------\n");
      mexPrintf("      Iteration no. %d\n", iter + 1);
      mexPrintf("      Inf-norm error = %.3e\n", static_cast<double>(max_res));
      mexPrintf("      2-norm error   = %.3e\n", static_cast<double>(sqrt(res2)));
      mexPrintf("      1-norm error   = %.3e\n", static_cast<double>(res1));
      mexPrintf("------------------------------------\n");
      mexEvalString("drawnow;");
    }
  if (cvg)
    return;
  else
    {
      if (stack_solve_algo == 5)
        Init_Gaussian_Elimination();
      else
        {
          x0_m = mxCreateDoubleMatrix(periods * size, 1, mxREAL);
          if (!x0_m)
            throw FatalException {"In Simulate_Newton_Two_Boundaries, can't allocate x0_m vector"};
          if (stack_solve_algo == 0 || stack_solve_algo == 4)
            tie(Ap, Ai, Ax, b)
                = Init_UMFPACK_Sparse_Two_Boundaries(x0_m, vector_table_conditional_local);
          else
            {
              b_m = mxCreateDoubleMatrix(periods * size, 1, mxREAL);
              if (!b_m)
                throw FatalException {
                    "In Simulate_Newton_Two_Boundaries, can't allocate b_m vector"};
              if (stack_solve_algo != 0 && stack_solve_algo != 4)
                {
                  A_m = mxCreateSparse(periods * size, periods * size, IM_i.size() * periods * 2,
                                       mxREAL);
                  if (!A_m)
                    throw FatalException {
                        "In Simulate_Newton_Two_Boundaries, can't allocate A_m matrix"};
                }
              Init_Matlab_Sparse_Two_Boundaries(A_m, b_m, x0_m);
            }
        }
      if (stack_solve_algo == 0 || stack_solve_algo == 4)
        {
          Solve_LU_UMFPack_Two_Boundaries(Ap, Ai, Ax, b, vector_table_conditional_local);
          mxDestroyArray(x0_m);
        }
      else if (stack_solve_algo == 2)
        Solve_Matlab_GMRES(A_m, b_m, true, x0_m);
      else if (stack_solve_algo == 3)
        Solve_Matlab_BiCGStab(A_m, b_m, true, x0_m, 1);
      else if (stack_solve_algo == 5)
        Solve_ByteCode_Symbolic_Sparse_GaussianElimination(symbolic);
    }
  using FloatSeconds = chrono::duration<double, chrono::seconds::period>;
  auto t2 {chrono::high_resolution_clock::now()};
  if (verbosity >= 1)
    {
      mexPrintf("(** %.2f seconds **)\n", FloatSeconds {t2 - t1}.count());
      mexEvalString("drawnow;");
    }
  if (!steady_state && stack_solve_algo == 4)
    {
      double ax = -0.1, bx = 1.1;

      auto [success, cx, fa, fb, fc] = mnbrak(ax, bx);
      if (!success)
        return;
      auto [success2, xmin] = golden(ax, bx, cx, 1e-1);
      if (!success2)
        return;
      slowc = xmin;
      if (verbosity >= 1)
        {
          auto t3 {chrono::high_resolution_clock::now()};
          mexPrintf("(** %.2f seconds **)\n", FloatSeconds {t3 - t2}.count());
          mexEvalString("drawnow;");
        }
    }
  if (tbreak_g == 0)
    tbreak_g = periods;
}

void
Interpreter::fixe_u()
{
  u_count = u_count_int * periods;
  u_count_alloc = 2 * u_count;
#ifdef DEBUG
  mexPrintf("fixe_u : alloc(%d double)\n", u_count_alloc);
#endif
  u = static_cast<double*>(mxMalloc(u_count_alloc * sizeof(double)));
  test_mxMalloc(u, __LINE__, __FILE__, __func__, u_count_alloc * sizeof(double));
#ifdef DEBUG
  mexPrintf("u=%d\n", u);
#endif
  fill_n(u, u_count_alloc, 0);
  u_count_init = u_count_int;
}
