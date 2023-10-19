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
#include <filesystem>
#include <numeric>
#include <cfenv>
#include <type_traits>
#include <chrono>

#include "Interpreter.hh"

constexpr double BIG = 1.0e+8, SMALL = 1.0e-5;

Interpreter::Interpreter(Evaluate &evaluator_arg, double *params_arg, double *y_arg, double *ya_arg, double *x_arg, double *steady_y_arg,
                         double *direction_arg, int y_size_arg,
                         int nb_row_x_arg, int periods_arg, int y_kmin_arg, int y_kmax_arg,
                         int maxit_arg_, double solve_tolf_arg, double markowitz_c_arg,
                         int minimal_solving_periods_arg, int stack_solve_algo_arg, int solve_algo_arg,
                         bool global_temporary_terms_arg, bool print_arg, mxArray *GlobalTemporaryTerms_arg,
                         bool steady_state_arg, bool block_decomposed_arg, int col_x_arg, int col_y_arg, const BasicSymbolTable &symbol_table_arg, int verbosity_arg) :
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
  T = nullptr;
  minimal_solving_periods = minimal_solving_periods_arg;
  stack_solve_algo = stack_solve_algo_arg;
  solve_algo = solve_algo_arg;
  global_temporary_terms = global_temporary_terms_arg;
  print = print_arg;
  col_x = col_x_arg;
  col_y = col_y_arg;
  GlobalTemporaryTerms = GlobalTemporaryTerms_arg;
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
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            compute_block_time(0, false, false);
          it_ = periods+y_kmin-1; // Do not leave it_ in inconsistent state
        }
      else
        {
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
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
      Per_y_ = it_*y_size;
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
                  y[Block_Contain[0].Variable + Per_y_] = ya - slowc * (r[0] / g1[0]);
                  if (fetestexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW))
                    {
                      res1 = numeric_limits<double>::quiet_NaN();
                      if (verbosity >= 1)
                        mexPrintf("      Singularity in block %d", block_num+1);
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
      y[Block_Contain[0].Variable + Per_y_] += -slowc * (rr / g1[0]);
      if (fetestexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW))
        {
          res1 = numeric_limits<double>::quiet_NaN();
          if (verbosity >= 1)
            mexPrintf("      Singularity in block %d", block_num+1);
        }
      iter++;
    }
  if (!cvg)
    throw FatalException{"In Solve Forward simple, convergence not achieved in block "
                         + to_string(block_num+1) + ", after " + to_string(iter) + " iterations"};
}

void
Interpreter::solve_simple_over_periods(bool forward)
{
  g1 = static_cast<double *>(mxMalloc(sizeof(double)));
  test_mxMalloc(g1, __LINE__, __FILE__, __func__, sizeof(double));
  r = static_cast<double *>(mxMalloc(sizeof(double)));
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
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            solve_simple_one_periods();
          it_= periods+y_kmin-1; // Do not leave it_ in inconsistent state
        }
      else
        {
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            solve_simple_one_periods();
          it_ = y_kmin; // Do not leave it_ in inconsistent state (see #1727)
        }
    }
  mxFree(g1);
  mxFree(r);
}

void
Interpreter::compute_complete_2b(bool no_derivatives, double *_res1, double *_res2, double *_max_res, int *_max_res_idx)
{
  res1 = 0;
  *_res1 = 0;
  *_res2 = 0;
  *_max_res = 0;
  for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
    {
      Per_u_ = (it_-y_kmin)*u_count_int;
      Per_y_ = it_*y_size;
      int shift = (it_-y_kmin) * size;
      compute_block_time(Per_u_, false, no_derivatives);
      if (!(isnan(res1) || isinf(res1)))
        for (int i = 0; i < size; i++)
          {
            double rr;
            rr = r[i];
            res[i+shift] = rr;
            if (max_res < fabs(rr))
              {
                *_max_res = fabs(rr);
                *_max_res_idx = i;
              }
            *_res2 += rr*rr;
            *_res1 += fabs(rr);
          }
      else
        return;
    }
  it_ = periods+y_kmin-1; // Do not leave it_ in inconsistent state
  return;
}

void
Interpreter::evaluate_a_block(bool initialization, bool single_block, const string &bin_base_name)
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
              residual[Block_Contain[j].Equation] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
        }
      else
        {
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              Per_y_ = it_*y_size;
              compute_block_time(0, true, false);
              if (single_block)
                for (int j = 0; j < size; j++)
                  residual[(it_-y_kmin)*size+j] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
              else
                for (int j = 0; j < size; j++)
                  residual[(it_-y_kmin)*y_size+Block_Contain[j].Equation] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
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
          if (!single_block)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              Per_y_ = it_*y_size;
              compute_block_time(0, true, false);
              if (!single_block)
                for (int j = 0; j < size; j++)
                  residual[(it_-y_kmin)*y_size+Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[(it_-y_kmin)*size+j] = r[j];
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
          if (!single_block)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              Per_y_ = it_*y_size;
              compute_block_time(0, true, false);
              if (!single_block)
                for (int j = 0; j < size; j++)
                  residual[(it_-y_kmin)*y_size+Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[(it_-y_kmin)*size+j] = r[j];
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
              residual[Block_Contain[j].Equation] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
        }
      else
        {
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              Per_y_ = it_*y_size;
              compute_block_time(0, true, false);
              if (single_block)
                for (int j = 0; j < size; j++)
                  residual[(it_-y_kmin)*size+j] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
              else
                for (int j = 0; j < size; j++)
                  residual[(it_-y_kmin)*y_size+Block_Contain[j].Equation] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
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
          if (!single_block)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              Per_y_ = it_*y_size;
              compute_block_time(0, true, false);
              if (!single_block)
                for (int j = 0; j < size; j++)
                  residual[(it_-y_kmin)*y_size+Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[(it_-y_kmin)*size+j] = r[j];
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
          if (!single_block)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              Per_y_ = it_*y_size;
              compute_block_time(0, true, false);
              if (!single_block)
                for (int j = 0; j < size; j++)
                  residual[(it_-y_kmin)*y_size+Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[(it_-y_kmin)*size+j] = r[j];
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
      for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
        {
          Per_u_ = (it_-y_kmin)*u_count_int;
          Per_y_ = it_*y_size;
          compute_block_time(Per_u_, true, false);
          if (!single_block)
            for (int j = 0; j < size; j++)
              residual[(it_-y_kmin)*y_size+Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[(it_-y_kmin)*size+j] = r[j];
        }
      mxFree(r);
      break;
    case BlockSimulationType::unknown:
      throw FatalException{"UNKNOWN block simulation type: impossible case"};
    }
}

int
Interpreter::simulate_a_block(const vector_table_conditional_local_type &vector_table_conditional_local, bool single_block, const string &bin_base_name)
{
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
        evaluate_a_block(true, single_block, bin_base_name);
      else
        {
          fixe_u(&u, u_count_int, u_count_int);
          Read_SparseMatrix(bin_base_name, size, 1, 0, 0, false, stack_solve_algo, solve_algo);
        }
      Per_u_ = 0;

      Simulate_Newton_One_Boundary(true);

      mxFree(u);
      mxFree(index_equa);
      mxFree(index_vara);
      fill_n(direction, y_size*col_y, 0);
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
          fixe_u(&u, u_count_int, u_count_int);
          Read_SparseMatrix(bin_base_name, size, 1, 0, 0, false, stack_solve_algo, solve_algo);
        }
      Per_u_ = 0;

      Simulate_Newton_One_Boundary(false);

      mxFree(index_equa);
      mxFree(index_vara);
      fill_n(direction, y_size*col_y, 0);
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
          while (!(cvg || (iter >= maxit_)))
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
              if (!(isnan(res1) || isinf(res1)))
                cvg = (max_res < solve_tolf);
              if (isnan(res1) || isinf(res1) || (stack_solve_algo == 4 && iter > 0))
                copy_n(y_save, y_size*(periods+y_kmax+y_kmin), y);
              u_count = u_count_saved;
              int prev_iter = iter;
              Simulate_Newton_Two_Boundaries(block_num, y_size, y_kmin, y_kmax, size, periods, cvg, minimal_solving_periods, stack_solve_algo, vector_table_conditional_local);
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

          cvg = false;
          Simulate_Newton_Two_Boundaries(block_num, y_size, y_kmin, y_kmax, size, periods, cvg, minimal_solving_periods, stack_solve_algo, vector_table_conditional_local);
          max_res = 0; max_res_idx = 0;
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
      fill_n(direction, y_size*col_y, 0);
      End_Solver();
      break;
    default:
      throw FatalException{"In simulate_a_block, Unknown type = " + to_string(static_cast<int>(type))};
      return ERROR_ON_EXIT;
    }
  return NO_ERROR_ON_EXIT;
}

void
Interpreter::check_for_controlled_exo_validity(int current_block, const vector<s_plan> &sconstrained_extended_path)
{
  vector<int> exogenous {evaluator.getCurrentBlockExogenous()};
  vector<int> endogenous {evaluator.getCurrentBlockEndogenous()};
  for (auto & it : sconstrained_extended_path)
    {
      if (find(endogenous.begin(), endogenous.end(), it.exo_num) != endogenous.end()
          && find(exogenous.begin(), exogenous.end(), it.var_num) == exogenous.end())
        throw FatalException{"\nThe conditional forecast involving as constrained variable "
            + symbol_table.getName(SymbolType::endogenous, it.exo_num)
            + " and as endogenized exogenous " + symbol_table.getName(SymbolType::exogenous, it.var_num)
            + " that do not appear in block=" + to_string(current_block+1)
            + ")\nYou should not use block in model options"};
      else if (find(endogenous.begin(), endogenous.end(), it.exo_num) != endogenous.end()
               && find(exogenous.begin(), exogenous.end(), it.var_num) != exogenous.end()
               && (type == BlockSimulationType::evaluateForward
                   || type == BlockSimulationType::evaluateBackward))
        throw FatalException{"\nThe conditional forecast cannot be implemented for the block="
            + to_string(current_block+1) + ") that has to be evaluated instead to be solved\nYou should not use block in model options"};
      else if (find(previous_block_exogenous.begin(), previous_block_exogenous.end(), it.var_num)
               != previous_block_exogenous.end())
        throw FatalException{"\nThe conditional forecast involves in the block "
            + to_string(current_block+1) + " the endogenized exogenous "
            + symbol_table.getName(SymbolType::exogenous, it.var_num)
            + " that appear also in a previous block\nYou should not use block in model options"};
    }
  for (auto it : exogenous)
    previous_block_exogenous.push_back(it);
}

pair<bool, vector<int>>
Interpreter::MainLoop(const string &bin_basename, bool evaluate, int block, bool constrained, const vector<s_plan> &sconstrained_extended_path, const vector_table_conditional_local_type &vector_table_conditional_local)
{
  initializeTemporaryTerms(global_temporary_terms);

  int nb_blocks {evaluator.get_block_number()};

  if (block >= nb_blocks)
    throw FatalException {"Interpreter::MainLoop: Input argument block = " + to_string(block+1)
        + " is greater than the number of blocks in the model ("
        + to_string(nb_blocks) + " see M_.block_structure" + (steady_state ? "_stat" : "") + ".block)"};

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
        residual = vector<double>(y_size*periods);
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
        check_for_controlled_exo_validity(current_block, sconstrained_extended_path);
      if (print)
        {
          if (steady_state)
            residual = vector<double>(size);
          else
            residual = vector<double>(size*periods);
          evaluator.printCurrentBlock();
        }
      else if (evaluate)
        {
#ifdef DEBUG
          mexPrintf("jacobian_block=mxCreateDoubleMatrix(%d, %d, mxREAL)\n", size, getCurrentBlockNbColJacob());
#endif
          jacobian_block[current_block] = mxCreateDoubleMatrix(size, evaluator.getCurrentBlockNbColJacob(), mxREAL);
          if (!steady_state)
            {
#ifdef DEBUG
              mexPrintf("allocates jacobian_exo_block( %d, %d, mxREAL)\n", size, evaluator.getCurrentBlockExoSize());
              mexPrintf("(0) Allocating Jacobian\n");
#endif

              jacobian_exo_block[current_block] = mxCreateDoubleMatrix(size, evaluator.getCurrentBlockExoSize(), mxREAL);
              jacobian_det_exo_block[current_block] = mxCreateDoubleMatrix(size, evaluator.getCurrentBlockExoDetSize(), mxREAL);
            }
          if (block >= 0)
            {
              if (steady_state)
                residual = vector<double>(size);
              else
                residual = vector<double>(size*periods);
            }
          evaluate_a_block(true, block >= 0, bin_basename);
        }
      else
        {
#ifdef DEBUG
          mexPrintf("endo in block %d, size=%d, type=%d, steady_state=%d, is_linear=%d, endo_nbr=%d, u_count_int=%d\n",
                    current_block+1, size, type, steady_state, is_linear, symbol_table_endo_nbr, u_count_int);
#endif
          bool result;
          if (sconstrained_extended_path.size())
            {
              jacobian_block[current_block] = mxCreateDoubleMatrix(size, evaluator.getCurrentBlockNbColJacob(), mxREAL);
              jacobian_exo_block[current_block] = mxCreateDoubleMatrix(size, evaluator.getCurrentBlockExoSize(), mxREAL);
              jacobian_det_exo_block[current_block] = mxCreateDoubleMatrix(size, evaluator.getCurrentBlockExoDetSize(), mxREAL);
              residual = vector<double>(size*periods);
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

pair<bool, vector<int>>
Interpreter::extended_path(const string &file_name, bool evaluate, int block, int nb_periods, const vector<s_plan> &sextended_path, const vector<s_plan> &sconstrained_extended_path, const vector<string> &dates, const table_conditional_global_type &table_conditional_global)
{
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
  int old_verbosity {verbosity};
  verbosity = 0;
  ostringstream res1;
  res1 << std::scientific << 2.54656875434865131;
  int real_max_length = res1.str().length();
  int date_length = dates[0].length();
  int table_length = 2 + date_length + 3 + endo_name_length_l + 3 + real_max_length + 3 + 3 + 2 + 6 + 2;
  string line;
  line.insert(line.begin(), table_length, '-');
  line.insert(line.length(), "\n");
  if (old_verbosity >= 1)
    {
      mexPrintf("\nExtended Path simulation:\n");
      mexPrintf("-------------------------\n");
      mexPrintf(line.c_str());
      string title = "|" + elastic("date", date_length+2, false) + "|" + elastic("variable", endo_name_length_l+2, false) + "|" + elastic("max. value", real_max_length+2, false) + "| iter. |" + elastic("cvg", 5, false) + "|\n";
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
          mexPrintf("|%s|", elastic(dates[t], date_length+2, false).c_str());
          mexEvalString("drawnow;");
        }
      for (const auto & it : sextended_path)
        x[y_kmin + (it.exo_num - 1) * nb_row_x] = it.value[t];

      vector_table_conditional_local.clear();
      if (auto it = table_conditional_global.find(t); it != table_conditional_global.end())
        vector_table_conditional_local = it->second;
      tie(r, blocks) = MainLoop(file_name, evaluate, block, true, sconstrained_extended_path, vector_table_conditional_local);
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
          mexPrintf("%s|%s| %4d  |  x  |\n", elastic(symbol_table.getName(SymbolType::endogenous, max_res_idx), endo_name_length_l+2, true).c_str(), elastic(res1.str(), real_max_length+2, false).c_str(), iter);
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
  if (T && !global_temporary_terms)
    mxFree(T);
  return {true, blocks};
}

pair<bool, vector<int>>
Interpreter::compute_blocks(const string &file_name, bool evaluate, int block)
{
  //The big loop on intructions
  vector<s_plan> s_plan_junk;
  vector_table_conditional_local_type vector_table_conditional_local_junk;

  auto [r, blocks] = MainLoop(file_name, evaluate, block, false, s_plan_junk, vector_table_conditional_local_junk);

  if (T && !global_temporary_terms)
    mxFree(T);
  return {true, blocks};
}

void
Interpreter::initializeTemporaryTerms(bool global_temporary_terms)
{
  int ntt { evaluator.getNumberOfTemporaryTerms() };

  if (steady_state)
    {
      if (T)
        mxFree(T);
      if (global_temporary_terms)
        {
          if (!GlobalTemporaryTerms)
            {
              mexPrintf("GlobalTemporaryTerms is nullptr\n");
              mexEvalString("drawnow;");
            }
          if (ntt != static_cast<int>(mxGetNumberOfElements(GlobalTemporaryTerms)))
            GlobalTemporaryTerms = mxCreateDoubleMatrix(ntt, 1, mxREAL);
          T = mxGetPr(GlobalTemporaryTerms);
        }
      else
        {
          T = static_cast<double *>(mxMalloc(ntt*sizeof(double)));
          test_mxMalloc(T, __LINE__, __FILE__, __func__, ntt*sizeof(double));
        }
    }
  else
    {
      if (T)
        mxFree(T);
      T = static_cast<double *>(mxMalloc(ntt*(periods+y_kmin+y_kmax)*sizeof(double)));
      test_mxMalloc(T, __LINE__, __FILE__, __func__, ntt*(periods+y_kmin+y_kmax)*sizeof(double));
    }
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

int
Interpreter::At_Row(int r, NonZeroElem **first) const
{
  *first = FNZE_R[r];
  return NbNZRow[r];
}

int
Interpreter::At_Pos(int r, int c, NonZeroElem **first) const
{
  *first = FNZE_R[r];
  while ((*first)->c_index != c)
    *first = (*first)->NZE_R_N;
  return NbNZRow[r];
}

int
Interpreter::At_Col(int c, NonZeroElem **first) const
{
  *first = FNZE_C[c];
  return NbNZCol[c];
}

int
Interpreter::At_Col(int c, int lag, NonZeroElem **first) const
{
  *first = FNZE_C[c];
  int i = 0;
  while ((*first)->lag_index != lag && *first)
    *first = (*first)->NZE_C_N;
  if (*first)
    {
      NonZeroElem *firsta = *first;
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
  return i;
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
Interpreter::Read_SparseMatrix(const string &file_name, int Size, int periods, int y_kmin, int y_kmax, bool two_boundaries, int stack_solve_algo, int solve_algo)
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
        throw FatalException{"In Read_SparseMatrix, " + binfile.string() + " cannot be opened"};
    }
  IM_i.clear();
  if (two_boundaries)
    {
      if (stack_solve_algo == 5)
        {
          for (int i = 0; i < u_count_init-Size; i++)
            {
              int val;
              SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char *>(&val), sizeof(val));
              IM_i[{ eq, var, lag }] = val;
            }
          for (int j = 0; j < Size; j++)
            IM_i[{ j, Size*(periods+y_kmax), 0 }] = j;
        }
      else if ((stack_solve_algo >= 0 && stack_solve_algo <= 4)
               || stack_solve_algo == 6)
        {
          for (int i = 0; i < u_count_init-Size; i++)
            {
              int val;
              SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char *>(&val), sizeof(val));
              IM_i[{ var - lag*Size, -lag, eq }] = val;
            }
          for (int j = 0; j < Size; j++)
            IM_i[{ Size*(periods+y_kmax), 0, j }] = j;
        }
      else
        throw FatalException{"Invalid value for solve_algo or stack_solve_algo"};
    }
  else
    {
      if ((stack_solve_algo == 5 && !steady_state) || (solve_algo == 5 && steady_state))
        {
          for (int i = 0; i < u_count_init; i++)
            {
              int val;
              SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char *>(&val), sizeof(val));
              IM_i[{ eq, var, lag }] = val;
            }
        }
      else if ((((stack_solve_algo >= 0 && stack_solve_algo <= 4)
                 || stack_solve_algo == 6) && !steady_state)
               || ((solve_algo >= 6 || solve_algo <= 8) && steady_state))
        {
          for (int i = 0; i < u_count_init; i++)
            {
              int val;
              SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char *>(&val), sizeof(val));
              IM_i[{ var - lag*Size, -lag, eq }] = val;
            }
        }
      else
        throw FatalException{"Invalid value for solve_algo or stack_solve_algo"};
    }
  index_vara = static_cast<int *>(mxMalloc(Size*(periods+y_kmin+y_kmax)*sizeof(int)));
  test_mxMalloc(index_vara, __LINE__, __FILE__, __func__, Size*(periods+y_kmin+y_kmax)*sizeof(int));
  for (int j = 0; j < Size; j++)
    SaveCode.read(reinterpret_cast<char *>(&index_vara[j]), sizeof(*index_vara));
  if (periods+y_kmin+y_kmax > 1)
    for (int i = 1; i < periods+y_kmin+y_kmax; i++)
      for (int j = 0; j < Size; j++)
        index_vara[j+Size*i] = index_vara[j+Size*(i-1)] + y_size;
  index_equa = static_cast<int *>(mxMalloc(Size*sizeof(int)));
  test_mxMalloc(index_equa, __LINE__, __FILE__, __func__, Size*sizeof(int));
  for (int j = 0; j < Size; j++)
    SaveCode.read(reinterpret_cast<char *>(&index_equa[j]), sizeof(*index_equa));
}

void
Interpreter::Simple_Init(int Size, const map<tuple<int, int, int>, int> &IM, bool &zero_solution)
{
  pivot = static_cast<int *>(mxMalloc(Size*sizeof(int)));
  test_mxMalloc(pivot, __LINE__, __FILE__, __func__, Size*sizeof(int));
  pivot_save = static_cast<int *>(mxMalloc(Size*sizeof(int)));
  test_mxMalloc(pivot_save, __LINE__, __FILE__, __func__, Size*sizeof(int));
  pivotk = static_cast<int *>(mxMalloc(Size*sizeof(int)));
  test_mxMalloc(pivotk, __LINE__, __FILE__, __func__, Size*sizeof(int));
  pivotv = static_cast<double *>(mxMalloc(Size*sizeof(double)));
  test_mxMalloc(pivotv, __LINE__, __FILE__, __func__, Size*sizeof(double));
  pivotva = static_cast<double *>(mxMalloc(Size*sizeof(double)));
  test_mxMalloc(pivotva, __LINE__, __FILE__, __func__, Size*sizeof(double));
  b = static_cast<int *>(mxMalloc(Size*sizeof(int)));
  test_mxMalloc(b, __LINE__, __FILE__, __func__, Size*sizeof(int));
  line_done = static_cast<bool *>(mxMalloc(Size*sizeof(bool)));
  test_mxMalloc(line_done, __LINE__, __FILE__, __func__, Size*sizeof(bool));

  mem_mngr.init_CHUNK_BLCK_SIZE(u_count);
  int i = Size*sizeof(NonZeroElem *);
  FNZE_R = static_cast<NonZeroElem **>(mxMalloc(i));
  test_mxMalloc(FNZE_R, __LINE__, __FILE__, __func__, i);
  FNZE_C = static_cast<NonZeroElem **>(mxMalloc(i));
  test_mxMalloc(FNZE_C, __LINE__, __FILE__, __func__, i);
  auto **temp_NZE_R = static_cast<NonZeroElem **>(mxMalloc(i));
  test_mxMalloc(temp_NZE_R, __LINE__, __FILE__, __func__, i);
  auto **temp_NZE_C = static_cast<NonZeroElem **>(mxMalloc(i));
  test_mxMalloc(temp_NZE_C, __LINE__, __FILE__, __func__, i);
  i = Size*sizeof(int);
  NbNZRow = static_cast<int *>(mxMalloc(i));
  test_mxMalloc(NbNZRow, __LINE__, __FILE__, __func__, i);
  NbNZCol = static_cast<int *>(mxMalloc(i));
  test_mxMalloc(NbNZCol, __LINE__, __FILE__, __func__, i);
  for (i = 0; i < Size; i++)
    {
      line_done[i] = false;
      FNZE_C[i] = nullptr;
      FNZE_R[i] = nullptr;
      temp_NZE_C[i] = nullptr;
      temp_NZE_R[i] = nullptr;
      NbNZRow[i] = 0;
      NbNZCol[i] = 0;
    }
  int u_count1 = Size;
  for (auto &[key, value] : IM)
    {
      auto &[eq, var, lag] = key;
      if (lag == 0) /*Build the index for sparse matrix containing the jacobian : u*/
        {
          NbNZRow[eq]++;
          NbNZCol[var]++;
          NonZeroElem *first = mem_mngr.mxMalloc_NZE();
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
  for (int i = 0; i < Size; i++)
    {
      b[i] = i;
      cum_abs_sum += fabs(u[i]);
    }
  if (cum_abs_sum < 1e-20)
    zero_solution = true;
  else
    zero_solution = false;

  mxFree(temp_NZE_R);
  mxFree(temp_NZE_C);
  u_count = u_count1;
}

void
Interpreter::Init_Matlab_Sparse_Simple(int Size, const map<tuple<int, int, int>, int> &IM, const mxArray *A_m, const mxArray *b_m, bool &zero_solution, const mxArray *x0_m) const
{
  double *b = mxGetPr(b_m);
  if (!b)
    throw FatalException{"In Init_Matlab_Sparse_Simple, can't retrieve b vector"};
  double *x0 = mxGetPr(x0_m);
  if (!x0)
    throw FatalException{"In Init_Matlab_Sparse_Simple, can't retrieve x0 vector"};
  mwIndex *Ai = mxGetIr(A_m);
  if (!Ai)
    throw FatalException{"In Init_Matlab_Sparse_Simple, can't allocate Ai index vector"};
  mwIndex *Aj = mxGetJc(A_m);
  if (!Aj)
    throw FatalException{"In Init_Matlab_Sparse_Simple, can't allocate Aj index vector"};
  double *A = mxGetPr(A_m);
  if (!A)
    throw FatalException{"In Init_Matlab_Sparse_Simple, can't retrieve A matrix"};

  for (int i = 0; i < y_size*(periods+y_kmin); i++)
    ya[i] = y[i];
#ifdef DEBUG
  unsigned int max_nze = mxGetNzmax(A_m);
#endif
  unsigned int NZE = 0;
  int last_var = 0;
  double cum_abs_sum = 0;
  for (int i = 0; i < Size; i++)
    {
      b[i] = u[i];
      cum_abs_sum += fabs(b[i]);
      x0[i] = y[i];
    }
  if (cum_abs_sum < 1e-20)
    zero_solution = true;
  else
    zero_solution = false;

  Aj[0] = 0;
  last_var = 0;
  for (auto &[key, index] : IM)
    {
      auto &[var, ignore, eq] = key;
      if (var != last_var)
        {
          Aj[1+last_var] = NZE;
          last_var = var;
        }
#ifdef DEBUG
      if (index < 0 || index >= u_count_alloc || index > Size + Size*Size)
        throw FatalException{"In Init_Matlab_Sparse_Simple, index (" + to_string(index)
                             + ") out of range for u vector max = "
                             + to_string(Size+Size*Size)
                             + " allocated = " + to_string(u_count_alloc)};
      if (NZE >= max_nze)
        throw FatalException{"In Init_Matlab_Sparse_Simple, exceeds the capacity of A_m sparse matrix"};
#endif
      A[NZE] = u[index];
      Ai[NZE] = eq;
      NZE++;
#ifdef DEBUG
      if (eq < 0 || eq >= Size)
        throw FatalException{"In Init_Matlab_Sparse_Simple, index (" + to_string(eq)
                             + ") out of range for b vector"};
      if (var < 0 || var >= Size)
        throw FatalException{"In Init_Matlab_Sparse_Simple, index (" + to_string(var)
                             + ") out of range for index_vara vector"};
      if (index_vara[var] < 0 || index_vara[var] >= y_size)
        throw FatalException{"In Init_Matlab_Sparse_Simple, index ("
                             + to_string(index_vara[var])
                             + ") out of range for y vector max=" + to_string(y_size)
                             +" (0)"};
#endif
    }
  Aj[Size] = NZE;
}

void
Interpreter::Init_UMFPACK_Sparse_Simple(int Size, const map<tuple<int, int, int>, int> &IM, SuiteSparse_long **Ap, SuiteSparse_long **Ai, double **Ax, double **b, bool &zero_solution, const mxArray *x0_m) const
{
  *b = static_cast<double *>(mxMalloc(Size * sizeof(double)));
  test_mxMalloc(*b, __LINE__, __FILE__, __func__, Size * sizeof(double));
  if (!(*b))
    throw FatalException{"In Init_UMFPACK_Sparse, can't retrieve b vector"};
  double *x0 = mxGetPr(x0_m);
  if (!x0)
    throw FatalException{"In Init_UMFPACK_Sparse_Simple, can't retrieve x0 vector"};
  *Ap = static_cast<SuiteSparse_long *>(mxMalloc((Size+1) * sizeof(SuiteSparse_long)));
  test_mxMalloc(*Ap, __LINE__, __FILE__, __func__, (Size+1) * sizeof(SuiteSparse_long));
  if (!(*Ap))
    throw FatalException{"In Init_UMFPACK_Sparse, can't allocate Ap index vector"};
  size_t prior_nz = IM.size();
  *Ai = static_cast<SuiteSparse_long *>(mxMalloc(prior_nz * sizeof(SuiteSparse_long)));
  test_mxMalloc(*Ai, __LINE__, __FILE__, __func__, prior_nz * sizeof(SuiteSparse_long));
  if (!(*Ai))
    throw FatalException{"In Init_UMFPACK_Sparse, can't allocate Ai index vector"};
  *Ax = static_cast<double *>(mxMalloc(prior_nz * sizeof(double)));
  test_mxMalloc(*Ax, __LINE__, __FILE__, __func__, prior_nz * sizeof(double));
  if (!(*Ax))
    throw FatalException{"In Init_UMFPACK_Sparse, can't retrieve Ax matrix"};
  for (int i = 0; i < Size; i++)
    {
      int eq = index_vara[i];
      ya[eq+it_*y_size] = y[eq+it_*y_size];
    }
#ifdef DEBUG
  unsigned int max_nze = prior_nz; //mxGetNzmax(A_m);
#endif
  unsigned int NZE = 0;
  int last_var = 0;
  double cum_abs_sum = 0;
  for (int i = 0; i < Size; i++)
    {
      (*b)[i] = u[i];
      cum_abs_sum += fabs((*b)[i]);
      x0[i] = y[i];
    }
  if (cum_abs_sum < 1e-20)
    zero_solution = true;
  else
    zero_solution = false;

  (*Ap)[0] = 0;
  last_var = 0;
  for (auto &[key, index] : IM)
    {
      auto &[var, ignore, eq] = key;
      if (var != last_var)
        {
          (*Ap)[1+last_var] = NZE;
          last_var = var;
        }
#ifdef DEBUG
      if (index < 0 || index >= u_count_alloc || index > Size + Size*Size)
        throw FatalException{"In Init_Matlab_Sparse_Simple, index (" + to_string(index)
                             + ") out of range for u vector max = "
                             + to_string(Size+Size*Size)
                             + " allocated = " + to_string(u_count_alloc)};
      if (NZE >= max_nze)
        throw FatalException{"In Init_Matlab_Sparse_Simple, exceeds the capacity of A_m sparse matrix"};
#endif
      (*Ax)[NZE] = u[index];
      (*Ai)[NZE] = eq;
      NZE++;
#ifdef DEBUG
      if (eq < 0 || eq >= Size)
        throw FatalException{"In Init_Matlab_Sparse_Simple, index (" + to_string(eq)
                             + ") out of range for b vector"};
      if (var < 0 || var >= Size)
        throw FatalException{"In Init_Matlab_Sparse_Simple, index (" + to_string(var)
                             + ") out of range for index_vara vector"};
      if (index_vara[var] < 0 || index_vara[var] >= y_size)
        throw FatalException{"In Init_Matlab_Sparse_Simple, index ("
                             + to_string(index_vara[var])
                             + ") out of range for y vector max=" + to_string(y_size)
                             + " (0)"};
#endif
    }
  (*Ap)[Size] = NZE;
}

int
Interpreter::find_exo_num(const vector<s_plan> &sconstrained_extended_path, int value)
{
  auto it = find_if(sconstrained_extended_path.begin(), sconstrained_extended_path.end(),
                    [=](auto v) { return v.exo_num == value; });
  if (it != sconstrained_extended_path.end())
    return it - sconstrained_extended_path.begin();
  else
    return -1;
}

int
Interpreter::find_int_date(const vector<pair<int, double>> &per_value, int value)
{
  auto it = find_if(per_value.begin(), per_value.end(), [=](auto v) { return v.first == value; });
  if (it != per_value.end())
    return it - per_value.begin();
  else
    return -1;
}

void
Interpreter::Init_UMFPACK_Sparse(int periods, int y_kmin, int y_kmax, int Size, const map<tuple<int, int, int>, int> &IM, SuiteSparse_long **Ap, SuiteSparse_long **Ai, double **Ax, double **b, const mxArray *x0_m, const vector_table_conditional_local_type &vector_table_conditional_local, int block_num) const
{
  int n = periods * Size;
  *b = static_cast<double *>(mxMalloc(n * sizeof(double)));
  if (!(*b))
    throw FatalException{"In Init_UMFPACK_Sparse, can't retrieve b vector"};
  double *x0 = mxGetPr(x0_m);
  if (!x0)
    throw FatalException{"In Init_UMFPACK_Sparse_Simple, can't retrieve x0 vector"};
  *Ap = static_cast<SuiteSparse_long *>(mxMalloc((n+1) * sizeof(SuiteSparse_long)));
  test_mxMalloc(*Ap, __LINE__, __FILE__, __func__, (n+1) * sizeof(SuiteSparse_long));
  if (!(*Ap))
    throw FatalException{"In Init_UMFPACK_Sparse, can't allocate Ap index vector"};
  size_t prior_nz = IM.size() * periods;
  *Ai = static_cast<SuiteSparse_long *>(mxMalloc(prior_nz * sizeof(SuiteSparse_long)));
  test_mxMalloc(*Ai, __LINE__, __FILE__, __func__, prior_nz * sizeof(SuiteSparse_long));
  if (!(*Ai))
    throw FatalException{"In Init_UMFPACK_Sparse, can't allocate Ai index vector"};
  *Ax = static_cast<double *>(mxMalloc(prior_nz * sizeof(double)));
  test_mxMalloc(*Ax, __LINE__, __FILE__, __func__, prior_nz * sizeof(double));
  if (!(*Ax))
    throw FatalException{"In Init_UMFPACK_Sparse, can't retrieve Ax matrix"};
  for (int i = 0; i < y_size*(periods+y_kmin); i++)
    ya[i] = y[i];
  unsigned int NZE = 0;
  int last_var = 0;
  for (int i = 0; i < periods*Size; i++)
    {
      (*b)[i] = 0;
      x0[i] = y[index_vara[Size*y_kmin+i]];
    }
  double *jacob_exo;
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
  (*Ap)[0] = 0;
  for (int t = 0; t < periods; t++)
    {
      last_var = -1;
      int var = 0;
      for (auto &[key, value] : IM)
        {
          var = get<0>(key);
          int eq = get<2>(key)+Size*t;
          int lag = -get<1>(key);
          int index = value + (t-lag) * u_count_init;
          if (var != last_var)
            {
              (*Ap)[1+last_var + t * Size] = NZE;
              last_var = var;
              if (var < Size*(periods+y_kmax))
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
              if (t == 0 && var < (periods+y_kmax)*Size
                  && lag == 0 && vector_table_conditional_local.size())
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
                          if (jacob_exo[k + row_x*flip_exo] != 0)
                            {
                              (*Ax)[NZE] = jacob_exo[k + row_x*flip_exo];
                              (*Ai)[NZE] = k;
                              NZE++;

#ifdef DEBUG
                              if (local_index < 0 || local_index >= Size * periods)
                                throw FatalException{"In Init_UMFPACK_Sparse, index ("
                                                     + to_string(local_index)
                                                     + ") out of range for b vector"};
                              if (k + row_x*flip_exo < 0 || k + row_x*flip_exo >= row_x * col_x)
                                throw FatalException{"In Init_UMFPACK_Sparse, index ("
                                                     + to_string(var+Size*(y_kmin+t+lag))
                                                     + ") out of range for jacob_exo vector"};
                              if (t+y_kmin+flip_exo*nb_row_x < 0
                                  || t+y_kmin+flip_exo*nb_row_x >= nb_row_x * this->col_x)
                                throw FatalException{"In Init_UMFPACK_Sparse, index ("
                                                     + to_string(index_vara[var+Size*(y_kmin+t+lag)])
                                                     + ") out of range for x vector max="
                                                     + to_string(nb_row_x * this->col_x)};
#endif
                              u[k] -= jacob_exo[k + row_x*flip_exo] * x[t+y_kmin+flip_exo*nb_row_x];
                            }
                        }
                    }
                }
            }

          if (var < (periods+y_kmax)*Size)
            {
              int ti_y_kmin = -min(t, y_kmin);
              int ti_y_kmax = min(periods-(t+1), y_kmax);
              int ti_new_y_kmax = min(t, y_kmax);
              int ti_new_y_kmin = -min(periods-(t+1), y_kmin);
              if (lag <= ti_new_y_kmax && lag >= ti_new_y_kmin) /*Build the index for sparse matrix containing the jacobian : u*/
                {
#ifdef DEBUG
                  if (NZE >= prior_nz)
                    throw FatalException{"In Init_UMFPACK_Sparse, exceeds the capacity of allocated sparse matrix"};
#endif
                  if (!fliped)
                    {
                      (*Ax)[NZE] = u[index];
                      (*Ai)[NZE] = eq - lag * Size;
                      NZE++;
                    }
                  else /*if (fliped)*/
                    {
#ifdef DEBUG
                      if (eq - lag * Size < 0 || eq  - lag * Size >= Size * periods)
                        throw FatalException{"In Init_UMFPACK_Sparse, index ("
                                             + to_string(eq  - lag * Size)
                                             + ") out of range for b vector"};
                      if (var+Size*(y_kmin+t) < 0
                          || var+Size*(y_kmin+t) >= Size*(periods+y_kmin+y_kmax))
                        throw FatalException{"In Init_UMFPACK_Sparse, index ("
                                             + to_string(var+Size*(y_kmin+t))
                                             + ") out of range for index_vara vector"};
                      if (index_vara[var+Size*(y_kmin+t)] < 0
                          || index_vara[var+Size*(y_kmin+t)] >= y_size*(periods+y_kmin+y_kmax))
                        throw FatalException{"In Init_UMFPACK_Sparse, index ("
                                             + to_string(index_vara[var+Size*(y_kmin+t)])
                                             + ") out of range for y vector max="
                                             + to_string(y_size*(periods+y_kmin+y_kmax))};
#endif
                      (*b)[eq - lag * Size] += u[index] * y[index_vara[var+Size*(y_kmin+t)]];
                    }

                }
              if (lag > ti_y_kmax || lag < ti_y_kmin)
                {
#ifdef DEBUG
                  if (eq < 0 || eq >= Size * periods)
                    throw FatalException{"In Init_UMFPACK_Sparse, index ("
                                         + to_string(eq)
                                         + ") out of range for b vector"};
                  if (var+Size*(y_kmin+t+lag) < 0
                      || var+Size*(y_kmin+t+lag) >= Size*(periods+y_kmin+y_kmax))
                    throw FatalException{"In Init_UMFPACK_Sparse, index ("
                                         + to_string(var+Size*(y_kmin+t+lag))
                                         + ") out of range for index_vara vector"};
                  if (index_vara[var+Size*(y_kmin+t+lag)] < 0
                      || index_vara[var+Size*(y_kmin+t+lag)] >= y_size*(periods+y_kmin+y_kmax))
                    throw FatalException{"In Init_UMFPACK_Sparse, index ("
                                         + to_string(index_vara[var+Size*(y_kmin+t+lag)])
                                         + ") out of range for y vector max="
                                         + to_string(y_size*(periods+y_kmin+y_kmax))};
#endif
                  (*b)[eq] += u[index+lag*u_count_init]*y[index_vara[var+Size*(y_kmin+t+lag)]];
                }
            }
          else /* ...and store it in the u vector*/
            {
#ifdef DEBUG
              if (index < 0 || index >= u_count_alloc)
                throw FatalException{"In Init_UMFPACK_Sparse, index (" + to_string(index)
                                     + ") out of range for u vector"};
              if (eq < 0 || eq >= (Size*periods))
                throw FatalException{"In Init_UMFPACK_Sparse, index (" + to_string(eq)
                                     + ") out of range for b vector"};
#endif
              (*b)[eq] += u[index];
            }
        }
    }
  (*Ap)[Size*periods] = NZE;
#ifdef DEBUG
  mexPrintf("*Ax = [");
  for (int i = 0; i < static_cast<int>(NZE); i++)
    mexPrintf("%f ", (*Ax)[i]);
  mexPrintf("]\n");

  mexPrintf("*Ap = [");
  for (int i = 0; i < n+1; i++)
    mexPrintf("%d ", (*Ap)[i]);
  mexPrintf("]\n");

  mexPrintf("*Ai = [");
  for (int i = 0; i < static_cast<int>(NZE); i++)
    mexPrintf("%d ", (*Ai)[i]);
  mexPrintf("]\n");
#endif
}

void
Interpreter::PrintM(int n, const double *Ax, const mwIndex *Ap, const mwIndex *Ai)
{
  int nnz = Ap[n];
  auto *A = static_cast<double *>(mxMalloc(n * n * sizeof(double)));
  test_mxMalloc(A, __LINE__, __FILE__, __func__, n * n * sizeof(double));
  fill_n(A, n * n, 0);
  int k = 0;
  for (int i = 0; i < n; i++)
    for (int j = Ap[i]; j < static_cast<int>(Ap[i + 1]); j++)
      {
        int row = Ai[j];
        A[row * n + i] = Ax[j];
        k++;
      }
  if (nnz != k)
    mexPrintf("Problem nnz(%d) != number of elements(%d)\n", nnz, k);
  mexPrintf("----------------------\n");
  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
        mexPrintf("%-6.3f ", A[i * n + j]);
      mexPrintf("\n");
    }
  mxFree(A);
}

void
Interpreter::Init_Matlab_Sparse(int periods, int y_kmin, int y_kmax, int Size, const map<tuple<int, int, int>, int> &IM, mxArray *A_m, mxArray *b_m, const mxArray *x0_m) const
{
  double *b = mxGetPr(b_m);

  if (!b)
    throw FatalException{"In Init_Matlab_Sparse, can't retrieve b vector"};
  double *x0 = mxGetPr(x0_m);
  if (!x0)
    throw FatalException{"In Init_Matlab_Sparse_Simple, can't retrieve x0 vector"};
  mwIndex *Aj = mxGetJc(A_m);
  if (!Aj)
    throw FatalException{"In Init_Matlab_Sparse, can't allocate Aj index vector"};
  mwIndex *Ai = mxGetIr(A_m);
  if (!Ai)
    throw FatalException{"In Init_Matlab_Sparse, can't allocate Ai index vector"};
  double *A = mxGetPr(A_m);
  if (!A)
    throw FatalException{"In Init_Matlab_Sparse, can't retrieve A matrix"};

  for (int i = 0; i < y_size*(periods+y_kmin); i++)
    ya[i] = y[i];
  unsigned int NZE = 0;
  int last_var = 0;
  for (int i = 0; i < periods*Size; i++)
    {
      b[i] = 0;
      x0[i] = y[index_vara[Size*y_kmin+i]];
    }
  Aj[0] = 0;
  for (int t = 0; t < periods; t++)
    {
      last_var = 0;
      for (auto &[key, value] : IM)
        {
          int var = get<0>(key);
          if (var != last_var)
            {
              Aj[1+last_var + t * Size] = NZE;
              last_var = var;
            }
          int eq = get<2>(key)+Size*t;
          int lag = -get<1>(key);
          int index = value + (t-lag)*u_count_init;
          if (var < (periods+y_kmax)*Size)
            {
              int ti_y_kmin = -min(t, y_kmin);
              int ti_y_kmax = min(periods-(t +1), y_kmax);
              int ti_new_y_kmax = min(t, y_kmax);
              int ti_new_y_kmin = -min(periods-(t+1), y_kmin);
              if (lag <= ti_new_y_kmax && lag >= ti_new_y_kmin) /*Build the index for sparse matrix containing the jacobian : u*/
                {
#ifdef DEBUG
                  if (index < 0 || index >= u_count_alloc || index > Size + Size*Size)
                    throw FatalException{"In Init_Matlab_Sparse, index (" + to_string(index)
                                         + ") out of range for u vector max = "
                                         + to_string(Size+Size*Size) + " allocated = "
                                         + to_string(u_count_alloc)};
                  if (NZE >= prior_nz)
                    throw FatalException{"In Init_Matlab_Sparse, exceeds the capacity of allocated sparse matrix"};
#endif
                  A[NZE] = u[index];
                  Ai[NZE] = eq - lag * Size;
                  NZE++;
                }
              if (lag > ti_y_kmax || lag < ti_y_kmin)
                {
#ifdef DEBUG
                  if (eq < 0 || eq >= Size * periods)
                    throw FatalException{"In Init_Matlab_Sparse, index (" + to_string(eq)
                                         + ") out of range for b vector"};
                  if (var+Size*(y_kmin+t+lag) < 0
                      || var+Size*(y_kmin+t+lag) >= Size*(periods+y_kmin+y_kmax))
                    throw FatalException{"In Init_Matlab_Sparse, index ("
                                         + to_string(var+Size*(y_kmin+t+lag))
                                         + ") out of range for index_vara vector"};
                  if (index_vara[var+Size*(y_kmin+t+lag)] < 0
                      || index_vara[var+Size*(y_kmin+t+lag)] >= y_size*(periods+y_kmin+y_kmax))
                    throw FatalException{"In Init_Matlab_Sparse, index ("
                                         + to_string(index_vara[var+Size*(y_kmin+t+lag)])
                                         + ") out of range for y vector max="
                                         + to_string(y_size*(periods+y_kmin+y_kmax))};
#endif
                  b[eq] += u[index+lag*u_count_init]*y[index_vara[var+Size*(y_kmin+t+lag)]];
                }
            }
          else /* ...and store it in the u vector*/
            {
#ifdef DEBUG
              if (index < 0 || index >= u_count_alloc)
                throw FatalException{"In Init_Matlab_Sparse, index (" + to_string(index)
                                     + ") out of range for u vector"};
              if (eq < 0 || eq >= (Size*periods))
                throw FatalException{"In Init_Matlab_Sparse, index (" + to_string(eq)
                                     + ") out of range for b vector"};
#endif
              b[eq] += u[index];
            }
        }
    }
  Aj[Size*periods] = NZE;
}

void
Interpreter::Init_GE(int periods, int y_kmin, int y_kmax, int Size, const map<tuple<int, int, int>, int> &IM)
{
  double tmp_b = 0.0;
  pivot = static_cast<int *>(mxMalloc(Size*periods*sizeof(int)));
  test_mxMalloc(pivot, __LINE__, __FILE__, __func__, Size*periods*sizeof(int));
  pivot_save = static_cast<int *>(mxMalloc(Size*periods*sizeof(int)));
  test_mxMalloc(pivot_save, __LINE__, __FILE__, __func__, Size*periods*sizeof(int));
  pivotk = static_cast<int *>(mxMalloc(Size*periods*sizeof(int)));
  test_mxMalloc(pivotk, __LINE__, __FILE__, __func__, Size*periods*sizeof(int));
  pivotv = static_cast<double *>(mxMalloc(Size*periods*sizeof(double)));
  test_mxMalloc(pivotv, __LINE__, __FILE__, __func__, Size*periods*sizeof(double));
  pivotva = static_cast<double *>(mxMalloc(Size*periods*sizeof(double)));
  test_mxMalloc(pivotva, __LINE__, __FILE__, __func__, Size*periods*sizeof(double));
  b = static_cast<int *>(mxMalloc(Size*periods*sizeof(int)));
  test_mxMalloc(b, __LINE__, __FILE__, __func__, Size*periods*sizeof(int));
  line_done = static_cast<bool *>(mxMalloc(Size*periods*sizeof(bool)));
  test_mxMalloc(line_done, __LINE__, __FILE__, __func__, Size*periods*sizeof(bool));
  mem_mngr.init_CHUNK_BLCK_SIZE(u_count);
  int i = (periods+y_kmax+1)*Size*sizeof(NonZeroElem *);
  FNZE_R = static_cast<NonZeroElem **>(mxMalloc(i));
  test_mxMalloc(FNZE_R, __LINE__, __FILE__, __func__, i);
  FNZE_C = static_cast<NonZeroElem **>(mxMalloc(i));
  test_mxMalloc(FNZE_C, __LINE__, __FILE__, __func__, i);
  auto **temp_NZE_R = static_cast<NonZeroElem **>(mxMalloc(i));
  test_mxMalloc(temp_NZE_R, __LINE__, __FILE__, __func__, i);
  auto **temp_NZE_C = static_cast<NonZeroElem **>(mxMalloc(i));
  test_mxMalloc(temp_NZE_C, __LINE__, __FILE__, __func__, i);
  i = (periods+y_kmax+1)*Size*sizeof(int);
  NbNZRow = static_cast<int *>(mxMalloc(i));
  test_mxMalloc(NbNZRow, __LINE__, __FILE__, __func__, i);
  NbNZCol = static_cast<int *>(mxMalloc(i));
  test_mxMalloc(NbNZCol, __LINE__, __FILE__, __func__, i);

  for (int i = 0; i < periods*Size; i++)
    {
      b[i] = 0;
      line_done[i] = false;
    }
  for (int i = 0; i < (periods+y_kmax+1)*Size; i++)
    {
      FNZE_C[i] = nullptr;
      FNZE_R[i] = nullptr;
      temp_NZE_C[i] = nullptr;
      temp_NZE_R[i] = nullptr;
      NbNZRow[i] = 0;
      NbNZCol[i] = 0;
    }
  int nnz = 0;
  //pragma omp parallel for ordered private(it4, ti_y_kmin, ti_y_kmax, eq, var, lag) schedule(dynamic)
  for (int t = 0; t < periods; t++)
    {
      int ti_y_kmin = -min(t, y_kmin);
      int ti_y_kmax = min(periods-(t+1), y_kmax);
      int eq = -1;
      //pragma omp ordered
      for (auto &[key, value] : IM)
        {
          int var = get<1>(key);
          if (eq != get<0>(key)+Size*t)
            tmp_b = 0;
          eq = get<0>(key)+Size*t;
          int lag = get<2>(key);
          if (var < (periods+y_kmax)*Size)
            {
              lag = get<2>(key);
              if (lag <= ti_y_kmax && lag >= ti_y_kmin) /*Build the index for sparse matrix containing the jacobian : u*/
                {
                  nnz++;
                  var += Size*t;
                  NbNZRow[eq]++;
                  NbNZCol[var]++;
                  NonZeroElem *first = mem_mngr.mxMalloc_NZE();
                  first->NZE_C_N = nullptr;
                  first->NZE_R_N = nullptr;
                  first->u_index = value+u_count_init*t;
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
              else /*Build the additive terms ooutside the simulation periods related to the first lags and the last leads...*/
                {
                  if (lag < ti_y_kmin)
                    tmp_b += u[value+u_count_init*t]*y[index_vara[var+Size*(y_kmin+t)]];
                  else
                    tmp_b += u[value+u_count_init*t]*y[index_vara[var+Size*(y_kmin+t)]];
                }
            }
          else /* ...and store it in the u vector*/
            {
              b[eq] = value+u_count_init*t;
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
          u_count_alloc += 5*u_count_alloc_save;
          u = static_cast<double *>(mxRealloc(u, u_count_alloc*sizeof(double)));
          if (!u)
            throw FatalException{"In Get_u, memory exhausted (realloc("
                                 + to_string(u_count_alloc*sizeof(double)) + "))"};
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
Interpreter::End_GE()
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
Interpreter::compare(int *save_op, int *save_opa, int *save_opaa, int beg_t, int periods, long nop4, int Size)
{
  long nop = nop4/2;
  double r = 0.0;
  bool OK = true;
  int *diff1 = static_cast<int *>(mxMalloc(nop*sizeof(int)));
  test_mxMalloc(diff1, __LINE__, __FILE__, __func__, nop*sizeof(int));
  int *diff2 = static_cast<int *>(mxMalloc(nop*sizeof(int)));
  test_mxMalloc(diff2, __LINE__, __FILE__, __func__, nop*sizeof(int));
  int max_save_ops_first = -1;
  long j = 0, i = 0;
  while (i < nop4 && OK)
    {
      t_save_op_s *save_op_s = reinterpret_cast<t_save_op_s *>(&save_op[i]);
      t_save_op_s *save_opa_s = reinterpret_cast<t_save_op_s *>(&save_opa[i]);
      t_save_op_s *save_opaa_s = reinterpret_cast<t_save_op_s *>(&save_opaa[i]);
      diff1[j] = save_op_s->first-save_opa_s->first;
      max_save_ops_first = max(max_save_ops_first, save_op_s->first+diff1[j]*(periods-beg_t));
      switch (save_op_s->operat)
        {
        case IFLD:
        case IFDIV:
          OK = (save_op_s->operat == save_opa_s->operat && save_opa_s->operat == save_opaa_s->operat
                && diff1[j] == (save_opa_s->first-save_opaa_s->first));
          i += 2;
          break;
        case IFLESS:
        case IFSUB:
          diff2[j] = save_op_s->second-save_opa_s->second;
          OK = (save_op_s->operat == save_opa_s->operat && save_opa_s->operat == save_opaa_s->operat
                && diff1[j] == (save_opa_s->first-save_opaa_s->first)
                && diff2[j] == (save_opa_s->second-save_opaa_s->second));
          i += 3;
          break;
        default:
          throw FatalException{"In compare, unknown operator = "
                               + to_string(save_op_s->operat)};
        }
      j++;
    }
  // the same pivot for all remaining periods
  if (OK)
    {
      for (int i = beg_t; i < periods; i++)
        for (int j = 0; j < Size; j++)
          pivot[i*Size+j] = pivot[(i-1)*Size+j]+Size;
      if (max_save_ops_first >= u_count_alloc)
        {
          u_count_alloc += max_save_ops_first;
          u = static_cast<double *>(mxRealloc(u, u_count_alloc*sizeof(double)));
          if (!u)
            throw FatalException{"In compare, memory exhausted (realloc("
                                 + to_string(u_count_alloc*sizeof(double)) + "))"};
        }
      for (int t = 1; t < periods-beg_t-y_kmax; t++)
        {
          int i = j = 0;
          while (i < nop4)
            {
              auto *save_op_s = reinterpret_cast<t_save_op_s *>(&save_op[i]);
              double *up = &u[save_op_s->first+t*diff1[j]];
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
                  *up -= u[save_op_s->second+t*diff2[j]]*r;;
                  i += 3;
                  break;
                case IFLESS:
                  *up = -u[save_op_s->second+t*diff2[j]]*r;
                  i += 3;
                  break;
                }
              j++;
            }
        }
      int t1 = max(1, periods-beg_t-y_kmax);
      int periods_beg_t = periods-beg_t;
      for (int t = t1; t < periods_beg_t; t++)
        {
          int i = j = 0;
          int gap = periods_beg_t-t;
          while (i < nop4)
            {
              if (auto *save_op_s = reinterpret_cast<t_save_op_s *>(&save_op[i]);
                  save_op_s->lag < gap)
                {
                  double *up = &u[save_op_s->first+t*diff1[j]];
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
                      *up -= u[save_op_s->second+t*diff2[j]]*r;
                      i += 3;
                      break;
                    case IFLESS:
                      *up = -u[save_op_s->second+t*diff2[j]]*r;
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
Interpreter::complete(int beg_t, int Size, int periods, int *b)
{
  double yy = 0.0;

  int size_of_save_code = (1+y_kmax)*Size*(Size+1+4)/2*4;
  int *save_code = static_cast<int *>(mxMalloc(size_of_save_code*sizeof(int)));
  test_mxMalloc(save_code, __LINE__, __FILE__, __func__, size_of_save_code*sizeof(int));
  int size_of_diff = (1+y_kmax)*Size*(Size+1+4);
  int *diff = static_cast<int *>(mxMalloc(size_of_diff*sizeof(int)));
  test_mxMalloc(diff, __LINE__, __FILE__, __func__, size_of_diff*sizeof(int));
  long cal_y = y_size*y_kmin;

  long i = (beg_t+1)*Size-1;
  long nop = 0;
  for (long j = i; j > i-Size; j--)
    {
      long pos = pivot[j];
      NonZeroElem *first;
      long nb_var = At_Row(pos, &first);
      first = first->NZE_R_N;
      nb_var--;
      save_code[nop] = IFLDZ;
      save_code[nop+1] = 0;
      save_code[nop+2] = 0;
      save_code[nop+3] = 0;
#ifdef DEBUG
      if ((nop+3) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
#endif
      nop += 4;
      for (long k = 0; k < nb_var; k++)
        {
          save_code[nop] = IFMUL;
          save_code[nop+1] = index_vara[first->c_index]+cal_y;
          save_code[nop+2] = first->u_index;
          save_code[nop+3] = first->lag_index;
#ifdef DEBUG
          if ((nop+3) >= size_of_save_code)
            mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
#endif
          nop += 4;
          first = first->NZE_R_N;
        }
      save_code[nop] = IFADD;
      save_code[nop+1] = b[pos];
      save_code[nop+2] = 0;
      save_code[nop+3] = 0;
#ifdef DEBUG
      if ((nop+3) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
#endif
      nop += 4;
      save_code[nop] = IFSTP;
      save_code[nop+1] = index_vara[j]+y_size*y_kmin;
      save_code[nop+2] = 0;
      save_code[nop+3] = 0;
#ifdef DEBUG
      if ((nop+2) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
#endif
      nop += 4;
    }
  i = beg_t*Size-1;
  long nop1 = 0, nopa = 0;
  for (long j = i; j > i-Size; j--)
    {
      long pos = pivot[j];
      NonZeroElem *first;
      long nb_var = At_Row(pos, &first);
      first = first->NZE_R_N;
      nb_var--;
      diff[nopa] = 0;
      diff[nopa+1] = 0;
      nopa += 2;
      nop1 += 4;
      for (long k = 0; k < nb_var; k++)
        {
          diff[nopa] = save_code[nop1+1]-(index_vara[first->c_index]+cal_y);
          diff[nopa+1] = save_code[nop1+2]-(first->u_index);
#ifdef DEBUG
          if ((nop1+2) >= size_of_save_code)
            mexPrintf("out of save_code[%d] (bound=%d)\n", nop1+2, size_of_save_code);
          if ((nopa+1) >= size_of_diff)
            mexPrintf("out of diff[%d] (bound=%d)\n", nopa+2, size_of_diff);
#endif
          nopa += 2;
          nop1 += 4;
          first = first->NZE_R_N;
        }
      diff[nopa] = save_code[nop1+1]-(b[pos]);
      diff[nopa+1] = 0;
#ifdef DEBUG
      if ((nop1+3) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop1+2, size_of_save_code);
      if ((nopa+1) >= size_of_diff)
        mexPrintf("out of diff[%d] (bound=%d)\n", nopa+2, size_of_diff);
#endif
      nopa += 2;
      nop1 += 4;
      diff[nopa] = save_code[nop1+1]-(index_vara[j]+y_size*y_kmin);
      diff[nopa+1] = 0;
#ifdef DEBUG
      if ((nop1+4) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop1+2, size_of_save_code);
      if ((nopa+1) >= size_of_diff)
        mexPrintf("out of diff[%d] (bound=%d)\n", nopa+2, size_of_diff);
#endif
      nopa += 2;
      nop1 += 4;
    }
  long max_var = (periods+y_kmin)*y_size;
  long min_var = y_kmin*y_size;
  for (int t = periods+y_kmin-1; t >= beg_t+y_kmin; t--)
    {
      int j = 0, k;
      int ti = t-y_kmin-beg_t;
      for (int i = 0; i < nop; i += 4)
        {
          switch (save_code[i])
            {
            case IFLDZ:
              yy = 0;
              break;
            case IFMUL:
              k = save_code[i+1]+ti*diff[j];
              if (k < max_var && k > min_var)
                yy += y[k]*u[save_code[i+2]+ti*diff[j+1]];
              break;
            case IFADD:
              yy = -(yy+u[save_code[i+1]+ti*diff[j]]);
              break;
            case IFSTP:
              k = save_code[i+1]+ti*diff[j];
              double err = yy - y[k];
              y[k] += slowc*(err);
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
Interpreter::bksub(int tbreak, int last_period, int Size, double slowc_l)
{
  for (int i = 0; i < y_size*(periods+y_kmin); i++)
    y[i] = ya[i];
  if (symbolic && tbreak)
    last_period = complete(tbreak, Size, periods, b);
  else
    last_period = periods;
  for (int t = last_period+y_kmin-1; t >= y_kmin; t--)
    {
      int ti = (t-y_kmin)*Size;
      int cal = y_kmin*Size;
      int cal_y = y_size*y_kmin;
      for (int i = ti-1; i >= ti-Size; i--)
        {
          int j = i+cal;
          int pos = pivot[i+Size];
          NonZeroElem *first;
          int nb_var = At_Row(pos, &first);
          first = first->NZE_R_N;
          nb_var--;
          int eq = index_vara[j]+y_size;
          double yy = 0;
          for (int k = 0; k < nb_var; k++)
            {
              yy += y[index_vara[first->c_index]+cal_y]*u[first->u_index];
              first = first->NZE_R_N;
            }
          yy = -(yy+y[eq]+u[b[pos]]);
          direction[eq] = yy;
          y[eq] += slowc_l*yy;
        }
    }
}

void
Interpreter::simple_bksub(int it_, int Size, double slowc_l)
{
  for (int i = 0; i < y_size; i++)
    y[i+it_*y_size] = ya[i+it_*y_size];
  for (int i = Size-1; i >= 0; i--)
    {
      int pos = pivot[i];
      NonZeroElem *first;
      int nb_var = At_Row(pos, &first);
      first = first->NZE_R_N;
      nb_var--;
      int eq = index_vara[i];
      double yy = 0;
      for (int k = 0; k < nb_var; k++)
        {
          yy += y[index_vara[first->c_index]+it_*y_size]*u[first->u_index];
          first = first->NZE_R_N;
        }
      yy = -(yy+y[eq+it_*y_size]+u[b[pos]]);
      direction[eq+it_*y_size] = yy;
      y[eq+it_*y_size] += slowc_l*yy;
    }
}

mxArray *
Interpreter::subtract_A_B(const mxArray *A_m, const mxArray *B_m)
{
  size_t n_A = mxGetN(A_m);
  size_t m_A = mxGetM(A_m);
  double *A_d = mxGetPr(A_m);
  size_t n_B = mxGetN(B_m);
  double *B_d = mxGetPr(B_m);
  mxArray *C_m = mxCreateDoubleMatrix(m_A, n_B, mxREAL);
  double *C_d = mxGetPr(C_m);
  for (int j = 0; j < static_cast<int>(n_A); j++)
    for (unsigned int i = 0; i < m_A; i++)
      {
        size_t index = j*m_A+i;
        C_d[index] = A_d[index] - B_d[index];
      }
  return C_m;
}

mxArray *
Interpreter::Sparse_subtract_SA_SB(const mxArray *A_m, const mxArray *B_m)
{
  size_t n_A = mxGetN(A_m);
  size_t m_A = mxGetM(A_m);
  mwIndex *A_i = mxGetIr(A_m);
  mwIndex *A_j = mxGetJc(A_m);
  size_t total_nze_A = A_j[n_A];
  double *A_d = mxGetPr(A_m);
  size_t n_B = mxGetN(B_m);
  mwIndex *B_i = mxGetIr(B_m);
  mwIndex *B_j = mxGetJc(B_m);
  size_t total_nze_B = B_j[n_B];
  double *B_d = mxGetPr(B_m);
  mxArray *C_m = mxCreateSparse(m_A, n_B, m_A*n_B, mxREAL);
  mwIndex *C_i = mxGetIr(C_m);
  mwIndex *C_j = mxGetJc(C_m);
  double *C_d = mxGetPr(C_m);
  unsigned int nze_B = 0, nze_C = 0, nze_A = 0;
  unsigned int A_col = 0, B_col = 0, C_col = 0;
  C_j[C_col] = 0;
  while (nze_A < total_nze_A || nze_B < total_nze_B)
    {
      while (nze_A >= static_cast<unsigned int>(A_j[A_col+1]) && (nze_A < total_nze_A))
        A_col++;
      size_t A_row = A_i[nze_A];
      while (nze_B >= static_cast<unsigned int>(B_j[B_col+1]) && (nze_B < total_nze_B))
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
              C_j[A_col+1] = nze_C++;
              C_col = A_col;
            }
          else if ((A_row < B_row && nze_A < total_nze_A) || nze_B == total_nze_B)
            {
              C_d[nze_C] = A_d[nze_A++];
              C_i[nze_C] = A_row;
              while (C_col < A_col)
                C_j[++C_col] = nze_C;
              C_j[A_col+1] = nze_C++;
              C_col = A_col;
            }
          else
            {
              C_d[nze_C] = -B_d[nze_B++];
              C_i[nze_C] = B_row;
              while (C_col < B_col)
                C_j[++C_col] = nze_C;
              C_j[B_col+1] = nze_C++;
              C_col = B_col;
            }
        }
      else if ((A_col < B_col && nze_A < total_nze_A) || nze_B == total_nze_B)
        {
          C_d[nze_C] = A_d[nze_A++];
          C_i[nze_C] = A_row;
          while (C_col < A_col)
            C_j[++C_col] = nze_C;
          C_j[A_col+1] = nze_C++;
          C_col = A_col;
        }
      else
        {
          C_d[nze_C] = -B_d[nze_B++];
          C_i[nze_C] = B_row;
          while (C_col < B_col)
            C_j[++C_col] = nze_C;
          C_j[B_col+1] = nze_C++;
          C_col = B_col;
        }
    }
  while (C_col < n_B)
    C_j[++C_col] = nze_C;
  mxSetNzmax(C_m, nze_C);
  return C_m;
}

mxArray *
Interpreter::mult_SAT_B(const mxArray *A_m, const mxArray *B_m)
{
  size_t n_A = mxGetN(A_m);
  mwIndex *A_i = mxGetIr(A_m);
  mwIndex *A_j = mxGetJc(A_m);
  double *A_d = mxGetPr(A_m);
  size_t n_B = mxGetN(B_m);
  double *B_d = mxGetPr(B_m);
  mxArray *C_m = mxCreateDoubleMatrix(n_A, n_B, mxREAL);
  double *C_d = mxGetPr(C_m);
  for (int j = 0; j < static_cast<int>(n_B); j++)
    for (unsigned int i = 0; i < n_A; i++)
      {
        double sum = 0;
        size_t nze_A = A_j[i];
        while (nze_A < static_cast<unsigned int>(A_j[i+1]))
          {
            size_t i_A = A_i[nze_A];
            sum += A_d[nze_A++] * B_d[i_A];
          }
        C_d[j*n_A+i] = sum;
      }
  return C_m;
}

mxArray *
Interpreter::Sparse_mult_SAT_B(const mxArray *A_m, const mxArray *B_m)
{
  size_t n_A = mxGetN(A_m);
  mwIndex *A_i = mxGetIr(A_m);
  mwIndex *A_j = mxGetJc(A_m);
  double *A_d = mxGetPr(A_m);
  size_t n_B = mxGetN(B_m);
  size_t m_B = mxGetM(B_m);
  double *B_d = mxGetPr(B_m);
  mxArray *C_m = mxCreateSparse(n_A, n_B, n_A*n_B, mxREAL);
  mwIndex *C_i = mxGetIr(C_m);
  mwIndex *C_j = mxGetJc(C_m);
  double *C_d = mxGetPr(C_m);
  unsigned int nze_C = 0;
  //unsigned int nze_A = 0;
  unsigned int C_col = 0;
  C_j[C_col] = 0;
  //#pragma omp parallel for
  for (unsigned int j = 0; j < n_B; j++)
    for (unsigned int i = 0; i < n_A; i++)
      {
        double sum = 0;
        size_t nze_A = A_j[i];
        while (nze_A < static_cast<unsigned int>(A_j[i+1]))
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

mxArray *
Interpreter::Sparse_mult_SAT_SB(const mxArray *A_m, const mxArray *B_m)
{
  size_t n_A = mxGetN(A_m);
  mwIndex *A_i = mxGetIr(A_m);
  mwIndex *A_j = mxGetJc(A_m);
  double *A_d = mxGetPr(A_m);
  size_t n_B = mxGetN(B_m);
  mwIndex *B_i = mxGetIr(B_m);
  mwIndex *B_j = mxGetJc(B_m);
  double *B_d = mxGetPr(B_m);
  mxArray *C_m = mxCreateSparse(n_A, n_B, n_A*n_B, mxREAL);
  mwIndex *C_i = mxGetIr(C_m);
  mwIndex *C_j = mxGetJc(C_m);
  double *C_d = mxGetPr(C_m);
  size_t nze_B = 0, nze_C = 0, nze_A = 0;
  unsigned int C_col = 0;
  C_j[C_col] = 0;
  for (unsigned int j = 0; j < n_B; j++)
    for (unsigned int i = 0; i < n_A; i++)
      {
        double sum = 0;
        nze_B = B_j[j];
        nze_A = A_j[i];
        while (nze_A < static_cast<unsigned int>(A_j[i+1]) && nze_B < static_cast<unsigned int>(B_j[j+1]))
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

mxArray *
Interpreter::Sparse_transpose(const mxArray *A_m)
{
  size_t n_A = mxGetN(A_m);
  size_t m_A = mxGetM(A_m);
  mwIndex *A_i = mxGetIr(A_m);
  mwIndex *A_j = mxGetJc(A_m);
  size_t total_nze_A = A_j[n_A];
  double *A_d = mxGetPr(A_m);
  mxArray *C_m = mxCreateSparse(n_A, m_A, total_nze_A, mxREAL);
  mwIndex *C_i = mxGetIr(C_m);
  mwIndex *C_j = mxGetJc(C_m);
  double *C_d = mxGetPr(C_m);
  unsigned int nze_C = 0, nze_A = 0;
  fill_n(C_j, m_A+1, 0);
  map<pair<mwIndex, unsigned int>, double> B2;
  for (unsigned int i = 0; i < n_A; i++)
    while (nze_A < static_cast<unsigned int>(A_j[i+1]))
      {
        C_j[A_i[nze_A]+1]++;
        B2[{ A_i[nze_A], i }] = A_d[nze_A];
        nze_A++;
      }
  for (unsigned int i = 0; i < m_A; i++)
    C_j[i+1] += C_j[i];
  for (auto &[key, val] : B2)
    {
      C_d[nze_C] = val;
      C_i[nze_C++] = key.second;
    }
  return C_m;
}

void
Interpreter::compute_block_time(int Per_u_, bool evaluate, bool no_derivatives)
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
      evaluator.evaluateBlock(it_, y, ya, y_size, x, nb_row_x, params, steady_y, u, Per_u_, T, periods+y_kmin+y_kmax, TEF, TEFD, TEFDD, r, g1, jacob, jacob_exo, jacob_exo_det, evaluate, no_derivatives);
    }
  catch (FloatingPointException &e)
    {
      res1 = numeric_limits<double>::quiet_NaN();
      if (verbosity >= 2)
        mexPrintf("%s\n      %s\n", e.message.c_str(), e.location.c_str());
    }
}

bool
Interpreter::compute_complete(bool no_derivatives, double &_res1, double &_res2, double &_max_res, int &_max_res_idx)
{
  bool result;
  res1 = 0;
  compute_block_time(0, false, no_derivatives);
  if (!(isnan(res1) || isinf(res1)))
    {
      _res1 = 0;
      _res2 = 0;
      _max_res = 0;
      for (int i = 0; i < size; i++)
        {
          double rr;
          rr = r[i];
          if (max_res < fabs(rr))
            {
              _max_res = fabs(rr);
              _max_res_idx = i;
            }
          _res2 += rr*rr;
          _res1 += fabs(rr);
        }
      result = true;
    }
  else
    result = false;
  return result;
}

bool
Interpreter::compute_complete(double lambda, double *crit)
{
  double res1_ = 0, res2_ = 0, max_res_ = 0;
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
      if (compute_complete(true, res1, res2, max_res, max_res_idx))
        res2_ = res2;
      else
        return false;
    }
  else
    {
      for (int it = y_kmin; it < periods+y_kmin; it++)
        for (int i = 0; i < size; i++)
          {
            int eq = index_vara[i];
            y[eq+it*y_size] = ya[eq+it*y_size] + lambda * direction[eq+it*y_size];
          }
      for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
        {
          Per_u_ = (it_-y_kmin)*u_count_int;
          Per_y_ = it_*y_size;
          if (compute_complete(true, res1, res2, max_res, max_res_idx))
            {
              res2_ += res2;
              res1_ += res1;
              if (max_res > max_res_)
                {
                  max_res = max_res_;
                  max_res_idx = max_res_idx_;
                }
            }
          else
            return false;
        }
      it_ = periods+y_kmin-1; // Do not leave it_ in inconsistent state
    }
  if (verbosity >= 2)
    mexPrintf("  lambda=%e, res2=%e\n", lambda, res2_);
  *crit = res2_/2;
  return true;
}

bool
Interpreter::mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc)
{
  constexpr double GOLD = 1.618034;
  constexpr double GLIMIT = 100.0;
  constexpr double TINY = 1.0e-20;

  auto sign = [](double a, double b) { return b >= 0.0 ? fabs(a) : -fabs(a); };

  if (verbosity >= 2)
    mexPrintf("bracketing *ax=%f, *bx=%f\n", *ax, *bx);
  if (!compute_complete(*ax, fa))
    return false;
  if (!compute_complete(*bx, fb))
    return false;
  if (*fb > *fa)
    {
      swap(*ax, *bx);
      swap(*fa, *fb);
    }
  *cx = (*bx)+GOLD*(*bx-*ax);
  if (!compute_complete(*cx, fc))
    return false;
  while (*fb > *fc)
    {
      double r = (*bx-*ax)*(*fb-*fc);
      double q = (*bx-*cx)*(*fb-*fa);
      double u = (*bx)-((*bx-*cx)*q-(*bx-*ax)*r)
        /(2.0*sign(fmax(fabs(q-r), TINY), q-r));
      double ulim = (*bx)+GLIMIT*(*cx-*bx);
      double fu;
      if ((*bx-u)*(u-*cx) > 0.0)
        {
          if (!compute_complete(u, &fu))
            return false;
          if (fu < *fc)
            {
              *ax = (*bx);
              *bx = u;
              *fa = (*fb);
              *fb = fu;
              return true;
            }
          else if (fu > *fb)
            {
              *cx = u;
              *fc = fu;
              return true;
            }
          u = (*cx)+GOLD*(*cx-*bx);
          if (!compute_complete(u, &fu))
            return false;
        }
      else if ((*cx-u)*(u-ulim) > 0.0)
        {
          if (!compute_complete(u, &fu))
            return false;
          if (fu < *fc)
            {
              *bx = *cx;
              *cx = u;
              u = *cx+GOLD*(*cx-*bx);
              *fb = *fc;
              *fc = fu;
              if (!compute_complete(u, &fu))
                return false;
            }
        }
      else if ((u-ulim)*(ulim-*cx) >= 0.0)
        {
          u = ulim;
          if (!compute_complete(u, &fu))
            return false;
        }
      else
        {
          u = (*cx)+GOLD*(*cx-*bx);
          if (!compute_complete(u, &fu))
            return false;
        }
      *ax = *bx;
      *bx = *cx;
      *cx = u;
      *fa = *fb;
      *fb = *fc;
      *fc = fu;
    }
  return true;
}

bool
Interpreter::golden(double ax, double bx, double cx, double tol, double solve_tolf, double *xmin)
{
  const double R = 0.61803399;
  const double C = (1.0-R);
  if (verbosity >= 2)
    mexPrintf("golden\n");
  int iter = 0, max_iter = 100;
  double f1, f2, x1, x2;
  double x0 = ax;
  double x3 = cx;
  if (fabs(cx-bx) > fabs(bx-ax))
    {
      x1 = bx;
      x2 = bx+C*(cx-bx);
    }
  else
    {
      x2 = bx;
      x1 = bx-C*(bx-ax);
    }
  if (!compute_complete(x1, &f1))
    return false;
  if (!compute_complete(x2, &f2))
    return false;
  while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2)) && f1 > solve_tolf && f2 > solve_tolf
         && iter < max_iter && abs(x1 - x2) > 1e-4)
    {
      if (f2 < f1)
        {
          x0 = x1;
          x1 = x2;
          x2 = R*x1+C*x3;
          f1 = f2;
          if (!compute_complete(x2, &f2))
            return false;
        }
      else
        {
          x3 = x2;
          x2 = x1;
          x1 = R*x2+C*x0;
          f2 = f1;
          if (!compute_complete(x1, &f1))
            return false;
        }
      iter++;
    }
  if (f1 < f2)
    {
      *xmin = x1;
      return true;
    }
  else
    {
      *xmin = x2;
      return true;
    }
}

void
Interpreter::Solve_Matlab_Relaxation(mxArray *A_m, mxArray *b_m, unsigned int Size, double slowc_l)
{
  double *b_m_d = mxGetPr(b_m);
  if (!b_m_d)
    throw FatalException{"In Solve_Matlab_Relaxation, can't retrieve b_m vector"};
  mwIndex *A_m_i = mxGetIr(A_m);
  if (!A_m_i)
    throw FatalException{"In Solve_Matlab_Relaxation, can't retrieve Ir vector of matrix A"};
  mwIndex *A_m_j = mxGetJc(A_m);
  if (!A_m_j)
    throw FatalException{"In Solve_Matlab_Relaxation, can't retrieve Jc vectior of matrix A"};
  double *A_m_d = mxGetPr(A_m);
  if (!A_m_d)
    throw FatalException{"In Solve_Matlab_Relaxation, can't retrieve double float data of matrix A"};

  /* Extract submatrices from the upper-left corner of A and subvectors at the
     beginning of b, so that the system looks like:
     ⎛B1 C1 …⎞   ⎛b1⎞
     ⎢A2 B2 …⎥   ⎢b2⎥
     ⎢   A3 …⎥ = ⎢ …⎥
     ⎢    … …⎥   ⎢ …⎥
     ⎝      …⎠   ⎝ …⎠
  */
  mxArray *B1 = mxCreateSparse(Size, Size, Size*Size, mxREAL);
  mwIndex *B1_i = mxGetIr(B1);
  mwIndex *B1_j = mxGetJc(B1);
  double *B1_d = mxGetPr(B1);
  mxArray *C1 = mxCreateSparse(Size, Size, Size*Size, mxREAL);
  mwIndex *C1_i = mxGetIr(C1);
  mwIndex *C1_j = mxGetJc(C1);
  double *C1_d = mxGetPr(C1);
  mxArray *A2 = mxCreateSparse(Size, Size, Size*Size, mxREAL);
  mwIndex *A2_i = mxGetIr(A2);
  mwIndex *A2_j = mxGetJc(A2);
  double *A2_d = mxGetPr(A2);
  mxArray *B2 = mxCreateSparse(Size, Size, Size*Size, mxREAL);
  mwIndex *B2_i = mxGetIr(B2);
  mwIndex *B2_j = mxGetJc(B2);
  double *B2_d = mxGetPr(B2);
  mxArray *A3 = mxCreateSparse(Size, Size, Size*Size, mxREAL);
  mwIndex *A3_i = mxGetIr(A3);
  mwIndex *A3_j = mxGetJc(A3);
  double *A3_d = mxGetPr(A3);

  mxArray *b1 = mxCreateDoubleMatrix(Size, 1, mxREAL);
  double *b1_d = mxGetPr(b1);
  mxArray *b2 = mxCreateDoubleMatrix(Size, 1, mxREAL);
  double *b2_d = mxGetPr(b2);

  unsigned int nze = 0; // Counter in nonzero elements of A
  unsigned int B1_nze = 0, C1_nze = 0, A2_nze = 0, B2_nze = 0, A3_nze = 0; // Same for submatrices

  for (size_t var = 0; var < 2*Size; var++) // Iterate over columns of A
    {
      if (var < Size)
        {
          b1_d[var] = b_m_d[var];

          B1_j[var] = B1_nze;
          A2_j[var] = A2_nze;
        }
      else
        {
          b2_d[var - Size] = b_m_d[var];

          C1_j[var - Size] = C1_nze;
          B2_j[var - Size] = B2_nze;
          A3_j[var - Size] = A3_nze;
        }

      while (static_cast<unsigned int>(A_m_j[var+1]) > nze)
        {
          size_t eq = A_m_i[nze];
          if (var < Size)
            {
              if (eq < Size)
                {
                  B1_i[B1_nze] = eq;
                  B1_d[B1_nze] = A_m_d[nze];
                  B1_nze++;
                }
              else // Here we know that eq < 2*Size, because of the structure of A
                {
                  A2_i[A2_nze] = eq - Size;
                  A2_d[A2_nze] = A_m_d[nze];
                  A2_nze++;
                }
            }
          else if (var < 2*Size)
            {
              if (eq < Size)
                {
                  C1_i[C1_nze] = eq;
                  C1_d[C1_nze] = A_m_d[nze];
                  C1_nze++;
                }
              else if (eq < 2*Size)
                {
                  B2_i[B2_nze] = eq - Size;
                  B2_d[B2_nze] = A_m_d[nze];
                  B2_nze++;
                }
              else // Here we know that eq < 3*Size, because of the structure of A
                {
                  A3_i[A3_nze] = eq - 2*Size;
                  A3_d[A3_nze] = A_m_d[nze];
                  A3_nze++;
                }
            }
          nze++;
        }
    }
  B1_j[Size] = B1_nze;
  C1_j[Size] = C1_nze;
  A2_j[Size] = A2_nze;
  B2_j[Size] = B2_nze;
  A3_j[Size] = A3_nze;

  vector<pair<mxArray *, mxArray *>> triangular_form;
  mxArray *d1 = nullptr;
  for (int t = 1; t <= periods; t++)
    {
      mxArray *B1_inv;
      mexCallMATLAB(1, &B1_inv, 1, &B1, "inv");

      // Compute subvector d1 of the triangularized system.
      mxArray *B1_inv_t = Sparse_transpose(B1_inv);
      mxDestroyArray(B1_inv);
      d1 = mult_SAT_B(B1_inv_t, b1);

      /* Compute block S1 of the triangularized system.
         Update B1, C1, B2, A2, A3, b1 and b2 for the next relaxation iteration.
         Save S1 and d1 for the subsequent backward iteration.
         E.g. at the end of the first iteration, the system will look like:
         ⎛ I S1  … …⎞   ⎛d1⎞
         ⎢   B1 C1 …⎥   ⎢b1⎥
         ⎢   A2 B2 …⎥ = ⎢b2⎥
         ⎢      A3 …⎥   ⎢ …⎥
         ⎝         …⎠   ⎝ …⎠
      */
      if (t < periods)
        {
          // Compute S1
          mxArray *S1 = Sparse_mult_SAT_SB(B1_inv_t, C1);

          // Update A2, B1, b1
          mxArray *A2_t = Sparse_transpose(A2);
          mxArray *tmp = Sparse_mult_SAT_SB(A2_t, S1);
          mxDestroyArray(B1);
          B1 = Sparse_subtract_SA_SB(B2, tmp);
          mxDestroyArray(tmp);

          tmp = mult_SAT_B(A2_t, d1);
          mxDestroyArray(A2_t);
          mxDestroyArray(b1);
          b1 = subtract_A_B(b2, tmp);
          mxDestroyArray(tmp);

          mxDestroyArray(A2);
          A2 = mxDuplicateArray(A3);

          // Save S1 and d1
          triangular_form.emplace_back(S1, d1);
        }

      mxDestroyArray(B1_inv_t);

      if (t < periods - 1)
        {
          // Update C1, B2, A3, b2
          C1_nze = B2_nze = A3_nze = 0;
          for (size_t var = (t+1)*Size; var < (t+2)*Size; var++)
            {
              b2_d[var - (t+1)*Size] = b_m_d[var];

              C1_j[var - (t+1)*Size] = C1_nze;
              B2_j[var - (t+1)*Size] = B2_nze;
              A3_j[var - (t+1)*Size] = A3_nze;

              while (static_cast<unsigned int>(A_m_j[var+1]) > nze)
                {
                  size_t eq = A_m_i[nze];
                  if (eq < (t+1) * Size)
                    {
                      C1_i[C1_nze] = eq - t*Size;
                      C1_d[C1_nze] = A_m_d[nze];
                      C1_nze++;
                    }
                  else if (eq < (t+2)*Size)
                    {
                      B2_i[B2_nze] = eq - (t+1)*Size;
                      B2_d[B2_nze] = A_m_d[nze];
                      B2_nze++;
                    }
                  else
                    {
                      A3_i[A3_nze] = eq - (t+2)*Size;
                      A3_d[A3_nze] = A_m_d[nze];
                      A3_nze++;
                    }
                  nze++;
                }
            }
          C1_j[Size] = C1_nze;
          B2_j[Size] = B2_nze;
          A3_j[Size] = A3_nze;
        }
    }

  // At this point, d1 contains the solution for the last period
  double *d1_d = mxGetPr(d1);
  for (unsigned i = 0; i < Size; i++)
    {
      int eq = index_vara[i+Size*(y_kmin+periods-1)];
      double yy = -(d1_d[i] + y[eq]);
      direction[eq] = yy;
      y[eq] += slowc_l * yy;
    }

  // Perform backward iteration to compute the solution for other periods
  for (int t = periods-2; t >= 0; t--)
    {
      auto [S1, d1_next] = triangular_form.back();
      triangular_form.pop_back();
      mxArray *S1_t = Sparse_transpose(S1);
      mxDestroyArray(S1);
      mxArray *tmp = mult_SAT_B(S1_t, d1);
      mxDestroyArray(S1_t);
      mxDestroyArray(d1);
      d1 = subtract_A_B(d1_next, tmp);
      d1_d = mxGetPr(d1);
      mxDestroyArray(d1_next);
      mxDestroyArray(tmp);
      for (unsigned i = 0; i < Size; i++)
        {
          int eq = index_vara[i+Size*(y_kmin+t)];
          double yy = -(d1_d[i] + y[eq]);
          direction[eq] = yy;
          y[eq] += slowc_l * yy;
        }
    }

  mxDestroyArray(B1);
  mxDestroyArray(C1);
  mxDestroyArray(A2);
  mxDestroyArray(B2);
  mxDestroyArray(A3);
  mxDestroyArray(b1);
  mxDestroyArray(b2);
  mxDestroyArray(A_m);
  mxDestroyArray(b_m);
  mxDestroyArray(d1);
}

void
Interpreter::End_Matlab_LU_UMFPack()
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

void
Interpreter::End_Solver()
{
  if (((stack_solve_algo == 0 || stack_solve_algo == 4) && !steady_state)
      || (solve_algo == 6 && steady_state))
    End_Matlab_LU_UMFPack();
}

void
Interpreter::Printfull_UMFPack(const SuiteSparse_long *Ap, const SuiteSparse_long *Ai, const double *Ax, const double *b, int n)
{
  double A[n*n];
  for (int i = 0; i < n*n; i++)
    A[i] = 0;
  int k = 0;
  for (int i = 0; i < n; i++)
    for (int j = Ap[i]; j < Ap[i+1]; j++)
      A[Ai[j] * n + i] = Ax[k++];
  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
        mexPrintf("%4.1f ", A[i*n+j]);
      mexPrintf("     %6.3f\n", b[i]);
    }
}

void
Interpreter::Print_UMFPack(const SuiteSparse_long *Ap, const SuiteSparse_long *Ai, const double *Ax, int n)
{
  int k = 0;
  for (int i = 0; i < n; i++)
    for (int j = Ap[i]; j < Ap[i+1]; j++)
      mexPrintf("(%d, %d)    %f\n", Ai[j]+1, i+1, Ax[k++]);
}

void
Interpreter::Solve_LU_UMFPack(SuiteSparse_long *Ap, SuiteSparse_long *Ai, double *Ax, double *b, int n, int Size, double slowc_l, bool is_two_boundaries, int it_, const vector_table_conditional_local_type &vector_table_conditional_local)
{
  SuiteSparse_long sys = 0;
  double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO], res[n];

  umfpack_dl_defaults(Control);
  SuiteSparse_long status = 0;
  if (iter == 0)
    {
      if (Symbolic)
        umfpack_dl_free_symbolic(&Symbolic);
      status = umfpack_dl_symbolic(n, n, Ap, Ai, Ax, &Symbolic, Control, Info);
      if (status != UMFPACK_OK)
        {
          umfpack_dl_report_info(Control, Info);
          umfpack_dl_report_status(Control, status);
          throw FatalException{"umfpack_dl_symbolic failed"};
        }
    }
  if (Numeric)
    umfpack_dl_free_numeric(&Numeric);
  status = umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);
  if (status != UMFPACK_OK)
    {
      umfpack_dl_report_info(Control, Info);
      umfpack_dl_report_status(Control, status);
      throw FatalException{"umfpack_dl_numeric failed"};
    }
  status = umfpack_dl_solve(sys, Ap, Ai, Ax, res, b, Numeric, Control, Info);
  if (status != UMFPACK_OK)
    {
      umfpack_dl_report_info(Control, Info);
      umfpack_dl_report_status(Control, status);
      throw FatalException{"umfpack_dl_solve failed"};
    }

  if (vector_table_conditional_local.size())
    {
      if (is_two_boundaries)
        for (int t = 0; t < n / Size; t++)
          if (t == 0)
            {
              for (int i = 0; i < Size; i++)
                {
                  bool fliped = vector_table_conditional_local[i].is_cond;
                  if (fliped)
                    {
                      int eq = index_vara[i+Size*(y_kmin)];
                      int flip_exo = vector_table_conditional_local[i].var_exo;
                      double yy = -(res[i] + x[y_kmin + flip_exo*nb_row_x]);
                      direction[eq] = 0;
                      x[flip_exo*nb_row_x + y_kmin] += slowc_l * yy;
                    }
                  else
                    {
                      int eq = index_vara[i+Size*(y_kmin)];
                      double yy = -(res[i] + y[eq]);
                      direction[eq] = yy;
                      y[eq] += slowc_l * yy;
                    }
                }
            }
          else
            {
              for (int i = 0; i < Size; i++)
                {
                  int eq = index_vara[i+Size*(t + y_kmin)];
                  double yy = -(res[i + Size * t] + y[eq]);
                  direction[eq] = yy;
                  y[eq] += slowc_l * yy;
                }
            }
      else
        for (int i = 0; i < n; i++)
          {
            int eq = index_vara[i];
            double yy = -(res[i] + y[eq+it_*y_size]);
            direction[eq] = yy;
            y[eq+it_*y_size] += slowc_l * yy;
          }
    }
  else
    {
      if (is_two_boundaries)
        for (int i = 0; i < n; i++)
          {
            int eq = index_vara[i+Size*y_kmin];
            double yy = -(res[i] + y[eq]);
            direction[eq] = yy;
            y[eq] += slowc_l * yy;
          }
      else
        for (int i = 0; i < n; i++)
          {
            int eq = index_vara[i];
            double yy = -(res[i] + y[eq+it_*y_size]);
            direction[eq] = yy;
            y[eq+it_*y_size] += slowc_l * yy;
          }
    }

  mxFree(Ap);
  mxFree(Ai);
  mxFree(Ax);
  mxFree(b);
}

void
Interpreter::Solve_LU_UMFPack(SuiteSparse_long *Ap, SuiteSparse_long *Ai, double *Ax, double *b, int n, int Size, double slowc_l, bool is_two_boundaries, int it_)
{
  SuiteSparse_long sys = 0;
  double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO], res[n];

  umfpack_dl_defaults(Control);
  SuiteSparse_long status = 0;
  if (iter == 0)
    {
      if (Symbolic)
        umfpack_dl_free_symbolic(&Symbolic);
      status = umfpack_dl_symbolic(n, n, Ap, Ai, Ax, &Symbolic, Control, Info);
      if (status != UMFPACK_OK)
        {
          umfpack_dl_report_info(Control, Info);
          umfpack_dl_report_status(Control, status);
          throw FatalException{"umfpack_dl_symbolic failed"};
        }
    }
  if (Numeric)
    umfpack_dl_free_numeric(&Numeric);
  status = umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);
  if (status != UMFPACK_OK)
    {
      umfpack_dl_report_info(Control, Info);
      umfpack_dl_report_status(Control, status);
      throw FatalException{"umfpack_dl_numeric failed"};
    }
  status = umfpack_dl_solve(sys, Ap, Ai, Ax, res, b, Numeric, Control, Info);
  if (status != UMFPACK_OK)
    {
      umfpack_dl_report_info(Control, Info);
      umfpack_dl_report_status(Control, status);
      throw FatalException{"umfpack_dl_solve failed"};
    }

  if (is_two_boundaries)
    for (int i = 0; i < n; i++)
      {
        int eq = index_vara[i+Size*y_kmin];
        double yy = -(res[i] + y[eq]);
        direction[eq] = yy;
        y[eq] += slowc_l * yy;
      }
  else
    for (int i = 0; i < n; i++)
      {
        int eq = index_vara[i];
        double yy = -(res[i] + y[eq+it_*y_size]);
        direction[eq] = yy;
        y[eq+it_*y_size] += slowc_l * yy;
      }
  mxFree(Ap);
  mxFree(Ai);
  mxFree(Ax);
  mxFree(b);
}

void
Interpreter::Solve_Matlab_GMRES(mxArray *A_m, mxArray *b_m, int Size, double slowc, int block, bool is_two_boundaries, int it_, mxArray *x0_m)
{
  size_t n = mxGetM(A_m);
  const char *field_names[] = {"droptol", "type"};
  mwSize dims[1] = { 1 };
  mxArray *Setup = mxCreateStructArray(1, dims, std::extent_v<decltype(field_names)>, field_names);
  mxSetFieldByNumber(Setup, 0, 0, mxCreateDoubleScalar(lu_inc_tol));
  mxSetFieldByNumber(Setup, 0, 1, mxCreateString("ilutp"));
  mxArray *lhs0[2];
  mxArray *rhs0[] = { A_m, Setup };
  if (mexCallMATLAB(std::extent_v<decltype(lhs0)>, lhs0, std::extent_v<decltype(rhs0)>, rhs0, "ilu"))
    throw FatalException("In GMRES, the incomplete LU decomposition (ilu) has failed");
  mxArray *L1 = lhs0[0];
  mxArray *U1 = lhs0[1];
  /*[za,flag1] = gmres(g1a,b,Blck_size,1e-6,Blck_size*periods,L1,U1);*/
  mxArray *rhs[] = { A_m, b_m, mxCreateDoubleScalar(Size), mxCreateDoubleScalar(1e-6),
    mxCreateDoubleScalar(static_cast<double>(n)), L1, U1, x0_m };
  mxArray *lhs[2];
  mexCallMATLAB(std::extent_v<decltype(lhs)>, lhs, std::extent_v<decltype(rhs)>, rhs, "gmres");
  mxArray *z = lhs[0];
  mxArray *flag = lhs[1];
  double *flag1 = mxGetPr(flag);
  mxDestroyArray(rhs0[1]);
  mxDestroyArray(rhs[2]);
  mxDestroyArray(rhs[3]);
  mxDestroyArray(rhs[4]);
  mxDestroyArray(rhs[5]);
  mxDestroyArray(rhs[6]);
  if (*flag1 > 0)
    {
      if (*flag1 == 1)
        mexWarnMsgTxt(("Error in bytecode: No convergence inside GMRES, in block "
                       + to_string(block+1)).c_str());
      else if (*flag1 == 2)
        mexWarnMsgTxt(("Error in bytecode: Preconditioner is ill-conditioned, in block "
                       + to_string(block+1)).c_str());
      else if (*flag1 == 3)
        mexWarnMsgTxt(("Error in bytecode: GMRES stagnated (Two consecutive iterates were the same.), in block "
                       + to_string(block+1)).c_str());
      lu_inc_tol /= 10;
    }
  else
    {
      double *res = mxGetPr(z);
      if (is_two_boundaries)
        for (int i = 0; i < static_cast<int>(n); i++)
          {
            int eq = index_vara[i+Size*y_kmin];
            double yy = -(res[i] + y[eq]);
            direction[eq] = yy;
            y[eq] += slowc * yy;
          }
      else
        for (int i = 0; i < static_cast<int>(n); i++)
          {
            int eq = index_vara[i];
            double yy = -(res[i] + y[eq+it_*y_size]);
            direction[eq] = yy;
            y[eq+it_*y_size] += slowc * yy;
          }
    }
  mxDestroyArray(A_m);
  mxDestroyArray(b_m);
  mxDestroyArray(x0_m);
  mxDestroyArray(z);
  mxDestroyArray(flag);
}

void
Interpreter::Solve_Matlab_BiCGStab(mxArray *A_m, mxArray *b_m, int Size, double slowc, int block, bool is_two_boundaries, int it_, mxArray *x0_m, int preconditioner)
{
  /* precond = 0  => Jacobi
     precond = 1  => Incomplet LU decomposition*/
  size_t n = mxGetM(A_m);
  mxArray *L1 = nullptr, *U1 = nullptr, *Diag = nullptr;

  if (preconditioner == 0)
    {
      mxArray *lhs0[1];
      mxArray *rhs0[] = { A_m, mxCreateDoubleScalar(0) };
      mexCallMATLAB(std::extent_v<decltype(lhs0)>, lhs0, std::extent_v<decltype(rhs0)>, rhs0, "spdiags");
      mxArray *tmp = lhs0[0];
      double *tmp_val = mxGetPr(tmp);
      Diag = mxCreateSparse(n, n, n, mxREAL);
      mwIndex *Diag_i = mxGetIr(Diag);
      mwIndex *Diag_j = mxGetJc(Diag);
      double *Diag_val = mxGetPr(Diag);
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
      const char *field_names[] = {"type", "droptol", "milu", "udiag", "thresh"};
      const int type = 0, droptol = 1, milu = 2, udiag = 3, thresh = 4;
      mwSize dims[1] = { static_cast<mwSize>(1) };
      mxArray *Setup = mxCreateStructArray(1, dims, std::extent_v<decltype(field_names)>, field_names);
      mxSetFieldByNumber(Setup, 0, type, mxCreateString("ilutp"));
      mxSetFieldByNumber(Setup, 0, droptol, mxCreateDoubleScalar(lu_inc_tol));
      mxSetFieldByNumber(Setup, 0, milu, mxCreateString("off"));
      mxSetFieldByNumber(Setup, 0, udiag, mxCreateDoubleScalar(0));
      mxSetFieldByNumber(Setup, 0, thresh, mxCreateDoubleScalar(1));
      mxArray *lhs0[2];
      mxArray *rhs0[] = { A_m, Setup };
      if (mexCallMATLAB(std::extent_v<decltype(lhs0)>, lhs0, std::extent_v<decltype(rhs0)>, rhs0, "ilu"))
        throw FatalException{"In BiCGStab, the incomplete LU decomposition (ilu) has failed"};
      L1 = lhs0[0];
      U1 = lhs0[1];
      mxDestroyArray(Setup);
    }
  double flags = 2;
  mxArray *z = nullptr;
  if (steady_state) /*Octave BicStab algorihtm involves a 0 division in case of a preconditionner equal to the LU decomposition of A matrix*/
    {
      mxArray *res = mult_SAT_B(Sparse_transpose(A_m), x0_m);
      double *resid = mxGetPr(res);
      double *b = mxGetPr(b_m);
      for (int i = 0; i < static_cast<int>(n); i++)
        resid[i] = b[i] - resid[i];
      mxArray *lhs[1];
      mxArray *rhs[] = { L1, res };
      mexCallMATLAB(std::extent_v<decltype(lhs)>, lhs, std::extent_v<decltype(rhs)>, rhs, "mldivide");
      mxArray *rhs2[] = { U1, lhs[0] };
      mexCallMATLAB(std::extent_v<decltype(lhs)>, lhs, std::extent_v<decltype(rhs2)>, rhs2, "mldivide");
      z = lhs[0];
      double *phat = mxGetPr(z);
      double *x0 = mxGetPr(x0_m);
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
          mxArray *rhs[] = { A_m, b_m, mxCreateDoubleScalar(1e-6),
            mxCreateDoubleScalar(static_cast<double>(n)), Diag };
          mxArray *lhs[2];
          mexCallMATLAB(std::extent_v<decltype(lhs)>, lhs, std::extent_v<decltype(rhs)>, rhs, "bicgstab");
          z = lhs[0];
          mxArray *flag = lhs[1];
          double *flag1 = mxGetPr(flag);
          flags = flag1[0];
          mxDestroyArray(flag);
          mxDestroyArray(rhs[2]);
          mxDestroyArray(rhs[3]);
          mxDestroyArray(rhs[4]);
        }
      else if (preconditioner == 1)
        {
          /*[za,flag1] = bicgstab(g1a,b,1e-6,Blck_size*periods,L1,U1);*/
          mxArray *rhs[] = { A_m, b_m, mxCreateDoubleScalar(1e-6),
            mxCreateDoubleScalar(static_cast<double>(n)), L1, U1, x0_m };
          mxArray *lhs[2];
          mexCallMATLAB(std::extent_v<decltype(lhs)>, lhs, std::extent_v<decltype(rhs)>, rhs, "bicgstab");
          z = lhs[0];
          mxArray *flag = lhs[1];
          double *flag1 = mxGetPr(flag);
          flags = flag1[0];
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
                       + to_string(block+1)).c_str());
      else if (flags == 2)
        mexWarnMsgTxt(("Error in bytecode: Preconditioner is ill-conditioned, in block "
                       + to_string(block+1)).c_str());
      else if (flags == 3)
        mexWarnMsgTxt(("Error in bytecode: BiCGStab stagnated (Two consecutive iterates were the same.), in block "
                       + to_string(block+1)).c_str());
      lu_inc_tol /= 10;
    }
  else
    {
      double *res = mxGetPr(z);
      if (is_two_boundaries)
        for (int i = 0; i < static_cast<int>(n); i++)
          {
            int eq = index_vara[i+Size*y_kmin];
            double yy = -(res[i] + y[eq]);
            direction[eq] = yy;
            y[eq] += slowc * yy;
          }
      else
        for (int i = 0; i < static_cast<int>(n); i++)
          {
            int eq = index_vara[i];
            double yy = -(res[i] + y[eq+it_*y_size]);
            direction[eq] = yy;
            y[eq+it_*y_size] += slowc * yy;
          }
    }
  mxDestroyArray(A_m);
  mxDestroyArray(b_m);
  mxDestroyArray(x0_m);
  mxDestroyArray(z);
}

void
Interpreter::Singular_display(int block, int Size)
{
  bool zero_solution;
  Simple_Init(Size, IM_i, zero_solution);
  mxArray *rhs[] = { mxCreateDoubleMatrix(Size, Size, mxREAL) };
  double *pind = mxGetPr(rhs[0]);
  for (int j = 0; j < Size * Size; j++)
    pind[j] = 0.0;
  for (int ii = 0; ii < Size; ii++)
    {
      NonZeroElem *first;
      int nb_eq = At_Col(ii, &first);
      for (int j = 0; j < nb_eq; j++)
        {
          int k = first->u_index;
          int jj = first->r_index;
          pind[ii * Size + jj] = u[k];
          first = first->NZE_C_N;
        }
    }
  mxArray *lhs[3];
  mexCallMATLAB(std::extent_v<decltype(lhs)>, lhs, std::extent_v<decltype(rhs)>, rhs, "svd");
  mxArray *SVD_u = lhs[0];
  mxArray *SVD_s = lhs[1];
  double *SVD_ps = mxGetPr(SVD_s);
  double *SVD_pu = mxGetPr(SVD_u);
  for (int i = 0; i < Size; i++)
    if (abs(SVD_ps[i * (1 + Size)]) < 1e-12)
      {
        mexPrintf(" The following equations form a linear combination:\n    ");
        double max_u = 0;
        for (int j = 0; j < Size; j++)
          if (abs(SVD_pu[j + i * Size]) > abs(max_u))
            max_u = SVD_pu[j + i * Size];
        vector<int> equ_list;
        for (int j = 0; j < Size; j++)
          {
            double rr = SVD_pu[j + i * Size] / max_u;
            if (rr < -1e-10)
              {
                equ_list.push_back(j);
                if (rr != -1)
                  mexPrintf(" - %3.2f*Dequ_%d_dy", abs(rr), j+1);
                else
                  mexPrintf(" - Dequ_%d_dy", j+1);
              }
            else if (rr > 1e-10)
              {
                equ_list.push_back(j);
                if (j > 0)
                  if (rr != 1)
                    mexPrintf(" + %3.2f*Dequ_%d_dy", rr, j+1);
                  else
                    mexPrintf(" + Dequ_%d_dy", j+1);
                else if (rr != 1)
                  mexPrintf(" %3.2f*Dequ_%d_dy", rr, j+1);
                else
                  mexPrintf(" Dequ_%d_dy", j+1);
              }
          }
        mexPrintf(" = 0\n");
      }
  mxDestroyArray(lhs[0]);
  mxDestroyArray(lhs[1]);
  mxDestroyArray(lhs[2]);
  if (block > 1)
    throw FatalException{"In Solve_ByteCode_Sparse_GaussianElimination, singular system in block "
                         + to_string(block+1)};
  else
    throw FatalException{"In Solve_ByteCode_Sparse_GaussianElimination, singular system"};
}

bool
Interpreter::Solve_ByteCode_Sparse_GaussianElimination(int Size, int blck, int it_)
{
  int pivj = 0, pivk = 0;
  NonZeroElem **bc = static_cast<NonZeroElem **>(mxMalloc(Size*sizeof(*bc)));
  test_mxMalloc(bc, __LINE__, __FILE__, __func__, Size*sizeof(*bc));
  double *piv_v = static_cast<double *>(mxMalloc(Size*sizeof(double)));
  test_mxMalloc(piv_v, __LINE__, __FILE__, __func__, Size*sizeof(double));
  int *pivj_v = static_cast<int *>(mxMalloc(Size*sizeof(int)));
  test_mxMalloc(pivj_v, __LINE__, __FILE__, __func__, Size*sizeof(int));
  int *pivk_v = static_cast<int *>(mxMalloc(Size*sizeof(int)));
  test_mxMalloc(pivk_v, __LINE__, __FILE__, __func__, Size*sizeof(int));
  int *NR = static_cast<int *>(mxMalloc(Size*sizeof(int)));
  test_mxMalloc(NR, __LINE__, __FILE__, __func__, Size*sizeof(int));

  for (int i = 0; i < Size; i++)
    {
      /*finding the max-pivot*/
      double piv = 0, piv_abs = 0;
      NonZeroElem *first;
      int nb_eq = At_Col(i, &first);
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
                  if (blck > 1)
                    mexPrintf("Error: singular system in Simulate_NG in block %d\n", blck+1);
                  else
                    mexPrintf("Error: singular system in Simulate_NG");
                }
              return true;
            }
          else
            {
              if (blck > 1)
                throw FatalException{"In Solve_ByteCode_Sparse_GaussianElimination, singular system in block "
                                     + to_string(blck+1)};
              else
                throw FatalException{"In Solve_ByteCode_Sparse_GaussianElimination, singular system"};
            }
        }
      double markovitz = 0, markovitz_max = -9e70;
      if (!one)
        for (int j = 0; j < l; j++)
          {
            if (N_max > 0 && NR[j] > 0)
              {
                if (fabs(piv_v[j]) > 0)
                  {
                    if (markowitz_c > 0)
                      markovitz = exp(log(fabs(piv_v[j])/piv_abs)
                                      -markowitz_c*log(static_cast<double>(NR[j])
                                                       /static_cast<double>(N_max)));
                    else
                      markovitz = fabs(piv_v[j])/piv_abs;
                    }
                  else
                    markovitz = 0;
                }
              else
                markovitz = fabs(piv_v[j])/piv_abs;
              if (markovitz > markovitz_max)
                {
                  piv = piv_v[j];
                  pivj = pivj_v[j]; //Line number
                  pivk = pivk_v[j]; //positi
                  markovitz_max = markovitz;
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
                      markovitz = exp(log(fabs(piv_v[j])/piv_abs)
                                      -markowitz_c*log(static_cast<double>(NR[j])
                                                       /static_cast<double>(N_max)));
                    else
                      markovitz = fabs(piv_v[j])/piv_abs;
                  }
                else
                  markovitz = 0;
              }
            else
              markovitz = fabs(piv_v[j])/piv_abs;
            if (NR[j] == 1)
              {
                piv = piv_v[j];
                pivj = pivj_v[j]; //Line number
                pivk = pivk_v[j]; //positi
                markovitz_max = markovitz;
              }
          }
      pivot[i] = pivj;
      pivotk[i] = pivk;
      pivotv[i] = piv;
      line_done[pivj] = true;

      /*divide all the non zeros elements of the line pivj by the max_pivot*/
      int nb_var = At_Row(pivj, &first);
      for (int j = 0; j < nb_var; j++)
        {
          u[first->u_index] /= piv;
          first = first->NZE_R_N;
        }
      u[b[pivj]] /= piv;
      /*subtract the elements on the non treated lines*/
      nb_eq = At_Col(i, &first);
      NonZeroElem *first_piva;
      int nb_var_piva = At_Row(pivj, &first_piva);
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
          NonZeroElem *first_piv = first_piva;
          NonZeroElem *first_sub;
          int nb_var_sub = At_Row(row, &first_sub);
          int l_sub = 0, l_piv = 0;
          int sub_c_index = first_sub->c_index, piv_c_index = first_piv->c_index;
          while (l_sub < nb_var_sub || l_piv < nb_var_piv)
            if (l_sub < nb_var_sub && (sub_c_index < piv_c_index || l_piv >= nb_var_piv))
              {
                first_sub = first_sub->NZE_R_N;
                if (first_sub)
                  sub_c_index = first_sub->c_index;
                else
                  sub_c_index = Size;
                l_sub++;
              }
            else if (sub_c_index > piv_c_index || l_sub >= nb_var_sub)
              {
                int tmp_u_count = Get_u();
                Insert(row, first_piv->c_index, tmp_u_count, 0);
                u[tmp_u_count] = -u[first_piv->u_index]*first_elem;
                first_piv = first_piv->NZE_R_N;
                if (first_piv)
                  piv_c_index = first_piv->c_index;
                else
                  piv_c_index = Size;
                l_piv++;
              }
            else
              {
                if (i == sub_c_index)
                  {
                    NonZeroElem *firsta = first;
                    NonZeroElem *first_suba = first_sub->NZE_R_N;
                    Delete(first_sub->r_index, first_sub->c_index);
                    first = firsta->NZE_C_N;
                    first_sub = first_suba;
                    if (first_sub)
                      sub_c_index = first_sub->c_index;
                    else
                      sub_c_index = Size;
                    l_sub++;
                    first_piv = first_piv->NZE_R_N;
                    if (first_piv)
                      piv_c_index = first_piv->c_index;
                    else
                      piv_c_index = Size;
                    l_piv++;
                  }
                else
                  {
                    u[first_sub->u_index] -= u[first_piv->u_index]*first_elem;
                    first_sub = first_sub->NZE_R_N;
                    if (first_sub)
                      sub_c_index = first_sub->c_index;
                    else
                      sub_c_index = Size;
                    l_sub++;
                    first_piv = first_piv->NZE_R_N;
                    if (first_piv)
                      piv_c_index = first_piv->c_index;
                    else
                      piv_c_index = Size;
                    l_piv++;
                  }
              }
          u[b[row]] -= u[b[pivj]]*first_elem;
        }
    }
  double slowc_lbx = slowc;
  for (int i = 0; i < y_size; i++)
    ya[i+it_*y_size] = y[i+it_*y_size];

  slowc_save = slowc;
  simple_bksub(it_, Size, slowc_lbx);
  End_GE();
  mxFree(piv_v);
  mxFree(pivj_v);
  mxFree(pivk_v);
  mxFree(NR);
  mxFree(bc);
  return false;
}

void
Interpreter::Solve_ByteCode_Symbolic_Sparse_GaussianElimination(int Size, bool symbolic, int Block_number)
{
  /*Triangularisation at each period of a block using a simple gaussian Elimination*/
  int *save_op = nullptr, *save_opa = nullptr, *save_opaa = nullptr;
  long int nop = 0, nopa = 0;
  bool record = false;

  int pivj = 0, pivk = 0;
  int tbreak = 0, last_period = periods;

  double *piv_v = static_cast<double *>(mxMalloc(Size*sizeof(double)));
  test_mxMalloc(piv_v, __LINE__, __FILE__, __func__, Size*sizeof(double));
  int *pivj_v = static_cast<int *>(mxMalloc(Size*sizeof(int)));
  test_mxMalloc(pivj_v, __LINE__, __FILE__, __func__, Size*sizeof(int));
  int *pivk_v = static_cast<int *>(mxMalloc(Size*sizeof(int)));
  test_mxMalloc(pivk_v, __LINE__, __FILE__, __func__, Size*sizeof(int));
  int *NR = static_cast<int *>(mxMalloc(Size*sizeof(int)));
  test_mxMalloc(NR, __LINE__, __FILE__, __func__, Size*sizeof(int));
  NonZeroElem **bc = static_cast<NonZeroElem **>(mxMalloc(Size*sizeof(NonZeroElem *)));
  test_mxMalloc(bc, __LINE__, __FILE__, __func__, Size*sizeof(NonZeroElem *));

  for (int t = 0; t < periods; t++)
    {
#ifdef MATLAB_MEX_FILE
      if (utIsInterruptPending())
        throw UserException{};
#endif

      if (record && symbolic)
        {
          save_op = static_cast<int *>(mxMalloc(nop*sizeof(int)));
          test_mxMalloc(save_op, __LINE__, __FILE__, __func__, nop*sizeof(int));
          nopa = nop;
        }
      nop = 0;
      Clear_u();
      int ti = t*Size;
      for (int i = ti; i < Size+ti; i++)
        {
          /*finding the max-pivot*/
          double piv = 0, piv_abs = 0;
          NonZeroElem *first;
          int nb_eq = At_Col(i, 0, &first);
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
                              markovitz = exp(log(fabs(piv_v[j])/piv_abs)
                                              -markowitz_c*log(static_cast<double>(NR[j])
                                                               /static_cast<double>(N_max)));
                            else
                              markovitz = fabs(piv_v[j])/piv_abs;
                          }
                        else
                          markovitz = 0;
                      }
                    else
                      markovitz = fabs(piv_v[j])/piv_abs;
                    if (markovitz > markovitz_max)
                      {
                        piv = piv_v[j];
                        pivj = pivj_v[j]; //Line number
                        pivk = pivk_v[j]; //positi
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
                              markovitz = exp(log(fabs(piv_v[j])/piv_abs)
                                              -markowitz_c*log(static_cast<double>(NR[j])
                                                               /static_cast<double>(N_max)));
                            else
                              markovitz = fabs(piv_v[j])/piv_abs;
                          }
                        else
                          markovitz = 0;
                      }
                    else
                      markovitz = fabs(piv_v[j])/piv_abs;
                    if (NR[j] == 1)
                      {
                        piv = piv_v[j];
                        pivj = pivj_v[j]; //Line number
                        pivk = pivk_v[j]; //positi
                        markovitz_max = markovitz;
                        NR_max = NR[j];
                      }
                  }
              if (fabs(piv) < eps && verbosity >= 1)
                mexPrintf("==> Error NR_max=%d, N_max=%d and piv=%f, piv_abs=%f, markovitz_max=%f\n", NR_max, N_max, piv, piv_abs, markovitz_max);
              if (NR_max == 0 && verbosity >= 1)
                mexPrintf("==> Error NR_max=0 and piv=%f, markovitz_max=%f\n", piv, markovitz_max);
              pivot[i] = pivj;
              pivot_save[i] = pivj;
              pivotk[i] = pivk;
              pivotv[i] = piv;
            }
          else
            {
              pivj = pivot[i-Size]+Size;
              pivot[i] = pivj;
              At_Pos(pivj, i, &first);
              pivk = first->u_index;
              piv = u[pivk];
              piv_abs = fabs(piv);
            }
          line_done[pivj] = true;

          if (record && symbolic)
            {
              if (nop+1 >= nopa)
                {
                  nopa = static_cast<long>(mem_increasing_factor*static_cast<double>(nopa));
                  save_op = static_cast<int *>(mxRealloc(save_op, nopa*sizeof(int)));
                }
              t_save_op_s *save_op_s = reinterpret_cast<t_save_op_s *>(&save_op[nop]);
              save_op_s->operat = IFLD;
              save_op_s->first = pivk;
              save_op_s->lag = 0;
              nop += 2;
              if (piv_abs < eps)
                {
                  if (Block_number > 1)
                    throw FatalException{"In Solve_ByteCode_Symbolic_Sparse_GaussianElimination, singular system in block "
                                         + to_string(Block_number+1)};
                  else
                    throw FatalException{"In Solve_ByteCode_Symbolic_Sparse_GaussianElimination, singular system"};
                }
              /*divide all the non zeros elements of the line pivj by the max_pivot*/
              int nb_var = At_Row(pivj, &first);
              for (int j = 0; j < nb_var; j++)
                {
                  u[first->u_index] /= piv;
                  if (nop+j*2+1 >= nopa)
                    {
                      nopa = static_cast<long>(mem_increasing_factor*static_cast<double>(nopa));
                      save_op = static_cast<int *>(mxRealloc(save_op, nopa*sizeof(int)));
                    }
                  save_op_s = reinterpret_cast<t_save_op_s *>(&save_op[nop+j*2]);
                  save_op_s->operat = IFDIV;
                  save_op_s->first = first->u_index;
                  save_op_s->lag = first->lag_index;
                  first = first->NZE_R_N;
                }
              nop += nb_var*2;
              u[b[pivj]] /= piv;
              if (nop+1 >= nopa)
                {
                  nopa = static_cast<long>(mem_increasing_factor*static_cast<double>(nopa));
                  save_op = static_cast<int *>(mxRealloc(save_op, nopa*sizeof(int)));
                }
              save_op_s = reinterpret_cast<t_save_op_s *>(&save_op[nop]);
              save_op_s->operat = IFDIV;
              save_op_s->first = b[pivj];
              save_op_s->lag = 0;
              nop += 2;
              // Subtract the elements on the non treated lines
              nb_eq = At_Col(i, &first);
              NonZeroElem *first_piva;
              int nb_var_piva = At_Row(pivj, &first_piva);

              int nb_eq_todo = 0;
              for (int j = 0; j < nb_eq && first; j++)
                {
                  if (!line_done[first->r_index])
                    bc[nb_eq_todo++] = first;
                  first = first->NZE_C_N;
                }
              for (int j = 0; j < nb_eq_todo; j++)
                {
                  t_save_op_s *save_op_s_l;
                  NonZeroElem *first = bc[j];
                  int row = first->r_index;
                  double first_elem = u[first->u_index];
                  if (nop+1 >= nopa)
                    {
                      nopa = static_cast<long>(mem_increasing_factor*static_cast<double>(nopa));
                      save_op = static_cast<int *>(mxRealloc(save_op, nopa*sizeof(int)));
                    }
                  save_op_s_l = reinterpret_cast<t_save_op_s *>(&save_op[nop]);
                  save_op_s_l->operat = IFLD;
                  save_op_s_l->first = first->u_index;
                  save_op_s_l->lag = abs(first->lag_index);
                  nop += 2;

                  int nb_var_piv = nb_var_piva;
                  NonZeroElem *first_piv = first_piva;
                  NonZeroElem *first_sub;
                  int nb_var_sub = At_Row(row, &first_sub);
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
                            sub_c_index = Size*periods;
                          l_sub++;
                        }
                      else if (sub_c_index > piv_c_index || l_sub >= nb_var_sub)
                        {
                          // There is an nonzero element at row pivot but not at the current row=> insert a negative element in the current row
                          int tmp_u_count = Get_u();
                          int lag = first_piv->c_index/Size-row/Size;
                          Insert(row, first_piv->c_index, tmp_u_count, lag);
                          u[tmp_u_count] = -u[first_piv->u_index]*first_elem;
                          if (nop+2 >= nopa)
                            {
                              nopa = static_cast<long>(mem_increasing_factor*static_cast<double>(nopa));
                              save_op = static_cast<int *>(mxRealloc(save_op, nopa*sizeof(int)));
                            }
                          save_op_s_l = reinterpret_cast<t_save_op_s *>(&save_op[nop]);
                          save_op_s_l->operat = IFLESS;
                          save_op_s_l->first = tmp_u_count;
                          save_op_s_l->second = first_piv->u_index;
                          save_op_s_l->lag = max(first_piv->lag_index, abs(tmp_lag));
                          nop += 3;
                          first_piv = first_piv->NZE_R_N;
                          if (first_piv)
                            piv_c_index = first_piv->c_index;
                          else
                            piv_c_index = Size*periods;
                          l_piv++;
                        }
                      else /*first_sub->c_index==first_piv->c_index*/
                        {
                          if (i == sub_c_index)
                            {
                              NonZeroElem *firsta = first;
                              NonZeroElem *first_suba = first_sub->NZE_R_N;
                              Delete(first_sub->r_index, first_sub->c_index);
                              first = firsta->NZE_C_N;
                              first_sub = first_suba;
                              if (first_sub)
                                sub_c_index = first_sub->c_index;
                              else
                                sub_c_index = Size*periods;
                              l_sub++;
                              first_piv = first_piv->NZE_R_N;
                              if (first_piv)
                                piv_c_index = first_piv->c_index;
                              else
                                piv_c_index = Size*periods;
                              l_piv++;
                            }
                          else
                            {
                              u[first_sub->u_index] -= u[first_piv->u_index]*first_elem;
                              if (nop+3 >= nopa)
                                {
                                  nopa = static_cast<long>(mem_increasing_factor*static_cast<double>(nopa));
                                  save_op = static_cast<int *>(mxRealloc(save_op, nopa*sizeof(int)));
                                }
                              save_op_s_l = reinterpret_cast<t_save_op_s *>(&save_op[nop]);
                              save_op_s_l->operat = IFSUB;
                              save_op_s_l->first = first_sub->u_index;
                              save_op_s_l->second = first_piv->u_index;
                              save_op_s_l->lag = max(abs(tmp_lag), first_piv->lag_index);
                              nop += 3;
                              first_sub = first_sub->NZE_R_N;
                              if (first_sub)
                                sub_c_index = first_sub->c_index;
                              else
                                sub_c_index = Size*periods;
                              l_sub++;
                              first_piv = first_piv->NZE_R_N;
                              if (first_piv)
                                piv_c_index = first_piv->c_index;
                              else
                                piv_c_index = Size*periods;
                              l_piv++;
                            }
                        }
                    }
                  u[b[row]] -= u[b[pivj]]*first_elem;

                  if (nop+3 >= nopa)
                    {
                      nopa = static_cast<long>(mem_increasing_factor*static_cast<double>(nopa));
                      save_op = static_cast<int *>(mxRealloc(save_op, nopa*sizeof(int)));
                    }
                  save_op_s_l = reinterpret_cast<t_save_op_s *>(&save_op[nop]);
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
                  if (Block_number > 1)
                    throw FatalException{"In Solve_ByteCode_Symbolic_Sparse_GaussianElimination, singular system in block "
                                         + to_string(Block_number+1)};
                  else
                    throw FatalException{"In Solve_ByteCode_Symbolic_Sparse_GaussianElimination, singular system"};
                }
              // Divide all the non zeros elements of the line pivj by the max_pivot
              int nb_var = At_Row(pivj, &first);
              for (int j = 0; j < nb_var; j++)
                {
                  u[first->u_index] /= piv;
                  first = first->NZE_R_N;
                }
              nop += nb_var*2;
              u[b[pivj]] /= piv;
              nop += 2;
              // Subtract the elements on the non treated lines
              nb_eq = At_Col(i, &first);
              NonZeroElem *first_piva;
              int nb_var_piva = At_Row(pivj, &first_piva);

              int nb_eq_todo = 0;
              for (int j = 0; j < nb_eq && first; j++)
                {
                  if (!line_done[first->r_index])
                    bc[nb_eq_todo++] = first;
                  first = first->NZE_C_N;
                }
              for (int j = 0; j < nb_eq_todo; j++)
                {
                  NonZeroElem *first = bc[j];
                  int row = first->r_index;
                  double first_elem = u[first->u_index];
                  nop += 2;
                  int nb_var_piv = nb_var_piva;
                  NonZeroElem *first_piv = first_piva;
                  NonZeroElem *first_sub;
                  int nb_var_sub = At_Row(row, &first_sub);
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
                            sub_c_index = Size*periods;
                          l_sub++;
                        }
                      else if (sub_c_index > piv_c_index || l_sub >= nb_var_sub)
                        {
                          /* There is an nonzero element at row pivot but not
                             at the current row ⇒ insert a negative element in the
                             current row */
                          int tmp_u_count = Get_u();
                          int lag = first_piv->c_index/Size-row/Size;
                          Insert(row, first_piv->c_index, tmp_u_count, lag);
                          u[tmp_u_count] = -u[first_piv->u_index]*first_elem;
                          nop += 3;
                          first_piv = first_piv->NZE_R_N;
                          if (first_piv)
                            piv_c_index = first_piv->c_index;
                          else
                            piv_c_index = Size*periods;
                          l_piv++;
                        }
                      else /*first_sub->c_index==first_piv->c_index*/
                        {
                          if (i == sub_c_index)
                            {
                              NonZeroElem *firsta = first;
                              NonZeroElem *first_suba = first_sub->NZE_R_N;
                              Delete(first_sub->r_index, first_sub->c_index);
                              first = firsta->NZE_C_N;
                              first_sub = first_suba;
                              if (first_sub)
                                sub_c_index = first_sub->c_index;
                              else
                                sub_c_index = Size*periods;
                              l_sub++;
                              first_piv = first_piv->NZE_R_N;
                              if (first_piv)
                                piv_c_index = first_piv->c_index;
                              else
                                piv_c_index = Size*periods;
                              l_piv++;
                            }
                          else
                            {
                              u[first_sub->u_index] -= u[first_piv->u_index]*first_elem;
                              nop += 3;
                              first_sub = first_sub->NZE_R_N;
                              if (first_sub)
                                sub_c_index = first_sub->c_index;
                              else
                                sub_c_index = Size*periods;
                              l_sub++;
                              first_piv = first_piv->NZE_R_N;
                              if (first_piv)
                                piv_c_index = first_piv->c_index;
                              else
                                piv_c_index = Size*periods;
                              l_piv++;
                            }
                        }
                    }
                  u[b[row]] -= u[b[pivj]]*first_elem;
                  nop += 3;
                }
            }
        }
      if (symbolic)
        {
          if (t > static_cast<int>(periods*0.35))
            {
              symbolic = false;
              mxFree(save_opaa);
              mxFree(save_opa);
              mxFree(save_op);
            }
          else if (record && nop == nop1)
            {
              if (t > static_cast<int>(periods*0.35))
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
                  if (compare(save_op, save_opa, save_opaa, t, periods, nop, Size))
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
  double slowc_lbx = slowc;
  for (int i = 0; i < y_size*(periods+y_kmin); i++)
    ya[i] = y[i];
  slowc_save = slowc;
  bksub(tbreak, last_period, Size, slowc_lbx);
  End_GE();
}

void
Interpreter::Check_and_Correct_Previous_Iteration(int y_size, int size)
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
              y[eq+it_*y_size] = ya[eq+it_*y_size] + slowc_save * direction[eq+it_*y_size];
            }
          compute_complete(true, res1, res2, max_res, max_res_idx);
        }

      while (res2 > g0 && slowc_save > 1e-1)
        {
          prev_slowc_save = slowc_save;
          slowc_save /= 1.5;
          for (int i = 0; i < size; i++)
            {
              int eq = index_vara[i];
              y[eq+it_*y_size] = ya[eq+it_*y_size] + slowc_save * direction[eq+it_*y_size];
            }
          compute_complete(true, res1, res2, max_res, max_res_idx);
        }
      if (verbosity >= 2)
        mexPrintf("Error: Simulation diverging, trying to correct it using slowc=%f\n", slowc_save);
      for (int i = 0; i < size; i++)
        {
          int eq = index_vara[i];
          y[eq+it_*y_size] = ya[eq+it_*y_size] + slowc_save * direction[eq+it_*y_size];
        }
      compute_complete(false, res1, res2, max_res, max_res_idx);
    }
  else
    for (int i = 0; i < size; i++)
      {
        int eq = index_vara[i];
        y[eq+it_*y_size] = ya[eq+it_*y_size] + slowc_save * direction[eq+it_*y_size];
      }
  slowc_save = slowc;
}

bool
Interpreter::Simulate_One_Boundary(int block_num, int y_size, int size)
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
            mexPrintf("-> variable %s (%d) at time %d = %f direction = %f\n", get_variable(SymbolType::endogenous, j).c_str(), j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
          else
            mexPrintf("   variable %s (%d) at time %d = %f direction = %f\n", get_variable(SymbolType::endogenous, j).c_str(), j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
        }
#endif
      if (steady_state)
        {
          if (verbosity >= 1)
            {
              if (iter == 0)
                mexPrintf(" the initial values of endogenous variables are too far from the solution.\nChange them!\n");
              else
                mexPrintf(" dynare cannot improve the simulation in block %d at time %d (variable %d)\n", block_num+1, it_+1, index_vara[max_res_idx]+1);
              mexEvalString("drawnow;");
            }
        }
      else
        {
          if (iter == 0)
            throw FatalException{"In Simulate_One_Boundary, The initial values of endogenous variables are too far from the solution. Change them!"};
          else
            throw FatalException{"In Simulate_One_Boundary, Dynare cannot improve the simulation in block "
                                 + to_string(block_num+1) + " at time " + to_string(it_+1)
                                 + " (variable " + to_string(index_vara[max_res_idx]+1)};
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
              mexPrintf(preconditioner_print_out("MODEL STEADY STATE: (method=GMRES)\n", preconditioner, true).c_str());
              break;
            case 8:
              mexPrintf(preconditioner_print_out("MODEL STEADY STATE: (method=BiCGStab)\n", preconditioner, true).c_str());
              break;
            }
        }

      mexPrintf("------------------------------------\n");
      mexPrintf("      Iteration no. %d\n", iter+1);
      mexPrintf("      Inf-norm error = %.3e\n", static_cast<double>(max_res));
      mexPrintf("      2-norm error   = %.3e\n", static_cast<double>(sqrt(res2)));
      mexPrintf("      1-norm error   = %.3e\n", static_cast<double>(res1));
      mexPrintf("------------------------------------\n");
    }
  bool zero_solution;

  if ((solve_algo == 5 && steady_state) || (stack_solve_algo == 5 && !steady_state))
    Simple_Init(size, IM_i, zero_solution);
  else
    {
      x0_m = mxCreateDoubleMatrix(size, 1, mxREAL);
      if (!x0_m)
        throw FatalException{"In Simulate_One_Boundary, can't allocate x0_m vector"};
      if (!((solve_algo == 6 && steady_state)
            || ((stack_solve_algo == 0 || stack_solve_algo == 1 || stack_solve_algo == 4
                 || stack_solve_algo == 6) && !steady_state)))
        {
          b_m = mxCreateDoubleMatrix(size, 1, mxREAL);
          if (!b_m)
            throw FatalException{"In Simulate_One_Boundary, can't allocate b_m vector"};
          A_m = mxCreateSparse(size, size, min(static_cast<int>(IM_i.size()*2), size * size), mxREAL);
          if (!A_m)
            throw FatalException{"In Simulate_One_Boundary, can't allocate A_m matrix"};
          Init_Matlab_Sparse_Simple(size, IM_i, A_m, b_m, zero_solution, x0_m);
        }
      else
        {
          Init_UMFPACK_Sparse_Simple(size, IM_i, &Ap, &Ai, &Ax, &b, zero_solution, x0_m);
          if (Ap_save[size] != Ap[size])
            {
              mxFree(Ai_save);
              mxFree(Ax_save);
              Ai_save = static_cast<SuiteSparse_long *>(mxMalloc(Ap[size] * sizeof(SuiteSparse_long)));
              test_mxMalloc(Ai_save, __LINE__, __FILE__, __func__, Ap[size] * sizeof(SuiteSparse_long));
              Ax_save = static_cast<double *>(mxMalloc(Ap[size] * sizeof(double)));
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
        double yy = -(y[eq+it_*y_size]);
        direction[eq] = yy;
        y[eq+it_*y_size] += slowc * yy;
      }
  else
    {
      if ((solve_algo == 5 && steady_state) || (stack_solve_algo == 5 && !steady_state))
        singular_system = Solve_ByteCode_Sparse_GaussianElimination(size, block_num, it_);
      else if ((solve_algo == 7 && steady_state) || (stack_solve_algo == 2 && !steady_state))
        Solve_Matlab_GMRES(A_m, b_m, size, slowc, block_num, false, it_, x0_m);
      else if ((solve_algo == 8 && steady_state) || (stack_solve_algo == 3 && !steady_state))
        Solve_Matlab_BiCGStab(A_m, b_m, size, slowc, block_num, false, it_, x0_m, preconditioner);
      else if ((solve_algo == 6 && steady_state) || ((stack_solve_algo == 0 || stack_solve_algo == 1 || stack_solve_algo == 4 || stack_solve_algo == 6) && !steady_state))
        {
          Solve_LU_UMFPack(Ap, Ai, Ax, b, size, size, slowc, false, it_);
          mxDestroyArray(x0_m);
        }
    }
  return singular_system;
}

bool
Interpreter::solve_linear(int block_num, int y_size, int size, int iter)
{
  bool cvg = false;
  compute_complete(false, res1, res2, max_res, max_res_idx);
  cvg = (max_res < solve_tolf);
  if (!cvg || isnan(res1) || isinf(res1))
    {
      if (iter)
        Check_and_Correct_Previous_Iteration(y_size, size);
      bool singular_system = Simulate_One_Boundary(block_num, y_size, size);
      if (singular_system && verbosity >= 1)
        Singular_display(block_num, size);
    }
  return cvg;
}

void
Interpreter::solve_non_linear(int block_num, int y_size, int size)
{
  max_res_idx = 0;
  bool cvg = false;
  iter = 0;
  glambda2 = g0 = very_big;
  while (!cvg && iter < maxit_)
    {
      cvg = solve_linear(block_num, y_size, size, iter);
      g0 = res2;
      iter++;
    }
  if (!cvg)
    {
      if (steady_state)
        throw FatalException{"In Solve Forward/Backward Complete, convergence not achieved in block "
                             + to_string(block_num+1) + ", after " + to_string(iter)
                             + " iterations"};
      else
        throw FatalException{"In Solve Forward/Backward Complete, convergence not achieved in block "
                             + to_string(block_num+1) + ", at time " + to_string(it_)
                             + ", after " + to_string(iter) + " iterations"};
    }
}

void
Interpreter::Simulate_Newton_One_Boundary(bool forward)
{
  g1 = static_cast<double *>(mxMalloc(size*size*sizeof(double)));
  test_mxMalloc(g1, __LINE__, __FILE__, __func__, size*size*sizeof(double));
  r = static_cast<double *>(mxMalloc(size*sizeof(double)));
  test_mxMalloc(r, __LINE__, __FILE__, __func__, size*sizeof(double));
  iter = 0;
  if ((solve_algo == 6 && steady_state)
      || ((stack_solve_algo == 0 || stack_solve_algo == 1 || stack_solve_algo == 4 || stack_solve_algo == 6) && !steady_state))
    {
      Ap_save = static_cast<SuiteSparse_long *>(mxMalloc((size + 1) * sizeof(SuiteSparse_long)));
      test_mxMalloc(Ap_save, __LINE__, __FILE__, __func__, (size + 1) * sizeof(SuiteSparse_long));
      Ap_save[size] = 0;
      Ai_save = static_cast<SuiteSparse_long *>(mxMalloc(1 * sizeof(SuiteSparse_long)));
      test_mxMalloc(Ai_save, __LINE__, __FILE__, __func__, 1 * sizeof(SuiteSparse_long));
      Ax_save = static_cast<double *>(mxMalloc(1 * sizeof(double)));
      test_mxMalloc(Ax_save, __LINE__, __FILE__, __func__, 1 * sizeof(double));
      b_save = static_cast<double *>(mxMalloc((size) * sizeof(SuiteSparse_long)));
      test_mxMalloc(b_save, __LINE__, __FILE__, __func__, (size) * sizeof(SuiteSparse_long));
    }
  if (steady_state)
    {
      it_ = 0;
      if (!is_linear)
        solve_non_linear(block_num, y_size, size);
      else
        solve_linear(block_num, y_size, size, 0);
    }
  else if (forward)
    {
      if (!is_linear)
        for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
          solve_non_linear(block_num, y_size, size);
      else
        for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
          solve_linear(block_num, y_size, size, 0);
    }
  else
    {
      if (!is_linear)
        for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
          solve_non_linear(block_num, y_size, size);
      else
        for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
          solve_linear(block_num, y_size, size, 0);
    }
  if ((solve_algo == 6 && steady_state)
      || ((stack_solve_algo == 0 || stack_solve_algo == 1 || stack_solve_algo == 4 || stack_solve_algo == 6) && !steady_state))
    {
      mxFree(Ap_save);
      mxFree(Ai_save);
      mxFree(Ax_save);
      mxFree(b_save);
    }
  mxFree(g1);
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
Interpreter::Simulate_Newton_Two_Boundaries(int blck, int y_size, int y_kmin, int y_kmax,int Size, int periods, bool cvg, int minimal_solving_periods, int stack_solve_algo, const vector_table_conditional_local_type &vector_table_conditional_local)
{
  double top = 0.5;
  double bottom = 0.1;
  int preconditioner = 2;
  if (start_compare == 0)
    start_compare = y_kmin;
  u_count_alloc_save = u_count_alloc;
  auto t1 { chrono::high_resolution_clock::now() };
  nop1 = 0;
  mxArray *b_m = nullptr, *A_m = nullptr, *x0_m = nullptr;
  double *Ax = nullptr, *b;
  SuiteSparse_long *Ap = nullptr, *Ai = nullptr;

  if (isnan(res1) || isinf(res1) || (res2 > 12*g0 && iter > 0))
    {
      if (iter == 0 || fabs(slowc_save) < 1e-8)
        {
          if (verbosity >= 2)
            mexPrintf("res1 = %f, res2 = %f g0 = %f iter = %d\n", res1, res2, g0, iter);
          for (int j = 0; j < y_size; j++)
            {
              bool select = false;
              for (int i = 0; i < Size; i++)
                if (j == index_vara[i])
                  {
                    select = true;
                    break;
                  }
              if (verbosity >= 2)
                {
                  if (select)
                    mexPrintf("-> variable %s (%d) at time %d = %f direction = %f\n", symbol_table.getName(SymbolType::endogenous, j).c_str(), j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
                  else
                    mexPrintf("   variable %s (%d) at time %d = %f direction = %f\n", symbol_table.getName(SymbolType::endogenous, j).c_str(), j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
                }
            }
          if (iter == 0)
            throw FatalException{"In Simulate_Newton_Two_Boundaries, the initial values of endogenous variables are too far from the solution. Change them!"};
          else
            throw FatalException{"In Simulate_Newton_Two_Boundaries, dynare cannot improve the simulation in block "
                                 + to_string(blck+1) + " at time " + to_string(it_+1)
                                 + " (variable " + to_string(index_vara[max_res_idx]+1)
                                 + " = " + to_string(max_res) + ")"};
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
              double a = (1/(slowc_save * slowc_save) * t1
                          - 1/(prev_slowc_save * prev_slowc_save) * t2)
                / (slowc_save - prev_slowc_save);
              double b = (-prev_slowc_save/(slowc_save * slowc_save) * t1
                          + slowc_save/(prev_slowc_save * prev_slowc_save) * t2)
                / (slowc_save - prev_slowc_save);
              prev_slowc_save = slowc_save;
              slowc_save = max(min(-b + sqrt(b*b - 3 * a * gp0) / (3 * a),
                                   top * slowc_save), bottom * slowc_save);
            }
          glambda2 = res2;
          try_at_iteration++;
          if (slowc_save <= bottom)
            {
              for (int i = 0; i < y_size*(periods+y_kmin); i++)
                y[i] = ya[i]+direction[i];
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
            mexPrintf("The model cannot be evaluated, trying to correct it using slowc=%f\n", slowc_save);
          else
            mexPrintf("Simulation diverging, trying to correct it using slowc=%f\n", slowc_save);
        }
      for (int i = 0; i < y_size*(periods+y_kmin); i++)
        y[i] = ya[i]+slowc_save*direction[i];
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
      if (res1/res1a-1 > -0.3 && symbolic && iter > 0)
        {
          if (restart > 2)
            {
              if (verbosity >= 2)
                mexPrintf("Divergence or slowdown occurred during simulation.\nIn the next iteration, pivoting method will be applied to all periods.\n");
              symbolic = false;
              alt_symbolic = true;
              markowitz_c_s = markowitz_c;
              markowitz_c = 0;
            }
          else
            {
              if (verbosity >= 2)
                mexPrintf("Divergence or slowdown occurred during simulation.\nIn the next iteration, pivoting method will be applied for a longer period.\n");
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
              mexPrintf("MODEL SIMULATION: (method=Sparse LU)\n");
              break;
            case 1:
            case 6:
              mexPrintf("MODEL SIMULATION: (method=LBJ)\n");
              break;
            case 2:
              mexPrintf(preconditioner_print_out("MODEL SIMULATION: (method=GMRES)\n", preconditioner, false).c_str());
              break;
            case 3:
              mexPrintf(preconditioner_print_out("MODEL SIMULATION: (method=BiCGStab)\n", preconditioner, false).c_str());
              break;
            case 4:
              mexPrintf("MODEL SIMULATION: (method=Sparse LU & optimal path length)\n");
              break;
            case 5:
              mexPrintf("MODEL SIMULATION: (method=Sparse Gaussian Elimination)\n");
              break;
            }
        }
      mexPrintf("------------------------------------\n");
      mexPrintf("      Iteration no. %d\n", iter+1);
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
        Init_GE(periods, y_kmin, y_kmax, Size, IM_i);
      else
        {
          x0_m = mxCreateDoubleMatrix(periods*Size, 1, mxREAL);
          if (!x0_m)
            throw FatalException{"In Simulate_Newton_Two_Boundaries, can't allocate x0_m vector"};
          if (stack_solve_algo == 0 || stack_solve_algo == 4)
            Init_UMFPACK_Sparse(periods, y_kmin, y_kmax, Size, IM_i, &Ap, &Ai, &Ax, &b, x0_m, vector_table_conditional_local, blck);
          else
            {
              b_m = mxCreateDoubleMatrix(periods*Size, 1, mxREAL);
              if (!b_m)
                throw FatalException{"In Simulate_Newton_Two_Boundaries, can't allocate b_m vector"};
              if (stack_solve_algo != 0 && stack_solve_algo != 4)
                {
                  A_m = mxCreateSparse(periods*Size, periods*Size, IM_i.size()* periods*2, mxREAL);
                  if (!A_m)
                    throw FatalException{"In Simulate_Newton_Two_Boundaries, can't allocate A_m matrix"};
                }
              Init_Matlab_Sparse(periods, y_kmin, y_kmax, Size, IM_i, A_m, b_m, x0_m);
            }
        }
      if (stack_solve_algo == 0 || stack_solve_algo == 4)
        {
          Solve_LU_UMFPack(Ap, Ai, Ax, b, Size * periods, Size, slowc, true, 0, vector_table_conditional_local);
          mxDestroyArray(x0_m);
        }
      else if (stack_solve_algo == 1 || stack_solve_algo == 6)
        {
          Solve_Matlab_Relaxation(A_m, b_m, Size, slowc);
          mxDestroyArray(x0_m);
        }
      else if (stack_solve_algo == 2)
        Solve_Matlab_GMRES(A_m, b_m, Size, slowc, blck, true, 0, x0_m);
      else if (stack_solve_algo == 3)
        Solve_Matlab_BiCGStab(A_m, b_m, Size, slowc, blck, true, 0, x0_m, 1);
      else if (stack_solve_algo == 5)
        Solve_ByteCode_Symbolic_Sparse_GaussianElimination(Size, symbolic, blck);
    }
  using FloatSeconds = chrono::duration<double, chrono::seconds::period>;
  auto t2 { chrono::high_resolution_clock::now() };
  if (verbosity >= 1)
    {
      mexPrintf("(** %.2f seconds **)\n", FloatSeconds{t2 - t1}.count());
      mexEvalString("drawnow;");
    }
  if (!steady_state && stack_solve_algo == 4)
    {
      double ax = -0.1, bx = 1.1, cx = 0.5, fa, fb, fc, xmin;

      if (!mnbrak(&ax, &bx, &cx, &fa, &fb, &fc))
        return;
      if (!golden(ax, bx, cx, 1e-1, solve_tolf, &xmin))
        return;
      slowc = xmin;
      if (verbosity >= 1)
        {
          auto t3 { chrono::high_resolution_clock::now() };
          mexPrintf("(** %.2f seconds **)\n", FloatSeconds{t3 - t2}.count());
          mexEvalString("drawnow;");
        }
    }
  if (tbreak_g == 0)
    tbreak_g = periods;
}

void
Interpreter::fixe_u(double **u, int u_count_int, int max_lag_plus_max_lead_plus_1)
{
  u_count = u_count_int * periods;
  u_count_alloc = 2*u_count;
#ifdef DEBUG
  mexPrintf("fixe_u : alloc(%d double)\n", u_count_alloc);
#endif
  *u = static_cast<double *>(mxMalloc(u_count_alloc*sizeof(double)));
  test_mxMalloc(*u, __LINE__, __FILE__, __func__, u_count_alloc*sizeof(double));
#ifdef DEBUG
  mexPrintf("*u=%d\n", *u);
#endif
  fill_n(*u, u_count_alloc, 0);
  u_count_init = max_lag_plus_max_lead_plus_1;
}
