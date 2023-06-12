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

#include "Interpreter.hh"

constexpr double BIG = 1.0e+8, SMALL = 1.0e-5;

Interpreter::Interpreter(Evaluate &evaluator_arg, double *params_arg, double *y_arg, double *ya_arg, double *x_arg, double *steady_y_arg,
                         double *direction_arg, size_t y_size_arg,
                         size_t nb_row_x_arg, int periods_arg, int y_kmin_arg, int y_kmax_arg,
                         int maxit_arg_, double solve_tolf_arg, int y_decal_arg, double markowitz_c_arg,
                         string &filename_arg, int minimal_solving_periods_arg, int stack_solve_algo_arg, int solve_algo_arg,
                         bool global_temporary_terms_arg, bool print_arg, bool print_error_arg, mxArray *GlobalTemporaryTerms_arg,
                         bool steady_state_arg, bool block_decomposed_arg, bool print_it_arg, int col_x_arg, int col_y_arg, const BasicSymbolTable &symbol_table_arg)
: dynSparseMatrix {evaluator_arg, y_size_arg, y_kmin_arg, y_kmax_arg, print_it_arg, steady_state_arg, block_decomposed_arg, periods_arg, minimal_solving_periods_arg, symbol_table_arg, print_error_arg}
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
  print_it = print_it_arg;
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
  while (!(cvg || iter > maxit_))
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
                  mexPrintf("Reducing the path length in Newton step slowc=%f\n", slowc);
                  feclearexcept(FE_ALL_EXCEPT);
                  y[Block_Contain[0].Variable + Per_y_] = ya - slowc * (r[0] / g1[0]);
                  if (fetestexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW))
                    {
                      res1 = numeric_limits<double>::quiet_NaN();
                      if (print_error)
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
          if (print_error)
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
          mexPrintf("endo in block %d, size=%d, type=%d, steady_state=%d, print_it=%d, is_linear=%d, endo_nbr=%d, u_count_int=%d\n",
                    current_block+1, size, type, steady_state, print_it, is_linear, symbol_table_endo_nbr, u_count_int);
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
  bool r;
  vector<int> blocks;
  for (int t = 0; t < nb_periods; t++)
    {
      previous_block_exogenous.clear();
      if (old_print_it)
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
