/*
 * Copyright © 2013-2022 Dynare Team
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
#include <cmath>
#include <limits>

#include "Evaluate.hh"

#ifdef MATLAB_MEX_FILE
extern "C" bool utIsInterruptPending();
#endif

Evaluate::Evaluate(int y_size_arg, int y_kmin_arg, int y_kmax_arg, bool print_it_arg, bool steady_state_arg, int periods_arg, int minimal_solving_periods_arg) :
  print_it(print_it_arg), minimal_solving_periods(minimal_solving_periods_arg)
{
  symbol_table_endo_nbr = 0;
  Block_List_Max_Lag = 0;
  Block_List_Max_Lead = 0;
  u_count_int = 0;
  block = -1;
  y_size = y_size_arg;
  y_kmin = y_kmin_arg;
  y_kmax = y_kmax_arg;
  periods = periods_arg;
  steady_state = steady_state_arg;
}

double
Evaluate::pow1(double a, double b)
{
  double r = pow(a, b);
  if (isnan(r) || isinf(r))
    {
      res1 = std::numeric_limits<double>::quiet_NaN();
      r = 0.0000000000000000000000001;
      if (print_error)
        throw PowExceptionHandling(a, b);
    }
  return r;
}

double
Evaluate::divide(double a, double b)
{
  double r = a / b;
  if (isnan(r) || isinf(r))
    {
      res1 = std::numeric_limits<double>::quiet_NaN();
      r = 1e70;
      if (print_error)
        throw DivideExceptionHandling(a, b);
    }
  return r;
}

double
Evaluate::log1(double a)
{
  double r = log(a);
  if (isnan(r) || isinf(r))
    {
      res1 = std::numeric_limits<double>::quiet_NaN();
      r = -1e70;
      if (print_error)
        throw LogExceptionHandling(a);
    }
  return r;
}

double
Evaluate::log10_1(double a)
{
  double r = log(a);
  if (isnan(r) || isinf(r))
    {
      res1 = std::numeric_limits<double>::quiet_NaN();
      r = -1e70;
      if (print_error)
        throw Log10ExceptionHandling(a);
    }
  return r;
}

void
Evaluate::compute_block_time(int Per_u_, bool evaluate, bool no_derivative)
{
  int var = 0, lag = 0, op;
  unsigned int eq, pos_col;
  ostringstream tmp_out;
  double v1, v2, v3;
  bool go_on = true;
  double ll;
  double rr;
  double *jacob = nullptr, *jacob_other_endo = nullptr, *jacob_exo = nullptr, *jacob_exo_det = nullptr;
  EQN_block = block_num;
  stack<double> Stack;
  ExternalFunctionType function_type = ExternalFunctionType::withoutDerivative;

#ifdef DEBUG
  mexPrintf("compute_block_time\n");
#endif
  if (evaluate)
    {
      jacob = mxGetPr(jacobian_block[block_num]);
      if (!steady_state)
        {
          jacob_other_endo = mxGetPr(jacobian_other_endo_block[block_num]);
          jacob_exo = mxGetPr(jacobian_exo_block[block_num]);
          jacob_exo_det = mxGetPr(jacobian_det_exo_block[block_num]);
        }
    }
#ifdef MATLAB_MEX_FILE
  if (utIsInterruptPending())
    throw UserExceptionHandling();
#endif

  while (go_on)
    {
      switch (it_code->first)
        {
        case Tags::FNUMEXPR:
#ifdef DEBUG
          mexPrintf("FNUMEXPR\n");
#endif
          it_code_expr = it_code;
          switch (static_cast<FNUMEXPR_ *>(it_code->second)->get_expression_type())
            {
            case ExpressionType::TemporaryTerm:
#ifdef DEBUG
              mexPrintf("TemporaryTerm\n");
#endif
              EQN_type = ExpressionType::TemporaryTerm;
              EQN_equation = static_cast<FNUMEXPR_ *>(it_code->second)->get_equation();
#ifdef DEBUG
              mexPrintf("EQN_equation=%d\n", EQN_equation); mexEvalString("drawnow;");
#endif
              break;
            case ExpressionType::ModelEquation:
#ifdef DEBUG
              mexPrintf("ModelEquation\n");
#endif
              EQN_type = ExpressionType::ModelEquation;
              EQN_equation = static_cast<FNUMEXPR_ *>(it_code->second)->get_equation();
              break;
            case ExpressionType::FirstEndoDerivative:
#ifdef DEBUG
              mexPrintf("FirstEndoDerivative\n");
#endif
              EQN_type = ExpressionType::FirstEndoDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(it_code->second)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable1();
              EQN_lag1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag1();
              break;
            case ExpressionType::FirstOtherEndoDerivative:
#ifdef DEBUG
              mexPrintf("FirstOtherEndoDerivative\n");
#endif
              EQN_type = ExpressionType::FirstOtherEndoDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(it_code->second)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable1();
              EQN_lag1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag1();
              break;
            case ExpressionType::FirstExoDerivative:
#ifdef DEBUG
              mexPrintf("FirstExoDerivative\n");
#endif
              EQN_type = ExpressionType::FirstExoDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(it_code->second)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable1();
              EQN_lag1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag1();
              break;
            case ExpressionType::FirstExodetDerivative:
#ifdef DEBUG
              mexPrintf("FirstExodetDerivative\n");
#endif
              EQN_type = ExpressionType::FirstExodetDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(it_code->second)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable1();
              EQN_lag1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag1();
              break;
            case ExpressionType::FirstParamDerivative:
#ifdef DEBUG
              mexPrintf("FirstParamDerivative\n");
#endif
              EQN_type = ExpressionType::FirstParamDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(it_code->second)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable1();
              break;
            case ExpressionType::SecondEndoDerivative:
#ifdef DEBUG
              mexPrintf("SecondEndoDerivative\n");
#endif
              EQN_type = ExpressionType::SecondEndoDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(it_code->second)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable1();
              EQN_lag1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag1();
              EQN_dvar2 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable2();
              EQN_lag2 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag2();
              break;
            case ExpressionType::SecondExoDerivative:
#ifdef DEBUG
              mexPrintf("SecondExoDerivative\n");
#endif
              EQN_type = ExpressionType::SecondExoDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(it_code->second)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable1();
              EQN_lag1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag1();
              EQN_dvar2 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable2();
              EQN_lag2 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag2();
              break;
            case ExpressionType::SecondExodetDerivative:
#ifdef DEBUG
              mexPrintf("SecondExodetDerivative\n");
#endif
              EQN_type = ExpressionType::SecondExodetDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(it_code->second)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable1();
              EQN_lag1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag1();
              EQN_dvar2 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable2();
              EQN_lag2 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag2();
              break;
            case ExpressionType::SecondParamDerivative:
#ifdef DEBUG
              mexPrintf("SecondParamDerivative\n");
#endif
              EQN_type = ExpressionType::SecondParamDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(it_code->second)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable1();
              EQN_dvar2 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable2();
              break;
            case ExpressionType::ThirdEndoDerivative:
#ifdef DEBUG
              mexPrintf("ThirdEndoDerivative\n");
#endif
              EQN_type = ExpressionType::ThirdEndoDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(it_code->second)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable1();
              EQN_lag1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag1();
              EQN_dvar2 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable2();
              EQN_lag2 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag2();
              EQN_dvar3 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable3();
              EQN_lag3 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag3();
              break;
            case ExpressionType::ThirdExoDerivative:
#ifdef DEBUG
              mexPrintf("ThirdExoDerivative\n");
#endif
              EQN_type = ExpressionType::ThirdExoDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(it_code->second)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable1();
              EQN_lag1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag1();
              EQN_dvar2 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable2();
              EQN_lag2 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag2();
              EQN_dvar3 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable3();
              EQN_lag3 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag3();
              break;
            case ExpressionType::ThirdExodetDerivative:
#ifdef DEBUG
              mexPrintf("ThirdExodetDerivative\n");
#endif
              EQN_type = ExpressionType::ThirdExodetDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(it_code->second)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable1();
              EQN_lag1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag1();
              EQN_dvar2 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable2();
              EQN_lag2 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag2();
              EQN_dvar3 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable3();
              EQN_lag3 = static_cast<FNUMEXPR_ *>(it_code->second)->get_lag3();
              break;
            case ExpressionType::ThirdParamDerivative:
#ifdef DEBUG
              mexPrintf("ThirdParamDerivative\n");
#endif
              EQN_type = ExpressionType::ThirdParamDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(it_code->second)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable1();
              EQN_dvar2 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable2();
              EQN_dvar3 = static_cast<FNUMEXPR_ *>(it_code->second)->get_dvariable3();
              break;
            }
          break;
        case Tags::FLDV:
          //load a variable in the processor
          switch (static_cast<SymbolType>(static_cast<FLDV_ *>(it_code->second)->get_type()))
            {
            case SymbolType::parameter:
              var = static_cast<FLDV_ *>(it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDV Param[var=%d]\n", var);
              tmp_out << " params[" << var << "](" << params[var] << ")";
#endif
              Stack.push(params[var]);
              break;
            case SymbolType::endogenous:
              var = static_cast<FLDV_ *>(it_code->second)->get_pos();
              lag = static_cast<FLDV_ *>(it_code->second)->get_lead_lag();
#ifdef DEBUG
              if (evaluate)
                mexPrintf("FLDV y[var=%d, lag=%d, it_=%d], y_size=%d evaluate=%d, ya[%d]=%f\n", var, lag, it_, y_size, evaluate, (it_+lag)*y_size+var, ya[(it_+lag)*y_size+var]);
              else
                mexPrintf("FLDV y[var=%d, lag=%d, it_=%d], y_size=%d evaluate=%d, y[%d]=%f\n", var, lag, it_, y_size, evaluate, (it_+lag)*y_size+var, y[(it_+lag)*y_size+var]);
#endif
              if (evaluate)
                Stack.push(ya[(it_+lag)*y_size+var]);
              else
                Stack.push(y[(it_+lag)*y_size+var]);
#ifdef DEBUG
              tmp_out << " y[" << it_+lag << ", " << var << "](" << y[(it_+lag)*y_size+var] << ")";
#endif
              break;
            case SymbolType::exogenous:
              var = static_cast<FLDV_ *>(it_code->second)->get_pos();
              lag = static_cast<FLDV_ *>(it_code->second)->get_lead_lag();
#ifdef DEBUG
              mexPrintf("FLDV x[var=%d, lag=%d, it_=%d], nb_row_x=%d evaluate=%d x[%d]=%f\n", var, lag, it_, nb_row_x, evaluate, it_+lag+var*nb_row_x, x[it_+lag+var*nb_row_x]);
              //tmp_out << " x[" << it_+lag << ", " << var << "](" << x[it_+lag+var*nb_row_x] << ")";
#endif
              Stack.push(x[it_+lag+var*nb_row_x]);
              break;
            case SymbolType::exogenousDet:
              var = static_cast<FLDV_ *>(it_code->second)->get_pos();
              lag = static_cast<FLDV_ *>(it_code->second)->get_lead_lag();
              Stack.push(x[it_+lag+var*nb_row_xd]);
              break;
            case SymbolType::modelLocalVariable:
#ifdef DEBUG
              mexPrintf("FLDV a local variable in Block %d Stack.size()=%d", block_num, Stack.size());
              mexPrintf(" value=%f\n", Stack.top());
#endif
              break;
            default:
              mexPrintf("FLDV: Unknown variable type\n");
            }
          break;
        case Tags::FLDSV:
          //load a variable in the processor
          switch (static_cast<SymbolType>(static_cast<FLDSV_ *>(it_code->second)->get_type()))
            {
            case SymbolType::parameter:
              var = static_cast<FLDSV_ *>(it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV Param[var=%d]=%f\n", var, params[var]);
              tmp_out << " params[" << var << "](" << params[var] << ")";
#endif
              Stack.push(params[var]);
              break;
            case SymbolType::endogenous:
              var = static_cast<FLDSV_ *>(it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV y[var=%d]=%f\n", var, ya[var]);
              tmp_out << " y[" << var << "](" << y[var] << ")";
#endif
              if (evaluate)
                Stack.push(ya[var]);
              else
                Stack.push(y[var]);
              break;
            case SymbolType::exogenous:
              var = static_cast<FLDSV_ *>(it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV x[var=%d]\n", var);
              tmp_out << " x[" << var << "](" << x[var] << ")";
#endif
              Stack.push(x[var]);
              break;
            case SymbolType::exogenousDet:
              var = static_cast<FLDSV_ *>(it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV xd[var=%d]\n", var);
#endif
              Stack.push(x[var]);
              break;
            case SymbolType::modelLocalVariable:
#ifdef DEBUG
              mexPrintf("FLDSV a local variable in Block %d Stack.size()=%d", block_num, Stack.size());
              mexPrintf(" value=%f\n", Stack.top());
#endif
              break;
            default:
              mexPrintf("FLDSV: Unknown variable type\n");
            }
          break;
        case Tags::FLDVS:
          //load a variable in the processor
          switch (static_cast<SymbolType>(static_cast<FLDVS_ *>(it_code->second)->get_type()))
            {
            case SymbolType::parameter:
              var = static_cast<FLDVS_ *>(it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("params[%d]\n", var);
#endif
              Stack.push(params[var]);
              break;
            case SymbolType::endogenous:
              var = static_cast<FLDVS_ *>(it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDVS steady_y[%d]\n", var);
#endif
              Stack.push(steady_y[var]);
              break;
            case SymbolType::exogenous:
              var = static_cast<FLDVS_ *>(it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDVS x[%d] \n", var);
#endif
              Stack.push(x[var]);
              break;
            case SymbolType::exogenousDet:
              var = static_cast<FLDVS_ *>(it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDVS xd[%d]\n", var);
#endif
              Stack.push(x[var]);
              break;
            case SymbolType::modelLocalVariable:
#ifdef DEBUG
              mexPrintf("FLDVS a local variable in Block %d Stack.size()=%d", block_num, Stack.size());
              mexPrintf(" value=%f\n", Stack.top());
#endif
              break;
            default:
              mexPrintf("FLDVS: Unknown variable type\n");
            }
          break;
        case Tags::FLDT:
          //load a temporary variable in the processor
          var = static_cast<FLDT_ *>(it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("FLDT T[it_=%d var=%d, y_kmin=%d, y_kmax=%d == %d]=>%f\n", it_, var, y_kmin, y_kmax, var*(periods+y_kmin+y_kmax)+it_, T[var*(periods+y_kmin+y_kmax)+it_]);
          tmp_out << " T[" << it_ << ", " << var << "](" << T[var*(periods+y_kmin+y_kmax)+it_] << ")";
#endif
          Stack.push(T[var*(periods+y_kmin+y_kmax)+it_]);
          break;
        case Tags::FLDST:
          //load a temporary variable in the processor
          var = static_cast<FLDST_ *>(it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("FLDST T[%d]", var);
#endif
          Stack.push(T[var]);
#ifdef DEBUG
          mexPrintf("=%f\n", T[var]);
          tmp_out << " T[" << var << "](" << T[var] << ")";
#endif
          break;
        case Tags::FLDU:
          //load u variable in the processor
          var = static_cast<FLDU_ *>(it_code->second)->get_pos();
          var += Per_u_;
#ifdef DEBUG
          mexPrintf("FLDU u[%d]\n", var);
          tmp_out << " u[" << var << "](" << u[var] << ")";
#endif
          Stack.push(u[var]);
          break;
        case Tags::FLDSU:
          //load u variable in the processor
          var = static_cast<FLDSU_ *>(it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("FLDSU u[%d]\n", var);
          tmp_out << " u[" << var << "](" << u[var] << ")";
#endif
          Stack.push(u[var]);
          break;
        case Tags::FLDR:
          //load u variable in the processor
          var = static_cast<FLDR_ *>(it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("FLDR r[%d]\n", var);
#endif
          Stack.push(r[var]);
          break;
        case Tags::FLDZ:
          //load 0 in the processor
#ifdef DEBUG
          mexPrintf("FLDZ\n");
#endif
          Stack.push(0.0);
#ifdef DEBUG
          tmp_out << " 0";
#endif
          break;
        case Tags::FLDC:
          //load a numerical constant in the processor
          ll = static_cast<FLDC_ *>(it_code->second)->get_value();
#ifdef DEBUG
          mexPrintf("FLDC = %f\n", ll);
          tmp_out << " " << ll;
#endif

          Stack.push(ll);
          break;
        case Tags::FSTPV:
          //load a variable in the processor
          switch (static_cast<SymbolType>(static_cast<FSTPV_ *>(it_code->second)->get_type()))
            {
            case SymbolType::parameter:
              var = static_cast<FSTPV_ *>(it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FSTPV params[%d]\n", var);
#endif
              params[var] = Stack.top();
              Stack.pop();
              break;
            case SymbolType::endogenous:
              var = static_cast<FSTPV_ *>(it_code->second)->get_pos();
              lag = static_cast<FSTPV_ *>(it_code->second)->get_lead_lag();
              y[(it_+lag)*y_size+var] = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" y[%d, %d](%f)=%s\n", it_+lag, var, y[(it_+lag)*y_size+var], tmp_out.str().c_str());
              tmp_out.str("");
#endif
              Stack.pop();
              break;
            case SymbolType::exogenous:
              var = static_cast<FSTPV_ *>(it_code->second)->get_pos();
              lag = static_cast<FSTPV_ *>(it_code->second)->get_lead_lag();
              x[it_+lag+var*nb_row_x] = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" x[%d, %d](%f)=%s\n", it_+lag, var, x[it_+lag+var*nb_row_x], tmp_out.str().c_str());
              tmp_out.str("");
#endif

              Stack.pop();
              break;
            case SymbolType::exogenousDet:
              var = static_cast<FSTPV_ *>(it_code->second)->get_pos();
              lag = static_cast<FSTPV_ *>(it_code->second)->get_lead_lag();
              x[it_+lag+var*nb_row_xd] = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" x[%d, %d](%f)=%s\n", it_+lag, var, x[it_+lag+var*nb_row_xd], tmp_out.str().c_str());
              tmp_out.str("");
#endif
              Stack.pop();
              break;
            default:
              mexPrintf("FSTPV: Unknown variable type\n");
            }
          break;
        case Tags::FSTPSV:
          //load a variable in the processor
          switch (static_cast<SymbolType>(static_cast<FSTPSV_ *>(it_code->second)->get_type()))
            {
            case SymbolType::parameter:
              var = static_cast<FSTPSV_ *>(it_code->second)->get_pos();
              params[var] = Stack.top();
              Stack.pop();
              break;
            case SymbolType::endogenous:
              var = static_cast<FSTPSV_ *>(it_code->second)->get_pos();
              y[var] = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" y[%d](%f)=%s\n", var, y[var], tmp_out.str().c_str());
              tmp_out.str("");
#endif
              Stack.pop();
              break;
            case SymbolType::exogenous:
            case SymbolType::exogenousDet:
              var = static_cast<FSTPSV_ *>(it_code->second)->get_pos();
              x[var] = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" x[%d, %d](%f)=%s\n", it_+lag, var, x[var], tmp_out.str().c_str());
              tmp_out.str("");
#endif
              Stack.pop();
              break;
            default:
              mexPrintf("FSTPSV: Unknown variable type\n");
            }
          break;
        case Tags::FSTPT:
          //store in a temporary variable from the processor
#ifdef DEBUG
          mexPrintf("FSTPT\n");
#endif
          var = static_cast<FSTPT_ *>(it_code->second)->get_pos();
          T[var*(periods+y_kmin+y_kmax)+it_] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" T[%d, %d](%f)=%s\n", it_, var, T[var*(periods+y_kmin+y_kmax)+it_], tmp_out.str().c_str());
          tmp_out.str("");
#endif

          Stack.pop();
          break;
        case Tags::FSTPST:
          //store in a temporary variable from the processor
#ifdef DEBUG
          mexPrintf("FSTPST\n");
#endif
          var = static_cast<FSTPST_ *>(it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("var=%d\n", var);
#endif
          T[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" T[%d](%f)=%s\n", var, T[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;
        case Tags::FSTPU:
          //store in u variable from the processor
          var = static_cast<FSTPU_ *>(it_code->second)->get_pos();
          var += Per_u_;
#ifdef DEBUG
          mexPrintf("FSTPU\n");
          mexPrintf("var=%d\n", var);
#endif
          u[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" u[%d](%f)=%s\n", var, u[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;
        case Tags::FSTPSU:
          //store in u variable from the processor
          var = static_cast<FSTPSU_ *>(it_code->second)->get_pos();
#ifdef DEBUG
          /*if (var >= u_count_alloc || var < 0)
            mexPrintf("Erreur var=%d\n", var);*/
#endif
          u[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" u[%d](%f)=%s\n", var, u[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;
        case Tags::FSTPR:
          //store in residual variable from the processor
          var = static_cast<FSTPR_ *>(it_code->second)->get_pos();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf("FSTPR r[%d]", var);
          tmp_out.str("");
#endif
          r[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf("(%f)=%s\n", r[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;
        case Tags::FSTPG:
          //store in derivative (g) variable from the processor
#ifdef DEBUG
          mexPrintf("FSTPG\n");
          mexEvalString("drawnow;");
#endif
          var = static_cast<FSTPG_ *>(it_code->second)->get_pos();
          g1[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" g1[%d](%f)=%s\n", var, g1[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;

        case Tags::FSTPG2:
          //store in the jacobian matrix
          rr = Stack.top();
          if (EQN_type != ExpressionType::FirstEndoDerivative)
            throw FatalExceptionHandling(" in compute_block_time, impossible case " + to_string(static_cast<int>(EQN_type)) + " not implement in static jacobian\n");
          eq = static_cast<FSTPG2_ *>(it_code->second)->get_row();
          var = static_cast<FSTPG2_ *>(it_code->second)->get_col();
#ifdef DEBUG
          mexPrintf("FSTPG2 eq=%d, var=%d\n", eq, var);
          mexEvalString("drawnow;");
#endif
          jacob[eq + size*var] = rr;
          break;
        case Tags::FSTPG3:
          //store in derivative (g) variable from the processor
#ifdef DEBUG
          mexPrintf("FSTPG3 Evaluate=%d\n", evaluate);
          mexEvalString("drawnow;");
          if (!evaluate)
            {
              mexPrintf("impossible case!! \n");
              mexEvalString("drawnow;");
            }

#endif
          rr = Stack.top();
          switch (EQN_type)
            {
            case ExpressionType::FirstEndoDerivative:
              eq = static_cast<FSTPG3_ *>(it_code->second)->get_row();
              var = static_cast<FSTPG3_ *>(it_code->second)->get_col();
              lag = static_cast<FSTPG3_ *>(it_code->second)->get_lag();
              pos_col = static_cast<FSTPG3_ *>(it_code->second)->get_col_pos();
#ifdef DEBUG
              mexPrintf("Endo eq=%d, pos_col=%d, size=%d, jacob=%x\n", eq, pos_col, size, jacob);
              mexPrintf("jacob=%x\n", jacob);
#endif
              jacob[eq + size*pos_col] = rr;
              break;
            case ExpressionType::FirstOtherEndoDerivative:
              //eq = static_cast<FSTPG3_ *>(it_code->second)->get_row();
              eq = EQN_equation;
              var = static_cast<FSTPG3_ *>(it_code->second)->get_col();
              lag = static_cast<FSTPG3_ *>(it_code->second)->get_lag();
              pos_col = static_cast<FSTPG3_ *>(it_code->second)->get_col_pos();
#ifdef DEBUG
              mexPrintf("other_endo eq=%d, pos_col=%d, size=%d\n", eq, pos_col, size);
              mexEvalString("drawnow;");
#endif
              jacob_other_endo[eq + size*pos_col] = rr;
              break;
            case ExpressionType::FirstExoDerivative:
              //eq = static_cast<FSTPG3_ *>(it_code->second)->get_row();
              eq = EQN_equation;
              var = static_cast<FSTPG3_ *>(it_code->second)->get_col();
              lag = static_cast<FSTPG3_ *>(it_code->second)->get_lag();
              pos_col = static_cast<FSTPG3_ *>(it_code->second)->get_col_pos();
#ifdef DEBUG
              mexPrintf("Exo eq=%d, pos_col=%d, size=%d\n", eq, pos_col, size);
              mexEvalString("drawnow;");
#endif
              jacob_exo[eq + size*pos_col] = rr;
              break;
            case ExpressionType::FirstExodetDerivative:
              //eq = static_cast<FSTPG3_ *>(it_code->second)->get_row();
              eq = EQN_equation;
              var = static_cast<FSTPG3_ *>(it_code->second)->get_col();
              lag = static_cast<FSTPG3_ *>(it_code->second)->get_lag();
              pos_col = static_cast<FSTPG3_ *>(it_code->second)->get_col_pos();
#ifdef DEBUG
              mexPrintf("Exo det eq=%d, pos_col=%d, size=%d\n", eq, pos_col, size);
              mexEvalString("drawnow;");
#endif

              jacob_exo_det[eq + size*pos_col] = rr;
              break;
            default:
              throw FatalExceptionHandling(" in compute_block_time, variable " + to_string(static_cast<int>(EQN_type)) + " not used yet\n");
            }
// #ifdef DEBUG
//           tmp_out << "=>";
//           mexPrintf(" g1[%d](%f)=%s\n", var, g1[var], tmp_out.str().c_str());
//           tmp_out.str("");
// #endif
          Stack.pop();
          break;

        case Tags::FBINARY:
          op = static_cast<FBINARY_ *>(it_code->second)->get_op_type();
#ifdef DEBUG
          mexPrintf("FBINARY, op=%d\n", op);
#endif
          v2 = Stack.top();
          Stack.pop();
          v1 = Stack.top();
          Stack.pop();
          switch (static_cast<BinaryOpcode>(op))
            {
            case BinaryOpcode::plus:
              Stack.push(v1 + v2);
#ifdef DEBUG
              tmp_out << " |" << v1 << "+" << v2 << "|";
#endif
              break;
            case BinaryOpcode::minus:
              Stack.push(v1 - v2);
#ifdef DEBUG
              tmp_out << " |" << v1 << "-" << v2 << "|";
#endif
              break;
            case BinaryOpcode::times:
              Stack.push(v1 * v2);
#ifdef DEBUG
              tmp_out << " |" << v1 << "*" << v2 << "|";
#endif
              break;
            case BinaryOpcode::divide:
              double tmp;
#ifdef DEBUG
              mexPrintf("v1=%f / v2=%f\n", v1, v2);
#endif
              try
                {
                  tmp = divide(v1, v2);
                }
              catch (FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n", fpeh.GetErrorMsg().c_str(), error_location(evaluate, steady_state, size, block_num, it_, Per_u_).c_str());
                  go_on = false;
                }
              Stack.push(tmp);
#ifdef DEBUG
              tmp_out << " |" << v1 << "/" << v2 << "|";
#endif
              break;
            case BinaryOpcode::less:
              Stack.push(static_cast<double>(v1 < v2));
#ifdef DEBUG
              mexPrintf("v1=%f v2=%f v1 < v2 = %f\n", v1, v2, static_cast<double>(v1 < v2));
#endif
              break;
            case BinaryOpcode::greater:
              Stack.push(static_cast<double>(v1 > v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << ">" << v2 << "|";
#endif
              break;
            case BinaryOpcode::lessEqual:
              Stack.push(static_cast<double>(v1 <= v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << "<=" << v2 << "|";
#endif
              break;
            case BinaryOpcode::greaterEqual:
              Stack.push(static_cast<double>(v1 >= v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << ">=" << v2 << "|";
#endif
              break;
            case BinaryOpcode::equalEqual:
              Stack.push(static_cast<double>(v1 == v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << "==" << v2 << "|";
#endif
              break;
            case BinaryOpcode::different:
              Stack.push(static_cast<double>(v1 != v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << "!=" << v2 << "|";
#endif
              break;
            case BinaryOpcode::power:
#ifdef DEBUG
              mexPrintf("pow\n");
#endif
              try
                {
                  tmp = pow1(v1, v2);
                }
              catch (FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n", fpeh.GetErrorMsg().c_str(), error_location(evaluate, steady_state, size, block_num, it_, Per_u_).c_str());
                  go_on = false;
                }
              Stack.push(tmp);

#ifdef DEBUG
              tmp_out << " |" << v1 << "^" << v2 << "|";
#endif
              break;
            case BinaryOpcode::powerDeriv:
              {
                int derivOrder = static_cast<int>(nearbyint(Stack.top()));
                Stack.pop();
                try
                  {
                    if (fabs(v1) < near_zero && v2 > 0
                        && derivOrder > v2
                        && fabs(v2-nearbyint(v2)) < near_zero)
                      Stack.push(0.0);
                    else
                      {
                        double dxp = pow1(v1, v2-derivOrder);
                        for (int i = 0; i < derivOrder; i++)
                          dxp *= v2--;
                        Stack.push(dxp);
                      }
                  }
                catch (FloatingPointExceptionHandling &fpeh)
                  {
                    mexPrintf("%s      %s\n", fpeh.GetErrorMsg().c_str(), error_location(evaluate, steady_state, size, block_num, it_, Per_u_).c_str());
                    go_on = false;
                  }
              }

#ifdef DEBUG
              tmp_out << " |PowerDeriv(" << v1 << ", " << v2 << ")|";
#endif
              break;
            case BinaryOpcode::max:
              Stack.push(max(v1, v2));
#ifdef DEBUG
              tmp_out << " |max(" << v1 << "," << v2 << ")|";
#endif
              break;
            case BinaryOpcode::min:
              Stack.push(min(v1, v2));
#ifdef DEBUG
              tmp_out << " |min(" << v1 << "," << v2 << ")|";
#endif
              break;
            case BinaryOpcode::equal:
              // Nothing to do
              break;
            default:
              {
                mexPrintf("Error\n");
                throw FatalExceptionHandling(" in compute_block_time, unknown binary operator "
                                             + to_string(op) + "\n");
              }
            }
          break;
        case Tags::FUNARY:
          op = static_cast<FUNARY_ *>(it_code->second)->get_op_type();
          v1 = Stack.top();
          Stack.pop();
#ifdef DEBUG
          mexPrintf("FUNARY, op=%d\n", op);
#endif
          switch (static_cast<UnaryOpcode>(op))
            {
            case UnaryOpcode::uminus:
              Stack.push(-v1);
#ifdef DEBUG
              tmp_out << " |-(" << v1 << ")|";
#endif

              break;
            case UnaryOpcode::exp:
              Stack.push(exp(v1));
#ifdef DEBUG
              tmp_out << " |exp(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::log:
              double tmp;
              try
                {
                  tmp = log1(v1);
                }
              catch (FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n", fpeh.GetErrorMsg().c_str(), error_location(evaluate, steady_state, size, block_num, it_, Per_u_).c_str());
                  go_on = false;
                }
              Stack.push(tmp);

#ifdef DEBUG
              tmp_out << " |log(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::log10:
              try
                {
                  tmp = log10_1(v1);
                }
              catch (FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n", fpeh.GetErrorMsg().c_str(), error_location(evaluate, steady_state, size, block_num, it_, Per_u_).c_str());
                  go_on = false;
                }
              Stack.push(tmp);
#ifdef DEBUG
              tmp_out << " |log10(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::cos:
              Stack.push(cos(v1));
#ifdef DEBUG
              tmp_out << " |cos(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::sin:
              Stack.push(sin(v1));
#ifdef DEBUG
              tmp_out << " |sin(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::tan:
              Stack.push(tan(v1));
#ifdef DEBUG
              tmp_out << " |tan(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::acos:
              Stack.push(acos(v1));
#ifdef DEBUG
              tmp_out << " |acos(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::asin:
              Stack.push(asin(v1));
#ifdef DEBUG
              tmp_out << " |asin(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::atan:
              Stack.push(atan(v1));
#ifdef DEBUG
              tmp_out << " |atan(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::cosh:
              Stack.push(cosh(v1));
#ifdef DEBUG
              tmp_out << " |cosh(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::sinh:
              Stack.push(sinh(v1));
#ifdef DEBUG
              tmp_out << " |sinh(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::tanh:
              Stack.push(tanh(v1));
#ifdef DEBUG
              tmp_out << " |tanh(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::acosh:
              Stack.push(acosh(v1));
#ifdef DEBUG
              tmp_out << " |acosh(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::asinh:
              Stack.push(asinh(v1));
#ifdef DEBUG
              tmp_out << " |asinh(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::atanh:
              Stack.push(atanh(v1));
#ifdef DEBUG
              tmp_out << " |atanh(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::sqrt:
              Stack.push(sqrt(v1));
#ifdef DEBUG
              tmp_out << " |sqrt(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::erf:
              Stack.push(erf(v1));
#ifdef DEBUG
              tmp_out << " |erf(" << v1 << ")|";

#endif
              break;
            case UnaryOpcode::erfc:
              Stack.push(erfc(v1));
#ifdef DEBUG
              tmp_out << " |erfc(" << v1 << ")|";

#endif
              break;
            default:
              {
                mexPrintf("Error\n");
                throw FatalExceptionHandling(" in compute_block_time, unknown unary operator " + to_string(op) + "\n");
              }
            }
          break;
        case Tags::FTRINARY:
          op = static_cast<FTRINARY_ *>(it_code->second)->get_op_type();
          v3 = Stack.top();
          Stack.pop();
          v2 = Stack.top();
          Stack.pop();
          v1 = Stack.top();
          Stack.pop();
          switch (static_cast<TrinaryOpcode>(op))
            {
            case TrinaryOpcode::normcdf:
              Stack.push(0.5*(1+erf((v1-v2)/v3/M_SQRT2)));
#ifdef DEBUG
              tmp_out << " |normcdf(" << v1 << ", " << v2 << ", " << v3 << ")|";
#endif
              break;
            case TrinaryOpcode::normpdf:
              Stack.push(1/(v3*sqrt(2*M_PI)*exp(pow((v1-v2)/v3, 2)/2)));
#ifdef DEBUG
              tmp_out << " |normpdf(" << v1 << ", " << v2 << ", " << v3 << ")|";
#endif
              break;
            default:
              {
                mexPrintf("Error\n");
                throw FatalExceptionHandling(" in compute_block_time, unknown trinary operator " + to_string(op) + "\n");
              }
            }
          break;

        case Tags::FPUSH:
          break;

        case Tags::FCALL:
          {
#ifdef DEBUG
            mexPrintf("------------------------------\n");
            mexPrintf("CALL "); mexEvalString("drawnow;");
#endif
            auto *fc = static_cast<FCALL_ *>(it_code->second);
            string function_name = fc->get_function_name();
#ifdef DEBUG
            mexPrintf("function_name=%s ", function_name.c_str()); mexEvalString("drawnow;");
#endif
            unsigned int nb_input_arguments = fc->get_nb_input_arguments();
#ifdef DEBUG
            mexPrintf("nb_input_arguments=%d ", nb_input_arguments); mexEvalString("drawnow;");
#endif
            unsigned int nb_output_arguments = fc->get_nb_output_arguments();
#ifdef DEBUG
            mexPrintf("nb_output_arguments=%d\n", nb_output_arguments); mexEvalString("drawnow;");
#endif

            mxArray *output_arguments[3];
            string arg_func_name = fc->get_arg_func_name();
#ifdef DEBUG
            mexPrintf("arg_func_name.length() = %d\n", arg_func_name.length());
            mexPrintf("arg_func_name.c_str() = %s\n", arg_func_name.c_str());
#endif
            unsigned int nb_add_input_arguments = fc->get_nb_add_input_arguments();
            function_type = fc->get_function_type();
#ifdef DEBUG
            mexPrintf("function_type=%d ExternalFunctionWithoutDerivative=%d\n", function_type, ExternalFunctionType::withoutDerivative);
            mexEvalString("drawnow;");
#endif
            mxArray **input_arguments;
            switch (function_type)
              {
              case ExternalFunctionType::withoutDerivative:
              case ExternalFunctionType::withFirstDerivative:
              case ExternalFunctionType::withFirstAndSecondDerivative:
                {
                  input_arguments = static_cast<mxArray **>(mxMalloc(nb_input_arguments * sizeof(mxArray *)));
                  test_mxMalloc(input_arguments, __LINE__, __FILE__, __func__, nb_input_arguments * sizeof(mxArray *));
#ifdef DEBUG
                  mexPrintf("Stack.size()=%d\n", Stack.size());
                  mexEvalString("drawnow;");
#endif
                  for (unsigned int i = 0; i < nb_input_arguments; i++)
                    {
                      mxArray *vv = mxCreateDoubleScalar(Stack.top());
                      input_arguments[nb_input_arguments - i - 1] = vv;
                      Stack.pop();
                    }
                  if (mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str()))
                    throw FatalExceptionHandling(" external function: " + function_name + " not found");

                  double *rr = mxGetPr(output_arguments[0]);
                  Stack.push(*rr);
                  if (function_type == ExternalFunctionType::withFirstDerivative || function_type == ExternalFunctionType::withFirstAndSecondDerivative)
                    {
                      unsigned int indx = fc->get_indx();
                      double *FD1 = mxGetPr(output_arguments[1]);
                      size_t rows = mxGetN(output_arguments[1]);
                      for (unsigned int i = 0; i < rows; i++)
                        TEFD[{ indx, i }] = FD1[i];
                    }
                  if (function_type == ExternalFunctionType::withFirstAndSecondDerivative)
                    {
                      unsigned int indx = fc->get_indx();
                      double *FD2 = mxGetPr(output_arguments[2]);
                      size_t rows = mxGetM(output_arguments[2]);
                      size_t cols = mxGetN(output_arguments[2]);
                      unsigned int k = 0;
                      for (unsigned int j = 0; j < cols; j++)
                        for (unsigned int i = 0; i < rows; i++)
                          TEFDD[{ indx, i, j }] = FD2[k++];
                    }
                }
                break;
              case ExternalFunctionType::numericalFirstDerivative:
                {
                  input_arguments = static_cast<mxArray **>(mxMalloc((nb_input_arguments+1+nb_add_input_arguments) * sizeof(mxArray *)));
                  test_mxMalloc(input_arguments, __LINE__, __FILE__, __func__, (nb_input_arguments+1+nb_add_input_arguments) * sizeof(mxArray *));
                  mxArray *vv = mxCreateString(arg_func_name.c_str());
                  input_arguments[0] = vv;
                  vv = mxCreateDoubleScalar(fc->get_row());
                  input_arguments[1] = vv;
                  vv = mxCreateCellMatrix(1, nb_add_input_arguments);
                  for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                    {
                      double rr = Stack.top();
#ifdef DEBUG
                      mexPrintf("i=%d rr = %f Stack.size()=%d\n", i, rr, Stack.size());
#endif
                      mxSetCell(vv, nb_add_input_arguments - (i+1), mxCreateDoubleScalar(rr));
                      Stack.pop();
                    }
                  input_arguments[nb_input_arguments+nb_add_input_arguments] = vv;
#ifdef DEBUG
                  mexCallMATLAB(0, nullptr, 1, &input_arguments[0], "disp");
                  mexCallMATLAB(0, nullptr, 1, &input_arguments[1], "disp");
                  mexCallMATLAB(0, nullptr, 1, &input_arguments[2], "celldisp");
                  mexPrintf("OK\n");
                  mexEvalString("drawnow;");
#endif
                  nb_input_arguments = 3;
                  if (mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str()))
                    throw FatalExceptionHandling(" external function: " + function_name + " not found");
                  double *rr = mxGetPr(output_arguments[0]);
#ifdef DEBUG
                  mexPrintf("*rr=%f\n", *rr);
#endif
                  Stack.push(*rr);
                }
                break;
              case ExternalFunctionType::firstDerivative:
                {
                  input_arguments = static_cast<mxArray **>(mxMalloc(nb_input_arguments * sizeof(mxArray *)));
                  test_mxMalloc(input_arguments, __LINE__, __FILE__, __func__, nb_input_arguments * sizeof(mxArray *));
                  for (unsigned int i = 0; i < nb_input_arguments; i++)
                    {
                      mxArray *vv = mxCreateDoubleScalar(Stack.top());
                      input_arguments[(nb_input_arguments - 1) - i] = vv;
                      Stack.pop();
                    }
                  if (mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str()))
                    throw FatalExceptionHandling(" external function: " + function_name + " not found");
                  unsigned int indx = fc->get_indx();
                  double *FD1 = mxGetPr(output_arguments[0]);
                  size_t rows = mxGetN(output_arguments[0]);
                  for (unsigned int i = 0; i < rows; i++)
                    TEFD[{ indx, i }] = FD1[i];
                }
                break;
              case ExternalFunctionType::numericalSecondDerivative:
                {
                  input_arguments = static_cast<mxArray **>(mxMalloc((nb_input_arguments+1+nb_add_input_arguments) * sizeof(mxArray *)));
                  test_mxMalloc(input_arguments, __LINE__, __FILE__, __func__, (nb_input_arguments+1+nb_add_input_arguments) * sizeof(mxArray *));
                  mxArray *vv = mxCreateString(arg_func_name.c_str());
                  input_arguments[0] = vv;
                  vv = mxCreateDoubleScalar(fc->get_row());
                  input_arguments[1] = vv;
                  vv = mxCreateDoubleScalar(fc->get_col());
                  input_arguments[2] = vv;
                  vv = mxCreateCellMatrix(1, nb_add_input_arguments);
                  for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                    {
                      double rr = Stack.top();
#ifdef DEBUG
                      mexPrintf("i=%d rr = %f\n", i, rr);
#endif
                      mxSetCell(vv, (nb_add_input_arguments - 1) - i, mxCreateDoubleScalar(rr));
                      Stack.pop();
                    }
                  input_arguments[nb_input_arguments+nb_add_input_arguments] = vv;
#ifdef DEBUG
                  mexCallMATLAB(0, nullptr, 1, &input_arguments[0], "disp");
                  mexCallMATLAB(0, nullptr, 1, &input_arguments[1], "disp");
                  mexCallMATLAB(0, nullptr, 1, &input_arguments[2], "celldisp");
                  mexPrintf("OK\n");
                  mexEvalString("drawnow;");
#endif
                  nb_input_arguments = 3;
                  if (mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str()))
                    throw FatalExceptionHandling(" external function: " + function_name + " not found");
                  double *rr = mxGetPr(output_arguments[0]);
                  Stack.push(*rr);
                }
                break;
              case ExternalFunctionType::secondDerivative:
                {
                  input_arguments = static_cast<mxArray **>(mxMalloc(nb_input_arguments * sizeof(mxArray *)));
                  test_mxMalloc(input_arguments, __LINE__, __FILE__, __func__, nb_input_arguments * sizeof(mxArray *));
                  for (unsigned int i = 0; i < nb_input_arguments; i++)
                    {
                      mxArray *vv = mxCreateDoubleScalar(Stack.top());
                      input_arguments[i] = vv;
                      Stack.pop();
                    }
                  if (mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str()))
                    throw FatalExceptionHandling(" external function: " + function_name + " not found");
                  unsigned int indx = fc->get_indx();
                  double *FD2 = mxGetPr(output_arguments[2]);
                  size_t rows = mxGetM(output_arguments[0]);
                  size_t cols = mxGetN(output_arguments[0]);
                  unsigned int k = 0;
                  for (unsigned int j = 0; j < cols; j++)
                    for (unsigned int i = 0; i < rows; i++)
                      TEFDD[{ indx, i, j }] = FD2[k++];
                }
                break;
              }
          }
          break;
        case Tags::FSTPTEF:
          var = static_cast<FSTPTEF_ *>(it_code->second)->get_number();
#ifdef DEBUG
          mexPrintf("FSTPTEF\n");
          mexPrintf("var=%d Stack.size()=%d\n", var, Stack.size());
#endif
          TEF[var-1] = Stack.top();
#ifdef DEBUG
          mexPrintf("FSTP TEF[var-1]=%f done\n", TEF[var-1]);
          mexEvalString("drawnow;");
#endif
          Stack.pop();
          break;
        case Tags::FLDTEF:
          var = static_cast<FLDTEF_ *>(it_code->second)->get_number();
#ifdef DEBUG
          mexPrintf("FLDTEF\n");
          mexPrintf("var=%d Stack.size()=%d\n", var, Stack.size());
          mexPrintf("FLD TEF[var-1]=%f done\n", TEF[var-1]);
          mexEvalString("drawnow;");
#endif
          Stack.push(TEF[var-1]);
          break;
        case Tags::FSTPTEFD:
          {
            unsigned int indx = static_cast<FSTPTEFD_ *>(it_code->second)->get_indx();
            unsigned int row = static_cast<FSTPTEFD_ *>(it_code->second)->get_row();
#ifdef DEBUG
            mexPrintf("FSTPTEFD\n");
            mexPrintf("indx=%d Stack.size()=%d\n", indx, Stack.size());
#endif
            if (function_type == ExternalFunctionType::numericalFirstDerivative)
              {
                TEFD[{ indx, row-1 }] = Stack.top();
#ifdef DEBUG
                mexPrintf("FSTP TEFD[{ indx, row }]=%f done\n", TEFD[{ indx, row-1 }]);
                mexEvalString("drawnow;");
#endif
                Stack.pop();
              }
          }

          break;
        case Tags::FLDTEFD:
          {
            unsigned int indx = static_cast<FLDTEFD_ *>(it_code->second)->get_indx();
            unsigned int row = static_cast<FLDTEFD_ *>(it_code->second)->get_row();
#ifdef DEBUG
            mexPrintf("FLDTEFD\n");
            mexPrintf("indx=%d row=%d Stack.size()=%d\n", indx, row, Stack.size());
            mexPrintf("FLD TEFD[{ indx, row }]=%f done\n", TEFD[{ indx, row-1 }]);
            mexEvalString("drawnow;");
#endif
            Stack.push(TEFD[{ indx, row-1 }]);
          }
          break;
        case Tags::FSTPTEFDD:
          {
            unsigned int indx = static_cast<FSTPTEFDD_ *>(it_code->second)->get_indx();
            unsigned int row = static_cast<FSTPTEFDD_ *>(it_code->second)->get_row();
            unsigned int col = static_cast<FSTPTEFDD_ *>(it_code->second)->get_col();
#ifdef DEBUG
            mexPrintf("FSTPTEFD\n");
            mexPrintf("indx=%d Stack.size()=%d\n", indx, Stack.size());
#endif
            if (function_type == ExternalFunctionType::numericalSecondDerivative)
              {
                TEFDD[{ indx, row-1, col-1 }] = Stack.top();
#ifdef DEBUG
                mexPrintf("FSTP TEFDD[{ indx, row, col }]=%f done\n", TEFDD[{ indx, row, col }]);
                mexEvalString("drawnow;");
#endif
                Stack.pop();
              }
          }

          break;
        case Tags::FLDTEFDD:
          {
            unsigned int indx = static_cast<FLDTEFDD_ *>(it_code->second)->get_indx();
            unsigned int row = static_cast<FLDTEFDD_ *>(it_code->second)->get_row();
            unsigned int col = static_cast<FSTPTEFDD_ *>(it_code->second)->get_col();
#ifdef DEBUG
            mexPrintf("FLDTEFD\n");
            mexPrintf("indx=%d Stack.size()=%d\n", indx, Stack.size());
            mexPrintf("FLD TEFD[{ indx, row, col }]=%f done\n", TEFDD[{ indx, row, col }]);
            mexEvalString("drawnow;");
#endif
            Stack.push(TEFDD[{ indx, row-1, col-1 }]);
          }
          break;
        case Tags::FCUML:
          v1 = Stack.top();
          Stack.pop();
          v2 = Stack.top();
          Stack.pop();
          Stack.push(v1+v2);
          break;
        case Tags::FENDBLOCK:
          //it's the block end
#ifdef DEBUG
          mexPrintf("FENDBLOCK\n");
#endif
          go_on = false;
          break;
        case Tags::FBEGINBLOCK:
          mexPrintf("Impossible case in Bytecode\n");
          break;
        case Tags::FENDEQU:
          if (no_derivative)
            go_on = false;
          break;
        case Tags::FJMPIFEVAL:
          if (evaluate)
            {
#ifdef DEBUG
              mexPrintf("FJMPIFEVAL length=%d\n", static_cast<FJMPIFEVAL_ *>(it_code->second)->get_pos());
              mexEvalString("drawnow;");
#endif
              it_code += static_cast<FJMPIFEVAL_ *>(it_code->second)->get_pos() /* - 1*/;
            }
          break;
        case Tags::FJMP:
#ifdef DEBUG
          mexPrintf("FJMP length=%d\n", static_cast<FJMP_ *>(it_code->second)->get_pos());
          mexEvalString("drawnow;");
#endif
          it_code += static_cast<FJMP_ *>(it_code->second)->get_pos() /*- 1 */;
          break;
        case Tags::FOK:
          op = static_cast<FOK_ *>(it_code->second)->get_arg();
          if (Stack.size() > 0)
            throw FatalExceptionHandling(" in compute_block_time, stack not empty\n");
          break;
        default:
          throw FatalExceptionHandling(" in compute_block_time, unknown opcode " + to_string(static_cast<int>(it_code->first)) + "\n");
        }
      it_code++;
    }
#ifdef DEBUG
  mexPrintf("==> end of compute_block_time Block = %d\n", block_num);
  mexEvalString("drawnow;");
#endif
}

void
Evaluate::evaluate_over_periods(bool forward)
{
  if (steady_state)
    compute_block_time(0, false, false);
  else
    {
      auto begining = it_code;
      if (forward)
        {
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              compute_block_time(0, false, false);
            }
          it_ = periods+y_kmin-1; // Do not leave it_ in inconsistent state
        }
      else
        {
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              compute_block_time(0, false, false);
            }
          it_ = y_kmin; // Do not leave it_ in inconsistent state (see #1727)
        }
    }
}

void
Evaluate::solve_simple_one_periods()
{
  bool cvg = false;
  int iter = 0;
  double ya;
  double slowc = 1;
  res1 = 0;
  while (!(cvg || iter > maxit_))
    {
      it_code = start_code;
      Per_y_ = it_*y_size;
      ya = y[Block_Contain[0].Variable + Per_y_];
      compute_block_time(0, false, false);
      if (!isfinite(res1))
        {
          res1 = std::numeric_limits<double>::quiet_NaN();
          while ((isinf(res1) || isnan(res1)) && (slowc > 1e-9))
            {
              it_code = start_code;
              compute_block_time(0, false, false);
              if (!isfinite(res1))
                {
                  slowc /= 1.5;
                  mexPrintf("Reducing the path length in Newton step slowc=%f\n", slowc);
                  y[Block_Contain[0].Variable + Per_y_] = ya - slowc * divide(r[0], g1[0]);
                }
            }
        }
      double rr;
      rr = r[0];
      cvg = (fabs(rr) < solve_tolf);
      if (cvg)
        continue;
      try
        {
          y[Block_Contain[0].Variable + Per_y_] += -slowc *divide(rr, g1[0]);
        }
      catch (FloatingPointExceptionHandling &fpeh)
        {
          mexPrintf("%s      \n", fpeh.GetErrorMsg().c_str());
          mexPrintf("      Singularity in block %d", block_num+1);
        }
      iter++;
    }
  if (!cvg)
    throw FatalExceptionHandling(" in Solve Forward simple, convergence not achieved in block "
                                 + to_string(block_num+1) + ", after " + to_string(iter) + " iterations\n");
}

void
Evaluate::solve_simple_over_periods(bool forward)
{
  g1 = static_cast<double *>(mxMalloc(sizeof(double)));
  test_mxMalloc(g1, __LINE__, __FILE__, __func__, sizeof(double));
  r = static_cast<double *>(mxMalloc(sizeof(double)));
  test_mxMalloc(r, __LINE__, __FILE__, __func__, sizeof(double));
  start_code = it_code;
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
Evaluate::set_block(int size_arg, int type_arg, string file_name_arg, string bin_base_name_arg, int block_num_arg,
                    bool is_linear_arg, int symbol_table_endo_nbr_arg, int Block_List_Max_Lag_arg, int Block_List_Max_Lead_arg, int u_count_int_arg, int block_arg)
{
  size = size_arg;
  type = type_arg;
  file_name = move(file_name_arg);
  bin_base_name = move(bin_base_name_arg);
  block_num = block_num_arg;
  is_linear = is_linear_arg;
  symbol_table_endo_nbr = symbol_table_endo_nbr_arg;
  Block_List_Max_Lag = Block_List_Max_Lag_arg;
  Block_List_Max_Lead = Block_List_Max_Lead_arg;
  u_count_int = u_count_int_arg;
  block = block_arg;
}

void
Evaluate::evaluate_complete(bool no_derivatives)
{
  it_code = start_code;
  compute_block_time(0, false, no_derivatives);
}

void
Evaluate::compute_complete_2b(bool no_derivatives, double *_res1, double *_res2, double *_max_res, int *_max_res_idx)
{
  res1 = 0;
  *_res1 = 0;
  *_res2 = 0;
  *_max_res = 0;
  for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
    {
      Per_u_ = (it_-y_kmin)*u_count_int;
      Per_y_ = it_*y_size;
      it_code = start_code;
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

bool
Evaluate::compute_complete(bool no_derivatives, double &_res1, double &_res2, double &_max_res, int &_max_res_idx)
{
  bool result;
  res1 = 0;
  it_code = start_code;
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
Evaluate::compute_complete(double lambda, double *crit)
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
  mexPrintf("  lambda=%e, res2=%e\n", lambda, res2_);
  *crit = res2_/2;
  return true;
}
