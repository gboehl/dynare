/*
 * Copyright © 2007-2022 Dynare Team
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

#ifndef ERROR_HANDLING
#define ERROR_HANDLING

#include <vector>
#include <utility>
#include <string>
#include <map>
#include <tuple>
#include <cstddef>
#include <sstream>
#include <iostream>
#include <stack>
#define _USE_MATH_DEFINES
#include <cmath>

#include "dynmex.h"

#define BYTE_CODE
#include "CodeInterpreter.hh"

using namespace std;

constexpr int NO_ERROR_ON_EXIT = 0, ERROR_ON_EXIT = 1;

using code_liste_type = vector<pair<Tags, void *>>;
using it_code_type = code_liste_type::const_iterator;

class GeneralExceptionHandling
{
  string ErrorMsg;
public:
  GeneralExceptionHandling(string ErrorMsg_arg) : ErrorMsg{move(ErrorMsg_arg)}
  {
  };
  inline string
  GetErrorMsg()
  {
    return ErrorMsg;
  }
  inline void
  completeErrorMsg(const string &ErrorMsg_arg)
  {
    ErrorMsg += ErrorMsg_arg;
  }
};

class FloatingPointExceptionHandling : public GeneralExceptionHandling
{
public:
  FloatingPointExceptionHandling(const string &value) : GeneralExceptionHandling("Floating point error in bytecode: " + value)
  {
  }
};

class LogExceptionHandling : public FloatingPointExceptionHandling
{
  double value;
public:
  LogExceptionHandling(double value_arg) : FloatingPointExceptionHandling("log(X)"),
                                           value(value_arg)
  {
    completeErrorMsg(" with X=" + to_string(value) + "\n");
  }
};

class Log10ExceptionHandling : public FloatingPointExceptionHandling
{
  double value;
public:
  Log10ExceptionHandling(double value_arg) : FloatingPointExceptionHandling("log10(X)"),
                                             value(value_arg)
  {
    completeErrorMsg(" with X=" + to_string(value) + "\n");
  }
};

class DivideExceptionHandling : public FloatingPointExceptionHandling
{
  double value1, value2;
public:
  DivideExceptionHandling(double value1_arg, double value2_arg) : FloatingPointExceptionHandling("a/X"),
                                                                  value1(value1_arg),
                                                                  value2(value2_arg)
  {
    completeErrorMsg(" with X=" + to_string(value2) + "\n");
  }
};

class PowExceptionHandling : public FloatingPointExceptionHandling
{
  double value1, value2;
public:
  PowExceptionHandling(double value1_arg, double value2_arg) : FloatingPointExceptionHandling("X^a"),
                                                               value1(value1_arg),
                                                               value2(value2_arg)
  {
    if (fabs(value1) > 1e-10)
      completeErrorMsg(" with X=" + to_string(value1) + "\n");
    else
      completeErrorMsg(" with X=" + to_string(value1) + " and a=" + to_string(value2) + "\n");
  };
};

class UserExceptionHandling : public GeneralExceptionHandling
{
  double value;
public:
  UserExceptionHandling() : GeneralExceptionHandling("Fatal error in bytecode:")
  {
    completeErrorMsg(" User break\n");
  };
};

class FatalExceptionHandling : public GeneralExceptionHandling
{
public:
  FatalExceptionHandling(const string &ErrorMsg_arg)
    : GeneralExceptionHandling("Fatal error in bytecode:")
  {
    completeErrorMsg(ErrorMsg_arg);
  };
  FatalExceptionHandling() : GeneralExceptionHandling("")
  {
  };
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
#ifdef MATLAB_MEX_FILE
extern "C" bool utIsInterruptPending();
#endif

class ErrorMsg
{
private:
  bool is_load_variable_list;

public:
  double *y, *ya;
  int y_size;
  double *T;
  int nb_row_xd, nb_row_x, col_x, col_y;
  int y_kmin, y_kmax, periods;
  double *x, *params;
  double *u;
  double *steady_y, *steady_x;
  double *g2, *g1, *r, *res;
  vector<s_plan> splan, spfplan;
  vector<mxArray *> jacobian_block, jacobian_other_endo_block, jacobian_exo_block, jacobian_det_exo_block;
  map<unsigned int, double> TEF;
  map<pair<unsigned int, unsigned int>, double> TEFD;
  map<tuple<unsigned int, unsigned int, unsigned int>, double> TEFDD;

  ExpressionType EQN_type;
  it_code_type it_code_expr;
  size_t endo_name_length; // Maximum length of endogenous names
  vector<string> P_endo_names, P_exo_names, P_param_names;
  unsigned int EQN_equation, EQN_block, EQN_block_number;
  unsigned int EQN_dvar1, EQN_dvar2, EQN_dvar3;
  vector<tuple<string, SymbolType, unsigned int>> Variable_list;

  inline
  ErrorMsg()
  {
    mxArray *M_ = mexGetVariable("global", "M_");
    if (!M_)
      mexErrMsgTxt("Can't find global variable M_");

    auto get_field_names = [&](const char *symbol_type)
    {
      vector<string> r;
      if (mxGetFieldNumber(M_, symbol_type) != -1)
        {
          auto M_field = mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, symbol_type));
          if (!mxIsCell(M_field))
            mexErrMsgTxt((string{"M_."} + symbol_type + " is not a cell array").c_str());
          for (size_t i = 0; i < mxGetNumberOfElements(M_field); i++)
            {
              const mxArray *cell_mx = mxGetCell(M_field, i);
              if (!(cell_mx && mxIsChar(cell_mx)))
                mexErrMsgTxt((string{"M_."} + symbol_type + " contains a cell which is not a character array").c_str());
              r.emplace_back(mxArrayToString(cell_mx));
            }
        }
      return r;
    };
    P_endo_names = get_field_names("endo_names");
    P_exo_names = get_field_names("exo_names");
    P_param_names = get_field_names("param_names");

    endo_name_length = 0;
    for (const auto &n : P_endo_names)
      endo_name_length = max(endo_name_length, n.size());

    is_load_variable_list = false;
  }

  inline string
  add_underscore_to_fpe(const string &str)
  {
    string temp;
    int pos1 = -1, pos2 = -1;
    string tmp_n(str.length(), ' ');
    string dollar{"$"}, pound{"£"}, tilde{"~"};
    for (const char & i : str)
      {
        if (dollar.compare(&i) != 0 && pound.compare(&i) != 0)
          temp += i;
        else
          {
            if (dollar.compare(&i) == 0)
              pos1 = static_cast<int>(temp.length());
            else
              pos2 = static_cast<int>(temp.length());
            if (pos1 >= 0 && pos2 >= 0)
              {
                tmp_n.erase(pos1, pos2-pos1+1);
                tmp_n.insert(pos1, pos2-pos1, tilde[0]);
                pos1 = pos2 = -1;
              }
          }
      }
    temp += "\n" + tmp_n;
    return temp;
  }

  inline void
  load_variable_list()
  {
    for (unsigned int variable_num = 0; variable_num < P_endo_names.size(); variable_num++)
      Variable_list.emplace_back(P_endo_names[variable_num], SymbolType::endogenous, variable_num);
    for (unsigned int variable_num = 0; variable_num < P_exo_names.size(); variable_num++)
      Variable_list.emplace_back(P_exo_names[variable_num], SymbolType::exogenous, variable_num);
  }

  inline int
  get_ID(const string &variable_name, SymbolType *variable_type)
  {
    if (!is_load_variable_list)
      {
        load_variable_list();
        is_load_variable_list = true;
      }
    size_t n = Variable_list.size();
    int i = 0;
    bool notfound = true;
    while (notfound && i < static_cast<int>(n))
      {
        if (variable_name == get<0>(Variable_list[i]))
          {
            notfound = false;
            *variable_type = get<1>(Variable_list[i]);
            return get<2>(Variable_list[i]);
          }
        i++;
      }
    return -1;
  }

  inline string
  get_variable(SymbolType variable_type, unsigned int variable_num) const
  {
    switch (variable_type)
      {
      case SymbolType::endogenous:
        if (variable_num < P_endo_names.size())
          return P_endo_names[variable_num];
        else
          mexPrintf("=> Unknown endogenous variable # %d", variable_num);
        break;
      case SymbolType::exogenous:
      case SymbolType::exogenousDet:
        if (variable_num < P_exo_names.size())
          return P_exo_names[variable_num];
        else
          mexPrintf("=> Unknown exogenous variable # %d", variable_num);
        break;
      case SymbolType::parameter:
        if (variable_num < P_param_names.size())
          return P_param_names[variable_num];
        else
          mexPrintf("=> Unknown parameter # %d", variable_num);
        break;
      default:
        break;
      }
    cerr << "ErrorHandling::get_variable: Internal error";
    exit(EXIT_FAILURE); // Silence GCC warning
  }

  inline string
  error_location(bool evaluate, bool steady_state, int size, int block_num, int it_, int Per_u_)
  {
    ostringstream Error_loc;
    if (!steady_state)
      switch (EQN_type)
        {
        case ExpressionType::TemporaryTerm:
          if (EQN_block_number > 1)
            Error_loc << "temporary term " << EQN_equation+1 << " in block " << EQN_block+1 << " at time " << it_;
          else
            Error_loc << "temporary term " << EQN_equation+1 << " at time " << it_;
          break;
        case ExpressionType::ModelEquation:
          if (EQN_block_number > 1)
            Error_loc << "equation " << EQN_equation+1 << " in block " << EQN_block+1 << " at time " << it_;
          else
            Error_loc << "equation " << EQN_equation+1 << " at time " << it_;
          break;
        case ExpressionType::FirstEndoDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to endogenous variable " << get_variable(SymbolType::endogenous, EQN_dvar1) << " at time " << it_;
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to endogenous variable " << get_variable(SymbolType::endogenous, EQN_dvar1) << " at time " << it_;
          break;
        case ExpressionType::FirstOtherEndoDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to other endogenous variable "  << get_variable(SymbolType::endogenous, EQN_dvar1) << " at time " << it_;
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to other endogenous variable " << get_variable(SymbolType::endogenous, EQN_dvar1) << " at time " << it_;
          break;
        case ExpressionType::FirstExoDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to exogenous variable "  << get_variable(SymbolType::endogenous, EQN_dvar1) << " at time " << it_;
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to exogenous variable " << get_variable(SymbolType::endogenous, EQN_dvar1) << " at time " << it_;
          break;
        case ExpressionType::FirstExodetDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to deterministic exogenous variable "  << get_variable(SymbolType::endogenous, EQN_dvar1) << " at time " << it_;
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to deterministic exogenous variable " << get_variable(SymbolType::endogenous, EQN_dvar1) << " at time " << it_;
          break;
        case ExpressionType::FirstParamDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to parameter "  << get_variable(SymbolType::endogenous, EQN_dvar1) << " at time " << it_;
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to parameter " << get_variable(SymbolType::endogenous, EQN_dvar1) << " at time " << it_;
          break;
        default:
          return "???";
        }
    else
      switch (EQN_type)
        {
        case ExpressionType::TemporaryTerm:
          if (EQN_block_number > 1)
            Error_loc << "temporary term " << EQN_equation+1 << " in block " << EQN_block+1;
          else
            Error_loc << "temporary term " << EQN_equation+1;
          break;
        case ExpressionType::ModelEquation:
          if (EQN_block_number > 1)
            Error_loc << "equation " << EQN_equation+1 << " in block " << EQN_block+1;
          else
            Error_loc << "equation " << EQN_equation+1;
          break;
        case ExpressionType::FirstEndoDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to endogenous variable "  << get_variable(SymbolType::endogenous, EQN_dvar1);
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to endogenous variable " << get_variable(SymbolType::endogenous, EQN_dvar1);
          break;
        case ExpressionType::FirstOtherEndoDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to other endogenous variable "  << get_variable(SymbolType::endogenous, EQN_dvar1);
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to other endogenous variable " << get_variable(SymbolType::endogenous, EQN_dvar1);
          break;
        case ExpressionType::FirstExoDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to exogenous variable "  << get_variable(SymbolType::endogenous, EQN_dvar1);
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to exogenous variable " << get_variable(SymbolType::endogenous, EQN_dvar1);
          break;
        case ExpressionType::FirstExodetDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to deterministic exogenous variable "  << get_variable(SymbolType::endogenous, EQN_dvar1);
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to deterministic exogenous variable " << get_variable(SymbolType::endogenous, EQN_dvar1);
          break;
        case ExpressionType::FirstParamDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to parameter "  << get_variable(SymbolType::endogenous, EQN_dvar1);
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to parameter " << get_variable(SymbolType::endogenous, EQN_dvar1);
          break;
        default:
          return ("???");
        }
    it_code_type it_code_ret;
    Error_loc << endl << add_underscore_to_fpe("      " + print_expression(it_code_expr, evaluate, size, block_num, steady_state, Per_u_, it_, it_code_ret, true));
    return Error_loc.str();
  }

  inline string
  print_expression(it_code_type it_code, bool evaluate, int size, int block_num, bool steady_state, int Per_u_, int it_, it_code_type &it_code_ret, bool compute) const
  {
    int var, lag = 0, op, eq;
    stack<string> Stack;
    stack<double> Stackf;
    ostringstream tmp_out, tmp_out2;
    string v1, v2, v3;
    double v1f, v2f, v3f = 0.0;
    bool go_on = true;
    double ll;
    ExpressionType equation_type = ExpressionType::TemporaryTerm;
    size_t found;
    double *jacob = nullptr, *jacob_other_endo = nullptr, *jacob_exo = nullptr, *jacob_exo_det = nullptr;
    ExternalFunctionType function_type = ExternalFunctionType::withoutDerivative;

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

    while (go_on)
      {
#ifdef MATLAB_MEX_FILE
        if (utIsInterruptPending())
          throw UserExceptionHandling();
#endif
        switch (it_code->first)
          {
          case Tags::FNUMEXPR:
            switch (static_cast<FNUMEXPR_ *>(it_code->second)->get_expression_type())
              {
              case ExpressionType::TemporaryTerm:
                equation_type = ExpressionType::TemporaryTerm;
                break;
              case ExpressionType::ModelEquation:
                equation_type = ExpressionType::ModelEquation;
                break;
              case ExpressionType::FirstEndoDerivative:
                equation_type = ExpressionType::FirstEndoDerivative;
                break;
              case ExpressionType::FirstOtherEndoDerivative:
                equation_type = ExpressionType::FirstOtherEndoDerivative;
                break;
              case ExpressionType::FirstExoDerivative:
                equation_type = ExpressionType::FirstExoDerivative;
                break;
              case ExpressionType::FirstExodetDerivative:
                equation_type = ExpressionType::FirstExodetDerivative;
                break;
              case ExpressionType::FirstParamDerivative:
                equation_type = ExpressionType::FirstParamDerivative;
                break;
              case ExpressionType::SecondEndoDerivative:
                equation_type = ExpressionType::SecondEndoDerivative;
                break;
              case ExpressionType::SecondExoDerivative:
                equation_type = ExpressionType::SecondExoDerivative;
                break;
              case ExpressionType::SecondExodetDerivative:
                equation_type = ExpressionType::SecondExodetDerivative;
                break;
              case ExpressionType::SecondParamDerivative:
                equation_type = ExpressionType::SecondExodetDerivative;
                break;
              case ExpressionType::ThirdEndoDerivative:
                equation_type = ExpressionType::ThirdEndoDerivative;
                break;
              case ExpressionType::ThirdExoDerivative:
                equation_type = ExpressionType::ThirdExoDerivative;
                break;
              case ExpressionType::ThirdExodetDerivative:
                equation_type = ExpressionType::ThirdExodetDerivative;
                break;
              case ExpressionType::ThirdParamDerivative:
                equation_type = ExpressionType::ThirdExodetDerivative;
                break;
              default:
                ostringstream tmp;
                tmp << " in print_expression, expression type " << static_cast<int>(static_cast<FNUMEXPR_ *>(it_code->second)->get_expression_type()) << " not implemented yet\n";
                throw FatalExceptionHandling(tmp.str());
              }
            break;
          case Tags::FLDV:
            //load a variable in the processor
            switch (static_cast<SymbolType>(static_cast<FLDV_ *>(it_code->second)->get_type()))
              {
              case SymbolType::parameter:
                var = static_cast<FLDV_ *>(it_code->second)->get_pos();
                Stack.push(get_variable(SymbolType::parameter, var));
                if (compute)
                  Stackf.push(params[var]);
                break;
              case SymbolType::endogenous:
                var = static_cast<FLDV_ *>(it_code->second)->get_pos();
                lag = static_cast<FLDV_ *>(it_code->second)->get_lead_lag();
                tmp_out.str("");
                if (lag > 0)
                  tmp_out << get_variable(SymbolType::endogenous, var) << "(+" << lag << ")";
                else if (lag < 0)
                  tmp_out << get_variable(SymbolType::endogenous, var) << "(" << lag << ")";
                else
                  tmp_out << get_variable(SymbolType::endogenous, var);
                Stack.push(tmp_out.str());
                if (compute)
                  {
                    if (evaluate)
                      Stackf.push(ya[(it_+lag)*y_size+var]);
                    else
                      Stackf.push(y[(it_+lag)*y_size+var]);
                  }
                break;
              case SymbolType::exogenous:
                var = static_cast<FLDV_ *>(it_code->second)->get_pos();
                lag = static_cast<FLDV_ *>(it_code->second)->get_lead_lag();
                tmp_out.str("");
                if (lag != 0)
                  tmp_out << get_variable(SymbolType::exogenous, var) << "(" << lag << ")";
                else
                  tmp_out << get_variable(SymbolType::exogenous, var);
                Stack.push(tmp_out.str());
                if (compute)
                  Stackf.push(x[it_+lag+var*nb_row_x]);
                break;
              case SymbolType::exogenousDet:
                var = static_cast<FLDV_ *>(it_code->second)->get_pos();
                lag = static_cast<FLDV_ *>(it_code->second)->get_lead_lag();
                tmp_out.str("");
                if (lag != 0)
                  tmp_out << get_variable(SymbolType::exogenousDet, var) << "(" << lag << ")";
                else
                  tmp_out << get_variable(SymbolType::exogenousDet, var);
                Stack.push(tmp_out.str());
                if (compute)
                  Stackf.push(x[it_+lag+var*nb_row_xd]);
                break;
              case SymbolType::modelLocalVariable:
                break;
              default:
                mexPrintf("FLDV: Unknown variable type\n");
              }
            break;
          case Tags::FLDSV:
          case Tags::FLDVS:
            //load a variable in the processor
            switch (static_cast<SymbolType>(static_cast<FLDSV_ *>(it_code->second)->get_type()))
              {
              case SymbolType::parameter:
                var = static_cast<FLDSV_ *>(it_code->second)->get_pos();
                Stack.push(get_variable(SymbolType::parameter, var));
                if (compute)
                  Stackf.push(params[var]);
                break;
              case SymbolType::endogenous:
                var = static_cast<FLDSV_ *>(it_code->second)->get_pos();
                Stack.push(get_variable(SymbolType::endogenous, var));
                if (compute)
                  {
                    if (it_code->first == Tags::FLDSV)
                      {
                        if (evaluate)
                          Stackf.push(ya[var]);
                        else
                          Stackf.push(y[var]);
                      }
                    else
                      Stackf.push(steady_y[var]);
                  }
                break;
              case SymbolType::exogenous:
                var = static_cast<FLDSV_ *>(it_code->second)->get_pos();
                Stack.push(get_variable(SymbolType::exogenous, var));
                if (compute)
                  Stackf.push(x[var]);
                break;
              case SymbolType::exogenousDet:
                var = static_cast<FLDSV_ *>(it_code->second)->get_pos();
                Stack.push(get_variable(SymbolType::exogenousDet, var));
                if (compute)
                  Stackf.push(x[var]);
                break;
              case SymbolType::modelLocalVariable:
                break;
              default:
                mexPrintf("FLDSV: Unknown variable type\n");
              }
            break;
          case Tags::FLDT:
            //load a temporary variable in the processor
            var = static_cast<FLDT_ *>(it_code->second)->get_pos();
            tmp_out.str("");
            tmp_out << "T" << var+1;
            Stack.push(tmp_out.str());
            if (compute)
              Stackf.push(T[var*(periods+y_kmin+y_kmax)+it_]);
            break;
          case Tags::FLDST:
            //load a temporary variable in the processor
            var = static_cast<FLDST_ *>(it_code->second)->get_pos();
            tmp_out.str("");
            tmp_out << "T" << var+1;
            Stack.push(tmp_out.str());
            if (compute)
              Stackf.push(T[var]);
            break;
          case Tags::FLDU:
            //load u variable in the processor
            var = static_cast<FLDU_ *>(it_code->second)->get_pos();
            tmp_out.str("");
            tmp_out << "u(" << var+1 << " + it_)";
            Stack.push(tmp_out.str());
            var += Per_u_;
            if (compute)
              Stackf.push(u[var]);
            break;
          case Tags::FLDSU:
            //load u variable in the processor
            var = static_cast<FLDSU_ *>(it_code->second)->get_pos();
            tmp_out.str("");
            tmp_out << "u(" << var+1 << ")";
            Stack.push(tmp_out.str());
            if (compute)
              Stackf.push(u[var]);
            break;
          case Tags::FLDR:
            var = static_cast<FLDR_ *>(it_code->second)->get_pos();
            tmp_out.str("");
            tmp_out << "residual(" << var+1 << ")";
            Stack.push(tmp_out.str());
            if (compute)
              Stackf.push(r[var]);
            break;
          case Tags::FLDZ:
            //load 0 in the processor
            tmp_out.str("");
            tmp_out << 0;
            Stack.push(tmp_out.str());
            if (compute)
              Stackf.push(0.0);
            break;
          case Tags::FLDC:
            //load a numerical constant in the processor
            ll = static_cast<FLDC_ *>(it_code->second)->get_value();
            tmp_out.str("");
            tmp_out << ll;
            Stack.push(tmp_out.str());
            if (compute)
              Stackf.push(ll);
            break;
          case Tags::FSTPV:
            //load a variable in the processor
            go_on = false;
            switch (static_cast<SymbolType>(static_cast<FSTPV_ *>(it_code->second)->get_type()))
              {
              case SymbolType::parameter:
                var = static_cast<FSTPV_ *>(it_code->second)->get_pos();
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(SymbolType::parameter, var) << " = " << tmp_out2.str();
                Stack.pop();
                if (compute)
                  {
                    params[var] = Stackf.top();
                    Stackf.pop();
                  }
                break;
              case SymbolType::endogenous:
                var = static_cast<FSTPV_ *>(it_code->second)->get_pos();
                lag = static_cast<FSTPV_ *>(it_code->second)->get_lead_lag();
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(SymbolType::endogenous, var);
                if (lag > 0)
                  tmp_out << "(+" << lag << ")";
                else if (lag < 0)
                  tmp_out << "(" << lag << ")";
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                if (compute)
                  {
                    y[(it_+lag)*y_size+var] = Stackf.top();
                    Stackf.pop();
                  }
                break;
              case SymbolType::exogenous:
                var = static_cast<FSTPV_ *>(it_code->second)->get_pos();
                lag = static_cast<FSTPV_ *>(it_code->second)->get_lead_lag();
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(SymbolType::exogenous, var);
                if (lag != 0)
                  tmp_out << "(" << lag << ")";
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                if (compute)
                  {
                    x[it_+lag+var*nb_row_x] = Stackf.top();
                    Stackf.pop();
                  }
                break;
              case SymbolType::exogenousDet:
                var = static_cast<FSTPV_ *>(it_code->second)->get_pos();
                lag = static_cast<FSTPV_ *>(it_code->second)->get_lead_lag();
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(SymbolType::exogenousDet, var);
                if (lag != 0)
                  tmp_out << "(" << lag << ")";
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                if (compute)
                  {
                    x[it_+lag+var*nb_row_xd] = Stackf.top();
                    Stackf.pop();
                  }
                break;
              default:
                mexPrintf("FSTPV: Unknown variable type\n");
              }
            break;
          case Tags::FSTPSV:
            go_on = false;
            //load a variable in the processor
            switch (static_cast<SymbolType>(static_cast<FSTPSV_ *>(it_code->second)->get_type()))
              {
              case SymbolType::parameter:
                var = static_cast<FSTPSV_ *>(it_code->second)->get_pos();
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(SymbolType::parameter, var);
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                if (compute)
                  {
                    params[var] = Stackf.top();
                    Stackf.pop();
                  }
                break;
              case SymbolType::endogenous:
                var = static_cast<FSTPSV_ *>(it_code->second)->get_pos();
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(SymbolType::endogenous, var);
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                if (compute)
                  {
                    y[var] = Stackf.top();
                    Stackf.pop();
                  }
                break;
              case SymbolType::exogenous:
              case SymbolType::exogenousDet:
                var = static_cast<FSTPSV_ *>(it_code->second)->get_pos();
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(SymbolType::exogenous, var);
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                if (compute)
                  {
                    x[var] = Stackf.top();
                    Stackf.pop();
                  }
                break;
              default:
                mexPrintf("FSTPSV: Unknown variable type\n");
              }
            break;
          case Tags::FSTPT:
            go_on = false;
            //store in a temporary variable from the processor
            var = static_cast<FSTPT_ *>(it_code->second)->get_pos();
            tmp_out.str("");
            tmp_out << "T" << var+1 << " = " << Stack.top();
            Stack.pop();
            if (compute)
              {
                T[var*(periods+y_kmin+y_kmax)+it_] = Stackf.top();
                Stackf.pop();
              }
            break;
          case Tags::FSTPST:
            go_on = false;
            //store in a temporary variable from the processor
            var = static_cast<FSTPST_ *>(it_code->second)->get_pos();
            tmp_out.str("");
            tmp_out << "T" << var+1 << " = " << Stack.top();
            Stack.pop();
            if (compute)
              {
                T[var] = Stackf.top();
                Stackf.pop();
              }
            break;
          case Tags::FSTPU:
            go_on = false;
            //store in u variable from the processor
            var = static_cast<FSTPU_ *>(it_code->second)->get_pos();
            tmp_out.str("");
            tmp_out << "u(" << var+1 << " + it_) = " << Stack.top();
            var += Per_u_;
            Stack.pop();
            if (compute)
              {
                u[var] = Stackf.top();
                Stackf.pop();
              }
            break;
          case Tags::FSTPSU:
            go_on = false;
            //store in u variable from the processor
            var = static_cast<FSTPSU_ *>(it_code->second)->get_pos();
            tmp_out.str("");
            tmp_out << "u(" << var+1 << ") = " << Stack.top();
            Stack.pop();
            if (compute)
              {
                u[var] = Stackf.top();
                Stackf.pop();
              }
            break;
          case Tags::FSTPR:
            go_on = false;
            //store in residual variable from the processor
            var = static_cast<FSTPR_ *>(it_code->second)->get_pos();
            tmp_out.str("");
            tmp_out << "residual(" << var+1 << ") = " << Stack.top();
            Stack.pop();
            if (compute)
              {
                r[var] = Stackf.top();
                Stackf.pop();
              }
            break;
          case Tags::FSTPG:
            go_on = false;
            //store in derivative (g) variable from the processor
            var = static_cast<FSTPG_ *>(it_code->second)->get_pos();
            tmp_out.str("");
            tmp_out << "g1[" << var+1 << "] = " << Stack.top();
            Stack.pop();
            if (compute)
              {
                g1[var] = Stackf.top();
                Stackf.pop();
              }
            break;
          case Tags::FSTPG2:
            go_on = false;
            //store in derivative (g) variable from the processor
            eq = static_cast<FSTPG2_ *>(it_code->second)->get_row();
            var = static_cast<FSTPG2_ *>(it_code->second)->get_col();
            tmp_out.str("");
            tmp_out << "/*jacob(" << eq << ", " << var << ")*/ jacob(" << eq+size*var+1 << ") = " << Stack.top();
            Stack.pop();
            if (compute)
              {
                jacob[eq + size*var] = Stackf.top();
                Stackf.pop();
              }
            break;
          case Tags::FSTPG3:
            //store in derivative (g) variable from the processor
            double r;
            unsigned int pos_col;
            go_on = false;
            if (compute)
              {
                r = Stackf.top();
                Stackf.pop();
              }
            eq = static_cast<FSTPG3_ *>(it_code->second)->get_row();
            var = static_cast<FSTPG3_ *>(it_code->second)->get_col();
            lag = static_cast<FSTPG3_ *>(it_code->second)->get_lag();
            pos_col = static_cast<FSTPG3_ *>(it_code->second)->get_col_pos();
            switch (equation_type)
              {
              case ExpressionType::FirstEndoDerivative:
                if (compute)
                  jacob[eq + size*pos_col] = r;
                tmp_out.str("");
                tmp_out << "/*jacob(" << eq << ", " << pos_col << " var= " << var << ")*/ jacob(" << eq+size*pos_col+1 << ") = " << Stack.top();
                Stack.pop();
                break;
              case ExpressionType::FirstOtherEndoDerivative:
                if (compute)
                  jacob_other_endo[eq + size*pos_col] = r;
                tmp_out.str("");
                tmp_out << "jacob_other_endo(" << eq+size*pos_col+1 << ") = " << Stack.top();
                Stack.pop();
                break;
              case ExpressionType::FirstExoDerivative:
                if (compute)
                  jacob_exo[eq + size*pos_col] = r;
                tmp_out.str("");
                tmp_out << "/*jacob_exo(" << eq << ", " << pos_col << " var=" << var << ")*/ jacob_exo(" << eq+size*pos_col+1 << ") = " << Stack.top();
                Stack.pop();
                break;
              case ExpressionType::FirstExodetDerivative:
                if (compute)
                  jacob_exo_det[eq + size*pos_col] = r;
                tmp_out.str("");
                tmp_out << "/*jacob_exo_det(" << eq << ", " << pos_col << " var=" << var << ")*/ jacob_exo_det(" << eq+size*pos_col+1 << ") = " << Stack.top();
                Stack.pop();
                break;
              default:
                ostringstream tmp;
                tmp << " in compute_block_time, variable " << static_cast<int>(EQN_type) << " not used yet\n";
                //throw FatalExceptionHandling(tmp.str());
                mexPrintf("%s", tmp.str().c_str());
              }
            break;
          case Tags::FBINARY:
            op = static_cast<FBINARY_ *>(it_code->second)->get_op_type();
            v2 = Stack.top();
            Stack.pop();
            v1 = Stack.top();
            Stack.pop();
            if (compute)
              {
                v2f = Stackf.top();
                Stackf.pop();
                v1f = Stackf.top();
                Stackf.pop();
              }
            switch (static_cast<BinaryOpcode>(op))
              {
              case BinaryOpcode::plus:
                if (compute)
                  Stackf.push(v1f + v2f);
                tmp_out.str("");
                tmp_out << v1 << " + " << v2;
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::minus:
                if (compute)
                  Stackf.push(v1f - v2f);
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " - ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::times:
                if (compute)
                  Stackf.push(v1f * v2f);
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " * ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::divide:
                if (compute)
                  {
                    r = v1f / v2f;
                    Stackf.push(r);
                  }
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                if (compute)
                  {
                    if (isinf(r))
                      tmp_out << "$";
                    tmp_out << " / ";
                    if (isinf(r))
                      tmp_out << "£";
                  }
                else
                  tmp_out << " / ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::less:
                if (compute)
                  Stackf.push(static_cast<double>(v1f < v2f));
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " < ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::greater:
                if (compute)
                  Stackf.push(static_cast<double>(v1f > v2f));
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " > ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::lessEqual:
                if (compute)
                  Stackf.push(static_cast<double>(v1f <= v2f));
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " <= ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::greaterEqual:
                if (compute)
                  Stackf.push(static_cast<double>(v1f >= v2f));
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " >= ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::equalEqual:
                if (compute)
                  Stackf.push(static_cast<double>(v1f == v2f));
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " == ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::different:
                if (compute)
                  Stackf.push(static_cast<double>(v1f != v2f));
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " != ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::power:
                if (compute)
                  {
                    r = pow(v1f, v2f);
                    Stackf.push(r);
                  }
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                if (compute)
                  {
                    if (isnan(r))
                      tmp_out << "$ ^ " << "£";
                    else
                      tmp_out << " ^ ";
                  }
                else
                  tmp_out << " ^ ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::powerDeriv:
                {
                  v3 = Stack.top();
                  Stack.pop();
                  if (compute)
                    {
                      int derivOrder = static_cast<int>(nearbyint(Stackf.top()));
                      Stackf.pop();
                      if (fabs(v1f) < near_zero && v2f > 0
                          && derivOrder > v2f
                          && fabs(v2f-nearbyint(v2f)) < near_zero)
                        {
                          r = 0.0;
                          Stackf.push(r);
                        }
                      else
                        {
                          double dxp = pow(v1f, v2f-derivOrder);
                          for (int i = 0; i < derivOrder; i++)
                            dxp *= v2f--;
                          Stackf.push(dxp);
                          r = dxp;
                        }
                    }
                  tmp_out.str("");
                  if (compute)
                    {
                      if (isnan(r))
                        tmp_out << "$ PowerDeriv " << "£";
                      else
                        tmp_out << "PowerDeriv";
                    }
                  else
                    tmp_out << "PowerDeriv";
                  tmp_out << "(" << v1 << ", " << v2 << ", " << v3 << ")";
                  Stack.push(tmp_out.str());
                }
                break;
              case BinaryOpcode::max:
                if (compute)
                  Stackf.push(max(v1f, v2f));
                tmp_out.str("");
                tmp_out << "max(" << v1 << ", " << v2 << ")";
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::min:
                if (compute)
                  Stackf.push(min(v1f, v2f));
                tmp_out.str("");
                tmp_out << "min(" << v1 << ", " << v2 << ")";
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::equal:
              default:
                mexPrintf("Error unknown Unary operator=%d\n", op); mexEvalString("drawnow;");
                ;
              }
            break;
          case Tags::FUNARY:
            op = static_cast<FUNARY_ *>(it_code->second)->get_op_type();
            v1 = Stack.top();
            Stack.pop();
            if (compute)
              {
                v1f = Stackf.top();
                Stackf.pop();
              }
            switch (static_cast<UnaryOpcode>(op))
              {
              case UnaryOpcode::uminus:
                if (compute)
                  Stackf.push(-v1f);
                tmp_out.str("");
                tmp_out << " -" << v1;
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::exp:
                if (compute)
                  Stackf.push(exp(v1f));
                tmp_out.str("");
                tmp_out << "exp(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::log:
                if (compute)
                  {
                    r = log(v1f);
                    Stackf.push(r);
                  }
                tmp_out.str("");
                if (compute)
                  {
                    if (isnan(r))
                      tmp_out << "$log" << "£" << "(" << v1 << ")";
                    else
                      tmp_out << "log(" << v1 << ")";
                  }
                else
                  tmp_out << "log(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::log10:
                if (compute)
                  {
                    r = log10(v1f);
                    Stackf.push(r);
                  }
                tmp_out.str("");
                if (compute)
                  {
                    if (isnan(r))
                      tmp_out << "$log10" << "£" << "(" << v1 << ")";
                    else
                      tmp_out << "log10(" << v1 << ")";
                  }
                else
                  tmp_out << "log10(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::cos:
                if (compute)
                  Stackf.push(cos(v1f));
                tmp_out.str("");
                tmp_out << "cos(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::sin:
                if (compute)
                  Stackf.push(sin(v1f));
                tmp_out.str("");
                tmp_out << "sin(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::tan:
                if (compute)
                  Stackf.push(tan(v1f));
                tmp_out.str("");
                tmp_out << "tan(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::acos:
                if (compute)
                  Stackf.push(acos(v1f));
                tmp_out.str("");
                tmp_out << "acos(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::asin:
                if (compute)
                  Stackf.push(asin(v1f));
                tmp_out.str("");
                tmp_out << "asin(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::atan:
                if (compute)
                  Stackf.push(atan(v1f));
                tmp_out.str("");
                tmp_out << "atan(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::cosh:
                if (compute)
                  Stackf.push(cosh(v1f));
                tmp_out.str("");
                tmp_out << "cosh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::sinh:
                if (compute)
                  Stackf.push(sinh(v1f));
                tmp_out.str("");
                tmp_out << "sinh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::tanh:
                if (compute)
                  Stackf.push(tanh(v1f));
                tmp_out.str("");
                tmp_out << "tanh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::acosh:
                if (compute)
                  Stackf.push(acosh(v1f));
                tmp_out.str("");
                tmp_out << "acosh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::asinh:
                if (compute)
                  Stackf.push(asinh(v1f));
                tmp_out.str("");
                tmp_out << "asinh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::atanh:
                if (compute)
                  Stackf.push(atanh(v1f));
                tmp_out.str("");
                tmp_out << "atanh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::sqrt:
                if (compute)
                  Stackf.push(sqrt(v1f));
                tmp_out.str("");
                tmp_out << "sqrt(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::erf:
                if (compute)
                  Stackf.push(erf(v1f));
                tmp_out.str("");
                tmp_out << "erf(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::erfc:
                if (compute)
                  Stackf.push(erfc(v1f));
                tmp_out.str("");
                tmp_out << "erfc(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              default:
                mexPrintf("Error unknown Binary operator=%d\n", op); mexEvalString("drawnow;");
                ;
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
            if (compute)
              {
                v3f = Stackf.top();
                Stackf.pop();
                v2f = Stackf.top();
                Stackf.pop();
                v1f = Stackf.top();
                Stackf.pop();
              }
            switch (static_cast<TrinaryOpcode>(op))
              {
              case TrinaryOpcode::normcdf:
                if (compute)
                  Stackf.push(0.5*(1+erf((v1f-v2f)/v3f/M_SQRT2)));
                tmp_out.str("");
                tmp_out << "normcdf(" << v1 << ", " << v2 << ", " << v3 << ")";
                Stack.push(tmp_out.str());
                break;
              case TrinaryOpcode::normpdf:
                if (compute)
                  Stackf.push(1/(v3f*sqrt(2*M_PI)*exp(pow((v1f-v2f)/v3f, 2)/2)));
                tmp_out.str("");
                tmp_out << "normpdf(" << v1 << ", " << v2 << ", " << v3 << ")";
                Stack.push(tmp_out.str());
                break;
              default:
                mexPrintf("Error unknown Trinary operator=%d\n", op); mexEvalString("drawnow;");
              }
            break;
          case Tags::FCALL:
            {
              auto *fc = static_cast<FCALL_ *>(it_code->second);
              string function_name = fc->get_function_name();
              unsigned int nb_input_arguments = fc->get_nb_input_arguments();
              unsigned int nb_output_arguments = fc->get_nb_output_arguments();

              mxArray *output_arguments[3];
              string arg_func_name = fc->get_arg_func_name();
              unsigned int nb_add_input_arguments = fc->get_nb_add_input_arguments();
              function_type = fc->get_function_type();
              mxArray **input_arguments;
              switch (function_type)
                {
                case ExternalFunctionType::withoutDerivative:
                case ExternalFunctionType::withFirstDerivative:
                case ExternalFunctionType::withFirstAndSecondDerivative:
                  {
                    if (compute)
                      {
                        input_arguments = static_cast<mxArray **>(mxMalloc(nb_input_arguments * sizeof(mxArray *)));
                        for (unsigned int i = 0; i < nb_input_arguments; i++)
                          {
                            mxArray *vv = mxCreateDoubleScalar(Stackf.top());
                            input_arguments[nb_input_arguments - i - 1] = vv;
                            Stackf.pop();
                          }
                        mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str());
                        double *rr = mxGetPr(output_arguments[0]);
                        Stackf.push(*rr);
                      }
                    tmp_out.str("");
                    tmp_out << function_name << "(";
                    vector<string> ss(nb_input_arguments);
                    for (unsigned int i = 0; i < nb_input_arguments; i++)
                      {
                        ss[nb_input_arguments-i-1] = Stack.top();
                        Stack.pop();
                      }
                    for (unsigned int i = 0; i < nb_input_arguments; i++)
                      {
                        tmp_out << ss[i];
                        if (i < nb_input_arguments - 1)
                          tmp_out << ", ";
                      }
                    tmp_out << ")";
                    Stack.push(tmp_out.str());
                  }
                  break;
                case ExternalFunctionType::numericalFirstDerivative:
                  {
                    if (compute)
                      {
                        input_arguments = static_cast<mxArray **>(mxMalloc((nb_input_arguments+1+nb_add_input_arguments) * sizeof(mxArray *)));
                        mxArray *vv = mxCreateString(arg_func_name.c_str());
                        input_arguments[0] = vv;
                        vv = mxCreateDoubleScalar(fc->get_row());
                        input_arguments[1] = vv;
                        vv = mxCreateCellMatrix(1, nb_add_input_arguments);
                        for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                          {
                            double rr = Stackf.top();
                            mxSetCell(vv, nb_add_input_arguments - (i+1), mxCreateDoubleScalar(rr));
                            Stackf.pop();
                          }
                        input_arguments[nb_input_arguments+nb_add_input_arguments] = vv;
                        nb_input_arguments = 3;
                        mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str());
                        double *rr = mxGetPr(output_arguments[0]);
                        Stackf.push(*rr);
                      }
                    tmp_out.str("");
                    tmp_out << function_name << "(";
                    tmp_out << arg_func_name.c_str() << ", " << fc->get_row() << ", {";
                    vector<string> ss(nb_input_arguments);
                    for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                      {
                        ss[nb_add_input_arguments-i-1] = Stack.top();
                        Stack.pop();
                      }
                    for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                      {
                        tmp_out << ss[i];
                        if (i < nb_add_input_arguments - 1)
                          tmp_out << ", ";
                      }
                    tmp_out << "})";
                    Stack.push(tmp_out.str());
                  }
                  break;
                case ExternalFunctionType::firstDerivative:
                  {
                    if (compute)
                      {
                        input_arguments = static_cast<mxArray **>(mxMalloc(nb_input_arguments * sizeof(mxArray *)));
                        for (unsigned int i = 0; i < nb_input_arguments; i++)
                          {
                            mxArray *vv = mxCreateDoubleScalar(Stackf.top());
                            input_arguments[(nb_input_arguments - 1) - i] = vv;
                            Stackf.pop();
                          }
                        mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str());
                      }
                    tmp_out.str("");
                    tmp_out << function_name << "(";
                    vector<string> ss(nb_input_arguments);
                    for (unsigned int i = 0; i < nb_input_arguments; i++)
                      {
                        ss[nb_input_arguments-i-1] = Stack.top();
                        Stack.pop();
                      }
                    for (unsigned int i = 0; i < nb_input_arguments; i++)
                      {
                        tmp_out << ss[i];
                        if (i < nb_input_arguments - 1)
                          tmp_out << ", ";
                      }
                    tmp_out << ")";
                    Stack.push(tmp_out.str());
                  }
                  break;
                case ExternalFunctionType::numericalSecondDerivative:
                  {
                    if (compute)
                      {
                        input_arguments = static_cast<mxArray **>(mxMalloc((nb_input_arguments+1+nb_add_input_arguments) * sizeof(mxArray *)));
                        mxArray *vv = mxCreateString(arg_func_name.c_str());
                        input_arguments[0] = vv;
                        vv = mxCreateDoubleScalar(fc->get_row());
                        input_arguments[1] = vv;
                        vv = mxCreateDoubleScalar(fc->get_col());
                        input_arguments[2] = vv;
                        vv = mxCreateCellMatrix(1, nb_add_input_arguments);
                        for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                          {
                            double rr = Stackf.top();
                            mxSetCell(vv, (nb_add_input_arguments - 1) - i, mxCreateDoubleScalar(rr));
                            Stackf.pop();
                          }
                        input_arguments[nb_input_arguments+nb_add_input_arguments] = vv;
                        nb_input_arguments = 3;
                        mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str());
                        double *rr = mxGetPr(output_arguments[0]);
                        Stackf.push(*rr);
                      }
                    tmp_out.str("");
                    tmp_out << function_name << "(";
                    tmp_out << arg_func_name.c_str() << ", " << fc->get_row() << ", " << fc->get_col() << ", {";
                    vector<string> ss(nb_input_arguments);
                    for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                      {
                        ss[nb_add_input_arguments-i-1] = Stack.top();
                        Stack.pop();
                      }
                    for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                      {
                        tmp_out << ss[i];
                        if (i < nb_add_input_arguments - 1)
                          tmp_out << ", ";
                      }
                    tmp_out << "})";
                    Stack.push(tmp_out.str());
                  }
                  break;
                case ExternalFunctionType::secondDerivative:
                  {
                    if (compute)
                      {
                        input_arguments = static_cast<mxArray **>(mxMalloc(nb_input_arguments * sizeof(mxArray *)));
                        for (unsigned int i = 0; i < nb_input_arguments; i++)
                          {
                            mxArray *vv = mxCreateDoubleScalar(Stackf.top());
                            input_arguments[i] = vv;
                            Stackf.pop();
                          }
                        mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str());
                      }
                    tmp_out.str("");
                    tmp_out << function_name << "(";
                    vector<string> ss(nb_input_arguments);
                    for (unsigned int i = 0; i < nb_input_arguments; i++)
                      {
                        ss[nb_input_arguments-i-1] = Stack.top();
                        Stack.pop();
                      }
                    for (unsigned int i = 0; i < nb_input_arguments; i++)
                      {
                        tmp_out << ss[i];
                        if (i < nb_input_arguments - 1)
                          tmp_out << ", ";
                      }
                    tmp_out << ")";
                    Stack.push(tmp_out.str());
                  }
                  break;
                }
              break;
            }
          case Tags::FSTPTEF:
            go_on = false;
            var = static_cast<FSTPTEF_ *>(it_code->second)->get_number();
            if (compute)
              {
                Stackf.pop();
              }
            tmp_out.str("");
            switch (function_type)
              {
              case ExternalFunctionType::withoutDerivative:
                tmp_out << "TEF(" << var << ") = " << Stack.top();
                break;
              case ExternalFunctionType::withFirstDerivative:
                tmp_out << "[TEF(" << var << "), TEFD(" << var << ") ]= " << Stack.top();
                break;
              case ExternalFunctionType::withFirstAndSecondDerivative:
                tmp_out << "[TEF(" << var << "), TEFD(" << var << "), TEFDD(" << var << ") ]= " << Stack.top();
                break;
              default:
                break;
              }
            Stack.pop();
            break;
          case Tags::FLDTEF:
            var = static_cast<FLDTEF_ *>(it_code->second)->get_number();
            if (compute)
              {
                auto it = TEF.find(var-1);
                Stackf.push(it->second);
              }
            tmp_out.str("");
            tmp_out << "TEF(" << var << ")";
            Stack.push(tmp_out.str());

            break;
          case Tags::FSTPTEFD:
            {
              go_on = false;
              unsigned int indx = static_cast<FSTPTEFD_ *>(it_code->second)->get_indx();
              unsigned int row = static_cast<FSTPTEFD_ *>(it_code->second)->get_row();
              if (compute)
                {
                  Stackf.pop();
                }
              tmp_out.str("");
              if (function_type == ExternalFunctionType::numericalFirstDerivative)
                tmp_out << "TEFD(" << indx << ", " << row << ") = " << Stack.top();
              else if (function_type == ExternalFunctionType::firstDerivative)
                tmp_out << "TEFD(" << indx << ") = " << Stack.top();
              Stack.pop();
            }
            break;
          case Tags::FLDTEFD:
            {
              unsigned int indx = static_cast<FLDTEFD_ *>(it_code->second)->get_indx();
              unsigned int row = static_cast<FLDTEFD_ *>(it_code->second)->get_row();
              if (compute)
                {
                  auto it = TEFD.find({ indx, row-1 });
                  Stackf.push(it->second);
                }
              tmp_out.str("");
              tmp_out << "TEFD(" << indx << ", " << row << ")";
              Stack.push(tmp_out.str());
            }
            break;
          case Tags::FSTPTEFDD:
            {
              go_on = false;
              unsigned int indx = static_cast<FSTPTEFDD_ *>(it_code->second)->get_indx();
              unsigned int row = static_cast<FSTPTEFDD_ *>(it_code->second)->get_row();
              unsigned int col = static_cast<FSTPTEFDD_ *>(it_code->second)->get_col();
              if (compute)
                {
                  Stackf.pop();
                }
              tmp_out.str("");
              if (function_type == ExternalFunctionType::numericalSecondDerivative)
                tmp_out << "TEFDD(" << indx << ", " << row << ", " << col << ") = " << Stack.top();
              else if (function_type == ExternalFunctionType::secondDerivative)
                tmp_out << "TEFDD(" << indx << ") = " << Stack.top();
              Stack.pop();
            }

            break;
          case Tags::FLDTEFDD:
            {
              unsigned int indx = static_cast<FLDTEFDD_ *>(it_code->second)->get_indx();
              unsigned int row = static_cast<FLDTEFDD_ *>(it_code->second)->get_row();
              unsigned int col = static_cast<FSTPTEFDD_ *>(it_code->second)->get_col();
              if (compute)
                {
                  auto it = TEFDD.find({ indx, row-1, col-1 });
                  Stackf.push(it->second);
                }
              tmp_out.str("");
              tmp_out << "TEFDD(" << indx << ", " << row << ", " << col << ")";
              Stack.push(tmp_out.str());
            }
            break;
          case Tags::FJMPIFEVAL:
            tmp_out.str("");
            tmp_out << "if (~evaluate)";
            go_on = false;
            break;
          case Tags::FJMP:
            tmp_out.str("");
            tmp_out << "else";
            go_on = false;
            break;
          case Tags::FCUML:
            if (compute)
              {
                v1f = Stackf.top();
                Stackf.pop();
                v2f = Stackf.top();
                Stackf.pop();
                Stackf.push(v1f+v2f);
              }
            v1 = Stack.top();
            Stack.pop();
            v2 = Stack.top();
            Stack.pop();
            tmp_out.str("");
            tmp_out << v1 << " + " << v2;
            Stack.push(tmp_out.str());
            break;
          case Tags::FENDBLOCK:
          case Tags::FENDEQU:
            go_on = false;
            break;
          case Tags::FOK:
            break;
          default:
            mexPrintf("Error it_code->first=%d unknown\n", it_code->first); mexEvalString("drawnow;");
            throw FatalExceptionHandling(" in print_expression, unknown opcode "
                                         + to_string(static_cast<int>(it_code->first))
                                         + "!! FENDEQU="
                                         + to_string(static_cast<int>(Tags::FENDEQU)) + "\n");
          }
        it_code++;
      }
    it_code_ret = it_code;
    return tmp_out.str();
  }

  static inline void
  test_mxMalloc(void *z, int line, const string &file, const string &func, int amount)
  {
    if (!z && amount > 0)
      throw FatalExceptionHandling(" mxMalloc: out of memory " + to_string(amount) + " bytes required at line " + to_string(line) + " in function " + func + " (file " + file);
  }
};

#endif
