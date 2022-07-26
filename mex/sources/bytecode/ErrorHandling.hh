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
#include <optional>
#define _USE_MATH_DEFINES
#include <cmath>

#include "dynmex.h"

#define BYTECODE_MEX
#include "Bytecode.hh"

using namespace std;

constexpr int NO_ERROR_ON_EXIT = 0, ERROR_ON_EXIT = 1;

using it_code_type = instructions_list_t::const_iterator;

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
protected:
  ExpressionType EQN_type;
  int EQN_equation, EQN_block, EQN_block_number, EQN_dvar1;
  size_t endo_name_length; // Maximum length of endogenous names
  vector<string> P_endo_names;
private:
  bool is_load_variable_list;
  vector<string> P_exo_names, P_param_names;
  vector<tuple<string, SymbolType, unsigned int>> Variable_list;

public:
  ErrorMsg() : is_load_variable_list {false}
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
            mexErrMsgTxt(("M_."s + symbol_type + " is not a cell array").c_str());
          for (size_t i = 0; i < mxGetNumberOfElements(M_field); i++)
            {
              const mxArray *cell_mx = mxGetCell(M_field, i);
              if (!(cell_mx && mxIsChar(cell_mx)))
                mexErrMsgTxt(("M_."s + symbol_type + " contains a cell which is not a character array").c_str());
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
  }

private:
  /* Given a string which possibly contains a floating-point exception
    (materialized by an operator between braces), returns a string spanning two
    lines, the first line containing the original string without the braces,
    the second line containing tildes (~) under the faulty operator. */
  string
  add_underscore_to_fpe(const string &str)
  {
    string temp;
    int pos1 = -1, pos2 = -1;
    string tmp_n(str.length(), ' ');
    for (const char & i : str)
      {
        if (i != '{' && i != '}')
          temp += i;
        else
          {
            if (i == '{')
              pos1 = static_cast<int>(temp.length());
            else
              pos2 = static_cast<int>(temp.length());
            if (pos1 >= 0 && pos2 >= 0)
              {
                tmp_n.erase(pos1, pos2-pos1+1);
                tmp_n.insert(pos1, pos2-pos1, '~');
                pos1 = pos2 = -1;
              }
          }
      }
    temp += "\n" + tmp_n;
    return temp;
  }

  void
  load_variable_list()
  {
    if (exchange(is_load_variable_list, true))
      return;
    for (size_t variable_num {0}; variable_num < P_endo_names.size(); variable_num++)
      Variable_list.emplace_back(P_endo_names[variable_num], SymbolType::endogenous, variable_num);
    for (size_t variable_num {0}; variable_num < P_exo_names.size(); variable_num++)
      Variable_list.emplace_back(P_exo_names[variable_num], SymbolType::exogenous, variable_num);
  }

public:
  int
  get_ID(const string &variable_name, SymbolType *variable_type)
  {
    load_variable_list();
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

protected:
  string
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
        if (variable_num < P_exo_names.size())
          return P_exo_names[variable_num];
        else
          mexPrintf("=> Unknown exogenous variable # %d", variable_num);
        break;
      case SymbolType::exogenousDet:
        mexErrMsgTxt("get_variable: exogenous deterministic not supported");
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

  string
  error_location(it_code_type expr_begin, it_code_type faulty_op, bool steady_state, int it_)
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
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to exogenous variable "  << get_variable(SymbolType::exogenous, EQN_dvar1) << " at time " << it_;
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to exogenous variable " << get_variable(SymbolType::exogenous, EQN_dvar1) << " at time " << it_;
          break;
        case ExpressionType::FirstExodetDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to deterministic exogenous variable "  << get_variable(SymbolType::exogenousDet, EQN_dvar1) << " at time " << it_;
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to deterministic exogenous variable " << get_variable(SymbolType::exogenousDet, EQN_dvar1) << " at time " << it_;
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
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to exogenous variable "  << get_variable(SymbolType::exogenous, EQN_dvar1);
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to exogenous variable " << get_variable(SymbolType::exogenous, EQN_dvar1);
          break;
        case ExpressionType::FirstExodetDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to deterministic exogenous variable "  << get_variable(SymbolType::exogenousDet, EQN_dvar1);
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to deterministic exogenous variable " << get_variable(SymbolType::exogenousDet, EQN_dvar1);
          break;
        default:
          return ("???");
        }
    auto [expr_str, it_code_ret] = print_expression(expr_begin, faulty_op);
    Error_loc << endl << add_underscore_to_fpe("      " + expr_str);
    return Error_loc.str();
  }

  /* Prints a bytecode expression in human readable form.
     If faulty_op is not default constructed, it should point to a tag withing
     the expression that created a floating point exception, in which case the
     corresponding mathematical operator will be printed within braces.
     The second output argument points to the tag past the expression. */
  pair<string, it_code_type>
  print_expression(const it_code_type &expr_begin, const optional<it_code_type> &faulty_op = nullopt) const
  {
    int var, lag{0}, eq;
    UnaryOpcode op1;
    BinaryOpcode op2;
    TrinaryOpcode op3;
    stack<string> Stack;
    ostringstream tmp_out, tmp_out2;
    string v1, v2, v3;
    bool go_on = true;
    double ll;
    ExpressionType equation_type = ExpressionType::TemporaryTerm;
    size_t found;
    ExternalFunctionCallType call_type{ExternalFunctionCallType::levelWithoutDerivative};

    it_code_type it_code {expr_begin};

    while (go_on)
      {
#ifdef MATLAB_MEX_FILE
        if (utIsInterruptPending())
          throw UserExceptionHandling();
#endif
        switch ((*it_code)->op_code)
          {
          case Tags::FNUMEXPR:
            switch (static_cast<FNUMEXPR_ *>(*it_code)->get_expression_type())
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
              default:
                ostringstream tmp;
                tmp << " in print_expression, expression type " << static_cast<int>(static_cast<FNUMEXPR_ *>(*it_code)->get_expression_type()) << " not implemented yet\n";
                throw FatalExceptionHandling(tmp.str());
              }
            break;
          case Tags::FLDV:
            //load a variable in the processor
            switch (static_cast<FLDV_ *>(*it_code)->get_type())
              {
              case SymbolType::parameter:
                var = static_cast<FLDV_ *>(*it_code)->get_pos();
                Stack.push(get_variable(SymbolType::parameter, var));
                break;
              case SymbolType::endogenous:
                var = static_cast<FLDV_ *>(*it_code)->get_pos();
                lag = static_cast<FLDV_ *>(*it_code)->get_lead_lag();
                tmp_out.str("");
                if (lag > 0)
                  tmp_out << get_variable(SymbolType::endogenous, var) << "(+" << lag << ")";
                else if (lag < 0)
                  tmp_out << get_variable(SymbolType::endogenous, var) << "(" << lag << ")";
                else
                  tmp_out << get_variable(SymbolType::endogenous, var);
                Stack.push(tmp_out.str());
                break;
              case SymbolType::exogenous:
                var = static_cast<FLDV_ *>(*it_code)->get_pos();
                lag = static_cast<FLDV_ *>(*it_code)->get_lead_lag();
                tmp_out.str("");
                if (lag != 0)
                  tmp_out << get_variable(SymbolType::exogenous, var) << "(" << lag << ")";
                else
                  tmp_out << get_variable(SymbolType::exogenous, var);
                Stack.push(tmp_out.str());
                break;
              case SymbolType::exogenousDet:
                mexErrMsgTxt("FLDV: exogenous deterministic not supported");
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
            switch (static_cast<FLDSV_ *>(*it_code)->get_type())
              {
              case SymbolType::parameter:
                var = static_cast<FLDSV_ *>(*it_code)->get_pos();
                Stack.push(get_variable(SymbolType::parameter, var));
                break;
              case SymbolType::endogenous:
                var = static_cast<FLDSV_ *>(*it_code)->get_pos();
                Stack.push(get_variable(SymbolType::endogenous, var));
                break;
              case SymbolType::exogenous:
                var = static_cast<FLDSV_ *>(*it_code)->get_pos();
                Stack.push(get_variable(SymbolType::exogenous, var));
                break;
              case SymbolType::exogenousDet:
                mexErrMsgTxt("FLDSV: exogenous deterministic not supported");
                break;
              case SymbolType::modelLocalVariable:
                break;
              default:
                mexPrintf("FLDSV: Unknown variable type\n");
              }
            break;
          case Tags::FLDT:
            //load a temporary variable in the processor
            var = static_cast<FLDT_ *>(*it_code)->get_pos();
            tmp_out.str("");
            tmp_out << "T" << var+1;
            Stack.push(tmp_out.str());
            break;
          case Tags::FLDST:
            //load a temporary variable in the processor
            var = static_cast<FLDST_ *>(*it_code)->get_pos();
            tmp_out.str("");
            tmp_out << "T" << var+1;
            Stack.push(tmp_out.str());
            break;
          case Tags::FLDU:
            //load u variable in the processor
            var = static_cast<FLDU_ *>(*it_code)->get_pos();
            tmp_out.str("");
            tmp_out << "u(" << var+1 << " + it_)";
            Stack.push(tmp_out.str());
            break;
          case Tags::FLDSU:
            //load u variable in the processor
            var = static_cast<FLDSU_ *>(*it_code)->get_pos();
            tmp_out.str("");
            tmp_out << "u(" << var+1 << ")";
            Stack.push(tmp_out.str());
            break;
          case Tags::FLDR:
            var = static_cast<FLDR_ *>(*it_code)->get_pos();
            tmp_out.str("");
            tmp_out << "residual(" << var+1 << ")";
            Stack.push(tmp_out.str());
            break;
          case Tags::FLDZ:
            //load 0 in the processor
            tmp_out.str("");
            tmp_out << 0;
            Stack.push(tmp_out.str());
            break;
          case Tags::FLDC:
            //load a numerical constant in the processor
            ll = static_cast<FLDC_ *>(*it_code)->get_value();
            tmp_out.str("");
            tmp_out << ll;
            Stack.push(tmp_out.str());
            break;
          case Tags::FSTPV:
            //load a variable in the processor
            go_on = false;
            switch (static_cast<FSTPV_ *>(*it_code)->get_type())
              {
              case SymbolType::parameter:
                var = static_cast<FSTPV_ *>(*it_code)->get_pos();
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(SymbolType::parameter, var) << " = " << tmp_out2.str();
                Stack.pop();
                break;
              case SymbolType::endogenous:
                var = static_cast<FSTPV_ *>(*it_code)->get_pos();
                lag = static_cast<FSTPV_ *>(*it_code)->get_lead_lag();
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
                break;
              case SymbolType::exogenous:
                var = static_cast<FSTPV_ *>(*it_code)->get_pos();
                lag = static_cast<FSTPV_ *>(*it_code)->get_lead_lag();
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(SymbolType::exogenous, var);
                if (lag != 0)
                  tmp_out << "(" << lag << ")";
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                break;
              case SymbolType::exogenousDet:
                mexErrMsgTxt("FSTPV: exogenous deterministic not supported");
                break;
              default:
                mexPrintf("FSTPV: Unknown variable type\n");
              }
            break;
          case Tags::FSTPSV:
            go_on = false;
            //load a variable in the processor
            switch (static_cast<FSTPSV_ *>(*it_code)->get_type())
              {
              case SymbolType::parameter:
                var = static_cast<FSTPSV_ *>(*it_code)->get_pos();
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(SymbolType::parameter, var);
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                break;
              case SymbolType::endogenous:
                var = static_cast<FSTPSV_ *>(*it_code)->get_pos();
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(SymbolType::endogenous, var);
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                break;
              case SymbolType::exogenous:
                var = static_cast<FSTPSV_ *>(*it_code)->get_pos();
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(SymbolType::exogenous, var);
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                break;
              case SymbolType::exogenousDet:
                mexErrMsgTxt("FSTPSV: exogenous deterministic not supported");
                break;
              default:
                mexPrintf("FSTPSV: Unknown variable type\n");
              }
            break;
          case Tags::FSTPT:
            go_on = false;
            //store in a temporary variable from the processor
            var = static_cast<FSTPT_ *>(*it_code)->get_pos();
            tmp_out.str("");
            tmp_out << "T" << var+1 << " = " << Stack.top();
            Stack.pop();
            break;
          case Tags::FSTPST:
            go_on = false;
            //store in a temporary variable from the processor
            var = static_cast<FSTPST_ *>(*it_code)->get_pos();
            tmp_out.str("");
            tmp_out << "T" << var+1 << " = " << Stack.top();
            Stack.pop();
            break;
          case Tags::FSTPU:
            go_on = false;
            //store in u variable from the processor
            var = static_cast<FSTPU_ *>(*it_code)->get_pos();
            tmp_out.str("");
            tmp_out << "u(" << var+1 << " + it_) = " << Stack.top();
            Stack.pop();
            break;
          case Tags::FSTPSU:
            go_on = false;
            //store in u variable from the processor
            var = static_cast<FSTPSU_ *>(*it_code)->get_pos();
            tmp_out.str("");
            tmp_out << "u(" << var+1 << ") = " << Stack.top();
            Stack.pop();
            break;
          case Tags::FSTPR:
            go_on = false;
            //store in residual variable from the processor
            var = static_cast<FSTPR_ *>(*it_code)->get_pos();
            tmp_out.str("");
            tmp_out << "residual(" << var+1 << ") = " << Stack.top();
            Stack.pop();
            break;
          case Tags::FSTPG:
            go_on = false;
            //store in derivative (g) variable from the processor
            var = static_cast<FSTPG_ *>(*it_code)->get_pos();
            tmp_out.str("");
            tmp_out << "g1[" << var+1 << "] = " << Stack.top();
            Stack.pop();
            break;
          case Tags::FSTPG2:
            go_on = false;
            //store in derivative (g) variable from the processor
            eq = static_cast<FSTPG2_ *>(*it_code)->get_row();
            var = static_cast<FSTPG2_ *>(*it_code)->get_col();
            tmp_out.str("");
            tmp_out << "jacob(" << eq+1 << ", " << var+1 << ") = " << Stack.top();
            Stack.pop();
            break;
          case Tags::FSTPG3:
            //store in derivative (g) variable from the processor
            go_on = false;
            int pos_col;
            eq = static_cast<FSTPG3_ *>(*it_code)->get_row();
            var = static_cast<FSTPG3_ *>(*it_code)->get_col();
            lag = static_cast<FSTPG3_ *>(*it_code)->get_lag();
            pos_col = static_cast<FSTPG3_ *>(*it_code)->get_col_pos();
            tmp_out.str("");
            switch (equation_type)
              {
              case ExpressionType::FirstEndoDerivative:
                tmp_out << "jacob";
                break;
              case ExpressionType::FirstOtherEndoDerivative:
                tmp_out << "jacob_other_endo";
                break;
              case ExpressionType::FirstExoDerivative:
                tmp_out << "jacob_exo";
                break;
              case ExpressionType::FirstExodetDerivative:
                tmp_out << "jacob_exo_det";
                break;
              default:
                ostringstream tmp;
                tmp << " in compute_block_time, variable " << static_cast<int>(EQN_type) << " not used yet\n";
                //throw FatalExceptionHandling(tmp.str());
                mexPrintf("%s", tmp.str().c_str());
              }
            tmp_out << "(" << eq+1 << ", " << pos_col+1 << " [var=" << var+1 << ", lag=" << lag << "]) = " << Stack.top();
            Stack.pop();
            break;
          case Tags::FBINARY:
            op2 = static_cast<FBINARY_ *>(*it_code)->get_op_type();
            v2 = Stack.top();
            Stack.pop();
            v1 = Stack.top();
            Stack.pop();
            switch (op2)
              {
              case BinaryOpcode::plus:
                tmp_out.str("");
                tmp_out << v1 << " + " << v2;
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::minus:
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
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                if (it_code == faulty_op)
                  tmp_out << "{ / }";
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
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                if (it_code == faulty_op)
                  tmp_out << "{ ^ }";
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
                  tmp_out.str("");
                  if (it_code == faulty_op)
                    tmp_out << "{PowerDeriv}";
                  else
                    tmp_out << "PowerDeriv";
                  tmp_out << "(" << v1 << ", " << v2 << ", " << v3 << ")";
                  Stack.push(tmp_out.str());
                }
                break;
              case BinaryOpcode::max:
                tmp_out.str("");
                tmp_out << "max(" << v1 << ", " << v2 << ")";
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::min:
                tmp_out.str("");
                tmp_out << "min(" << v1 << ", " << v2 << ")";
                Stack.push(tmp_out.str());
                break;
              case BinaryOpcode::equal:
              default:
                mexPrintf("Error unknown binary operator=%d\n", static_cast<int>(op2)); mexEvalString("drawnow;");
                ;
              }
            break;
          case Tags::FUNARY:
            op1 = static_cast<FUNARY_ *>(*it_code)->get_op_type();
            v1 = Stack.top();
            Stack.pop();
            switch (op1)
              {
              case UnaryOpcode::uminus:
                tmp_out.str("");
                tmp_out << " -" << v1;
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::exp:
                tmp_out.str("");
                tmp_out << "exp(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::log:
                tmp_out.str("");
                if (it_code == faulty_op)
                  tmp_out << "{log}";
                else
                  tmp_out << "log";
                tmp_out << "(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::log10:
                tmp_out.str("");
                if (it_code == faulty_op)
                  tmp_out << "{log10}";
                else
                  tmp_out << "log10";
                tmp_out << "(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::cos:
                tmp_out.str("");
                tmp_out << "cos(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::sin:
                tmp_out.str("");
                tmp_out << "sin(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::tan:
                tmp_out.str("");
                tmp_out << "tan(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::acos:
                tmp_out.str("");
                tmp_out << "acos(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::asin:
                tmp_out.str("");
                tmp_out << "asin(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::atan:
                tmp_out.str("");
                tmp_out << "atan(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::cosh:
                tmp_out.str("");
                tmp_out << "cosh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::sinh:
                tmp_out.str("");
                tmp_out << "sinh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::tanh:
                tmp_out.str("");
                tmp_out << "tanh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::acosh:
                tmp_out.str("");
                tmp_out << "acosh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::asinh:
                tmp_out.str("");
                tmp_out << "asinh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::atanh:
                tmp_out.str("");
                tmp_out << "atanh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::sqrt:
                tmp_out.str("");
                tmp_out << "sqrt(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::erf:
                tmp_out.str("");
                tmp_out << "erf(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case UnaryOpcode::erfc:
                tmp_out.str("");
                tmp_out << "erfc(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              default:
                mexPrintf("Error unknown unary operator=%d\n", static_cast<int>(op1)); mexEvalString("drawnow;");
                ;
              }
            break;
          case Tags::FTRINARY:
            op3 = static_cast<FTRINARY_ *>(*it_code)->get_op_type();
            v3 = Stack.top();
            Stack.pop();
            v2 = Stack.top();
            Stack.pop();
            v1 = Stack.top();
            Stack.pop();
            switch (op3)
              {
              case TrinaryOpcode::normcdf:
                tmp_out.str("");
                tmp_out << "normcdf(" << v1 << ", " << v2 << ", " << v3 << ")";
                Stack.push(tmp_out.str());
                break;
              case TrinaryOpcode::normpdf:
                tmp_out.str("");
                tmp_out << "normpdf(" << v1 << ", " << v2 << ", " << v3 << ")";
                Stack.push(tmp_out.str());
                break;
              default:
                mexPrintf("Error unknown trinary operator=%d\n", static_cast<int>(op3)); mexEvalString("drawnow;");
              }
            break;
          case Tags::FCALL:
            {
              auto *fc = static_cast<FCALL_ *>(*it_code);
              string function_name = fc->get_function_name();
              int nb_input_arguments{fc->get_nb_input_arguments()};

              string arg_func_name = fc->get_arg_func_name();
              int nb_add_input_arguments{fc->get_nb_add_input_arguments()};
              call_type = fc->get_call_type();
              switch (call_type)
                {
                case ExternalFunctionCallType::levelWithoutDerivative:
                case ExternalFunctionCallType::levelWithFirstDerivative:
                case ExternalFunctionCallType::levelWithFirstAndSecondDerivative:
                  {
                    tmp_out.str("");
                    tmp_out << function_name << "(";
                    vector<string> ss(nb_input_arguments);
                    for (int i{0}; i < nb_input_arguments; i++)
                      {
                        ss[nb_input_arguments-i-1] = Stack.top();
                        Stack.pop();
                      }
                    for (int i{0}; i < nb_input_arguments; i++)
                      {
                        tmp_out << ss[i];
                        if (i < nb_input_arguments - 1)
                          tmp_out << ", ";
                      }
                    tmp_out << ")";
                    Stack.push(tmp_out.str());
                  }
                  break;
                case ExternalFunctionCallType::numericalFirstDerivative:
                  {
                    tmp_out.str("");
                    tmp_out << function_name << "(";
                    tmp_out << arg_func_name.c_str() << ", " << fc->get_row() << ", {";
                    vector<string> ss(nb_input_arguments);
                    for (int i{0}; i < nb_add_input_arguments; i++)
                      {
                        ss[nb_add_input_arguments-i-1] = Stack.top();
                        Stack.pop();
                      }
                    for (int i{0}; i < nb_add_input_arguments; i++)
                      {
                        tmp_out << ss[i];
                        if (i < nb_add_input_arguments - 1)
                          tmp_out << ", ";
                      }
                    tmp_out << "})";
                    Stack.push(tmp_out.str());
                  }
                  break;
                case ExternalFunctionCallType::separatelyProvidedFirstDerivative:
                  {
                    tmp_out.str("");
                    tmp_out << function_name << "(";
                    vector<string> ss(nb_input_arguments);
                    for (int i{0}; i < nb_input_arguments; i++)
                      {
                        ss[nb_input_arguments-i-1] = Stack.top();
                        Stack.pop();
                      }
                    for (int i{0}; i < nb_input_arguments; i++)
                      {
                        tmp_out << ss[i];
                        if (i < nb_input_arguments - 1)
                          tmp_out << ", ";
                      }
                    tmp_out << ")";
                    Stack.push(tmp_out.str());
                  }
                  break;
                case ExternalFunctionCallType::numericalSecondDerivative:
                  {
                    tmp_out.str("");
                    tmp_out << function_name << "(";
                    tmp_out << arg_func_name.c_str() << ", " << fc->get_row() << ", " << fc->get_col() << ", {";
                    vector<string> ss(nb_input_arguments);
                    for (int i{0}; i < nb_add_input_arguments; i++)
                      {
                        ss[nb_add_input_arguments-i-1] = Stack.top();
                        Stack.pop();
                      }
                    for (int i{0}; i < nb_add_input_arguments; i++)
                      {
                        tmp_out << ss[i];
                        if (i < nb_add_input_arguments - 1)
                          tmp_out << ", ";
                      }
                    tmp_out << "})";
                    Stack.push(tmp_out.str());
                  }
                  break;
                case ExternalFunctionCallType::separatelyProvidedSecondDerivative:
                  {
                    tmp_out.str("");
                    tmp_out << function_name << "(";
                    vector<string> ss(nb_input_arguments);
                    for (int i{0}; i < nb_input_arguments; i++)
                      {
                        ss[nb_input_arguments-i-1] = Stack.top();
                        Stack.pop();
                      }
                    for (int i{0}; i < nb_input_arguments; i++)
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
            var = static_cast<FSTPTEF_ *>(*it_code)->get_number();
            tmp_out.str("");
            switch (call_type)
              {
              case ExternalFunctionCallType::levelWithoutDerivative:
                tmp_out << "TEF(" << var << ") = " << Stack.top();
                break;
              case ExternalFunctionCallType::levelWithFirstDerivative:
                tmp_out << "[TEF(" << var << "), TEFD(" << var << ") ]= " << Stack.top();
                break;
              case ExternalFunctionCallType::levelWithFirstAndSecondDerivative:
                tmp_out << "[TEF(" << var << "), TEFD(" << var << "), TEFDD(" << var << ") ]= " << Stack.top();
                break;
              default:
                break;
              }
            Stack.pop();
            break;
          case Tags::FLDTEF:
            var = static_cast<FLDTEF_ *>(*it_code)->get_number();
            tmp_out.str("");
            tmp_out << "TEF(" << var << ")";
            Stack.push(tmp_out.str());

            break;
          case Tags::FSTPTEFD:
            {
              go_on = false;
              int indx{static_cast<FSTPTEFD_ *>(*it_code)->get_indx()};
              int row{static_cast<FSTPTEFD_ *>(*it_code)->get_row()};
              tmp_out.str("");
              if (call_type == ExternalFunctionCallType::numericalFirstDerivative)
                tmp_out << "TEFD(" << indx << ", " << row << ") = " << Stack.top();
              else if (call_type == ExternalFunctionCallType::separatelyProvidedFirstDerivative)
                tmp_out << "TEFD(" << indx << ") = " << Stack.top();
              Stack.pop();
            }
            break;
          case Tags::FLDTEFD:
            {
              int indx{static_cast<FLDTEFD_ *>(*it_code)->get_indx()};
              int row{static_cast<FLDTEFD_ *>(*it_code)->get_row()};
              tmp_out.str("");
              tmp_out << "TEFD(" << indx << ", " << row << ")";
              Stack.push(tmp_out.str());
            }
            break;
          case Tags::FSTPTEFDD:
            {
              go_on = false;
              int indx{static_cast<FSTPTEFDD_ *>(*it_code)->get_indx()};
              int row{static_cast<FSTPTEFDD_ *>(*it_code)->get_row()};
              int col{static_cast<FSTPTEFDD_ *>(*it_code)->get_col()};
              tmp_out.str("");
              if (call_type == ExternalFunctionCallType::numericalSecondDerivative)
                tmp_out << "TEFDD(" << indx << ", " << row << ", " << col << ") = " << Stack.top();
              else if (call_type == ExternalFunctionCallType::separatelyProvidedSecondDerivative)
                tmp_out << "TEFDD(" << indx << ") = " << Stack.top();
              Stack.pop();
            }

            break;
          case Tags::FLDTEFDD:
            {
              int indx{static_cast<FLDTEFDD_ *>(*it_code)->get_indx()};
              int row{static_cast<FLDTEFDD_ *>(*it_code)->get_row()};
              int col{static_cast<FSTPTEFDD_ *>(*it_code)->get_col()};
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
          case Tags::FENDBLOCK:
          case Tags::FENDEQU:
            go_on = false;
            break;
          default:
            mexPrintf("Error tag=%d unknown\n", (*it_code)->op_code); mexEvalString("drawnow;");
            throw FatalExceptionHandling(" in print_expression, unknown opcode "
                                         + to_string(static_cast<int>((*it_code)->op_code))
                                         + "!! FENDEQU="
                                         + to_string(static_cast<int>(Tags::FENDEQU)) + "\n");
          }
        it_code++;
      }
    return { tmp_out.str(), it_code };
  }
  
public:
  static void
  test_mxMalloc(void *z, int line, const string &file, const string &func, int amount)
  {
    if (!z && amount > 0)
      throw FatalExceptionHandling(" mxMalloc: out of memory " + to_string(amount) + " bytes required at line " + to_string(line) + " in function " + func + " (file " + file);
  }
};

#endif
