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
    switch (EQN_type)
      {
      case ExpressionType::TemporaryTerm:
        Error_loc << "temporary term";
        break;
      case ExpressionType::ModelEquation:
        Error_loc << "equation";
        break;
      case ExpressionType::FirstEndoDerivative:
      case ExpressionType::FirstOtherEndoDerivative:
      case ExpressionType::FirstExoDerivative:
      case ExpressionType::FirstExodetDerivative:
        Error_loc << "first order derivative of equation";
        break;
      }
    Error_loc << " " << EQN_equation+1;
    if (EQN_block_number > 1)
      Error_loc << " in block " << EQN_block+1;
    switch (EQN_type)
      {
      case ExpressionType::TemporaryTerm:
      case ExpressionType::ModelEquation:
        break;
      case ExpressionType::FirstEndoDerivative:
        Error_loc << " with respect to endogenous variable " << get_variable(SymbolType::endogenous, EQN_dvar1);
        break;
      case ExpressionType::FirstOtherEndoDerivative:
        Error_loc << " with respect to other endogenous variable " << get_variable(SymbolType::endogenous, EQN_dvar1);
        break;
      case ExpressionType::FirstExoDerivative:
        Error_loc << " with respect to exogenous variable " << get_variable(SymbolType::exogenous, EQN_dvar1);
        break;
      case ExpressionType::FirstExodetDerivative:
        Error_loc << " with respect to deterministic exogenous variable " << get_variable(SymbolType::exogenousDet, EQN_dvar1);
        break;
      }
    if (!steady_state)
      Error_loc << " at time " << it_;

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
    /* First element is output string, 2nd element is precedence of last
       operator, 3rd element is opcode if the last operator was a binary
       operator.
       Use same precedence values as in the preprocessor for MATLAB or JSON output:
       – 0: = (assignment)
       – 1: == !=
       – 2: < > <= >=
       – 3: + -
       – 4: * /
       – 5: ^ powerDeriv
       – 100: everything else */
    stack<tuple<string, int, optional<BinaryOpcode>>> Stack;
    bool go_on {true};
    ExpressionType equation_type {ExpressionType::ModelEquation};
    ExternalFunctionCallType call_type {ExternalFunctionCallType::levelWithoutDerivative};

    auto lag_to_string = [](int l) -> string
    {
      if (l > 0)
        return "(+" + to_string(l) + ")";
      else if (l < 0)
        return "(" + to_string(l) + ")";
      else
        return {};
    };

    // Helper for FST* tags
    auto assign_lhs = [&](const string &lhs)
    {
      auto &[str, prec, opcode] {Stack.top()};
      str.insert(0, lhs + " = ");
      prec = 0;
      opcode = BinaryOpcode::equal;
      go_on = false;
    };

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
                throw FatalExceptionHandling{" in print_expression, expression type "
                                             + to_string(static_cast<int>(static_cast<FNUMEXPR_ *>(*it_code)->get_expression_type()))
                                             + " not implemented yet\n"};
              }
            break;
          case Tags::FLDV:
            {
              int var {static_cast<FLDV_ *>(*it_code)->get_pos()};
              int lag {static_cast<FLDV_ *>(*it_code)->get_lead_lag()};
              switch (SymbolType type {static_cast<FLDV_ *>(*it_code)->get_type()};
                      type)
                {
                case SymbolType::parameter:
                case SymbolType::endogenous:
                case SymbolType::exogenous:
                case SymbolType::exogenousDet:
                  Stack.emplace(get_variable(type, var) + lag_to_string(lag), 100, nullopt);
                  break;
                default:
                  throw FatalExceptionHandling{"FLDV: Unknown variable type\n"};
              }
            }
            break;
          case Tags::FLDSV:
            {
              int var {static_cast<FLDSV_ *>(*it_code)->get_pos()};
              switch (SymbolType type {static_cast<FLDSV_ *>(*it_code)->get_type()};
                      type)
                {
                case SymbolType::parameter:
                case SymbolType::endogenous:
                case SymbolType::exogenous:
                case SymbolType::exogenousDet:
                  Stack.emplace(get_variable(type, var), 100, nullopt);
                  break;
                default:
                  throw FatalExceptionHandling{"FLDSV: Unknown variable type\n"};
                }
            }
            break;
          case Tags::FLDVS:
            {
              int var {static_cast<FLDVS_ *>(*it_code)->get_pos()};
              switch (SymbolType type {static_cast<FLDVS_ *>(*it_code)->get_type()};
                      type)
                {
                case SymbolType::parameter:
                case SymbolType::endogenous:
                case SymbolType::exogenous:
                case SymbolType::exogenousDet:
                  Stack.emplace(get_variable(type, var), 100, nullopt);
                  break;
                default:
                  throw FatalExceptionHandling{"FLDVS: Unknown variable type\n"};
                }
            }
            break;
          case Tags::FLDT:
            Stack.emplace("T" + to_string(static_cast<FLDT_ *>(*it_code)->get_pos() + 1), 100, nullopt);
            break;
          case Tags::FLDST:
            Stack.emplace("T" + to_string(static_cast<FLDST_ *>(*it_code)->get_pos() + 1), 100, nullopt);
            break;
          case Tags::FLDU:
            Stack.emplace("u(" + to_string(static_cast<FLDU_ *>(*it_code)->get_pos() + 1) + " + it_)", 100, nullopt);
            break;
          case Tags::FLDSU:
            Stack.emplace("u(" + to_string(static_cast<FLDSU_ *>(*it_code)->get_pos() + 1) + ")", 100, nullopt);
            break;
          case Tags::FLDR:
            Stack.emplace("residual(" + to_string(static_cast<FLDR_ *>(*it_code)->get_pos() + 1) + ")", 100, nullopt);
            break;
          case Tags::FLDZ:
            Stack.emplace("0", 100, nullopt);
            break;
          case Tags::FLDC:
            Stack.emplace(to_string(static_cast<FLDC_ *>(*it_code)->get_value()), 100, nullopt);
            break;
          case Tags::FSTPV:
            {
              int var {static_cast<FSTPV_ *>(*it_code)->get_pos()};
              int lag {static_cast<FSTPV_ *>(*it_code)->get_lead_lag()};
              switch (SymbolType type {static_cast<FSTPV_ *>(*it_code)->get_type()};
                      type)
                {
                case SymbolType::parameter:
                case SymbolType::endogenous:
                case SymbolType::exogenous:
                case SymbolType::exogenousDet:
                  assign_lhs(get_variable(type, var) + lag_to_string(lag));
                  break;
                default:
                  throw FatalExceptionHandling{"FSTPV: Unknown variable type\n"};
                }
            }
            break;
          case Tags::FSTPSV:
            {
              int var {static_cast<FSTPSV_ *>(*it_code)->get_pos()};
              switch (SymbolType type {static_cast<FSTPSV_ *>(*it_code)->get_type()};
                      type)
                {
                case SymbolType::parameter:
                case SymbolType::endogenous:
                case SymbolType::exogenous:
                case SymbolType::exogenousDet:
                  assign_lhs(get_variable(type, var));
                  break;
                default:
                  throw FatalExceptionHandling{"FSTPSV: Unknown variable type\n"};
                }
            }
            break;
          case Tags::FSTPT:
            assign_lhs("T" + to_string(static_cast<FSTPT_ *>(*it_code)->get_pos() + 1));
            break;
          case Tags::FSTPST:
            assign_lhs("T" + to_string(static_cast<FSTPST_ *>(*it_code)->get_pos() + 1));
            break;
          case Tags::FSTPU:
            assign_lhs("u(" + to_string(static_cast<FSTPU_ *>(*it_code)->get_pos() + 1) + " + it_)");
            break;
          case Tags::FSTPSU:
            assign_lhs("u(" + to_string(static_cast<FSTPSU_ *>(*it_code)->get_pos() + 1) + ")");
            break;
          case Tags::FSTPR:
            assign_lhs("residual(" + to_string(static_cast<FSTPR_ *>(*it_code)->get_pos() + 1) + ")");
            break;
          case Tags::FSTPG:
            assign_lhs("g1(" + to_string(static_cast<FSTPG_ *>(*it_code)->get_pos() + 1) + ")");
            break;
          case Tags::FSTPG2:
            {
              int eq {static_cast<FSTPG2_ *>(*it_code)->get_row()};
              int var {static_cast<FSTPG2_ *>(*it_code)->get_col()};
              assign_lhs("jacob(" + to_string(eq+1) + ", " + to_string(var+1) + ")");
            }
            break;
          case Tags::FSTPG3:
            {
              int eq {static_cast<FSTPG3_ *>(*it_code)->get_row()};
              int var {static_cast<FSTPG3_ *>(*it_code)->get_col()};
              int lag {static_cast<FSTPG3_ *>(*it_code)->get_lag()};
              int col_pos {static_cast<FSTPG3_ *>(*it_code)->get_col_pos()};
              string matrix_name { [&]
              {
                switch (equation_type)
                  {
                  case ExpressionType::FirstEndoDerivative:
                    return "jacob";
                  case ExpressionType::FirstOtherEndoDerivative:
                    return "jacob_other_endo";
                  case ExpressionType::FirstExoDerivative:
                    return "jacob_exo";
                  case ExpressionType::FirstExodetDerivative:
                    return "jacob_exo_det";
                  default:
                    throw FatalExceptionHandling{" unknown equation type " + to_string(static_cast<int>(equation_type)) + "\n"};
                  }
              }() };

              assign_lhs(matrix_name + "(" + to_string(eq+1) + ", " + to_string(col_pos+1)
                         + " [var=" + to_string(var+1) + ", lag=" + to_string(lag) + "])");
            }
            break;
          case Tags::FUNARY:
            {
              UnaryOpcode op {static_cast<FUNARY_ *>(*it_code)->get_op_type()};
              auto [arg, prec_arg, op_arg] {Stack.top()};
              Stack.pop();

              ostringstream s;

              // Always put parenthesis around uminus nodes
              if (op == UnaryOpcode::uminus)
                s << "(";

              s << [&]
              {
                switch (op)
                  {
                  case UnaryOpcode::uminus:
                    return "-";
                  case UnaryOpcode::exp:
                    return "exp";
                  case UnaryOpcode::log:
                    return it_code == faulty_op ? "{log}" : "log";
                  case UnaryOpcode::log10:
                    return it_code == faulty_op ? "{log10}" : "log10";
                  case UnaryOpcode::cos:
                    return "cos";
                  case UnaryOpcode::sin:
                    return "sin";
                  case UnaryOpcode::tan:
                    return "tan";
                  case UnaryOpcode::acos:
                    return "acos";
                  case UnaryOpcode::asin:
                    return "asin";
                  case UnaryOpcode::atan:
                    return "atan";
                  case UnaryOpcode::cosh:
                    return "cosh";
                  case UnaryOpcode::sinh:
                    return "sinh";
                  case UnaryOpcode::tanh:
                    return "tanh";
                  case UnaryOpcode::acosh:
                    return "acosh";
                  case UnaryOpcode::asinh:
                    return "asinh";
                  case UnaryOpcode::atanh:
                    return "atanh";
                  case UnaryOpcode::sqrt:
                    return "sqrt";
                  case UnaryOpcode::cbrt:
                    return "cbrt";
                  case UnaryOpcode::erf:
                    return "erf";
                  case UnaryOpcode::erfc:
                    return "erfc";
                  case UnaryOpcode::abs:
                    return "abs";
                  case UnaryOpcode::sign:
                    return "sign";
                  case UnaryOpcode::steadyState:
                    return "steady_state";
                  case UnaryOpcode::steadyStateParamDeriv:
                  case UnaryOpcode::steadyStateParam2ndDeriv:
                    throw FatalExceptionHandling{"Unexpected derivative of steady_state operator"};
                  case UnaryOpcode::expectation:
                    throw FatalExceptionHandling{"Unexpected expectation operator"};
                  case UnaryOpcode::diff:
                    return "diff";
                  case UnaryOpcode::adl:
                    return "adl";
                  }
                throw FatalExceptionHandling{"Unknown opcode"};
              }();

              /* Print argument. Enclose it with parentheses if:
                 - current opcode is not uminus, or
                 - current opcode is uminus and argument has lowest precedence */
              bool close_parenthesis {false};
              if (op != UnaryOpcode::uminus
                  || (op == UnaryOpcode::uminus && prec_arg < 100))
                {
                  s << "(";
                  close_parenthesis = true;
                }
              s << arg;
              if (close_parenthesis)
                s << ")";

              // Close parenthesis for uminus
              if (op == UnaryOpcode::uminus)
                s << ")";

              Stack.emplace(s.str(), 100, nullopt);
            }
            break;
          case Tags::FBINARY:
            {
              BinaryOpcode op {static_cast<FBINARY_ *>(*it_code)->get_op_type()};
              auto [arg2, prec_arg2, op_arg2] {Stack.top()};
              Stack.pop();
              auto [arg1, prec_arg1, op_arg1] {Stack.top()};
              Stack.pop();

              if (op == BinaryOpcode::max)
                Stack.emplace("max(" + arg1 + ", " + arg2 + ")", 100, nullopt);
              else if (op == BinaryOpcode::min)
                Stack.emplace("min(" + arg1 + ", " + arg2 + ")", 100, nullopt);
              else if (op == BinaryOpcode::powerDeriv)
                {
                  auto [arg3, prec_arg3, op_arg3] {Stack.top()};
                  Stack.pop();
                  Stack.emplace((it_code == faulty_op ? "{PowerDeriv}(" : "PowerDeriv(") + arg1
                                + ", " + arg2 + ", " + arg2 + ")", 100, nullopt);
                }
              else
                {
                  ostringstream s;
                  int prec { [&]
                  {
                    switch (op)
                      {
                      case BinaryOpcode::equal:
                        return 0;
                      case BinaryOpcode::equalEqual:
                      case BinaryOpcode::different:
                        return 1;
                      case BinaryOpcode::lessEqual:
                      case BinaryOpcode::greaterEqual:
                      case BinaryOpcode::less:
                      case BinaryOpcode::greater:
                        return 2;
                      case BinaryOpcode::plus:
                      case BinaryOpcode::minus:
                        return 3;
                      case BinaryOpcode::times:
                      case BinaryOpcode::divide:
                        return 4;
                      case BinaryOpcode::power:
                      case BinaryOpcode::powerDeriv:
                        return 5;
                      case BinaryOpcode::min:
                      case BinaryOpcode::max:
                        return 100;
                      }
                    throw FatalExceptionHandling{"Unknown opcode"};
                  }() };

                  /* Print left argument. If left argument has a lower
                     precedence, or if current and left argument are both power
                     operators, add parenthesis around left argument. */
                  bool close_parenthesis {false};
                  if (prec_arg1 < prec
                      || (op == BinaryOpcode::power && op_arg1 == BinaryOpcode::power))
                    {
                      s << "(";
                      close_parenthesis = true;
                    }
                  s << arg1;
                  if (close_parenthesis)
                    s << ")";

                  s << [&]
                  {
                    switch (op)
                    {
                    case BinaryOpcode::plus:
                      return "+";
                    case BinaryOpcode::minus:
                      return "-";
                    case BinaryOpcode::times:
                      return "*";
                    case BinaryOpcode::divide:
                      return it_code == faulty_op ? "{ / }" : "/";
                    case BinaryOpcode::less:
                      return " < ";
                    case BinaryOpcode::greater:
                      return " > ";
                    case BinaryOpcode::lessEqual:
                      return " <= ";
                    case BinaryOpcode::greaterEqual:
                      return " >= ";
                    case BinaryOpcode::equalEqual:
                      return " == ";
                    case BinaryOpcode::different:
                      return " != ";
                    case BinaryOpcode::power:
                      return it_code == faulty_op ? "{ ^ }" : "^";
                    case BinaryOpcode::powerDeriv:
                    case BinaryOpcode::max:
                    case BinaryOpcode::min:
                      throw FatalExceptionHandling{"Should not arrive here"};
                    case BinaryOpcode::equal:
                      return " = ";
                    }
                    throw FatalExceptionHandling{"Unknown opcode"};
                  }();

                  /* Print right argument. Add parenthesis around right argument if:
                     - its precedence is lower than that of the current node
                     - it is a power operator and current operator is also a power operator
                     - it has same precedence as current operator and current operator is
                     either a minus or a divide */
                  close_parenthesis = false;
                  if (prec_arg2 < prec
                      || (op == BinaryOpcode::power && op_arg2 == BinaryOpcode::power)
                      || (prec_arg2 == prec && (op == BinaryOpcode::minus || op == BinaryOpcode::divide)))
                    {
                      s << "(";
                      close_parenthesis = true;
                    }
                  s << arg2;
                  if (close_parenthesis)
                    s << ")";

                  Stack.emplace(s.str(), prec, op);
                }
            }
            break;
          case Tags::FTRINARY:
            {
              TrinaryOpcode op {static_cast<FTRINARY_ *>(*it_code)->get_op_type()};
              auto [arg3, prec_arg3, op_arg3] {Stack.top()};
              Stack.pop();
              auto [arg2, prec_arg2, op_arg2] {Stack.top()};
              Stack.pop();
              auto [arg1, prec_arg1, op_arg1] {Stack.top()};
              Stack.pop();

              string opname { [&]
              {
                switch (op)
                  {
                  case TrinaryOpcode::normcdf:
                    return "normcdf";
                  case TrinaryOpcode::normpdf:
                    return "normpdf";
                  }
                throw FatalExceptionHandling{"Unknown opcode"};
              }() };

              Stack.emplace(opname + "(" + arg1 + ", " + arg2 + ", " + arg3 + ")", 100, nullopt);
            }
            break;
          case Tags::FCALL:
            {
              auto *fc = static_cast<FCALL_ *>(*it_code);
              string function_name {fc->get_function_name()};
              int nb_input_arguments{fc->get_nb_input_arguments()};
              int nb_add_input_arguments{fc->get_nb_add_input_arguments()};
              string arg_func_name {fc->get_arg_func_name()};
              ostringstream s;

              auto print_args = [&](int nargs)
              {
                vector<string> ss(nargs);
                for (int i {0}; i < nargs; i++)
                  {
                    ss[nargs-i-1] = get<0>(Stack.top());
                    Stack.pop();
                  }
                for (int i {0}; i < nargs; i++)
                  {
                    if (i > 0)
                      s << ", ";
                    s << ss[i];
                  }
              };

              call_type = fc->get_call_type();
              switch (call_type)
                {
                case ExternalFunctionCallType::levelWithoutDerivative:
                case ExternalFunctionCallType::levelWithFirstDerivative:
                case ExternalFunctionCallType::levelWithFirstAndSecondDerivative:
                case ExternalFunctionCallType::separatelyProvidedFirstDerivative:
                case ExternalFunctionCallType::separatelyProvidedSecondDerivative:
                  s << function_name << "(";
                  print_args(nb_input_arguments);
                  s << ")";
                  break;
                case ExternalFunctionCallType::numericalFirstDerivative:
                  s << function_name << "(" << arg_func_name << ", " << fc->get_row()+1 << ", {";
                  print_args(nb_add_input_arguments);
                  s << "})";
                  break;
                case ExternalFunctionCallType::numericalSecondDerivative:
                  s << function_name << "(" << arg_func_name << ", " << fc->get_row()+1 << ", " << fc->get_col()+1 << ", {";
                  print_args(nb_add_input_arguments);
                  s << "})";
                  break;
                }
              Stack.emplace(s.str(), 100, nullopt);
            }
            break;
          case Tags::FSTPTEF:
            {
              int indx {static_cast<FSTPTEF_ *>(*it_code)->get_number()};
              switch (call_type)
              {
              case ExternalFunctionCallType::levelWithoutDerivative:
                assign_lhs("TEF(" + to_string(indx+1) + ")");
                break;
              case ExternalFunctionCallType::levelWithFirstDerivative:
                assign_lhs("[TEF(" + to_string(indx+1) + "), TEFD("  + to_string(indx+1) + ") ]");
                break;
              case ExternalFunctionCallType::levelWithFirstAndSecondDerivative:
                assign_lhs("[TEF(" + to_string(indx+1) + "), TEFD("  + to_string(indx+1) + "), TEFDD("  + to_string(indx+1) + ") ]");
                break;
              default:
                throw FatalExceptionHandling{"Unexpected external function call type"};
              }
            }
            break;
          case Tags::FLDTEF:
            Stack.emplace("TEF(" + to_string(static_cast<FLDTEF_ *>(*it_code)->get_number()+1) + ")", 100, nullopt);
            break;
          case Tags::FSTPTEFD:
            {
              int indx {static_cast<FSTPTEFD_ *>(*it_code)->get_indx()};
              int row {static_cast<FSTPTEFD_ *>(*it_code)->get_row()};
              if (call_type == ExternalFunctionCallType::numericalFirstDerivative)
                assign_lhs("TEFD(" + to_string(indx+1) + ", " + to_string(row+1) + ")");
              else if (call_type == ExternalFunctionCallType::separatelyProvidedFirstDerivative)
                assign_lhs("TEFD(" + to_string(indx+1) + ")");
            }
            break;
          case Tags::FLDTEFD:
            {
              int indx {static_cast<FLDTEFD_ *>(*it_code)->get_indx()};
              int row {static_cast<FLDTEFD_ *>(*it_code)->get_row()};
              Stack.emplace("TEFD(" + to_string(indx+1) + ", " + to_string(row+1) + ")", 100, nullopt);
            }
            break;
          case Tags::FSTPTEFDD:
            {
              int indx {static_cast<FSTPTEFDD_ *>(*it_code)->get_indx()};
              int row {static_cast<FSTPTEFDD_ *>(*it_code)->get_row()};
              int col {static_cast<FSTPTEFDD_ *>(*it_code)->get_col()};
              if (call_type == ExternalFunctionCallType::numericalSecondDerivative)
                assign_lhs("TEFDD(" + to_string(indx+1) + ", " + to_string(row+1) + ", " + to_string(col+1) + ")");
              else if (call_type == ExternalFunctionCallType::separatelyProvidedSecondDerivative)
                assign_lhs("TEFDD(" + to_string(indx+1) + ")");
            }
            break;
          case Tags::FLDTEFDD:
            {
              int indx {static_cast<FLDTEFDD_ *>(*it_code)->get_indx()};
              int row {static_cast<FLDTEFDD_ *>(*it_code)->get_row()};
              int col {static_cast<FSTPTEFDD_ *>(*it_code)->get_col()};
              Stack.emplace("TEFDD(" + to_string(indx+1) + ", " + to_string(row+1) + ", " + to_string(col+1) + ")", 100, nullopt);
            }
            break;
          case Tags::FJMPIFEVAL:
            Stack.emplace("if (~evaluate)", 100, nullopt);
            go_on = false;
            break;
          case Tags::FJMP:
            Stack.emplace("else", 100, nullopt);
            go_on = false;
            break;
          case Tags::FENDBLOCK:
          case Tags::FENDEQU:
            go_on = false;
            break;
          default:
            throw FatalExceptionHandling(" in print_expression, unknown opcode "
                                         + to_string(static_cast<int>((*it_code)->op_code)) + "\n");
          }
        it_code++;
      }
    return { get<0>(Stack.top()), it_code };
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
