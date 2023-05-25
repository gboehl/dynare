/*
 * Copyright © 2013-2023 Dynare Team
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
#include <stack>
#include <cfenv>
#include <numbers>

#include <dynmex.h>

#include "Evaluate.hh"
#include "CommonEnums.hh"
#include "ErrorHandling.hh"

Evaluate::Evaluate(const filesystem::path &codfile, bool steady_state_arg, const BasicSymbolTable &symbol_table_arg) :
  symbol_table {symbol_table_arg},
  steady_state {steady_state_arg}
{
  ifstream CompiledCode {codfile, ios::in | ios::binary | ios::ate};
  if (!CompiledCode.is_open())
    throw FatalException {codfile.string() + " cannot be opened"};
  auto Code_Size {CompiledCode.tellg()};
  raw_bytecode = make_unique<char[]>(Code_Size);
  auto code {raw_bytecode.get()};
  CompiledCode.seekg(0);
  CompiledCode.read(code, Code_Size);
  CompiledCode.close();

  bool done {false};
  while (!done)
    {
      BytecodeInstruction *instr {reinterpret_cast<BytecodeInstruction *>(code)};
      switch (*reinterpret_cast<Tags *>(code))
        {
        case Tags::FLDZ:
# ifdef DEBUGL
          mexPrintf("FLDZ\n");
# endif
          code += sizeof(FLDZ_);
          break;
        case Tags::FEND:
# ifdef DEBUGL
          mexPrintf("FEND\n");
# endif
          code += sizeof(FEND_);
          done = true;
          break;
        case Tags::FENDBLOCK:
# ifdef DEBUGL
          mexPrintf("FENDBLOCK\n");
# endif
          code += sizeof(FENDBLOCK_);
          break;
        case Tags::FENDEQU:
# ifdef DEBUGL
          mexPrintf("FENDEQU\n");
# endif
          code += sizeof(FENDEQU_);
          break;
        case Tags::FDIMT:
# ifdef DEBUGL
          mexPrintf("FDIMT\n");
# endif
          code += sizeof(FDIMT_);
          break;
        case Tags::FDIMST:
# ifdef DEBUGL
          mexPrintf("FDIMST\n");
# endif
          code += sizeof(FDIMST_);
          break;
        case Tags::FNUMEXPR:
# ifdef DEBUGL
          mexPrintf("FNUMEXPR\n");
# endif
          code += sizeof(FNUMEXPR_);
          break;
        case Tags::FLDC:
# ifdef DEBUGL
          mexPrintf("FLDC\n");
# endif
          code += sizeof(FLDC_);
          break;
        case Tags::FLDU:
# ifdef DEBUGL
          mexPrintf("FLDU\n");
# endif
          code += sizeof(FLDU_);
          break;
        case Tags::FLDSU:
# ifdef DEBUGL
          mexPrintf("FLDSU\n");
# endif
          code += sizeof(FLDSU_);
          break;
        case Tags::FLDR:
# ifdef DEBUGL
          mexPrintf("FLDR\n");
# endif
          code += sizeof(FLDR_);
          break;
        case Tags::FLDT:
# ifdef DEBUGL
          mexPrintf("FLDT\n");
# endif
          code += sizeof(FLDT_);
          break;
        case Tags::FLDST:
# ifdef DEBUGL
          mexPrintf("FLDST\n");
# endif
          code += sizeof(FLDST_);
          break;
        case Tags::FSTPT:
# ifdef DEBUGL
          mexPrintf("FSTPT\n");
# endif
          code += sizeof(FSTPT_);
          break;
        case Tags::FSTPST:
# ifdef DEBUGL
          mexPrintf("FSTPST\n");
# endif
          code += sizeof(FSTPST_);
          break;
        case Tags::FSTPR:
# ifdef DEBUGL
          mexPrintf("FSTPR\n");
# endif
          code += sizeof(FSTPR_);
          break;
        case Tags::FSTPU:
# ifdef DEBUGL
          mexPrintf("FSTPU\n");
# endif
          code += sizeof(FSTPU_);
          break;
        case Tags::FSTPSU:
# ifdef DEBUGL
          mexPrintf("FSTPSU\n");
# endif
          code += sizeof(FSTPSU_);
          break;
        case Tags::FSTPG:
# ifdef DEBUGL
          mexPrintf("FSTPG\n");
# endif
          code += sizeof(FSTPG_);
          break;
        case Tags::FSTPG2:
# ifdef DEBUGL
          mexPrintf("FSTPG2\n");
# endif
          code += sizeof(FSTPG2_);
          break;
        case Tags::FSTPG3:
# ifdef DEBUGL
          mexPrintf("FSTPG3\n");
# endif
          code += sizeof(FSTPG3_);
          break;
        case Tags::FUNARY:
# ifdef DEBUGL
          mexPrintf("FUNARY\n");
# endif
          code += sizeof(FUNARY_);
          break;
        case Tags::FBINARY:
# ifdef DEBUGL
          mexPrintf("FBINARY\n");
# endif
          code += sizeof(FBINARY_);
          break;
        case Tags::FTRINARY:
# ifdef DEBUGL
          mexPrintf("FTRINARY\n");
# endif
          code += sizeof(FTRINARY_);
          break;
        case Tags::FLDVS:
# ifdef DEBUGL
          mexPrintf("FLDVS\n");
# endif
          code += sizeof(FLDVS_);
          break;
        case Tags::FLDSV:
# ifdef DEBUGL
          mexPrintf("FLDSV\n");
# endif
          code += sizeof(FLDSV_);
          break;
        case Tags::FSTPSV:
# ifdef DEBUGL
          mexPrintf("FSTPSV\n");
# endif
          code += sizeof(FSTPSV_);
          break;
        case Tags::FLDV:
# ifdef DEBUGL
          mexPrintf("FLDV\n");
# endif
          code += sizeof(FLDV_);
          break;
        case Tags::FSTPV:
# ifdef DEBUGL
          mexPrintf("FSTPV\n");
# endif
          code += sizeof(FSTPV_);
          break;
        case Tags::FBEGINBLOCK:
# ifdef DEBUGL
          mexPrintf("FBEGINBLOCK\n");
# endif
          deserialized_special_instrs.push_back(make_unique<FBEGINBLOCK_>(code));
          begin_block.push_back(instructions_list.size());
          nb_blocks++;
          instr = deserialized_special_instrs.back().get();
          break;
        case Tags::FJMPIFEVAL:
# ifdef DEBUGL
          mexPrintf("FJMPIFEVAL\n");
# endif
          code += sizeof(FJMPIFEVAL_);
          break;
        case Tags::FJMP:
# ifdef DEBUGL
          mexPrintf("FJMP\n");
# endif
          code += sizeof(FJMP_);
          break;
        case Tags::FCALL:
# ifdef DEBUGL
          mexPrintf("FCALL\n");
# endif
          deserialized_special_instrs.push_back(make_unique<FCALL_>(code));
          instr = deserialized_special_instrs.back().get();
          break;
        case Tags::FLDTEF:
# ifdef DEBUGL
          mexPrintf("FLDTEF\n");
# endif
          code += sizeof(FLDTEF_);
          break;
        case Tags::FSTPTEF:
# ifdef DEBUGL
          mexPrintf("FSTPTEF\n");
# endif
          code += sizeof(FSTPTEF_);
          break;
        case Tags::FLDTEFD:
# ifdef DEBUGL
          mexPrintf("FLDTEFD\n");
# endif
          code += sizeof(FLDTEFD_);
          break;
        case Tags::FSTPTEFD:
# ifdef DEBUGL
          mexPrintf("FSTPTEFD\n");
# endif
          code += sizeof(FSTPTEFD_);
          break;
        case Tags::FLDTEFDD:
# ifdef DEBUGL
          mexPrintf("FLDTEFDD\n");
# endif
          code += sizeof(FLDTEFDD_);
          break;
        case Tags::FSTPTEFDD:
# ifdef DEBUGL
          mexPrintf("FSTPTEFDD\n");
# endif
          code += sizeof(FSTPTEFDD_);
          break;
        default:
          throw FatalException {"Unknown tag value=" + to_string(static_cast<int>(*reinterpret_cast<Tags *>(code)))};
        }
      instructions_list.push_back(instr);
    }
}

string
Evaluate::error_location(it_code_type expr_begin, it_code_type faulty_op, int it_) const
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
    case ExpressionType::FirstExoDerivative:
    case ExpressionType::FirstExodetDerivative:
      Error_loc << "first order derivative of equation";
      break;
    }
  Error_loc << " " << EQN_equation+1;
  if (nb_blocks > 1)
    Error_loc << " in block " << block_num+1;
  switch (EQN_type)
    {
    case ExpressionType::TemporaryTerm:
    case ExpressionType::ModelEquation:
      break;
    case ExpressionType::FirstEndoDerivative:
      Error_loc << " with respect to endogenous variable " << symbol_table.getName(SymbolType::endogenous, EQN_dvar1);
      break;
    case ExpressionType::FirstExoDerivative:
      Error_loc << " with respect to exogenous variable " << symbol_table.getName(SymbolType::exogenous, EQN_dvar1);
      break;
    case ExpressionType::FirstExodetDerivative:
      Error_loc << " with respect to deterministic exogenous variable " << symbol_table.getName(SymbolType::exogenousDet, EQN_dvar1);
      break;
    }
  if (!steady_state)
    Error_loc << " at time " << it_;

  /* Given a string which possibly contains a floating-point exception
     (materialized by an operator between braces), returns a string spanning two
     lines, the first line containing the original string without the braces,
     the second line containing tildes (~) under the faulty operator. */
  auto add_underscore_to_fpe = [](const string &str)
  {
    string line1;
    optional<size_t> pos1, pos2;
    string line2(str.length(), ' ');
    for (char i : str)
      {
        if (i != '{' && i != '}')
          line1 += i;
        else
          {
            if (i == '{')
              pos1 = line1.length();
            else
              pos2 = line1.length();
            if (pos1 && pos2)
              {
                line2.replace(*pos1, *pos2-*pos1, *pos2-*pos1, '~');
                pos1.reset();
                pos2.reset();
              }
          }
      }
    return line1 + "\n" + line2;
  };

  auto [expr_str, it_code_ret] = print_expression(expr_begin, faulty_op);
  Error_loc << endl << add_underscore_to_fpe("      " + expr_str);
  return Error_loc.str();
}

pair<string, Evaluate::it_code_type>
Evaluate::print_expression(const Evaluate::it_code_type &expr_begin, const optional<it_code_type> &faulty_op) const
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
        throw UserException{};
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
            case ExpressionType::FirstExoDerivative:
              equation_type = ExpressionType::FirstExoDerivative;
              break;
            case ExpressionType::FirstExodetDerivative:
              equation_type = ExpressionType::FirstExodetDerivative;
              break;
            default:
              throw FatalException{"In print_expression, expression type "
                                   + to_string(static_cast<int>(static_cast<FNUMEXPR_ *>(*it_code)->get_expression_type()))
                                   + " not implemented yet"};
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
                Stack.emplace(symbol_table.getName(type, var) + lag_to_string(lag), 100, nullopt);
                break;
              default:
                throw FatalException{"FLDV: Unknown variable type"};
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
                Stack.emplace(symbol_table.getName(type, var), 100, nullopt);
                break;
              default:
                throw FatalException{"FLDSV: Unknown variable type"};
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
                Stack.emplace(symbol_table.getName(type, var), 100, nullopt);
                break;
              default:
                throw FatalException{"FLDVS: Unknown variable type"};
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
          {
            // We don’t use std::to_string(), because it uses fixed formatting
            ostringstream s;
            s << defaultfloat << static_cast<FLDC_ *>(*it_code)->get_value();
            Stack.emplace(s.str(), 100, nullopt);
          }
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
                assign_lhs(symbol_table.getName(type, var) + lag_to_string(lag));
                break;
              default:
                throw FatalException{"FSTPV: Unknown variable type"};
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
                assign_lhs(symbol_table.getName(type, var));
                break;
              default:
                throw FatalException{"FSTPSV: Unknown variable type"};
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
                case ExpressionType::FirstExoDerivative:
                  return "jacob_exo";
                case ExpressionType::FirstExodetDerivative:
                  return "jacob_exo_det";
                default:
                  throw FatalException{"Unknown equation type " + to_string(static_cast<int>(equation_type))};
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
                  return it_code == faulty_op ? "{tan}" : "tan";
                case UnaryOpcode::acos:
                  return it_code == faulty_op ? "{acos}" : "acos";
                case UnaryOpcode::asin:
                  return it_code == faulty_op ? "{asin}" : "asin";
                case UnaryOpcode::atan:
                  return "atan";
                case UnaryOpcode::cosh:
                  return it_code == faulty_op ? "{cosh}" : "cosh";
                case UnaryOpcode::sinh:
                  return it_code == faulty_op ? "{sinh}" : "sinh";
                case UnaryOpcode::tanh:
                  return "tanh";
                case UnaryOpcode::acosh:
                  return it_code == faulty_op ? "{acosh}" : "acosh";
                case UnaryOpcode::asinh:
                  return "asinh";
                case UnaryOpcode::atanh:
                  return it_code == faulty_op ? "{atanh}" : "atanh";
                case UnaryOpcode::sqrt:
                  return it_code == faulty_op ? "{sqrt}" : "sqrt";
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
                  throw FatalException{"Unexpected derivative of steady_state operator"};
                case UnaryOpcode::expectation:
                  throw FatalException{"Unexpected expectation operator"};
                case UnaryOpcode::diff:
                  return "diff";
                case UnaryOpcode::adl:
                  return "adl";
                }
              throw FatalException{"Unknown opcode"};
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
                              + ", " + arg2 + ", " + arg3 + ")", 100, nullopt);
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
                  throw FatalException{"Unknown opcode"};
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
                    throw FatalException{"Should not arrive here"};
                  case BinaryOpcode::equal:
                    return " = ";
                  }
                  throw FatalException{"Unknown opcode"};
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
              throw FatalException{"Unknown opcode"};
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
              throw FatalException{"Unexpected external function call type"};
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
          throw FatalException{"In print_expression, unknown opcode "
                               + to_string(static_cast<int>((*it_code)->op_code))};
        }
      it_code++;
    }
  return { get<0>(Stack.top()), it_code };
}

void
Evaluate::evaluateBlock(int it_, double *__restrict__ y, const double *__restrict__ ya, int y_size, double *__restrict__ x, int nb_row_x, double *__restrict__ params, const double *__restrict__ steady_y, double *__restrict__ u, int Per_u_, double *__restrict__ T, int T_nrows, map<int, double> &TEF, map<pair<int, int>, double> &TEFD, map<tuple<int, int, int>, double> &TEFDD, double *__restrict__ r, double *__restrict__ g1, double *__restrict__ jacob, double *__restrict__ jacob_exo, double *__restrict__ jacob_exo_det, bool evaluate, bool no_derivatives)
{
  auto it_code { currentBlockBeginning() };
  int var{0}, lag{0};
  UnaryOpcode op1;
  BinaryOpcode op2;
  TrinaryOpcode op3;
  unsigned int eq, pos_col;
#ifdef DEBUG
  ostringstream tmp_out;
#endif
  double v1, v2, v3;
  bool go_on = true;
  double ll;
  double rr;
  stack<double> Stack;
  ExternalFunctionCallType call_type{ExternalFunctionCallType::levelWithoutDerivative};
  it_code_type it_code_expr;

#ifdef MATLAB_MEX_FILE
  if (utIsInterruptPending())
    throw UserException{};
#endif

  while (go_on)
    {
      switch ((*it_code)->op_code)
        {
        case Tags::FNUMEXPR:
#ifdef DEBUG
          mexPrintf("FNUMEXPR\n");
#endif
          it_code_expr = it_code;
          switch (static_cast<FNUMEXPR_ *>(*it_code)->get_expression_type())
            {
            case ExpressionType::TemporaryTerm:
#ifdef DEBUG
              mexPrintf("TemporaryTerm\n");
#endif
              EQN_type = ExpressionType::TemporaryTerm;
              EQN_equation = static_cast<FNUMEXPR_ *>(*it_code)->get_equation();
#ifdef DEBUG
              mexPrintf("EQN_equation=%d\n", EQN_equation); mexEvalString("drawnow;");
#endif
              break;
            case ExpressionType::ModelEquation:
#ifdef DEBUG
              mexPrintf("ModelEquation\n");
#endif
              EQN_type = ExpressionType::ModelEquation;
              EQN_equation = static_cast<FNUMEXPR_ *>(*it_code)->get_equation();
              break;
            case ExpressionType::FirstEndoDerivative:
#ifdef DEBUG
              mexPrintf("FirstEndoDerivative\n");
#endif
              EQN_type = ExpressionType::FirstEndoDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(*it_code)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(*it_code)->get_dvariable1();
              EQN_lag1 = static_cast<FNUMEXPR_ *>(*it_code)->get_lag1();
              break;
            case ExpressionType::FirstExoDerivative:
#ifdef DEBUG
              mexPrintf("FirstExoDerivative\n");
#endif
              EQN_type = ExpressionType::FirstExoDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(*it_code)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(*it_code)->get_dvariable1();
              EQN_lag1 = static_cast<FNUMEXPR_ *>(*it_code)->get_lag1();
              break;
            case ExpressionType::FirstExodetDerivative:
#ifdef DEBUG
              mexPrintf("FirstExodetDerivative\n");
#endif
              EQN_type = ExpressionType::FirstExodetDerivative;
              EQN_equation = static_cast<FNUMEXPR_ *>(*it_code)->get_equation();
              EQN_dvar1 = static_cast<FNUMEXPR_ *>(*it_code)->get_dvariable1();
              EQN_lag1 = static_cast<FNUMEXPR_ *>(*it_code)->get_lag1();
              break;
            }
          break;
        case Tags::FLDV:
          //load a variable in the processor
          switch (static_cast<FLDV_ *>(*it_code)->get_type())
            {
            case SymbolType::parameter:
              var = static_cast<FLDV_ *>(*it_code)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDV Param[var=%d]\n", var);
              tmp_out << " params[" << var << "](" << params[var] << ")";
#endif
              Stack.push(params[var]);
              break;
            case SymbolType::endogenous:
              var = static_cast<FLDV_ *>(*it_code)->get_pos();
              lag = static_cast<FLDV_ *>(*it_code)->get_lead_lag();
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
              var = static_cast<FLDV_ *>(*it_code)->get_pos();
              lag = static_cast<FLDV_ *>(*it_code)->get_lead_lag();
#ifdef DEBUG
              mexPrintf("FLDV x[var=%d, lag=%d, it_=%d], nb_row_x=%d evaluate=%d x[%d]=%f\n", var, lag, it_, nb_row_x, evaluate, it_+lag+var*nb_row_x, x[it_+lag+var*nb_row_x]);
              //tmp_out << " x[" << it_+lag << ", " << var << "](" << x[it_+lag+var*nb_row_x] << ")";
#endif
              Stack.push(x[it_+lag+var*nb_row_x]);
              break;
            case SymbolType::exogenousDet:
              throw FatalException{"FLDV: exogenous deterministic not supported"};
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
          switch (static_cast<FLDSV_ *>(*it_code)->get_type())
            {
            case SymbolType::parameter:
              var = static_cast<FLDSV_ *>(*it_code)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV Param[var=%d]=%f\n", var, params[var]);
              tmp_out << " params[" << var << "](" << params[var] << ")";
#endif
              Stack.push(params[var]);
              break;
            case SymbolType::endogenous:
              var = static_cast<FLDSV_ *>(*it_code)->get_pos();
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
              var = static_cast<FLDSV_ *>(*it_code)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV x[var=%d]\n", var);
              tmp_out << " x[" << var << "](" << x[var] << ")";
#endif
              Stack.push(x[var]);
              break;
            case SymbolType::exogenousDet:
              throw FatalException{"FLDSV: exogenous deterministic not supported"};
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
          switch (static_cast<FLDVS_ *>(*it_code)->get_type())
            {
            case SymbolType::parameter:
              var = static_cast<FLDVS_ *>(*it_code)->get_pos();
#ifdef DEBUG
              mexPrintf("params[%d]\n", var);
#endif
              Stack.push(params[var]);
              break;
            case SymbolType::endogenous:
              var = static_cast<FLDVS_ *>(*it_code)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDVS steady_y[%d]\n", var);
#endif
              Stack.push(steady_y[var]);
              break;
            case SymbolType::exogenous:
              var = static_cast<FLDVS_ *>(*it_code)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDVS x[%d] \n", var);
#endif
              Stack.push(x[var]);
              break;
            case SymbolType::exogenousDet:
              throw FatalException{"FLDVS: exogenous deterministic not supported"};
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
          var = static_cast<FLDT_ *>(*it_code)->get_pos();
#ifdef DEBUG
          mexPrintf("FLDT T[it_=%d var=%d, y_kmin=%d, y_kmax=%d == %d]=>%f\n", it_, var, y_kmin, y_kmax, var*T_nrows+it_, T[var*T_nrows+it_]);
          tmp_out << " T[" << it_ << ", " << var << "](" << T[var*T_nrows+it_] << ")";
#endif
          Stack.push(T[var*T_nrows+it_]);
          break;
        case Tags::FLDST:
          //load a temporary variable in the processor
          var = static_cast<FLDST_ *>(*it_code)->get_pos();
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
          var = static_cast<FLDU_ *>(*it_code)->get_pos();
          var += Per_u_;
#ifdef DEBUG
          mexPrintf("FLDU u[%d]\n", var);
          tmp_out << " u[" << var << "](" << u[var] << ")";
#endif
          Stack.push(u[var]);
          break;
        case Tags::FLDSU:
          //load u variable in the processor
          var = static_cast<FLDSU_ *>(*it_code)->get_pos();
#ifdef DEBUG
          mexPrintf("FLDSU u[%d]\n", var);
          tmp_out << " u[" << var << "](" << u[var] << ")";
#endif
          Stack.push(u[var]);
          break;
        case Tags::FLDR:
          //load u variable in the processor
          var = static_cast<FLDR_ *>(*it_code)->get_pos();
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
          ll = static_cast<FLDC_ *>(*it_code)->get_value();
#ifdef DEBUG
          mexPrintf("FLDC = %f\n", ll);
          tmp_out << " " << ll;
#endif

          Stack.push(ll);
          break;
        case Tags::FSTPV:
          //load a variable in the processor
          switch (static_cast<FSTPV_ *>(*it_code)->get_type())
            {
            case SymbolType::parameter:
              var = static_cast<FSTPV_ *>(*it_code)->get_pos();
#ifdef DEBUG
              mexPrintf("FSTPV params[%d]\n", var);
#endif
              params[var] = Stack.top();
              Stack.pop();
              break;
            case SymbolType::endogenous:
              var = static_cast<FSTPV_ *>(*it_code)->get_pos();
              lag = static_cast<FSTPV_ *>(*it_code)->get_lead_lag();
              y[(it_+lag)*y_size+var] = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" y[%d, %d](%f)=%s\n", it_+lag, var, y[(it_+lag)*y_size+var], tmp_out.str().c_str());
              tmp_out.str("");
#endif
              Stack.pop();
              break;
            case SymbolType::exogenous:
              var = static_cast<FSTPV_ *>(*it_code)->get_pos();
              lag = static_cast<FSTPV_ *>(*it_code)->get_lead_lag();
              x[it_+lag+var*nb_row_x] = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" x[%d, %d](%f)=%s\n", it_+lag, var, x[it_+lag+var*nb_row_x], tmp_out.str().c_str());
              tmp_out.str("");
#endif

              Stack.pop();
              break;
            case SymbolType::exogenousDet:
              throw FatalException{"FSTPV: exogenous deterministic not supported"};
              break;
            default:
              mexPrintf("FSTPV: Unknown variable type\n");
            }
          break;
        case Tags::FSTPSV:
          //load a variable in the processor
          switch (static_cast<FSTPSV_ *>(*it_code)->get_type())
            {
            case SymbolType::parameter:
              var = static_cast<FSTPSV_ *>(*it_code)->get_pos();
              params[var] = Stack.top();
              Stack.pop();
              break;
            case SymbolType::endogenous:
              var = static_cast<FSTPSV_ *>(*it_code)->get_pos();
              y[var] = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" y[%d](%f)=%s\n", var, y[var], tmp_out.str().c_str());
              tmp_out.str("");
#endif
              Stack.pop();
              break;
            case SymbolType::exogenous:
              var = static_cast<FSTPSV_ *>(*it_code)->get_pos();
              x[var] = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" x[%d, %d](%f)=%s\n", it_+lag, var, x[var], tmp_out.str().c_str());
              tmp_out.str("");
#endif
              Stack.pop();
              break;
            case SymbolType::exogenousDet:
              throw FatalException{"FSTPSV: exogenous deterministic not supported"};
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
          var = static_cast<FSTPT_ *>(*it_code)->get_pos();
          T[var*T_nrows+it_] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" T[%d, %d](%f)=%s\n", it_, var, T[var*T_nrows+it_], tmp_out.str().c_str());
          tmp_out.str("");
#endif

          Stack.pop();
          break;
        case Tags::FSTPST:
          //store in a temporary variable from the processor
#ifdef DEBUG
          mexPrintf("FSTPST\n");
#endif
          var = static_cast<FSTPST_ *>(*it_code)->get_pos();
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
          var = static_cast<FSTPU_ *>(*it_code)->get_pos();
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
          var = static_cast<FSTPSU_ *>(*it_code)->get_pos();
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
          var = static_cast<FSTPR_ *>(*it_code)->get_pos();
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
          var = static_cast<FSTPG_ *>(*it_code)->get_pos();
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
            throw FatalException{"In compute_block_time, impossible case " + to_string(static_cast<int>(EQN_type))
                                 + " not implement in static jacobian"};
          eq = static_cast<FSTPG2_ *>(*it_code)->get_row();
          var = static_cast<FSTPG2_ *>(*it_code)->get_col();
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
              eq = static_cast<FSTPG3_ *>(*it_code)->get_row();
              var = static_cast<FSTPG3_ *>(*it_code)->get_col();
              lag = static_cast<FSTPG3_ *>(*it_code)->get_lag();
              pos_col = static_cast<FSTPG3_ *>(*it_code)->get_col_pos();
#ifdef DEBUG
              mexPrintf("Endo eq=%d, pos_col=%d, size=%d, jacob=%x\n", eq, pos_col, size, jacob);
              mexPrintf("jacob=%x\n", jacob);
#endif
              jacob[eq + size*pos_col] = rr;
              break;
            case ExpressionType::FirstExoDerivative:
              //eq = static_cast<FSTPG3_ *>(*it_code)->get_row();
              eq = EQN_equation;
              var = static_cast<FSTPG3_ *>(*it_code)->get_col();
              lag = static_cast<FSTPG3_ *>(*it_code)->get_lag();
              pos_col = static_cast<FSTPG3_ *>(*it_code)->get_col_pos();
#ifdef DEBUG
              mexPrintf("Exo eq=%d, pos_col=%d, size=%d\n", eq, pos_col, size);
              mexEvalString("drawnow;");
#endif
              jacob_exo[eq + size*pos_col] = rr;
              break;
            case ExpressionType::FirstExodetDerivative:
              //eq = static_cast<FSTPG3_ *>(*it_code)->get_row();
              eq = EQN_equation;
              var = static_cast<FSTPG3_ *>(*it_code)->get_col();
              lag = static_cast<FSTPG3_ *>(*it_code)->get_lag();
              pos_col = static_cast<FSTPG3_ *>(*it_code)->get_col_pos();
#ifdef DEBUG
              mexPrintf("Exo det eq=%d, pos_col=%d, size=%d\n", eq, pos_col, size);
              mexEvalString("drawnow;");
#endif

              jacob_exo_det[eq + size*pos_col] = rr;
              break;
            default:
              throw FatalException{"In compute_block_time, variable " + to_string(static_cast<int>(EQN_type)) + " not used yet"};
            }
// #ifdef DEBUG
//           tmp_out << "=>";
//           mexPrintf(" g1[%d](%f)=%s\n", var, g1[var], tmp_out.str().c_str());
//           tmp_out.str("");
// #endif
          Stack.pop();
          break;

        case Tags::FBINARY:
          op2 = static_cast<FBINARY_ *>(*it_code)->get_op_type();
#ifdef DEBUG
          mexPrintf("FBINARY, op=%d\n", static_cast<int>(op2));
#endif
          v2 = Stack.top();
          Stack.pop();
          v1 = Stack.top();
          Stack.pop();
          switch (op2)
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
              {
                feclearexcept(FE_ALL_EXCEPT);
                double tmp {v1 / v2};
                if (fetestexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW))
                  throw DivideException{v1, v2, error_location(it_code_expr, it_code, it_)};
                Stack.push(tmp);
              }
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
              {
                feclearexcept(FE_ALL_EXCEPT);
                double tmp {pow(v1, v2)};
                if (fetestexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW))
                  throw PowException{v1, v2, error_location(it_code_expr, it_code, it_)};
                Stack.push(tmp);
              }
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
                    if (fabs(v1) < power_deriv_near_zero && v2 > 0
                        && derivOrder > v2
                        && fabs(v2-nearbyint(v2)) < power_deriv_near_zero)
                      Stack.push(0.0);
                    else
                      {
                        feclearexcept(FE_ALL_EXCEPT);
                        double dxp {pow(v1, v2-derivOrder)};
                        if (fetestexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW))
                          throw PowException{v1, v2-derivOrder, error_location(it_code_expr, it_code, it_)};
                        for (int i = 0; i < derivOrder; i++)
                          dxp *= v2--;
                        Stack.push(dxp);
                      }
                  }
                catch (FloatingPointException &fpeh)
                  {
                    mexPrintf("%s\n      %s\n", fpeh.message.c_str(), error_location(it_code_expr, it_code, it_).c_str());
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
            }
          break;
        case Tags::FUNARY:
          op1 = static_cast<FUNARY_ *>(*it_code)->get_op_type();
          v1 = Stack.top();
          Stack.pop();
#ifdef DEBUG
          mexPrintf("FUNARY, op=%d\n", static_cast<int>(op1));
#endif
          switch (op1)
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
              {
                feclearexcept(FE_ALL_EXCEPT);
                double tmp {log(v1)};
                if (fetestexcept(FE_DIVBYZERO | FE_INVALID))
                  throw UnaryOpException{"log", v1, error_location(it_code_expr, it_code, it_)};
                Stack.push(tmp);
              }
#ifdef DEBUG
              tmp_out << " |log(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::log10:
              {
                feclearexcept(FE_ALL_EXCEPT);
                double tmp {log10(v1)};
                if (fetestexcept(FE_DIVBYZERO | FE_INVALID))
                  throw UnaryOpException{"log10", v1, error_location(it_code_expr, it_code, it_)};
                Stack.push(tmp);
              }
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
              {
                feclearexcept(FE_ALL_EXCEPT);
                double tmp {tan(v1)};
                if (fetestexcept(FE_OVERFLOW))
                  throw UnaryOpException{"tan", v1, error_location(it_code_expr, it_code, it_)};
                Stack.push(tmp);
              }
#ifdef DEBUG
              tmp_out << " |tan(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::acos:
              {
                feclearexcept(FE_ALL_EXCEPT);
                double tmp {acos(v1)};
                if (fetestexcept(FE_INVALID))
                  throw UnaryOpException{"acos", v1, error_location(it_code_expr, it_code, it_)};
                Stack.push(tmp);
              }
#ifdef DEBUG
              tmp_out << " |acos(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::asin:
              {
                feclearexcept(FE_ALL_EXCEPT);
                double tmp {asin(v1)};
                if (fetestexcept(FE_INVALID))
                  throw UnaryOpException{"asin", v1, error_location(it_code_expr, it_code, it_)};
                Stack.push(tmp);
              }
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
              {
                feclearexcept(FE_ALL_EXCEPT);
                double tmp {cosh(v1)};
                if (fetestexcept(FE_OVERFLOW))
                  throw UnaryOpException{"cosh", v1, error_location(it_code_expr, it_code, it_)};
                Stack.push(tmp);
              }
#ifdef DEBUG
              tmp_out << " |cosh(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::sinh:
              {
                feclearexcept(FE_ALL_EXCEPT);
                double tmp {sinh(v1)};
                if (fetestexcept(FE_OVERFLOW))
                  throw UnaryOpException{"sinh", v1, error_location(it_code_expr, it_code, it_)};
                Stack.push(tmp);
              }
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
              {
                feclearexcept(FE_ALL_EXCEPT);
                double tmp {acosh(v1)};
                if (fetestexcept(FE_INVALID))
                  throw UnaryOpException{"acosh", v1, error_location(it_code_expr, it_code, it_)};
                Stack.push(tmp);
              }
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
              {
                feclearexcept(FE_ALL_EXCEPT);
                double tmp {atanh(v1)};
                if (fetestexcept(FE_INVALID | FE_DIVBYZERO))
                  throw UnaryOpException{"atanh", v1, error_location(it_code_expr, it_code, it_)};
                Stack.push(tmp);
              }
#ifdef DEBUG
              tmp_out << " |atanh(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::sqrt:
              {
                feclearexcept(FE_ALL_EXCEPT);
                double tmp {sqrt(v1)};
                if (fetestexcept(FE_INVALID))
                  throw UnaryOpException{"sqrt", v1, error_location(it_code_expr, it_code, it_)};
                Stack.push(tmp);
              }
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
            case UnaryOpcode::cbrt:
              Stack.push(cbrt(v1));
#ifdef DEBUG
              tmp_out << " |cbrt(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::abs:
              Stack.push(abs(v1));
#ifdef DEBUG
              tmp_out << " |abs(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::sign:
              Stack.push((v1 > 0) ? 1 : ((v1 < 0) ? -1 : 0));
#ifdef DEBUG
              tmp_out << " |sign(" << v1 << ")|";
#endif
              break;
            case UnaryOpcode::steadyState:
              throw FatalException {"Internal error: operator steady_state should not appear"};
            case UnaryOpcode::steadyStateParamDeriv:
              throw FatalException {"Internal error: 1st derivative w.r.t. parameters of operator steady_state should not appear"};
            case UnaryOpcode::steadyStateParam2ndDeriv:
              throw FatalException {"Internal error: 2nd derivative w.r.t. parameters of operator steady_state should not appear"};
            case UnaryOpcode::expectation:
              throw FatalException {"Internal error: operator expectation should not appear"};
            case UnaryOpcode::diff:
              throw FatalException {"Internal error: operator diff should not appear"};
            case UnaryOpcode::adl:
              throw FatalException {"Internal error: operator adl should not appear"};
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
              Stack.push(0.5*(1+erf((v1-v2)/v3/numbers::sqrt2)));
#ifdef DEBUG
              tmp_out << " |normcdf(" << v1 << ", " << v2 << ", " << v3 << ")|";
#endif
              break;
            case TrinaryOpcode::normpdf:
              Stack.push(1/(v3*sqrt(2*numbers::pi)*exp(pow((v1-v2)/v3, 2)/2)));
#ifdef DEBUG
              tmp_out << " |normpdf(" << v1 << ", " << v2 << ", " << v3 << ")|";
#endif
              break;
            }
          break;

        case Tags::FCALL:
          {
#ifdef DEBUG
            mexPrintf("------------------------------\n");
            mexPrintf("CALL "); mexEvalString("drawnow;");
#endif
            auto *fc = static_cast<FCALL_ *>(*it_code);
            string function_name = fc->get_function_name();
#ifdef DEBUG
            mexPrintf("function_name=%s ", function_name.c_str()); mexEvalString("drawnow;");
#endif
            int nb_input_arguments{fc->get_nb_input_arguments()};
#ifdef DEBUG
            mexPrintf("nb_input_arguments=%d ", nb_input_arguments); mexEvalString("drawnow;");
#endif
            int nb_output_arguments{fc->get_nb_output_arguments()};
#ifdef DEBUG
            mexPrintf("nb_output_arguments=%d\n", nb_output_arguments); mexEvalString("drawnow;");
#endif

            mxArray *output_arguments[3];
            string arg_func_name = fc->get_arg_func_name();
#ifdef DEBUG
            mexPrintf("arg_func_name.length() = %d\n", arg_func_name.length());
            mexPrintf("arg_func_name.c_str() = %s\n", arg_func_name.c_str());
#endif
            int nb_add_input_arguments{fc->get_nb_add_input_arguments()};
            call_type = fc->get_call_type();
#ifdef DEBUG
            mexPrintf("call_type=%d ExternalFunctionCallTypeWithoutDerivative=%d\n", call_type, ExternalFunctionCallType::levelWithoutDerivative);
            mexEvalString("drawnow;");
#endif
            mxArray **input_arguments;
            switch (call_type)
              {
              case ExternalFunctionCallType::levelWithoutDerivative:
              case ExternalFunctionCallType::levelWithFirstDerivative:
              case ExternalFunctionCallType::levelWithFirstAndSecondDerivative:
                {
                  input_arguments = static_cast<mxArray **>(mxMalloc(nb_input_arguments * sizeof(mxArray *)));
                  test_mxMalloc(input_arguments, __LINE__, __FILE__, __func__, nb_input_arguments * sizeof(mxArray *));
#ifdef DEBUG
                  mexPrintf("Stack.size()=%d\n", Stack.size());
                  mexEvalString("drawnow;");
#endif
                  for (int i{0}; i < nb_input_arguments; i++)
                    {
                      mxArray *vv = mxCreateDoubleScalar(Stack.top());
                      input_arguments[nb_input_arguments - i - 1] = vv;
                      Stack.pop();
                    }
                  if (mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str()))
                    throw FatalException{"External function: " + function_name + " not found"};

                  double *rr = mxGetPr(output_arguments[0]);
                  Stack.push(*rr);
                  if (call_type == ExternalFunctionCallType::levelWithFirstDerivative
                      || call_type == ExternalFunctionCallType::levelWithFirstAndSecondDerivative)
                    {
                      int indx{fc->get_indx()};
                      double *FD1 = mxGetPr(output_arguments[1]);
                      size_t rows = mxGetN(output_arguments[1]);
                      for (int i{0}; i < static_cast<int>(rows); i++)
                        TEFD[{ indx, i }] = FD1[i];
                    }
                  if (call_type == ExternalFunctionCallType::levelWithFirstAndSecondDerivative)
                    {
                      int indx{fc->get_indx()};
                      double *FD2 = mxGetPr(output_arguments[2]);
                      size_t rows = mxGetM(output_arguments[2]);
                      size_t cols = mxGetN(output_arguments[2]);
                      int k{0};
                      for (int j{0}; j < static_cast<int>(cols); j++)
                        for (int i{0}; i < static_cast<int>(rows); i++)
                          TEFDD[{ indx, i, j }] = FD2[k++];
                    }
                }
                break;
              case ExternalFunctionCallType::numericalFirstDerivative:
                {
                  input_arguments = static_cast<mxArray **>(mxMalloc((nb_input_arguments+1+nb_add_input_arguments) * sizeof(mxArray *)));
                  test_mxMalloc(input_arguments, __LINE__, __FILE__, __func__, (nb_input_arguments+1+nb_add_input_arguments) * sizeof(mxArray *));
                  mxArray *vv = mxCreateString(arg_func_name.c_str());
                  input_arguments[0] = vv;
                  vv = mxCreateDoubleScalar(fc->get_row());
                  input_arguments[1] = vv;
                  vv = mxCreateCellMatrix(1, nb_add_input_arguments);
                  for (int i = 0; i < nb_add_input_arguments; i++)
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
                    throw FatalException{"External function: " + function_name + " not found"};
                  double *rr = mxGetPr(output_arguments[0]);
#ifdef DEBUG
                  mexPrintf("*rr=%f\n", *rr);
#endif
                  Stack.push(*rr);
                }
                break;
              case ExternalFunctionCallType::separatelyProvidedFirstDerivative:
                {
                  input_arguments = static_cast<mxArray **>(mxMalloc(nb_input_arguments * sizeof(mxArray *)));
                  test_mxMalloc(input_arguments, __LINE__, __FILE__, __func__, nb_input_arguments * sizeof(mxArray *));
                  for (int i{0}; i < nb_input_arguments; i++)
                    {
                      mxArray *vv = mxCreateDoubleScalar(Stack.top());
                      input_arguments[(nb_input_arguments - 1) - i] = vv;
                      Stack.pop();
                    }
                  if (mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str()))
                    throw FatalException{"External function: " + function_name + " not found"};
                  int indx{fc->get_indx()};
                  double *FD1 = mxGetPr(output_arguments[0]);
                  size_t rows = mxGetN(output_arguments[0]);
                  for (int i{0}; i < static_cast<int>(rows); i++)
                    TEFD[{ indx, i }] = FD1[i];
                }
                break;
              case ExternalFunctionCallType::numericalSecondDerivative:
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
                  for (int i{0}; i < nb_add_input_arguments; i++)
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
                    throw FatalException{"External function: " + function_name + " not found"};
                  double *rr = mxGetPr(output_arguments[0]);
                  Stack.push(*rr);
                }
                break;
              case ExternalFunctionCallType::separatelyProvidedSecondDerivative:
                {
                  input_arguments = static_cast<mxArray **>(mxMalloc(nb_input_arguments * sizeof(mxArray *)));
                  test_mxMalloc(input_arguments, __LINE__, __FILE__, __func__, nb_input_arguments * sizeof(mxArray *));
                  for (int i{0}; i < nb_input_arguments; i++)
                    {
                      mxArray *vv = mxCreateDoubleScalar(Stack.top());
                      input_arguments[i] = vv;
                      Stack.pop();
                    }
                  if (mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str()))
                    throw FatalException{"External function: " + function_name + " not found"};
                  int indx{fc->get_indx()};
                  double *FD2 = mxGetPr(output_arguments[2]);
                  size_t rows = mxGetM(output_arguments[0]);
                  size_t cols = mxGetN(output_arguments[0]);
                  int k{0};
                  for (int j{0}; j < static_cast<int>(cols); j++)
                    for (int i{0}; i < static_cast<int>(rows); i++)
                      TEFDD[{ indx, i, j }] = FD2[k++];
                }
                break;
              }
          }
          break;
        case Tags::FSTPTEF:
          var = static_cast<FSTPTEF_ *>(*it_code)->get_number();
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
          var = static_cast<FLDTEF_ *>(*it_code)->get_number();
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
            int indx{static_cast<FSTPTEFD_ *>(*it_code)->get_indx()};
            int row{static_cast<FSTPTEFD_ *>(*it_code)->get_row()};
#ifdef DEBUG
            mexPrintf("FSTPTEFD\n");
            mexPrintf("indx=%d Stack.size()=%d\n", indx, Stack.size());
#endif
            if (call_type == ExternalFunctionCallType::numericalFirstDerivative)
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
            int indx{static_cast<FLDTEFD_ *>(*it_code)->get_indx()};
            int row{static_cast<FLDTEFD_ *>(*it_code)->get_row()};
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
            int indx{static_cast<FSTPTEFDD_ *>(*it_code)->get_indx()};
            int row{static_cast<FSTPTEFDD_ *>(*it_code)->get_row()};
            int col{static_cast<FSTPTEFDD_ *>(*it_code)->get_col()};
#ifdef DEBUG
            mexPrintf("FSTPTEFD\n");
            mexPrintf("indx=%d Stack.size()=%d\n", indx, Stack.size());
#endif
            if (call_type == ExternalFunctionCallType::numericalSecondDerivative)
              {
                TEFDD[{ indx, row-1, col-1 }] = Stack.top();
#ifdef DEBUG
                mexPrintf("FSTP TEFDD[{ indx, row-1, col-1 }]=%f done\n", TEFDD[{ indx, row-1, col-1 }]);
                mexEvalString("drawnow;");
#endif
                Stack.pop();
              }
          }

          break;
        case Tags::FLDTEFDD:
          {
            int indx{static_cast<FLDTEFDD_ *>(*it_code)->get_indx()};
            int row{static_cast<FLDTEFDD_ *>(*it_code)->get_row()};
            int col{static_cast<FSTPTEFDD_ *>(*it_code)->get_col()};
#ifdef DEBUG
            mexPrintf("FLDTEFD\n");
            mexPrintf("indx=%d Stack.size()=%d\n", indx, Stack.size());
            mexPrintf("FLD TEFD[{ indx, row-1, col-1 }]=%f done\n", TEFDD[{ indx, row-1, col-1 }]);
            mexEvalString("drawnow;");
#endif
            Stack.push(TEFDD[{ indx, row-1, col-1 }]);
          }
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
          if (no_derivatives)
            go_on = false;
          break;
        case Tags::FJMPIFEVAL:
          if (evaluate)
            {
#ifdef DEBUG
              mexPrintf("FJMPIFEVAL length=%d\n", static_cast<FJMPIFEVAL_ *>(*it_code)->get_pos());
              mexEvalString("drawnow;");
#endif
              it_code += static_cast<FJMPIFEVAL_ *>(*it_code)->get_pos() /* - 1*/;
            }
          break;
        case Tags::FJMP:
#ifdef DEBUG
          mexPrintf("FJMP length=%d\n", static_cast<FJMP_ *>(*it_code)->get_pos());
          mexEvalString("drawnow;");
#endif
          it_code += static_cast<FJMP_ *>(*it_code)->get_pos() /*- 1 */;
          break;
        default:
          throw FatalException{"In compute_block_time, unknown opcode " + to_string(static_cast<int>((*it_code)->op_code))};
        }
      it_code++;
    }
#ifdef DEBUG
  mexPrintf("==> end of compute_block_time Block = %d\n", block_num);
  mexEvalString("drawnow;");
#endif
}

void
Evaluate::gotoBlock(int block)
{
  block_num = block;

  auto *fb {currentBlockTag()};
  if (fb->op_code != Tags::FBEGINBLOCK)
    throw FatalException {"Evaluate::gotoBlock: internal inconsistency"};

  size = fb->get_size();
}

void
Evaluate::printCurrentBlock()
{
  auto it_code { currentBlockBeginning() };
  mexPrintf("\nBlock %d\n", block_num+1);
  mexPrintf("----------\n");
  bool go_on {true};
  bool space {false};
  while (go_on)
    {
      if ((*it_code)->op_code == Tags::FENDBLOCK)
        go_on = false;
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

int
Evaluate::getNumberOfTemporaryTerms() const
{
  BytecodeInstruction *instr {instructions_list.front()};
  if (steady_state)
    {
      if (instr->op_code == Tags::FDIMST)
        return reinterpret_cast<FDIMST_ *>(instr)->get_size();
      else
        throw FatalException {"Evaluate::getNumberOfTemporaryTerms: static .cod file does not begin with FDIMST!"};
    }
  else
    {
      if (instr->op_code == Tags::FDIMT)
        return reinterpret_cast<FDIMT_ *>(instr)->get_size();
      else
        throw FatalException {"Evaluate::getNumberOfTemporaryTerms: dynamic .cod file does not begin with FDIMT!"};
    }
}
