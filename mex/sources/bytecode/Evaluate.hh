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

#ifndef _EVALUATE_HH
#define _EVALUATE_HH

#include <vector>
#include <string>
#include <map>
#include <optional>

#include "Bytecode.hh"
#include "ErrorHandling.hh"
#include "BasicSymbolTable.hh"

#include <dynmex.h>

using instructions_list_t = vector<BytecodeInstruction *>;
using it_code_type = instructions_list_t::const_iterator;

class CodeLoad
{
private:
  char *code;
  int nb_blocks;
  vector<size_t> begin_block;
public:
  int
  get_block_number() const
  {
    return nb_blocks;
  }

  size_t
  get_begin_block(int block) const
  {
    return begin_block[block];
  }
  instructions_list_t
  get_op_code(const filesystem::path &codfile)
  {
    instructions_list_t tags_liste;
    ifstream CompiledCode;
    streamoff Code_Size;
    CompiledCode.open(codfile, ios::in | ios::binary| ios::ate);
    if (!CompiledCode.is_open())
      return tags_liste;
    Code_Size = CompiledCode.tellg();
    CompiledCode.seekg(ios::beg);
    code = static_cast<char *>(mxMalloc(Code_Size));
    CompiledCode.seekg(0);
    CompiledCode.read(reinterpret_cast<char *>(code), Code_Size);
    CompiledCode.close();
    nb_blocks = 0;
    bool done = false;
    int instruction = 0;
    while (!done)
      {
        BytecodeInstruction *instr {reinterpret_cast<BytecodeInstruction *>(code)};
        switch (*reinterpret_cast<Tags *>(code))
          {
          case Tags::FLDZ:
# ifdef DEBUGL
            mexPrintf("FLDZ = %d size = %d\n", Tags::FLDZ, sizeof(FLDZ_));
# endif
            tags_liste.push_back(instr);
            code += sizeof(FLDZ_);
            break;
          case Tags::FEND:
# ifdef DEBUGL
            mexPrintf("FEND\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FEND_);
            done = true;
            break;
          case Tags::FENDBLOCK:
# ifdef DEBUGL
            mexPrintf("FENDBLOCK\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FENDBLOCK_);
            break;
          case Tags::FENDEQU:
# ifdef DEBUGL
            mexPrintf("FENDEQU\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FENDEQU_);
            break;
          case Tags::FDIMT:
# ifdef DEBUGL
            mexPrintf("FDIMT = %d size = %d\n", Tags::FDIMT, sizeof(FDIMT_));
# endif
            tags_liste.push_back(instr);
            code += sizeof(FDIMT_);
            break;
          case Tags::FDIMST:
# ifdef DEBUGL
            mexPrintf("FDIMST\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FDIMST_);
            break;
          case Tags::FNUMEXPR:
# ifdef DEBUGL
            mexPrintf("FNUMEXPR\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FNUMEXPR_);
            break;
          case Tags::FLDC:
# ifdef DEBUGL
            mexPrintf("FLDC\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FLDC_);
            break;
          case Tags::FLDU:
# ifdef DEBUGL
            mexPrintf("FLDU\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FLDU_);
            break;
          case Tags::FLDSU:
# ifdef DEBUGL
            mexPrintf("FLDSU\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FLDSU_);
            break;
          case Tags::FLDR:
# ifdef DEBUGL
            mexPrintf("FLDR\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FLDR_);
            break;
          case Tags::FLDT:
# ifdef DEBUGL
            mexPrintf("FLDT\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FLDT_);
            break;
          case Tags::FLDST:
# ifdef DEBUGL
            mexPrintf("FLDST\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FLDST_);
            break;
          case Tags::FSTPT:
# ifdef DEBUGL
            mexPrintf("FSTPT = %d size = %d\n", Tags::FSTPT, sizeof(FSTPT_));
# endif
            tags_liste.push_back(instr);
            code += sizeof(FSTPT_);
            break;
          case Tags::FSTPST:
# ifdef DEBUGL
            mexPrintf("FSTPST\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FSTPST_);
            break;
          case Tags::FSTPR:
# ifdef DEBUGL
            mexPrintf("FSTPR\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FSTPR_);
            break;
          case Tags::FSTPU:
# ifdef DEBUGL
            mexPrintf("FSTPU\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FSTPU_);
            break;
          case Tags::FSTPSU:
# ifdef DEBUGL
            mexPrintf("FSTPSU\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FSTPSU_);
            break;
          case Tags::FSTPG:
# ifdef DEBUGL
            mexPrintf("FSTPG\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FSTPG_);
            break;
          case Tags::FSTPG2:
# ifdef DEBUGL
            mexPrintf("FSTPG2\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FSTPG2_);
            break;
          case Tags::FSTPG3:
# ifdef DEBUGL
            mexPrintf("FSTPG3\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FSTPG3_);
            break;
          case Tags::FUNARY:
# ifdef DEBUGL
            mexPrintf("FUNARY\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FUNARY_);
            break;
          case Tags::FBINARY:
# ifdef DEBUGL
            mexPrintf("FBINARY\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FBINARY_);
            break;
          case Tags::FTRINARY:
# ifdef DEBUGL
            mexPrintf("FTRINARY\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FTRINARY_);
            break;
          case Tags::FLDVS:
# ifdef DEBUGL
            mexPrintf("FLDVS\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FLDVS_);
            break;
          case Tags::FLDSV:
# ifdef DEBUGL
            mexPrintf("FLDSV\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FLDSV_);
            break;
          case Tags::FSTPSV:
# ifdef DEBUGL
            mexPrintf("FSTPSV\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FSTPSV_);
            break;
          case Tags::FLDV:
# ifdef DEBUGL
            mexPrintf("FLDV\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FLDV_);
            break;
          case Tags::FSTPV:
# ifdef DEBUGL
            mexPrintf("FSTPV\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FSTPV_);
            break;
          case Tags::FBEGINBLOCK:
# ifdef DEBUGL
            mexPrintf("FBEGINBLOCK\n");
# endif
            {
              auto *fbegin_block = new FBEGINBLOCK_{code};
              begin_block.push_back(tags_liste.size());
              tags_liste.push_back(fbegin_block);
              nb_blocks++;
            }
            break;
          case Tags::FJMPIFEVAL:
# ifdef DEBUGL
            mexPrintf("FJMPIFEVAL\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FJMPIFEVAL_);
            break;
          case Tags::FJMP:
# ifdef DEBUGL
            mexPrintf("FJMP\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FJMP_);
            break;
          case Tags::FCALL:
            {
# ifdef DEBUGL
              mexPrintf("FCALL\n");
# endif
              auto *fcall = new FCALL_{code};
              tags_liste.push_back(fcall);
# ifdef DEBUGL
              mexPrintf("FCALL finish\n"); mexEvalString("drawnow;");
              mexPrintf("-- *code=%d\n", *code); mexEvalString("drawnow;");
# endif
            }
            break;
          case Tags::FLDTEF:
# ifdef DEBUGL
            mexPrintf("FLDTEF\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FLDTEF_);
            break;
          case Tags::FSTPTEF:
# ifdef DEBUGL
            mexPrintf("FSTPTEF\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FSTPTEF_);
            break;
          case Tags::FLDTEFD:
# ifdef DEBUGL
            mexPrintf("FLDTEFD\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FLDTEFD_);
            break;
          case Tags::FSTPTEFD:
# ifdef DEBUGL
            mexPrintf("FSTPTEFD\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FSTPTEFD_);
            break;
          case Tags::FLDTEFDD:
# ifdef DEBUGL
            mexPrintf("FLDTEFDD\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FLDTEFDD_);
            break;
          case Tags::FSTPTEFDD:
# ifdef DEBUGL
            mexPrintf("FSTPTEFDD\n");
# endif
            tags_liste.push_back(instr);
            code += sizeof(FSTPTEFDD_);
            break;
          default:
            mexPrintf("Unknown Tag value=%d code=%x\n", *code, code);
            done = true;
          }
        instruction++;
      }
    return tags_liste;
  }
};

class Evaluate
{
private:
  ExpressionType EQN_type;
  int EQN_equation, EQN_block, EQN_dvar1;
  int EQN_lag1, EQN_lag2, EQN_lag3;
  map<int, double> TEF;
  map<pair<int, int>, double> TEFD;
  map<tuple<int, int, int>, double> TEFDD;

  string error_location(it_code_type expr_begin, it_code_type faulty_op, bool steady_state, int it_) const;

protected:
  BasicSymbolTable &symbol_table;
  int EQN_block_number;
  double *y, *ya;
  int y_size;
  double *T;
  int nb_row_x;
  int y_kmin, y_kmax, periods;
  double *x, *params;
  double *u;
  double *steady_y;
  double *g1, *r, *res;
  vector<mxArray *> jacobian_block, jacobian_exo_block, jacobian_det_exo_block;
  mxArray *GlobalTemporaryTerms;
  it_code_type start_code, end_code;
  double pow1(double a, double b);
  double divide(double a, double b);
  double log1(double a);
  double log10_1(double a);
  void evaluate_over_periods(bool forward);
  void solve_simple_one_periods();
  void solve_simple_over_periods(bool forward);
  void compute_block_time(int Per_u_, bool evaluate, bool no_derivatives);
  instructions_list_t code_liste;
  it_code_type it_code;
  int Per_u_, Per_y_;
  int it_;
  int maxit_;
  double *direction;
  double solve_tolf;
  bool print_error;
  double res1, res2, max_res;
  int max_res_idx;
  vector<Block_contain_type> Block_Contain;

  int size;
  int *index_vara;

  BlockSimulationType type;
  int block_num, symbol_table_endo_nbr, u_count_int, block;
  string file_name, bin_base_name;
  bool is_linear;

  bool steady_state;

  /* Prints a bytecode expression in human readable form.
     If faulty_op is not default constructed, it should point to a tag withing
     the expression that created a floating point exception, in which case the
     corresponding mathematical operator will be printed within braces.
     The second output argument points to the tag past the expression. */
  pair<string, it_code_type> print_expression(const it_code_type &expr_begin, const optional<it_code_type> &faulty_op = nullopt) const;

public:
  Evaluate(int y_size_arg, int y_kmin_arg, int y_kmax_arg, bool steady_state_arg, int periods_arg, BasicSymbolTable &symbol_table_arg);
  void set_block(int size_arg, BlockSimulationType type_arg, string file_name_arg, string bin_base_name_arg, int block_num_arg,
                 bool is_linear_arg, int symbol_table_endo_nbr_arg, int u_count_int_arg, int block_arg);
  void evaluate_complete(bool no_derivatives);
  bool compute_complete(bool no_derivatives, double &res1, double &res2, double &max_res, int &max_res_idx);
  void compute_complete_2b(bool no_derivatives, double *_res1, double *_res2, double *_max_res, int *_max_res_idx);

  bool compute_complete(double lambda, double *crit);
};

#endif // _EVALUATE_HH
