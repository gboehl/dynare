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
#include <memory>
#include <filesystem>

#include "Bytecode.hh"
#include "BasicSymbolTable.hh"

class Evaluate
{
private:
  using instructions_list_t = vector<BytecodeInstruction *>;
  using it_code_type = instructions_list_t::const_iterator;

  // Memory copy of the contents of the .cod file
  unique_ptr<char[]> raw_bytecode;

  /* Owns read instructions that have their specialized deserializing
     constructors (and are thus not part of the “code” memory block) */
  vector<unique_ptr<BytecodeInstruction>> deserialized_special_instrs;

  /* List of deserialized instructions
     Those are either pointers inside “raw_bytecode” or “deserialized_special_instrs” */
  instructions_list_t instructions_list;

   // Number of blocks in the model
  int nb_blocks {0};

  // Index of beginnings of blocks within instructions_list
  vector<size_t> begin_block;

  ExpressionType EQN_type;
  int EQN_equation, EQN_block, EQN_dvar1;
  int EQN_lag1, EQN_lag2, EQN_lag3;
  map<int, double> TEF;
  map<pair<int, int>, double> TEFD;
  map<tuple<int, int, int>, double> TEFDD;

  string error_location(it_code_type expr_begin, it_code_type faulty_op, bool steady_state, int it_) const;

  FBEGINBLOCK_ *
  currentBlockTag() const
  {
    return reinterpret_cast<FBEGINBLOCK_ *>(instructions_list[begin_block[block_num]]);
  }

  // Returns iterator to first instruction in the current block (after FBEGINBLOCK)
  it_code_type
  currentBlockBeginning()
  {
    return instructions_list.begin() + begin_block[block_num] + 1;
  }

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
  void evaluateBlock(int Per_u_, bool evaluate, bool no_derivatives);
  int it_;

  int block_num; // Index of the current block
  int size; // Size of the current block
  BlockSimulationType type;
  bool is_linear;
  int symbol_table_endo_nbr, u_count_int;
  vector<Block_contain_type> Block_Contain;

  bool steady_state;

  /* Prints a bytecode expression in human readable form.
     If faulty_op is not default constructed, it should point to a tag within
     the expression that created a floating point exception, in which case the
     corresponding mathematical operator will be printed within braces.
     The second output argument points to the tag past the expression. */
  pair<string, it_code_type> print_expression(const it_code_type &expr_begin, const optional<it_code_type> &faulty_op = nullopt) const;

  // Prints current block
  void printCurrentBlock();

  void gotoBlock(int block);

  void initializeTemporaryTerms(bool global_temporary_terms);

  auto
  getCurrentBlockExogenous() const
  {
    return currentBlockTag()->get_exogenous();
  }
  auto
  getCurrentBlockEndogenous() const
  {
    return currentBlockTag()->get_endogenous();
  }
  auto
  getCurrentBlockNbColJacob() const
  {
    return currentBlockTag()->get_nb_col_jacob();
  }
  auto
  getCurrentBlockExoSize() const
  {
    return currentBlockTag()->get_exo_size();
  }
  auto
  getCurrentBlockExoDetSize() const
  {
    return currentBlockTag()->get_det_exo_size();
  }

public:
  Evaluate(int y_size_arg, int y_kmin_arg, int y_kmax_arg, bool steady_state_arg, int periods_arg, BasicSymbolTable &symbol_table_arg);
  // TODO: integrate into the constructor
  void loadCodeFile(const filesystem::path &codfile);

  int
  get_block_number() const
  {
    return nb_blocks;
  }
};

#endif // _EVALUATE_HH
