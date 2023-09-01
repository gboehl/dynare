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
#include <deque>

#include "Bytecode.hh"
#include "BasicSymbolTable.hh"

class Evaluate
{
private:
  using instructions_list_t = vector<BytecodeInstruction *>;
  using it_code_type = instructions_list_t::const_iterator;

  const BasicSymbolTable &symbol_table;
  const bool steady_state; // Whether this is a static or dynamic .cod file

  // Memory copy of the contents of the .cod file
  unique_ptr<char[]> raw_bytecode;

  /* Owns read instructions that have their specialized deserializing
     constructors (and are thus not part of the “code” memory block). We use
     std::deque for storing them, because that class guarantees the stability
     of iterators, and thus of pointers to elements; we store such pointers in
     the “instructions_list” data member. */
  deque<FBEGINBLOCK_> deserialized_fbeginblock;
  deque<FCALL_> deserialized_fcall;

  /* List of deserialized instructions
     Those are either pointers inside “raw_bytecode” or “deserialized_{fbeginblock,fcall}” */
  instructions_list_t instructions_list;

   // Number of blocks in the model
  int nb_blocks {0};

  // Index of beginnings of blocks within instructions_list
  vector<size_t> begin_block;

  int block_num; // Index of the current block
  int size; // Size of the current block

  ExpressionType EQN_type;
  int EQN_equation, EQN_dvar1;
  int EQN_lag1, EQN_lag2, EQN_lag3;

  string error_location(it_code_type expr_begin, it_code_type faulty_op, int it_) const;

  /* Prints a bytecode expression in human readable form.
     If faulty_op is not default constructed, it should point to a tag within
     the expression that created a floating point exception, in which case the
     corresponding mathematical operator will be printed within braces.
     The second output argument points to the tag past the expression. */
  pair<string, it_code_type> print_expression(const it_code_type &expr_begin, const optional<it_code_type> &faulty_op = nullopt) const;

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

public:
  Evaluate(const filesystem::path &codfile, bool steady_state_arg, const BasicSymbolTable &symbol_table_arg);

  void evaluateBlock(int it_, double *__restrict__ y, const double *__restrict__ ya, int y_size, double *__restrict__ x, int nb_row_x, double *__restrict__ params, const double *__restrict__ steady_y, double *__restrict__ u, int Per_u_, double *__restrict__ T, int T_nrows, map<int, double> &TEF, map<pair<int, int>, double> &TEFD, map<tuple<int, int, int>, double> &TEFDD, double *__restrict__ r, double *__restrict__ g1, double *__restrict__ jacob, double *__restrict__ jacob_exo, double *__restrict__ jacob_exo_det, bool evaluate, bool no_derivatives);

  // Prints current block
  void printCurrentBlock();

  void gotoBlock(int block);

  int getNumberOfTemporaryTerms() const;

  auto
  getCurrentBlockSize() const
  {
    return currentBlockTag()->get_size();
  }
  auto
  getCurrentBlockType() const
  {
    return currentBlockTag()->get_type();
  }
  auto
  isCurrentBlockLinear() const
  {
    return currentBlockTag()->get_is_linear();
  }
  auto
  getCurrentBlockEquationsAndVariables() const
  {
    return currentBlockTag()->get_Block_Contain();
  }
  auto
  getCurrentBlockUCount() const
  {
    return currentBlockTag()->get_u_count_int();
  }
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

  int
  get_block_number() const
  {
    return nb_blocks;
  }
};

#endif // _EVALUATE_HH
