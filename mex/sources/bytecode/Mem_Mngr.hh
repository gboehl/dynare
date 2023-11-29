/*
 * Copyright © 2007-2021 Dynare Team
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

#ifndef _MEM_MNGR_HH
#define _MEM_MNGR_HH

#include <fstream>
#include <string>
#include <vector>

#include "ErrorHandling.hh"

using namespace std;

struct NonZeroElem
{
  int u_index;
  int r_index, c_index, lag_index;
  NonZeroElem *NZE_R_N, *NZE_C_N;
};

using v_NonZeroElem = vector<NonZeroElem*>;

class Mem_Mngr
{
public:
  void init_Mem();
  void mxFree_NZE(void* pos);
  NonZeroElem* mxMalloc_NZE();
  void init_CHUNK_BLCK_SIZE(int u_count);
  void Free_All();
  Mem_Mngr();
  void fixe_file_name(string filename_arg);
  bool swp_f;

private:
  v_NonZeroElem Chunk_Stack;
  unsigned int CHUNK_SIZE, CHUNK_BLCK_SIZE, Nb_CHUNK;
  unsigned int CHUNK_heap_pos;
  NonZeroElem** NZE_Mem_add;
  NonZeroElem* NZE_Mem;
  vector<NonZeroElem*> NZE_Mem_Allocated;
  int swp_f_b;
  fstream SaveCode_swp;
  string filename_mem;
};

#endif // _MEM_MNGR_HH
