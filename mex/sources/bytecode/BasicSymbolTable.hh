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

#ifndef BASIC_SYMBOL_TABLE_HH
#define BASIC_SYMBOL_TABLE_HH

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "dynmex.h"

#include "Bytecode.hh"

using namespace std;

class BasicSymbolTable
{
public:
  BasicSymbolTable(const mxArray* M_);
  [[nodiscard]] string getName(SymbolType type, int tsid) const;
  [[nodiscard]] pair<SymbolType, int> getIDAndType(const string& name) const;
  [[nodiscard]] size_t maxEndoNameLength() const;

private:
  vector<string> endo_names, param_names, exo_names, exo_det_names;
  map<string, pair<SymbolType, int>> name_to_id_and_type;
};

#endif
