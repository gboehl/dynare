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

#ifndef _BASIC_SYMBOL_TABLE_HH
#define _BASIC_SYMBOL_TABLE_HH

#include <string>
#include <vector>
#include <map>
#include <utility>

#include "Bytecode.hh"

using namespace std;

class BasicSymbolTable
{
public:
  BasicSymbolTable();
  string getName(SymbolType type, int tsid) const;
  pair<SymbolType, int> getIDAndType(const string &name) const;
  size_t maxEndoNameLength() const;
private:
  vector<string> endo_names, param_names, exo_names, exo_det_names;
  map<string, pair<SymbolType, int>> name_to_id_and_type;
};

#endif // _BASIC_SYMBOL_TABLE_HH
