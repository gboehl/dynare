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

#include "BasicSymbolTable.hh"

#include "dynmex.h"

BasicSymbolTable::BasicSymbolTable()
{
  mxArray *M_ = mexGetVariable("global", "M_");
  if (!M_)
    mexErrMsgTxt("Can't find global variable M_");

  auto get_field_names = [&](const char *field_name, SymbolType type)
  {
    vector<string> r;
    if (mxGetFieldNumber(M_, field_name) != -1)
      {
        auto M_field = mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, field_name));
        if (!mxIsCell(M_field))
          mexErrMsgTxt(("M_."s + field_name + " is not a cell array").c_str());
        for (size_t i {0}; i < mxGetNumberOfElements(M_field); i++)
          {
            const mxArray *cell_mx = mxGetCell(M_field, i);
            if (!(cell_mx && mxIsChar(cell_mx)))
              mexErrMsgTxt(("M_."s + field_name + " contains a cell which is not a character array").c_str());
            r.emplace_back(mxArrayToString(cell_mx));
          }
      }

    // Fill the reverse map
    for (size_t i {0}; i < r.size(); i++)
      name_to_id_and_type.emplace(r[i], pair{type, i});

    return r;
  };
  endo_names = get_field_names("endo_names", SymbolType::endogenous);
  param_names = get_field_names("param_names", SymbolType::parameter);
  exo_names = get_field_names("exo_names", SymbolType::exogenous);
  exo_det_names = get_field_names("exo_det_names", SymbolType::exogenousDet);
}

string
BasicSymbolTable::getName(SymbolType type, int tsid) const
{
  try
    {
      switch (type)
        {
        case SymbolType::endogenous:
          return endo_names.at(tsid);
        case SymbolType::exogenous:
          return exo_names.at(tsid);
        case SymbolType::exogenousDet:
          return exo_det_names.at(tsid);
        case SymbolType::parameter:
          return param_names.at(tsid);
        default:
          mexErrMsgTxt(("Unsupported symbol type: " + to_string(static_cast<int>(type))).c_str());
        }
    }
  catch (out_of_range &)
    {
      mexErrMsgTxt(("Unknown symbol with ID " + to_string(tsid) + " and type " + to_string(static_cast<int>(type))).c_str());
    }
  exit(EXIT_FAILURE); // Silence GCC warning
}

pair<SymbolType, int>
BasicSymbolTable::getIDAndType(const string &name) const
{
  try
    {
      return name_to_id_and_type.at(name);
    }
  catch (out_of_range &)
    {
      mexErrMsgTxt(("Unknown symbol: " + name).c_str());
    }
  exit(EXIT_FAILURE); // Silence GCC warning
}

size_t
BasicSymbolTable::maxEndoNameLength() const
{
  size_t r {0};
  for (const auto &n : endo_names)
    r = max(r, n.size());
  return r;
}
