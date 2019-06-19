/*
 * Copyright © 2006 Ondra Kamenik
 * Copyright © 2019 Dynare Team
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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "static_atoms.hh"
#include "utils/cc/exception.hh"

using namespace ogp;

void
StaticAtoms::import_atoms(const DynamicAtoms &da, OperationTree &otree, Tintintmap &tmap)
{
  Constants::import_constants(da, otree, tmap);

  for (int i = 0; i < da.get_name_storage().num(); i++)
    {
      const string &name = da.get_name_storage().get_name(i);
      register_name(name);
      int tnew = otree.add_nulary();
      assign(name, tnew);
      if (da.is_referenced(name))
        {
          const DynamicAtoms::Tlagmap &lmap = da.lagmap(name);
          for (const auto &it : lmap)
            {
              int told = it.second;
              tmap.emplace(told, tnew);
            }
        }
    }
}

int
StaticAtoms::check(const string &name) const
{
  if (DynamicAtoms::is_string_constant(name))
    return Constants::check(name);
  else
    return check_variable(name);
}

int
StaticAtoms::index(const string &name) const
{
  auto it = vars.find(name);
  if (it == vars.end())
    return -1;
  else
    return it->second;
}

void
StaticAtoms::assign(const string &name, int t)
{
  if (DynamicAtoms::is_string_constant(name))
    {
      double val = std::stod(name);
      add_constant(t, val);
    }
  else
    {
      varnames.insert(name);
      vars.emplace(name, t);
      indices.emplace(t, name);
    }
}

vector<int>
StaticAtoms::variables() const
{
  vector<int> res;
  for (const auto &var : vars)
    res.push_back(var.second);
  return res;
}

void
StaticAtoms::register_name(string name)
{
  varnames.insert(name);
  varorder.push_back(std::move(name));
}

void
StaticAtoms::print() const
{
  std::cout << "constants:\n";
  Constants::print();
  std::cout << "variable names:\n";
  varnames.print();
  std::cout << "map to tree indices:\n";
  for (auto var : vars)
    std::cout << var.first << u8"\t→\t" << var.second << "\n";
}
