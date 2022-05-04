/*
 * Copyright © 2006 Ondra Kamenik
 * Copyright © 2019-2022 Dynare Team
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

#include "utils/cc/exception.hh"

#include "static_fine_atoms.hh"
#include "parser_exception.hh"

using namespace ogp;

void
StaticFineAtoms::import_atoms(const FineAtoms &fa, OperationTree &otree, Tintintmap &tmap)
{
  StaticAtoms::import_atoms(fa, otree, tmap);

  // we just need to put parameters, endovars, and exovars to
  // respective vectors, the names are already in the storage

  // parameters
  auto &fa_params = fa.get_params();
  for (const auto &fa_param : fa_params)
    register_param(fa_param);

  // endogenous
  auto &fa_endovars = fa.get_endovars();
  for (const auto &fa_endovar : fa_endovars)
    register_endo(fa_endovar);

  // exogenous
  auto &fa_exovars = fa.get_exovars();
  for (const auto &fa_exovar : fa_exovars)
    register_exo(fa_exovar);

  parsing_finished();
}

void
StaticFineAtoms::import_atoms(const FineAtoms &fa, OperationTree &otree, Tintintmap &tmap,
                              const char *dummy)
{
  StaticAtoms::import_atoms(fa, otree, tmap);

  // we just need to put parameters, endovars, and exovars to
  // respective vectors, the names are already in the storage

  // parameters
  auto &fa_params = fa.get_params();
  for (const auto &fa_param : fa_params)
    register_param(fa_param);

  // endogenous
  auto &fa_endovars = fa.get_endovars();
  for (unsigned int i = 0; i < fa_endovars.size(); i++)
    register_endo(fa_endovars[fa.y2outer_endo()[i]]);

  // exogenous
  auto &fa_exovars = fa.get_exovars();
  for (unsigned int i = 0; i < fa_exovars.size(); i++)
    register_exo(fa_exovars[fa.y2outer_exo()[i]]);

  parsing_finished();
}

int
StaticFineAtoms::check_variable(const string &name) const
{
  if (!varnames.query(name))
    throw ParserException(string("Variable <")+name+"> not declared.", 0);
  return index(name);
}

void
StaticFineAtoms::parsing_finished()
{
  // build der_atoms, and endo_atoms_map and exo_atoms_map
  der_atoms.clear();
  endo_atoms_map.clear();
  exo_atoms_map.clear();

  // go through all endo and exo insert tree indices, ignore names
  // whose tree index is -1 (those which are not referenced)
  for (auto &endovar : endovars)
    {
      int t = index(endovar);
      if (t != -1)
        {
          endo_atoms_map.push_back(der_atoms.size());
          der_atoms.push_back(t);
        }
    }
  for (auto &exovar : exovars)
    {
      int t = index(exovar);
      if (t != -1)
        {
          exo_atoms_map.push_back(der_atoms.size());
          der_atoms.push_back(t);
        }
    }
}

int
StaticFineAtoms::name2outer_param(const string &name) const
{
  auto it = param_outer_map.find(name);
  if (it == param_outer_map.end())
    throw ogu::Exception(__FILE__, __LINE__,
                         "Name is not a parameter in StaticFineAtoms::name2outer_param");
  return it->second;
}

int
StaticFineAtoms::name2outer_endo(const string &name) const
{
  auto it = endo_outer_map.find(name);
  if (it == endo_outer_map.end())
    throw ogu::Exception(__FILE__, __LINE__,
                         "Name is not an endogenous variable in StaticFineAtoms::name2outer_endo");
  return it->second;
}

int
StaticFineAtoms::name2outer_exo(const string &name) const
{
  auto it = exo_outer_map.find(name);
  if (it == exo_outer_map.end())
    throw ogu::Exception(__FILE__, __LINE__,
                         "Name is not an exogenous variable in StaticFineAtoms::name2outer_exo");
  return it->second;
}

void
StaticFineAtoms::register_uniq_endo(string name)
{
  if (varnames.query(name))
    throw ogp::ParserException(string("Endogenous variable <")+name+"> is not unique.", 0);
  varnames.insert(name);
  register_endo(std::move(name));
}

void
StaticFineAtoms::register_uniq_exo(string name)
{
  if (varnames.query(name))
    throw ogp::ParserException(string("Exogenous variable <")+name+"> is not unique.", 0);
  varnames.insert(name);
  register_exo(std::move(name));
}

void
StaticFineAtoms::register_uniq_param(string name)
{
  if (varnames.query(name))
    throw ogp::ParserException(string("Parameter <")+name+"> is not unique.", 0);
  varnames.insert(name);
  register_param(std::move(name));
}

void
StaticFineAtoms::print() const
{
  StaticAtoms::print();
  std::cout << "endo atoms map:\n";
  for (unsigned int i = 0; i < endo_atoms_map.size(); i++)
    std::cout << i << " → " << endo_atoms_map[i] << "\n";
  std::cout << "exo atoms map:\n";
  for (unsigned int i = 0; i < exo_atoms_map.size(); i++)
    std::cout << i << " → " << exo_atoms_map[i] << "\n";
  std::cout << "der atoms:\n";
  for (unsigned int i = 0; i < der_atoms.size(); i++)
    std::cout << i << "\t" << der_atoms[i] << "\n";
}

void
StaticFineAtoms::register_endo(string name)
{
  if (!varnames.query(name))
    throw ogp::ParserException(string("Endogenous variable <")
                               +name+"> not found in storage.", 0);
  endovars.push_back(name);
  endo_outer_map.emplace(std::move(name), endovars.size()-1);
}

void
StaticFineAtoms::register_exo(string name)
{
  if (!varnames.query(name))
    throw ogp::ParserException(string("Exogenous variable <")
                               +name+"> not found in storage.", 0);
  exovars.push_back(name);
  exo_outer_map.emplace(std::move(name), exovars.size()-1);
}

void
StaticFineAtoms::register_param(string name)
{
  if (!varnames.query(name))
    throw ogp::ParserException(string("Parameter <")+name+"> not found in storage.", 0);
  params.push_back(name);
  param_outer_map.emplace(std::move(name), params.size()-1);
}
