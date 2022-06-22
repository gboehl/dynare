/*
 * Copyright © 2005 Ondra Kamenik
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

#include "parser_exception.hh"
#include "fine_atoms.hh"

using namespace ogp;

AllvarOuterOrdering::AllvarOuterOrdering(const vector<string> &allvar_outer,
                                         const FineAtoms &a)
  : atoms(a), allvar(),
    endo2all(a.get_endovars().size(), -1),
    exo2all(a.get_exovars().size(), -1)
{
  // fill in the allvar from allvar_outer
  for (const auto &s : allvar_outer)
    {
      if (atoms.varnames.query(s))
        allvar.push_back(s);
      else
        throw ogu::Exception(__FILE__, __LINE__,
                             "Variable " + s + " is not a declared symbol in AllvarOuterOrdering constructor");
    }

  // fill in endo2all and exo2all
  for (unsigned int i = 0; i < allvar.size(); i++)
    {
      auto it = atoms.endo_outer_map.find(allvar[i]);
      if (it != atoms.endo_outer_map.end())
        endo2all[it->second] = i;
      else
        {
          it = atoms.exo_outer_map.find(allvar[i]);
          if (it != atoms.exo_outer_map.end())
            exo2all[it->second] = i;
          else
            throw ogu::Exception(__FILE__, __LINE__,
                                 "Name " + allvar[i] + " is neither endogenous nor exogenous variable in AllvarOuterOrdering constructor");
        }
    }

  // check whether everything has been filled
  unsigned int iendo = 0;
  while (iendo < endo2all.size() && endo2all[iendo] != -1)
    iendo++;
  unsigned int iexo = 0;
  while (iexo < exo2all.size() && exo2all[iexo] != -1)
    iexo++;
  if (iendo < endo2all.size())
    throw ogu::Exception(__FILE__, __LINE__,
                         "Endogenous variable " + atoms.get_endovars()[iendo]
                         +" not found in outer all ordering in AllvarOuterOrdering constructor");
  if (iexo < exo2all.size())
    throw ogu::Exception(__FILE__, __LINE__,
                         "Exogenous variable " + atoms.get_exovars()[iexo]
                         +" not found in outer all ordering in AllvarOuterOrdering constructor");
}

AllvarOuterOrdering::AllvarOuterOrdering(const AllvarOuterOrdering &avo,
                                         const FineAtoms &a)
  : atoms(a), allvar(avo.allvar),
    endo2all(avo.endo2all),
    exo2all(avo.exo2all)
{
}

FineAtoms::FineAtoms(const FineAtoms &fa)
  : DynamicAtoms(fa), params(), endovars(), exovars(),
    der_atoms(fa.der_atoms),
    endo_atoms_map(fa.endo_atoms_map),
    exo_atoms_map(fa.exo_atoms_map)
{
  // fill in params
  for (const auto &param : fa.params)
    {
      if (!varnames.query(param))
        throw ogu::Exception(__FILE__, __LINE__,
                             "Parameter " + param + " does not exist in FineAtoms copy cosntructor");
      params.push_back(param);
      param_outer_map.emplace(param, params.size()-1);
    }
  // fill in endovars
  for (const auto &endovar : fa.endovars)
    {
      if (!varnames.query(endovar))
        throw ogu::Exception(__FILE__, __LINE__,
                             "Endo variable " + endovar + " does not exist in FineAtoms copy constructor");
      endovars.push_back(endovar);
      endo_outer_map.emplace(endovar, endovars.size()-1);
    }
  // fill in exovars
  for (const auto &exovar : fa.exovars)
    {
      if (!varnames.query(exovar))
        throw ogu::Exception(__FILE__, __LINE__,
                             "Exo variable " + exovar + " does not exist in FineAtoms copy cosntructor");
      exovars.push_back(exovar);
      exo_outer_map.emplace(exovar, exovars.size()-1);
    }

  if (fa.endo_order)
    endo_order = fa.endo_order->clone(endovars, *this);

  if (fa.exo_order)
    exo_order = fa.exo_order->clone(exovars, *this);

  if (fa.allvar_order)
    allvar_order = std::make_unique<AllvarOuterOrdering>(*(fa.allvar_order), *this);
}

int
FineAtoms::check_variable(const string &name) const
{
  string str;
  int ll;
  parse_variable(name, str, ll);
  if (varnames.query(str))
    return DynamicAtoms::check_variable(name);
  else
    {
      throw ParserException("Variable <"+str+"> not declared.", 0);
      return -1;
    }
}

int
FineAtoms::num_exo_periods() const
{
  int mlead, mlag;
  exovarspan(mlead, mlag);
  return mlead-mlag+1;
}

void
FineAtoms::parsing_finished(VarOrdering::ord_type ot)
{
  make_internal_orderings(ot);

  // by default, concatenate outer endo and outer exo and make it as
  // allvar outer:
  vector<string> allvar_tmp;
  allvar_tmp.insert(allvar_tmp.end(), endovars.begin(), endovars.end());
  allvar_tmp.insert(allvar_tmp.end(), exovars.begin(), exovars.end());

  allvar_order = std::make_unique<AllvarOuterOrdering>(allvar_tmp, *this);
}

void
FineAtoms::parsing_finished(VarOrdering::ord_type ot,
                            const vector<string> &allvar)
{
  make_internal_orderings(ot);
  allvar_order = std::make_unique<AllvarOuterOrdering>(allvar, *this);
}

const vector<string> &
FineAtoms::get_allvar() const
{
  if (!allvar_order)
    throw ogu::Exception(__FILE__, __LINE__,
                         "FineAtoms::get_allvars called before parsing_finished");

  return allvar_order->get_allvar();
}

const vector<int> &
FineAtoms::outer_endo2all() const
{
  if (!allvar_order)
    throw ogu::Exception(__FILE__, __LINE__,
                         "FineAtoms::outer_endo2all called before parsing_finished");

  return allvar_order->get_endo2all();
}

const vector<int> &
FineAtoms::outer_exo2all() const
{
  if (!allvar_order)
    throw ogu::Exception(__FILE__, __LINE__,
                         "FineAtoms::outer_exo2all called before parsing_finished");

  return allvar_order->get_exo2all();
}

vector<int>
FineAtoms::variables() const
{
  if (endo_order)
    return der_atoms;
  else
    {
      throw ogu::Exception(__FILE__, __LINE__,
                           "FineAtoms::variables called before parsing_finished");
      return {};
    }
}

int
FineAtoms::nstat() const
{
  if (endo_order)
    return endo_order->nstat();
  else
    {
      throw ogu::Exception(__FILE__, __LINE__,
                           "FineAtoms::nstat called before parsing_finished");
      return -1;
    }
}

int
FineAtoms::npred() const
{
  if (endo_order)
    return endo_order->npred();
  else
    {
      throw ogu::Exception(__FILE__, __LINE__,
                           "FineAtoms::npred called before parsing_finished");
      return -1;
    }
}

int
FineAtoms::nboth() const
{
  if (endo_order)
    return endo_order->nboth();
  else
    {
      throw ogu::Exception(__FILE__, __LINE__,
                           "FineAtoms::nboth called before parsing_finished");
      return -1;
    }
}

int
FineAtoms::nforw() const
{
  if (endo_order)
    return endo_order->nforw();
  else
    {
      throw ogu::Exception(__FILE__, __LINE__,
                           "FineAtoms::nforw called before parsing_finished");
      return -1;
    }
}

int
FineAtoms::get_pos_of_endo(int t) const
{
  if (endo_order)
    return endo_order->get_pos_of(t);
  else
    {
      throw ogu::Exception(__FILE__, __LINE__,
                           "FineAtoms::get_pos_of_endo called before parsing_finished");
      return -1;
    }
}

int
FineAtoms::get_pos_of_exo(int t) const
{
  if (exo_order)
    return exo_order->get_pos_of(t);
  else
    {
      throw ogu::Exception(__FILE__, __LINE__,
                           "FineAtoms::get_pos_of_exo called before parsing_finished");
      return -1;
    }
}

int
FineAtoms::get_pos_of_all(int t) const
{
  if (endo_order && exo_order)
    {
      if (endo_order->check(t))
        return endo_order->get_pos_of(t);
      else if (exo_order->check(t))
        return endo_order->length() + exo_order->get_pos_of(t);
      else
        {
          throw ogu::Exception(__FILE__, __LINE__,
                               "Atom is not endo nor exo in FineAtoms::get_pos_of_all");
          return -1;
        }
    }
  else
    {
      throw ogu::Exception(__FILE__, __LINE__,
                           "FineAtoms::get_pos_of_exo called before parsing_finished");
      return -1;
    }
}

const vector<int> &
FineAtoms::y2outer_endo() const
{
  if (!endo_order)
    throw ogu::Exception(__FILE__, __LINE__,
                         "FineAtoms::y2outer_endo called before parsing_finished");
  return endo_order->get_y2outer();
}

const vector<int> &
FineAtoms::outer2y_endo() const
{
  if (!endo_order)
    throw ogu::Exception(__FILE__, __LINE__,
                         "FineAtoms::outer2y_endo called before parsing_finished");
  return endo_order->get_outer2y();
}

const vector<int> &
FineAtoms::y2outer_exo() const
{
  if (!exo_order)
    throw ogu::Exception(__FILE__, __LINE__,
                         "FineAtoms::y2outer_endo called before parsing_finished");
  return exo_order->get_y2outer();
}

const vector<int> &
FineAtoms::outer2y_exo() const
{
  if (!exo_order)
    throw ogu::Exception(__FILE__, __LINE__,
                         "FineAtoms::outer2y_exo called before parsing_finished");
  return exo_order->get_outer2y();
}

const vector<int> &
FineAtoms::get_endo_atoms_map() const
{
  if (!endo_order)
    throw ogu::Exception(__FILE__, __LINE__,
                         "FineAtoms::get_endo_atoms_map called before parsing_finished");
  return endo_atoms_map;
}

const vector<int> &
FineAtoms::get_exo_atoms_map() const
{
  if (!exo_order)
    throw ogu::Exception(__FILE__, __LINE__,
                         "FineAtoms::get_exo_atoms_map called before parsing_finished");
  return exo_atoms_map;
}

int
FineAtoms::name2outer_param(const string &name) const
{
  auto it = param_outer_map.find(name);
  if (it == param_outer_map.end())
    throw ogu::Exception(__FILE__, __LINE__,
                         "Name is not a parameter in FineAtoms::name2outer_param");
  return it->second;
}

int
FineAtoms::name2outer_endo(const string &name) const
{
  auto it = endo_outer_map.find(name);
  if (it == endo_outer_map.end())
    throw ogu::Exception(__FILE__, __LINE__,
                         "Name is not an endogenous variable in FineAtoms::name2outer_endo");
  return it->second;
}

int
FineAtoms::name2outer_exo(const string &name) const
{
  auto it = exo_outer_map.find(name);
  if (it == exo_outer_map.end())
    throw ogu::Exception(__FILE__, __LINE__,
                         "Name is not an exogenous variable in FineAtoms::name2outer_exo");
  return it->second;
}

int
FineAtoms::name2outer_allvar(const string &name) const
{
  if (!allvar_order)
    throw ogu::Exception(__FILE__, __LINE__,
                         "FineAtoms::name2outer_allvar called beore parsing_finished");

  auto it = endo_outer_map.find(name);
  if (it != endo_outer_map.end())
    return allvar_order->get_endo2all()[it->second];
  else
    {
      it = exo_outer_map.find(name);
      if (it != exo_outer_map.end())
        return allvar_order->get_exo2all()[it->second];
    }

  throw ogu::Exception(__FILE__, __LINE__,
                       "Name " + name + " is neither endo nor exo variable in FineAtoms::name2outer_allvar");
  return -1;
}

void
FineAtoms::register_uniq_endo(string name)
{
  if (varnames.query(name))
    throw ogp::ParserException("Endogenous variable <"+name+"> is not unique.", 0);
  varnames.insert(name);
  endovars.push_back(name);
  endo_outer_map.emplace(std::move(name), endovars.size()-1);
}

void
FineAtoms::register_uniq_exo(string name)
{
  if (varnames.query(name))
    throw ogp::ParserException("Exogenous variable <"+name+"> is not unique.", 0);
  varnames.insert(name);
  exovars.push_back(name);
  exo_outer_map.emplace(std::move(name), exovars.size()-1);
}

void
FineAtoms::register_uniq_param(string name)
{
  if (varnames.query(name))
    throw ogp::ParserException("Parameter <"+name+"> is not unique.", 0);
  varnames.insert(name);
  params.push_back(name);
  param_outer_map.emplace(std::move(name), params.size()-1);
}

void
FineAtoms::make_internal_orderings(VarOrdering::ord_type ot)
{
  bool endo_ordering_done = false;
  bool exo_ordering_done = false;

  order_type = ot;

  int mlead, mlag;
  endovarspan(mlead, mlag);
  if (mlag >= -1 && mlead <= 1)
    {
      // make endo ordering
      if (ot == VarOrdering::pbspbfbf)
        endo_order = std::make_unique<EndoVarOrdering1>(endovars, *this);
      else
        endo_order = std::make_unique<EndoVarOrdering2>(endovars, *this);
      endo_order->do_ordering();
      endo_ordering_done = true;
    }

  exovarspan(mlead, mlag);
  if (mlag == 0 && mlead == 0)
    {
      // make exo ordering
      exo_order = std::make_unique<ExoVarOrdering>(exovars, *this);
      exo_order->do_ordering();
      exo_ordering_done = true;
    }

  if (endo_ordering_done && exo_ordering_done)
    {
      // concatenate der atoms from endo_order and exo_order
      der_atoms.clear();
      der_atoms.insert(der_atoms.end(),
                       endo_order->get_der_atoms().begin(),
                       endo_order->get_der_atoms().end());
      der_atoms.insert(der_atoms.end(),
                       exo_order->get_der_atoms().begin(),
                       exo_order->get_der_atoms().end());

      // create endo_atoms_map; der_atoms is a concatenation, so it is easy
      int endo_atoms = endo_order->get_der_atoms().size();
      endo_atoms_map.clear();
      for (int i = 0; i < endo_atoms; i++)
        endo_atoms_map.push_back(i);
      // create exo_atoms_map
      int exo_atoms = exo_order->get_der_atoms().size();
      exo_atoms_map.clear();
      for (int i = 0; i < exo_atoms; i++)
        exo_atoms_map.push_back(endo_atoms + i);
    }
}

void
FineAtoms::print() const
{
  DynamicAtoms::print();
  if (endo_order)
    {
      std::cout << "Endo ordering:\n";
      endo_order->print();
    }
  else
    std::cout << "Endo ordering not created.\n";

  if (exo_order)
    {
      std::cout << "Exo ordering:\n";
      exo_order->print();
    }
  else
    std::cout << "Exo ordering not created.\n";

  std::cout << "endo atoms map:\n";
  for (unsigned int i = 0; i < endo_atoms_map.size(); i++)
    std::cout << i << " → " << endo_atoms_map[i] << "\n";
  std::cout << "exo atoms map:\n";
  for (unsigned int i = 0; i < exo_atoms_map.size(); i++)
    std::cout << i << " → " << exo_atoms_map[i] << "\n";
}
