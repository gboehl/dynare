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

#include "parser/cc/parser_exception.hh"
#include "utils/cc/exception.hh"

#include "dynare_atoms.hh"

#include <string>
#include <cmath>
#include <limits>
#include <sstream>
#include <iomanip>

using namespace ogdyn;
using std::string;

void
DynareStaticAtoms::register_name(string name)
{
  if (varnames.query(name))
    throw ogp::ParserException("The name "+name+" is not unique.", 0);
  StaticAtoms::register_name(std::move(name));
}

int
DynareStaticAtoms::check_variable(const string &name) const
{
  if (!varnames.query(name))
    throw ogp::ParserException("Unknown name <"+name+">", 0);
  auto it = vars.find(name);
  if (it == vars.end())
    return -1;
  else
    return it->second;
}

void
DynareDynamicAtoms::parse_variable(const string &in, std::string &out, int &ll) const
{
  ll = 0;
  auto left = in.find_first_of("({");
  if (left != string::npos)
    {
      out = in.substr(0, left);
      left++;
      auto right = in.find_first_of(")}", left);
      if (string::npos == right)
        throw ogp::ParserException("Syntax error when parsing Dynare atom <"+in+">.", 0);
      ll = std::stoi(in.substr(left, right-left));
    }
  else
    out = in;
}

void
DynareDynamicAtoms::register_uniq_endo(string name)
{
  FineAtoms::register_uniq_endo(name);
  atom_type.emplace(std::move(name), atype::endovar);
}

void
DynareDynamicAtoms::register_uniq_exo(string name)
{
  FineAtoms::register_uniq_exo(name);
  atom_type.emplace(std::move(name), atype::exovar);
}

void
DynareDynamicAtoms::register_uniq_param(string name)
{
  FineAtoms::register_uniq_param(name);
  atom_type.emplace(std::move(name), atype::param);
}

bool
DynareDynamicAtoms::is_type(const string &name, atype tp) const
{
  auto it = atom_type.find(name);
  if (it != atom_type.end() && it->second == tp)
    return true;
  else
    return false;
}

void
DynareDynamicAtoms::print() const
{
  SAtoms::print();
  std::cout << "Name types:\n";
  for (auto it : atom_type)
    std::cout << "name=" << it.first << " type="
              << (it.second == atype::endovar ? "endovar" : it.second == atype::exovar ? "exovar" : "param")
              << '\n';
}

std::string
DynareDynamicAtoms::convert(int t) const
{
  if (t < ogp::OperationTree::num_constants)
    {
      throw ogu::Exception(__FILE__, __LINE__,
                           "Tree index is a built-in constant in DynareDynamicAtoms::convert");
      return {};
    }
  if (is_constant(t))
    {
      double v = get_constant_value(t);
      std::ostringstream buf;
      buf << std::setprecision(std::numeric_limits<double>::max_digits10)
          << v;
      return buf.str();
    }

  const string &s = name(t);
  if (is_type(s, atype::endovar))
    {
      int ll = lead(t);
      if (ll)
        return s + '(' + std::to_string(ll) + ')';
    }

  return s;
}

void
DynareAtomValues::setValues(ogp::EvalTree &et) const
{
  // set constants
  atoms.setValues(et);

  // set parameteres
  for (unsigned int i = 0; i < atoms.get_params().size(); i++)
    if (atoms.is_referenced(atoms.get_params()[i]))
      {
        const ogp::DynamicAtoms::Tlagmap &lmap = atoms.lagmap(atoms.get_params()[i]);
        for (auto it : lmap)
          {
            int t = it.second;
            et.set_nulary(t, paramvals[i]);
          }
      }

  // set endogenous
  for (unsigned int outer_i = 0; outer_i < atoms.get_endovars().size(); outer_i++)
    if (atoms.is_referenced(atoms.get_endovars()[outer_i]))
      {
        const ogp::DynamicAtoms::Tlagmap &lmap = atoms.lagmap(atoms.get_endovars()[outer_i]);
        for (auto it : lmap)
          {
            int ll = it.first;
            int t = it.second;
            int i = atoms.outer2y_endo()[outer_i];
            if (ll == -1)
              et.set_nulary(t, yym[i-atoms.nstat()]);
            else if (ll == 0)
              et.set_nulary(t, yy[i]);
            else
              et.set_nulary(t, yyp[i-atoms.nstat()-atoms.npred()]);
          }
      }

  // set exogenous
  for (unsigned int outer_i = 0; outer_i < atoms.get_exovars().size(); outer_i++)
    if (atoms.is_referenced(atoms.get_exovars()[outer_i]))
      {
        const ogp::DynamicAtoms::Tlagmap &lmap = atoms.lagmap(atoms.get_exovars()[outer_i]);
        for (auto it : lmap)
          {
            int ll = it.first;
            if (ll == 0) // this is always true because of checks
              {
                int t = it.second;
                int i = atoms.outer2y_exo()[outer_i];
                et.set_nulary(t, xx[i]);
              }
          }
      }
}

void
DynareStaticSteadyAtomValues::setValues(ogp::EvalTree &et) const
{
  // set constants
  atoms_static.setValues(et);

  // set parameters
  for (auto name : atoms_static.get_params())
    {
      int t = atoms_static.index(name);
      if (t != -1)
        {
          int idyn = atoms.name2outer_param(name);
          et.set_nulary(t, paramvals[idyn]);
        }
    }

  // set endogenous
  for (auto name : atoms_static.get_endovars())
    {
      int t = atoms_static.index(name);
      if (t != -1)
        {
          int idyn = atoms.outer2y_endo()[atoms.name2outer_endo(name)];
          et.set_nulary(t, yy[idyn]);
        }
    }

  // set exogenous
  for (auto name : atoms_static.get_exovars())
    {
      int t = atoms_static.index(name);
      if (t != -1)
        et.set_nulary(t, 0.0);
    }
}

DynareSteadySubstitutions::DynareSteadySubstitutions(const ogp::FineAtoms &a,
                                                     const ogp::OperationTree &tree,
                                                     const Tsubstmap &subst,
                                                     const Vector &pvals, Vector &yy)
  : atoms(a), y(yy)
{
  // fill the vector of left and right hand sides
  for (auto it : subst)
    {
      left_hand_sides.push_back(it.first);
      right_hand_sides.push_back(it.second);
    }

  // evaluate right hand sides
  DynareSteadyAtomValues dsav(atoms, pvals, y);
  ogp::FormulaCustomEvaluator fe(tree, right_hand_sides);
  fe.eval(dsav, *this);
}

void
DynareSteadySubstitutions::load(int i, double res)
{
  const string &name = left_hand_sides[i];
  int iouter = atoms.name2outer_endo(name);
  int iy = atoms.outer2y_endo()[iouter];
  if (!std::isfinite(y[iy]))
    y[iy] = res;
}

DynareStaticSteadySubstitutions::
DynareStaticSteadySubstitutions(const ogp::FineAtoms &a, const ogp::StaticFineAtoms &sa,
                                const ogp::OperationTree &tree,
                                const Tsubstmap &subst,
                                const Vector &pvals, Vector &yy)
  : atoms(a), atoms_static(sa), y(yy)
{
  // fill the vector of left and right hand sides
  for (const auto &it : subst)
    {
      left_hand_sides.push_back(it.first);
      right_hand_sides.push_back(it.second);
    }

  // evaluate right hand sides
  DynareStaticSteadyAtomValues dsav(atoms, atoms_static, pvals, y);
  ogp::FormulaCustomEvaluator fe(tree, right_hand_sides);
  fe.eval(dsav, *this);
}

void
DynareStaticSteadySubstitutions::load(int i, double res)
{
  const string &name = left_hand_sides[i];
  int iouter = atoms.name2outer_endo(name);
  int iy = atoms.outer2y_endo()[iouter];
  if (!std::isfinite(y[iy]))
    y[iy] = res;
}
