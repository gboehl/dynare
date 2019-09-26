/*
 * Copyright © 2005 Ondra Kamenik
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

#include "utils/cc/exception.hh"
#include "dynamic_atoms.hh"

using namespace ogp;

void
NameStorage::insert(string name)
{
  if (!query(name))
    {
      name_store.push_back(name);
      name_set.insert(std::move(name));
    }
}

void
NameStorage::print() const
{
  for (auto i : name_store)
    std::cout << i << '\n';
}

void
Constants::import_constants(const Constants &c, OperationTree &otree, Tintintmap &tmap)
{
  for (auto it : c.cmap)
    {
      int told = it.first;
      int tnew = otree.add_nulary();
      tmap.emplace(told, tnew);
      add_constant(tnew, it.second);
    }
}

void
Constants::setValues(EvalTree &et) const
{
  for (const auto & it : cmap)
    et.set_nulary(it.first, it.second);
}

void
Constants::add_constant(int t, double val)
{
  cmap.emplace(t, val);
  cinvmap.emplace(val, t);
}

bool
Constants::is_constant(int t) const
{
  if (t < OperationTree::num_constants)
    return true;
  auto it = cmap.find(t);
  return (it != cmap.end());
}

double
Constants::get_constant_value(int t) const
{
  auto it = cmap.find(t);
  if (it != cmap.end())
    return it->second;
  else
    throw ogu::Exception(__FILE__, __LINE__,
                         "Tree index is not constant in Constants::get_constant_value");
}

int
Constants::check(const string &str) const
{
  double d = std::stod(str);
  auto it = cinvmap.find(d);
  if (it != cinvmap.end())
    return it->second;
  else
    return -1;
}

void
Constants::print() const
{
  for (const auto &it : cmap)
    std::cout << "$" << it.first << ":  " << it.second << "\n";
}

int
DynamicAtoms::check(const string &name) const
{
  if (is_string_constant(name))
    return Constants::check(name);

  return check_variable(name);
}

int
DynamicAtoms::check_variable(const string &name) const
{
  string str;
  int ll;
  parse_variable(name, str, ll);
  auto it = vars.find(str);

  if (it != vars.end())
    {
      const Tlagmap &lmap = it->second;
      auto itt = lmap.find(ll);
      if (itt != lmap.end())
        return itt->second;
    }
  return -1;
}

void
DynamicAtoms::assign(const string &name, int t)
{
  if (is_string_constant(name))
    assign_constant(name, t);
  else
    assign_variable(name, t);
}

void
DynamicAtoms::assign_constant(const string &name, int t)
{
  double val = std::stod(name);
  add_constant(t, val);
}

// parse the name and then call assing_variable(varname, ll, t)

void
DynamicAtoms::assign_variable(const string &name, int t)
{
  int ll;
  string str;
  parse_variable(name, str, ll);
  // here str is just name without lead/lag
  varnames.insert(str);

  assign_variable(str, ll, t);
}

void
DynamicAtoms::assign_variable(const string &varname, int ll, int t)
{
  if (indices.end() != indices.find(t))
    throw ogu::Exception(__FILE__, __LINE__,
                         "Attempt to assign already allocated tree index");

  auto it = vars.find(varname);
  if (it != vars.end())
    {
      Tlagmap &lmap = it->second;
      if (lmap.end() != lmap.find(ll))
        throw ogu::Exception(__FILE__, __LINE__,
                             "Attempt to assign already allocated variable");
      lmap.emplace(ll, t);
    }
  else
    {
      Tlagmap lmap;
      lmap.emplace(ll, t);
      vars.emplace(varname, lmap);
    }
  indices.emplace(t, varname);

  nv++;
  minlag = std::min(ll, minlag);
  maxlead = std::max(ll, maxlead);
}

void
DynamicAtoms::unassign_variable(const string &varname, int ll, int t)
{
  auto it = vars.find(varname);
  if (it != vars.end())
    {
      Tlagmap &lmap = it->second;
      auto itt = lmap.find(ll);
      if (itt != lmap.end())
        {
          if (itt->second == t)
            {
              // erase it from the lagmap; if it becomes empty,
              // erase the lagmap from varmap
              lmap.erase(itt);
              if (lmap.size() == 0)
                vars.erase(it);
              // erase it from the indices
              auto ittt = indices.find(t);
              if (ittt != indices.end())
                indices.erase(ittt);

              nv--;
              if (ll == minlag || ll == maxlead)
                update_minmaxll();
            }
          else
            throw ogu::Exception(__FILE__, __LINE__,
                                 "Tree index inconsistent in DynamicAtoms::unassign_variable");
        }
      else
        throw ogu::Exception(__FILE__, __LINE__,
                             "Lead/lag of the variable not found in DynamicAtoms::unassign_variable");
    }
  else
    throw ogu::Exception(__FILE__, __LINE__,
                         "Variable not found in DynamicAtoms::unassign_variable");
}

void
DynamicAtoms::update_minmaxll()
{
  minlag = std::numeric_limits<int>::max();
  maxlead = std::numeric_limits<int>::min();
  for (const auto &it : vars)
    {
      const Tlagmap &lmap = it.second;
      for (auto itt : lmap)
        {
          int ll = itt.first;
          minlag = std::min(ll, minlag);
          maxlead = std::max(ll, maxlead);
        }
    }
}

vector<int>
DynamicAtoms::variables() const
{
  vector<int> res;
  for (const auto & var : vars)
    {
      const Tlagmap &lmap = var.second;
      for (auto itt : lmap)
        res.push_back(itt.second);
    }
  return res;
}

void
DynamicAtoms::varspan(int t, int &mlead, int &mlag) const
{
  auto it = indices.find(t);
  if (indices.end() == it)
    {
      mlead = std::numeric_limits<int>::min();
      mlag = std::numeric_limits<int>::max();
      return;
    }
  varspan(it->second, mlead, mlag);
}

void
DynamicAtoms::varspan(const string &name, int &mlead, int &mlag) const
{
  auto it = vars.find(name);
  if (vars.end() == it)
    {
      mlead = std::numeric_limits<int>::min();
      mlag = std::numeric_limits<int>::max();
      return;
    }
  const Tlagmap &lmap = it->second;
  auto beg = lmap.begin();
  auto end = lmap.rbegin();
  mlag = beg->first;
  mlead = end->first;
}

void
DynamicAtoms::varspan(const vector<string> &names, int &mlead, int &mlag) const
{
  mlead = std::numeric_limits<int>::min();
  mlag = std::numeric_limits<int>::max();
  for (const auto &name : names)
    {
      int lag, lead;
      varspan(name, lead, lag);
      mlead = std::max(lead, mlead);
      mlag = std::min(lag, mlag);
    }
}

bool
DynamicAtoms::is_named_atom(int t) const
{
  return indices.end() != indices.find(t);
}

int
DynamicAtoms::index(const string &name, int ll) const
{
  auto it = vars.find(name);
  if (vars.end() != it)
    {
      const Tlagmap &lmap = it->second;
      auto itt = lmap.find(ll);
      if (lmap.end() != itt)
        return itt->second;
    }
  return -1;
}

bool
DynamicAtoms::is_referenced(const string &name) const
{
  return vars.find(name) != vars.end();
}

const DynamicAtoms::Tlagmap &
DynamicAtoms::lagmap(const string &name) const
{
  auto it = vars.find(name);
  if (vars.end() == it)
    throw ogu::Exception(__FILE__, __LINE__,
                         std::string("Couldn't find the name ")
                         + name + " in DynamicAtoms::lagmap");
  return it->second;
}

const string &
DynamicAtoms::name(int t) const
{
  auto it = indices.find(t);
  if (indices.end() == it)
    throw ogu::Exception(__FILE__, __LINE__,
                         "Couldn't find tree index in DynamicAtoms::name");
  return it->second;
}

int
DynamicAtoms::lead(int t) const
{
  const string &nam = name(t);
  const Tlagmap &lmap = lagmap(nam);
  auto it = lmap.begin();
  while (it != lmap.end() && it->second != t)
    ++it;
  if (lmap.end() == it)
    throw ogu::Exception(__FILE__, __LINE__,
                         "Couldn't find the three index in DynamicAtoms::lead");
  return it->first;
}

void
DynamicAtoms::print() const
{
  std::cout << "names:\n";
  varnames.print();
  std::cout << "constants:\n";
  Constants::print();
  std::cout << "variables:\n";
  for (const auto & var : vars)
    {
      const Tlagmap &lmap = var.second;
      for (auto itt : lmap)
        std::cout << "$" << itt.second << ": " << var.first << "(" << itt.first << ")\n";
    }
  std::cout << "indices:\n";
  for (auto indice : indices)
    std::cout << "t=" << indice.first << u8" ⇒ " << indice.second << "\n";
}

/** Note that the str has been parsed by the lexicographic
 * analyzer. It can be either a variable or a double. So it is easy to
 * recognize it by the first character. */
bool
DynamicAtoms::is_string_constant(const string &str)
{
  return str[0] == '.' || str[0] == '-' || (str[0] >= '0' && str[0] <= '9');
}

VarOrdering::VarOrdering(const VarOrdering &vo, const vector<string> &vnames,
                         const DynamicAtoms &a)
  : n_stat(vo.n_stat), n_pred(vo.n_pred), n_both(vo.n_both), n_forw(vo.n_forw),
    der_atoms(vo.der_atoms), positions(vo.positions),
    outer2y(vo.outer2y), y2outer(vo.y2outer), varnames(vnames), atoms(a)
{
}

bool
VarOrdering::check(int t) const
{
  return positions.find(t) != positions.end();
}

int
VarOrdering::get_pos_of(int t) const
{
  auto it = positions.find(t);
  if (it != positions.end())
    return it->second;
  else
    throw ogu::Exception(__FILE__, __LINE__,
                         "Couldn't find the tree index in VarOrdering::get_pos_of");
}

void
VarOrdering::do_general(ord_type ordering)
{
  // auxiliary vectors for setting der_atoms and map
  vector<int> pred_minus;
  vector<int> both_minus;
  vector<int> stat;
  vector<int> pred_pad;
  vector<int> both_pad;
  vector<int> forw_pad;
  vector<int> both_plus;
  vector<int> forw_plus;

  // auxiliary vectors for setting y2outer and outer2y
  vector<int> y2o_stat;
  vector<int> y2o_pred;
  vector<int> y2o_both;
  vector<int> y2o_forw;

  for (unsigned int i = 0; i < varnames.size(); i++)
    {
      const string &ss = varnames[i];
      int lead;
      int lag;
      atoms.varspan(ss, lead, lag);
      if (lag == 0 && lead == 0)
        {
          stat.push_back(atoms.index(ss, 0));
          y2o_stat.push_back(i);
        }
      else if (lag == -1 && lead < 1)
        {
          pred_minus.push_back(atoms.index(ss, -1));
          pred_pad.push_back(atoms.index(ss, 0));
          y2o_pred.push_back(i);
        }
      else if (lag > -1 && lead == 1)
        {
          forw_pad.push_back(atoms.index(ss, 0));
          forw_plus.push_back(atoms.index(ss, 1));
          y2o_forw.push_back(i);
        }
      else if (lag == -1 && lead == 1)
        {
          both_minus.push_back(atoms.index(ss, -1));
          both_pad.push_back(atoms.index(ss, 0));
          both_plus.push_back(atoms.index(ss, 1));
          y2o_both.push_back(i);
        }
      else
        throw ogu::Exception(__FILE__, __LINE__,
                             "A wrong lag/lead of a variable in VarOrdering::do_pbspbfbf");
    }

  // here we fill ords according to ordering
  vector<int> *ords[8];
  if (ordering == pbspbfbf)
    {
      ords[0] = &pred_minus;
      ords[1] = &both_minus;
      ords[2] = &stat;
      ords[3] = &pred_pad;
      ords[4] = &both_pad;
      ords[5] = &forw_pad;
      ords[6] = &both_plus;
      ords[7] = &forw_plus;
    }
  else if (ordering == bfspbfpb)
    {
      ords[0] = &both_plus;
      ords[1] = &forw_plus;
      ords[2] = &stat;
      ords[3] = &pred_pad;
      ords[4] = &both_pad;
      ords[5] = &forw_pad;
      ords[6] = &pred_minus;
      ords[7] = &both_minus;
    }
  else // BEWARE: when implementing a new ordering, check also the code below setting y2outer
    throw ogu::Exception(__FILE__, __LINE__,
                         "Ordering not implemented in VarOrdering::do_general");

  // make der_atoms and positions
  int off = 0;
  for (auto & ord : ords)
    for (unsigned int j = 0; j < ord->size(); j++, off++)
      if ((*ord)[j] != -1)
        {
          der_atoms.push_back((*ord)[j]);
          positions.emplace((*ord)[j], off);
        }

  // set integer constants
  n_stat = stat.size();
  n_pred = pred_pad.size();
  n_both = both_pad.size();
  n_forw = forw_pad.size();

  // make y2outer mapping
  y2outer.insert(y2outer.end(), y2o_stat.begin(), y2o_stat.end());
  y2outer.insert(y2outer.end(), y2o_pred.begin(), y2o_pred.end());
  y2outer.insert(y2outer.end(), y2o_both.begin(), y2o_both.end());
  y2outer.insert(y2outer.end(), y2o_forw.begin(), y2o_forw.end());
  // make outer2y mapping
  outer2y.resize(y2outer.size(), -1);
  for (unsigned int i = 0; i < y2outer.size(); i++)
    outer2y[y2outer[i]] = i;
}

void
VarOrdering::do_increasing_time()
{
  // get maxlead and minlag of the variables
  int mlag, mlead;
  atoms.varspan(varnames, mlead, mlag);
  // setup the matrix of tree indices, if there is no occurrence,
  // the index is set to -1
  vector<int> ll_init(varnames.size(), -1);
  vector<vector<int>> tree_ind(mlead-mlag+1, ll_init);
  for (unsigned int iv = 0; iv < varnames.size(); iv++)
    {
      try
        {
          const DynamicAtoms::Tlagmap &lmap = atoms.lagmap(varnames[iv]);
          for (auto it : lmap)
            {
              int ll = it.first;
              int t = it.second;
              tree_ind[ll-mlag][iv] = t;
            }
        }
      catch (const ogu::Exception &e)
        {
          // ignore the error of not found variable in the tree
        }
    }

  // setup der_atoms and positions
  for (int ll = mlag; ll <= mlead; ll++)
    for (unsigned int iv = 0; iv < varnames.size(); iv++)
      {
        int t = tree_ind[ll-mlag][iv];
        if (t != -1)
          {
            der_atoms.push_back(t);
            int pos = (ll-mlag)*varnames.size() + iv;
            positions.emplace(t, pos);
          }
      }

  // set outer2y and y2outer to identities
  for (unsigned int iv = 0; iv < varnames.size(); iv++)
    {
      outer2y.push_back(iv);
      y2outer.push_back(iv);
    }

  // set n_stat, n_pred, n_both, and n_forw
  for (auto varname : varnames)
    {
      int mmlag, mmlead;
      atoms.varspan(varname, mmlead, mmlag);
      if (mmlead == 0 && mmlag == 0)
        n_stat++;
      else if (mmlead <= 0 && mmlag < 0)
        n_pred++;
      else if (mmlead > 0 && mmlag >= 0)
        n_forw++;
      else if (mmlead > 0 && mmlag < 0)
        n_both++;
      else if (mmlead < mmlag)
        // variable does not occur in the tree, cound as static
        n_stat++;
      else
        throw ogu::Exception(__FILE__, __LINE__,
                             "A wrong lag/lead of a variable in VarOrdering::do_increasing_time");
    }
}

void
VarOrdering::print() const
{
  std::cout << "nstat=" << n_stat << ", npred=" << n_pred << ", nboth=" << n_both
            << ", nforw=" << n_forw << "\n"
            << "der_atoms:\n";
  for (int der_atom : der_atoms)
    std::cout << " " << der_atom;
  std::cout << "\nmap:\n";
  for (auto position : positions)
    std::cout << " [" << position.first << u8"→" << position.second << "]";
  std::cout << "\ny2outer:\n";
  for (int i : y2outer)
    std::cout << " " <<  i;
  std::cout << "\nouter2y:\n";
  for (int i : outer2y)
    std::cout << " " << i;
  std::cout << "\n";
}
