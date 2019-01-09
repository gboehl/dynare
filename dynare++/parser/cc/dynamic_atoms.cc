// Copyright (C) 2005, Ondra Kamenik

// $Id: dynamic_atoms.cpp 1362 2007-07-10 11:50:18Z kamenik $

#include "utils/cc/exception.hh"
#include "dynamic_atoms.hh"

#include <cstring>
#include <climits>

using namespace ogp;

NameStorage::NameStorage(const NameStorage &stor)
{
  for (auto i : stor.name_store)
    {
      auto *str = new char[strlen(i)+1];
      strcpy(str, i);
      name_store.push_back(str);
      name_set.insert(str);
    }
}

NameStorage::~NameStorage()
{
  while (name_store.size() > 0)
    {
      delete [] name_store.back();
      name_store.pop_back();
    }
}

const char *
NameStorage::query(const char *name) const
{
  auto it = name_set.find(name);
  if (it == name_set.end())
    return NULL;
  else
    return (*it);
}

const char *
NameStorage::insert(const char *name)
{
  auto it = name_set.find(name);
  if (it == name_set.end())
    {
      auto *str = new char[strlen(name)+1];
      strcpy(str, name);
      name_store.push_back(str);
      name_set.insert(str);
      return str;
    }
  else
    {
      return (*it);
    }
}

void
NameStorage::print() const
{
  for (auto i : name_store)
    printf("%s\n", i);
}

void
Constants::import_constants(const Constants &c, OperationTree &otree, Tintintmap &tmap)
{
  for (auto it : c.cmap)
    {
      int told = it.first;
      int tnew = otree.add_nulary();
      tmap.insert(Tintintmap::value_type(told, tnew));
      add_constant(tnew, it.second);
    }
}

void
Constants::setValues(EvalTree &et) const
{
  Tconstantmap::const_iterator it;
  for (it = cmap.begin(); it != cmap.end(); ++it)
    et.set_nulary((*it).first, (*it).second);
}

void
Constants::add_constant(int t, double val)
{
  cmap.insert(Tconstantmap::value_type(t, val));
  cinvmap.insert(Tconstantinvmap::value_type(val, t));
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
    return (*it).second;
  else
    {
      throw ogu::Exception(__FILE__, __LINE__,
                           "Tree index is not constant in Constants::get_constant_value");
      return 0;
    }
}

int
Constants::check(const char *str) const
{
  double d;
  sscanf(str, "%lf", &d);
  auto it = cinvmap.find(d);
  if (it != cinvmap.end())
    return (*it).second;
  else
    return -1;
}

void
Constants::print() const
{
  Tconstantmap::const_iterator it;
  for (it = cmap.begin(); it != cmap.end(); ++it)
    printf("$%d:  %8.4g\n", (*it).first, (*it).second);
}

DynamicAtoms::DynamicAtoms()
  : nv(0), minlag(INT_MAX), maxlead(INT_MIN)
{
}

DynamicAtoms::DynamicAtoms(const DynamicAtoms &da)
  : Constants(da),
    varnames(da.varnames), vars(), indices(),
    nv(da.nv), minlag(da.minlag), maxlead(da.maxlead)
{
  // copy vars
  for (const auto & var : da.vars)
    vars.insert(Tvarmap::value_type(varnames.query(var.first),
                                    var.second));
  // copy indices
  for (auto indice : da.indices)
    indices.insert(Tindexmap::value_type(indice.first,
                                         varnames.query(indice.second)));
}

int
DynamicAtoms::check(const char *name) const
{
  if (is_string_constant(name))
    return Constants::check(name);

  return check_variable(name);
}

int
DynamicAtoms::check_variable(const char *name) const
{
  string str;
  int ll;
  parse_variable(name, str, ll);
  auto it = vars.find(str.c_str());

  if (it != vars.end())
    {
      const Tlagmap &lmap = (*it).second;
      auto itt = lmap.find(ll);
      if (itt != lmap.end())
        return (*itt).second;
    }
  return -1;
}

void
DynamicAtoms::assign(const char *name, int t)
{
  if (is_string_constant(name))
    assign_constant(name, t);
  else
    assign_variable(name, t);
}

void
DynamicAtoms::assign_constant(const char *name, int t)
{
  double val;
  sscanf(name, "%lf", &val);
  add_constant(t, val);
}

// parse the name and then call assing_variable(varname, ll, t)

void
DynamicAtoms::assign_variable(const char *name, int t)
{
  int ll;
  string str;
  parse_variable(name, str, ll);
  // here str is just name without lead/lag
  const char *ss = varnames.insert(str.c_str());

  assign_variable(ss, ll, t);
}

void
DynamicAtoms::assign_variable(const char *varname, int ll, int t)
{
  if (indices.end() != indices.find(t))
    throw ogu::Exception(__FILE__, __LINE__,
                         "Attempt to assign already allocated tree index");

  auto it = vars.find(varname);
  if (it != vars.end())
    {
      Tlagmap &lmap = (*it).second;
      if (lmap.end() != lmap.find(ll))
        throw ogu::Exception(__FILE__, __LINE__,
                             "Attempt to assign already allocated variable");
      lmap.insert(Tlagmap::value_type(ll, t));
    }
  else
    {
      Tlagmap lmap;
      lmap.insert(Tlagmap::value_type(ll, t));
      vars.insert(Tvarmap::value_type(varname, lmap));
    }
  indices.insert(Tindexmap::value_type(t, varname));

  nv++;
  if (ll < minlag)
    minlag = ll;
  if (ll > maxlead)
    maxlead = ll;
}

void
DynamicAtoms::unassign_variable(const char *varname, int ll, int t)
{
  auto it = vars.find(varname);
  if (it != vars.end())
    {
      Tlagmap &lmap = (*it).second;
      auto itt = lmap.find(ll);
      if (itt != lmap.end())
        {
          if ((*itt).second == t)
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
  minlag = INT_MAX;
  maxlead = INT_MIN;
  for (Tvarmap::const_iterator it = vars.begin(); it != vars.end(); ++it)
    {
      const Tlagmap &lmap = (*it).second;
      for (auto itt : lmap)
        {
          int ll = itt.first;
          if (ll < minlag)
            minlag = ll;
          if (ll > maxlead)
            maxlead = ll;
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
      mlead = INT_MIN;
      mlag = INT_MAX;
      return;
    }
  varspan((*it).second, mlead, mlag);
}

void
DynamicAtoms::varspan(const char *name, int &mlead, int &mlag) const
{
  auto it = vars.find(name);
  if (vars.end() == it)
    {
      mlead = INT_MIN;
      mlag = INT_MAX;
      return;
    }
  const Tlagmap &lmap = (*it).second;
  auto beg = lmap.begin();
  auto end = lmap.rbegin();
  mlag = (*beg).first;
  mlead = (*end).first;
}

void
DynamicAtoms::varspan(const vector<const char *> &names, int &mlead, int &mlag) const
{
  mlead = INT_MIN;
  mlag = INT_MAX;
  for (auto name : names)
    {
      int lag, lead;
      varspan(name, lead, lag);
      if (lead > mlead)
        mlead = lead;
      if (lag < mlag)
        mlag = lag;
    }
}

bool
DynamicAtoms::is_named_atom(int t) const
{
  return (indices.end() != indices.find(t));
}

int
DynamicAtoms::index(const char *name, int ll) const
{
  auto it = vars.find(name);
  if (vars.end() != it)
    {
      const Tlagmap &lmap = (*it).second;
      auto itt = lmap.find(ll);
      if (lmap.end() != itt)
        return (*itt).second;
    }
  return -1;
}

bool
DynamicAtoms::is_referenced(const char *name) const
{
  auto it = vars.find(name);
  return it != vars.end();
}

const DynamicAtoms::Tlagmap &
DynamicAtoms::lagmap(const char *name) const
{
  auto it = vars.find(name);
  if (vars.end() == it)
    throw ogu::Exception(__FILE__, __LINE__,
                         std::string("Couldn't find the name ")
                         + name + " in DynamicAtoms::lagmap");
  return (*it).second;
}

const char *
DynamicAtoms::name(int t) const
{
  auto it = indices.find(t);
  if (indices.end() == it)
    throw ogu::Exception(__FILE__, __LINE__,
                         "Couldn't find tree index in DynamicAtoms::name");
  return (*it).second;
}

int
DynamicAtoms::lead(int t) const
{
  const char *nam = name(t);
  const Tlagmap &lmap = lagmap(nam);
  auto it = lmap.begin();
  while (it != lmap.end() && (*it).second != t)
    ++it;
  if (lmap.end() == it)
    throw ogu::Exception(__FILE__, __LINE__,
                         "Couldn't find the three index in DynamicAtoms::lead");
  return (*it).first;
}

void
DynamicAtoms::print() const
{
  printf("names:\n");
  varnames.print();
  printf("constants:\n");
  Constants::print();
  printf("variables:\n");
  for (const auto & var : vars)
    {
      const Tlagmap &lmap = var.second;
      for (auto itt : lmap)
        printf("$%d: %s(%d)\n", itt.second, var.first, itt.first);
    }
  printf("indices:\n");
  for (auto indice : indices)
    printf("t=%d ==> %s\n", indice.first, indice.second);
}

/** Note that the str has been parsed by the lexicographic
 * analyzer. It can be either a variable or a double. So it is easy to
 * recognize it by the first character. */
bool
DynamicAtoms::is_string_constant(const char *str)
{
  return str[0] == '.' || str[0] == '-' || (str[0] >= '0' && str[0] <= '9');
}

VarOrdering::VarOrdering(const VarOrdering &vo, const vector<const char *> &vnames,
                         const DynamicAtoms &a)
  : n_stat(vo.n_stat), n_pred(vo.n_pred), n_both(vo.n_both), n_forw(vo.n_forw),
    der_atoms(vo.der_atoms), positions(vo.positions),
    outer2y(vo.outer2y), y2outer(vo.y2outer), varnames(vnames), atoms(a)
{
}

bool
VarOrdering::check(int t) const
{
  auto it = positions.find(t);
  return it != positions.end();
}

int
VarOrdering::get_pos_of(int t) const
{
  auto it = positions.find(t);
  if (it != positions.end())
    {
      return (*it).second;
    }
  else
    {
      throw ogu::Exception(__FILE__, __LINE__,
                           "Couldn't find the tree index in VarOrdering::get_pos_of");
      return -1;
    }
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
      const char *ss = varnames[i];
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
        {
          throw ogu::Exception(__FILE__, __LINE__,
                               "A wrong lag/lead of a variable in VarOrdering::do_pbspbfbf");
        }
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
  else // BEWARE: when implementing a new ordering, check also a
    {// code below setting y2outer
      throw ogu::Exception(__FILE__, __LINE__,
                           "Ordering not implemented in VarOrdering::do_general");
    }

  // make der_atoms and positions
  int off = 0;
  for (auto & ord : ords)
    for (unsigned int j = 0; j < ord->size(); j++, off++)
      if ((*ord)[j] != -1)
        {
          der_atoms.push_back((*ord)[j]);
          positions.insert(std::pair<int, int>((*ord)[j], off));
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
  vector<vector<int> > tree_ind(mlead-mlag+1, ll_init);
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
            positions.insert(map<int, int>::value_type(t, pos));
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
        {
          n_stat++;
        }
      else if (mmlead <= 0 && mmlag < 0)
        {
          n_pred++;
        }
      else if (mmlead > 0 && mmlag >= 0)
        {
          n_forw++;
        }
      else if (mmlead > 0 && mmlag < 0)
        {
          n_both++;
        }
      else if (mmlead < mmlag)
        {
          // variable does not occur in the tree, cound as static
          n_stat++;
        }
      else
        {
          throw ogu::Exception(__FILE__, __LINE__,
                               "A wrong lag/lead of a variable in VarOrdering::do_increasing_time");
        }
    }
}

void
VarOrdering::print() const
{
  printf("nstat=%d, npred=%d, nboth=%d, nforw=%d\n", n_stat, n_pred, n_both, n_forw);
  printf("der_atoms:\n");
  for (int der_atom : der_atoms)
    printf(" %d", der_atom);
  printf("\nmap:\n");
  for (auto position : positions)
    printf(" [%d->%d]", position.first, position.second);
  printf("\ny2outer:\n");
  for (int i : y2outer)
    printf(" %d", i);
  printf("\nouter2y:\n");
  for (int i : outer2y)
    printf(" %d", i);
  printf("\n");
}

// Local Variables:
// mode:C++
// End:
