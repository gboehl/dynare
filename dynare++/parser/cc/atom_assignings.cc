// Copyright © 2006, Ondra Kamenik

// $Id: atom_assignings.cpp 92 2007-04-19 11:38:21Z ondra $

#include "atom_assignings.hh"
#include "location.hh"
#include "parser_exception.hh"

#include "utils/cc/exception.hh"

#include <limits>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace ogp;

AtomAssignings::AtomAssignings(const AtomAssignings &aa, ogp::StaticAtoms &a)
  : atoms(a), expr(aa.expr, atoms), left_names(aa.left_names),
    lname2expr(aa.lname2expr), order(aa.order)
{
}

/** A global symbol for passing info to the AtomAssignings from
 * asgn_parse(). */
AtomAssignings *aparser;

/** The declaration of functions defined in asgn_ll.cc and asgn_tab.cc
 * generated from assign.lex assign.y */
void *asgn__scan_string(const char *);
void asgn__destroy_buffer(void *);
void asgn_parse();
extern location_type asgn_lloc;

void
AtomAssignings::parse(const string &stream)
{
  asgn_lloc.off = 0;
  asgn_lloc.ll = 0;
  void *p = asgn__scan_string(stream.c_str());
  aparser = this;
  asgn_parse();
  asgn__destroy_buffer(p);
}

void
AtomAssignings::error(string mes)
{
  throw ParserException(std::move(mes), asgn_lloc.off);
}

void
AtomAssignings::add_assignment_to_double(string name, double val)
{
  // if left hand side is a registered atom, insert it to tree
  int t;
  try
    {
      if (atoms.check(name))
        t = expr.add_nulary(name);
      else
        t = -1;
    }
  catch (const ParserException &e)
    {
      t = -1;
    }
  // register left hand side in order
  order.push_back(t);

  // add the double to the tree
  std::ostringstream buf;
  buf << std::setprecision(std::numeric_limits<double>::max_digits10)
      << val;
  try
    {
      expr.parse(buf.str());
    }
  catch (const ParserException &e)
    {
      // should never happen
      throw ParserException(string("Error parsing double ")+buf.str()+": "+e.message(), 0);
    }

  // register name of the left hand side and put to lname2expr
  left_names.insert(name);
  lname2expr.emplace(std::move(name), order.size()-1);
}

void
AtomAssignings::add_assignment(int asgn_off, const string &str, int name_len,
                               int right_off, int right_len)
{
  // the order of doing things here is important: since the
  // FormulaParser requires that all references from the i-th tree
  // refere to trees with index lass than i, so to capture also a
  // nulary term for the left hand side, it must be inserted to the
  // expression tree before the expression is parsed.

  // find the name in the atoms
  string name = str.substr(0, name_len);

  // if left hand side is a registered atom, insert it to tree
  int t;
  try
    {
      t = atoms.check(name);
      if (t == -1)
        t = expr.add_nulary(name);
    }
  catch (const ParserException &e)
    {
      atoms.register_name(name);
      t = expr.add_nulary(name);
    }
  // register left hand side in order
  order.push_back(t);

  // parse expression on the right
  try
    {
      expr.parse(str.substr(right_off, right_len));
    }
  catch (const ParserException &e)
    {
      throw ParserException(e, asgn_off+right_off);
    }

  // register name of the left hand side and put to lname2expr
  left_names.insert(name);
  if (lname2expr.find(name) != lname2expr.end())
    {
      // Prevent the occurrence of #415
      std::cerr << "Changing the value of " << name << " is not supported. Aborting." << std::endl;
      exit(EXIT_FAILURE);
    }
  lname2expr[name] = order.size()-1;
}

void
AtomAssignings::apply_subst(const AtomSubstitutions::Toldnamemap &mm)
{
  // go through all old variables and see what are their derived new
  // variables
  for (const auto & it : mm)
    {
      const string &oldname = it.first;
      const AtomSubstitutions::Tshiftnameset &sset = it.second;
      if (!sset.empty())
        {
          int told = atoms.index(oldname);
          if (told < 0 && !atoms.get_name_storage().query(oldname))
            atoms.register_name(oldname);
          if (told == -1)
            told = expr.add_nulary(oldname);
          // at least one substitution here, so make an expression
          expr.add_formula(told);
          // say that this expression is not assigned to any atom
          order.push_back(-1);
          // now go through all new names derived from the old name and
          // reference to the newly added formula
          for (const auto & itt : sset)
            {
              const string &newname = itt.first;
              left_names.insert(newname);
              lname2expr.emplace(newname, expr.nformulas()-1);
            }
        }
    }
}

void
AtomAssignings::print() const
{
  std::cout << "Atom Assignings\nExpressions:\n";
  expr.print();
  std::cout << "Left names:\n";
  for (auto it : lname2expr)
    std::cout << it.first << u8" ⇒ " << expr.formula(it.second) << " (t=" << order[it.second] << ")\n";
}

void
AtomAsgnEvaluator::setValues(EvalTree &et) const
{
  // set values of constants
  aa.atoms.setValues(et);

  // set values of variables to NaN or to user set values
  double nan = std::numeric_limits<double>::quiet_NaN();
  for (int i = 0; i < aa.atoms.nvar(); i++)
    {
      const string &ss = aa.atoms.name(i);
      int t = aa.atoms.index(ss);
      if (t >= 0)
        {
          auto it = user_values.find(t);
          if (it == user_values.end())
            et.set_nulary(t, nan);
          else
            et.set_nulary(t, it->second);
        }
    }
}

void
AtomAsgnEvaluator::set_user_value(const string &name, double val)
{
  int t = aa.atoms.index(name);
  if (t >= 0)
    {
      auto it = user_values.find(t);
      if (it == user_values.end())
        user_values.emplace(t, val);
      else
        it->second = val;
    }
}

void
AtomAsgnEvaluator::load(int i, double res)
{
  // set the value
  operator[](i) = res;
  // if i-th expression is atom, set its value to this EvalTree
  int t = aa.order[i];
  if (t >= 0)
    etree.set_nulary(t, res);
}

double
AtomAsgnEvaluator::get_value(const string &name) const
{
  auto it = aa.lname2expr.find(name);
  if (it == aa.lname2expr.end())
    return std::numeric_limits<double>::quiet_NaN();
  else
    return operator[](it->second);
}
