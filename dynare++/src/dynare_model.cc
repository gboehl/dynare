/*
 * Copyright © 2006-2011 Ondra Kamenik
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

#include "parser/cc/parser_exception.hh"
#include "parser/cc/location.hh"
#include "utils/cc/exception.hh"
#include "dynare_model.hh"
#include "dynare_exception.hh"
#include "planner_builder.hh"
#include "forw_subst_builder.hh"

#include <string>
#include <cmath>
#include <limits>
#include <ostream>
#include <memory>
#include <algorithm>
#include <iomanip>

using namespace ogdyn;

ParsedMatrix::ParsedMatrix(const ogp::MatrixParser &mp)
  : TwoDMatrix(mp.nrows(), mp.ncols())
{
  zeros();
  for (ogp::MPIterator it = mp.begin(); it != mp.end(); ++it)
    get(it.row(), it.col()) = *it;
}

DynareModel::DynareModel()
  : atoms(), eqs(atoms)
{
}

DynareModel::DynareModel(const DynareModel &dm)
  : atoms(dm.atoms), eqs(dm.eqs, atoms), order(dm.order),
    t_plobjective(dm.t_plobjective),
    t_pldiscount(dm.t_pldiscount)
{
  if (dm.param_vals)
    param_vals = std::make_unique<Vector>(const_cast<const Vector &>(*dm.param_vals));
  if (dm.init_vals)
    init_vals = std::make_unique<Vector>(const_cast<const Vector &>(*dm.init_vals));
  if (dm.vcov_mat)
    vcov_mat = std::make_unique<TwoDMatrix>(const_cast<const TwoDMatrix &>(*dm.vcov_mat));
  if (dm.old_atoms)
    old_atoms = std::make_unique<DynareDynamicAtoms>(static_cast<const DynareDynamicAtoms &>(*dm.old_atoms));
  if (dm.atom_substs)
    atom_substs = std::make_unique<ogp::AtomSubstitutions>(*dm.atom_substs, *old_atoms, atoms);
  if (dm.pbuilder)
    pbuilder = std::make_unique<PlannerBuilder>(*dm.pbuilder, *this);
  if (dm.fbuilder)
    fbuilder = std::make_unique<ForwSubstBuilder>(*dm.fbuilder, *this);
}

const PlannerInfo *
DynareModel::get_planner_info() const
{
  if (pbuilder)
    return &(pbuilder->get_info());
  return nullptr;
}

const ForwSubstInfo *
DynareModel::get_forw_subst_info() const
{
  if (fbuilder)
    return &(fbuilder->get_info());
  return nullptr;
}

const ogp::SubstInfo *
DynareModel::get_subst_info() const
{
  if (atom_substs)
    return &(atom_substs->get_info());
  return nullptr;
}

void
DynareModel::setInitOuter(const Vector &x)
{
  if (x.length() != atoms.ny())
    throw DynareException(__FILE__, __LINE__,
                          "Wrong length of vector in DynareModel::setInitOuter");
  for (int i = 0; i < atoms.ny(); i++)
    (*init_vals)[i] = x[atoms.y2outer_endo()[i]];
}

void
DynareModel::print() const
{
  std::cout << "all atoms:\n";
  atoms.print();
  std::cout << "formulas:\n";
  DebugOperationFormatter dof(*this);
  for (int i = 0; i < eqs.nformulas(); i++)
    {
      int tf = eqs.formula(i);
      std::cout << "formula " << tf << "%d:\n";
      eqs.getTree().print_operation_tree(tf, std::cout, dof);
    }
}

void
DynareModel::dump_model(std::ostream &os) const
{
  // endogenous variable declaration
  os << "var";
  for (auto i : atoms.get_endovars())
    os << " " << i;
  os << ";\n\n";

  // exogenous variables
  os << "varexo";
  for (auto i : atoms.get_exovars())
    os << " " << i;
  os << ";\n\n";

  // parameters
  os << "parameters";
  for (auto i : atoms.get_params())
    os << " " << i;
  os << ";\n\n";

  // parameter values
  os.precision(16);
  for (int i = 0; i < static_cast<int>(atoms.get_params().size()); i++)
    os << atoms.get_params()[i] << "=" << getParams()[i] << ";\n";
  os << "\n\n";

  // model section
  ogp::OperationStringConvertor osc(atoms, getParser().getTree());
  os << "model;\n";
  for (int i = 0; i < getParser().nformulas(); i++)
    {
      os << "// Equation " << i << "\n0 = ";
      int t = getParser().formula(i);
      os << osc.convert(getParser().getTree().operation(t), t);
      os << ";\n";
    }
  os << "end;\n";

  // initval as steady state
  os << "initval;\n";
  for (int i = 0; i < static_cast<int>(atoms.get_endovars().size()); i++)
    os << atoms.get_endovars()[atoms.y2outer_endo()[i]] << "=" << getInit()[i] << ";\n";
  os << "end;\n";
}

void
DynareModel::add_name(std::string name, int flag)
{
  if (flag == 1)
    // endogenous
    atoms.register_uniq_endo(name);
  else if (flag == 2)
    // exogenous
    atoms.register_uniq_exo(name);
  else if (flag == 3)
    // parameter
    atoms.register_uniq_param(name);
  else
    throw DynareException(__FILE__, __LINE__, "Unrecognized flag value.");
}

void
DynareModel::check_model() const
{
  if (order == -1)
    throw DynareException(__FILE__, __LINE__,
                          "Order of approximation not set in DynareModel::check_model");

  if (atoms.ny() != eqs.nformulas())
    throw DynareException(__FILE__, __LINE__, "Model has " + std::to_string(eqs.nformulas())
                          + " equations for " + std::to_string(atoms.ny()) + " endogenous variables");

  /* check whether all nulary terms of all formulas in eqs are either constant
     or assigned to a name */
  for (int i = 0; i < eqs.nformulas(); i++)
    {
      int ft = eqs.formula(i);
      const unordered_set<int> &nuls = eqs.nulary_of_term(ft);
      for (int nul : nuls)
        if (!atoms.is_constant(nul) && !atoms.is_named_atom(nul))
          throw DynareException(__FILE__, __LINE__,
                                "Dangling nulary term found, internal error.");
    }

  int mlag, mlead;
  atoms.exovarspan(mlead, mlag);
  if (atoms.nexo() > 0 && (mlead != 0 || mlag != 0))
    throw DynareException(__FILE__, __LINE__,
                          "The model contains occurrences of lagged/leaded exogenous variables");

  atoms.endovarspan(mlead, mlag);
  if (mlead > 1 || mlag < -1)
    throw DynareException(__FILE__, __LINE__,
                          "The model contains occurrences of too lagged/leaded endogenous variables");

  // check the dimension of vcov matrix
  if (getAtoms().nexo() != getVcov().nrows())
    throw DynareException(__FILE__, __LINE__,
                          "Dimension of VCOV matrix does not correspond to the shocks");
}

int
DynareModel::variable_shift(int t, int tshift)
{
  const string &name = atoms.name(t);
  if (atoms.is_type(name, DynareDynamicAtoms::atype::param)
      || atoms.is_constant(t))
    throw DynareException(__FILE__, __LINE__,
                          "The tree index is not a variable in DynareModel::variable_shift");
  int ll = atoms.lead(t) + tshift;
  int res = atoms.index(name, ll);
  if (res == -1)
    res = eqs.add_nulary(name + '(' + std::to_string(ll) + ')');
  return res;
}

void
DynareModel::variable_shift_map(const unordered_set<int> &a_set, int tshift,
                                map<int, int> &s_map)
{
  s_map.clear();
  for (int t : a_set)
    // make shift map only for non-constants and non-parameters
    if (!atoms.is_constant(t))
      {
        const string &name = atoms.name(t);
        if (atoms.is_type(name, DynareDynamicAtoms::atype::endovar)
            || atoms.is_type(name, DynareDynamicAtoms::atype::exovar))
          {
            int tt = variable_shift(t, tshift);
            s_map.emplace(t, tt);
          }
      }
}

void
DynareModel::termspan(int t, int &mlead, int &mlag) const
{
  mlead = std::numeric_limits<int>::min();
  mlag = std::numeric_limits<int>::max();
  const unordered_set<int> &nul_terms = eqs.nulary_of_term(t);
  for (int nul_term : nul_terms)
    {
      if (!atoms.is_constant(nul_term)
          && (atoms.is_type(atoms.name(nul_term), DynareDynamicAtoms::atype::endovar)
              || atoms.is_type(atoms.name(nul_term), DynareDynamicAtoms::atype::exovar)))
        {
          int ll = atoms.lead(nul_term);
          mlag = std::min(ll, mlag);
          mlead = std::max(ll, mlead);
        }
    }
}

bool
DynareModel::is_constant_term(int t) const
{
  const unordered_set<int> &nul_terms = eqs.nulary_of_term(t);
  for (int nul_term : nul_terms)
    if (!atoms.is_constant(nul_term)
        && !atoms.is_type(atoms.name(nul_term), DynareDynamicAtoms::atype::param))
      return false;
  return true;
}

unordered_set<int>
DynareModel::get_nonlinear_subterms(int t) const
{
  NLSelector nls(*this);
  return eqs.getTree().select_terms(t, nls);
}

void
DynareModel::substitute_atom_for_term(const string &name, int ll, int t)
{
  /* if the term t is itself a named atom (parameter, exo, endo), then we have
     to unassign it first */
  if (atoms.is_named_atom(t))
    atoms.unassign_variable(atoms.name(t), atoms.lead(t), t);
  /* assign allocated tree index for the term now to name(ll) */
  atoms.assign_variable(name, ll, t);
  // make operation t nulary in operation tree
  eqs.nularify(t);
}

void
DynareModel::final_job()
{
  if (t_plobjective != -1 && t_pldiscount != -1)
    {
      /* at this moment include all equations and all variables; in future we
         will exclude purely exogenous processes; TODO */
      PlannerBuilder::Tvarset vset;
      for (int i = 0; i < atoms.ny(); i++)
        vset.insert(atoms.get_endovars()[i]);
      PlannerBuilder::Teqset eset;
      for (int i = 0; i < eqs.nformulas(); i++)
        eset.push_back(i);

      // construct the planner builder, this adds a lot of stuff to the model
      pbuilder = std::make_unique<PlannerBuilder>(*this, vset, eset);
    }

  // construct ForwSubstBuilder
  fbuilder = std::make_unique<ForwSubstBuilder>(*this);

  // call parsing_finished (this will define an outer ordering of all variables)
  atoms.parsing_finished(ogp::VarOrdering::bfspbfpb);
  // make a copy of atoms and name it old_atoms
  old_atoms = std::make_unique<DynareDynamicAtoms>(atoms);
  // construct empty substitutions from old_atoms to atoms
  atom_substs = std::make_unique<ogp::AtomSubstitutions>(*old_atoms, atoms);
  /* do the actual substitution, it will also call parsing_finished for atoms
     which creates internal orderings */
  atoms.substituteAllLagsAndExo1Leads(eqs, *atom_substs);
}

extern ogp::location_type dynglob_lloc;

DynareParser::DynareParser(const string &stream, int ord)
  : DynareModel(),
    pa_atoms(), paramset(pa_atoms),
    ia_atoms(), initval(ia_atoms), vcov(),
    model_beg(0), model_end(-1),
    paramset_beg(0), paramset_end(-1),
    initval_beg(0), initval_end(-1),
    vcov_beg(0), vcov_end(-1),
    order_beg(0), order_end(-1),
    plobjective_beg(0), plobjective_end(-1),
    pldiscount_beg(0), pldiscount_end(-1)
{
  // global parse
  try
    {
      parse_glob(stream);
    }
  catch (const ogp::ParserException &e)
    {
      throw ogp::ParserException(e, dynglob_lloc.off);
    }
  // setting parameters parse
  try
    {
      if (paramset_end > paramset_beg)
        paramset.parse(stream.substr(paramset_beg, paramset_end-paramset_beg));
    }
  catch (const ogp::ParserException &e)
    {
      throw ogp::ParserException(e, paramset_beg);
    }
  // model parse
  try
    {
      if (model_end > model_beg)
        eqs.parse(stream.substr(model_beg, model_end-model_beg));
      else
        throw ogp::ParserException("Model section not found.", 0);
    }
  catch (const ogp::ParserException &e)
    {
      throw ogp::ParserException(e, model_beg);
    }
  // initval setting parse
  try
    {
      if (initval_end > initval_beg)
        initval.parse(stream.substr(initval_beg, initval_end-initval_beg));
    }
  catch (const ogp::ParserException &e)
    {
      throw ogp::ParserException(e, initval_beg);
    }
  // vcov parse
  try
    {
      if (vcov_end > vcov_beg)
        vcov.parse(stream.substr(vcov_beg, vcov_end-vcov_beg));
    }
  catch (const ogp::ParserException &e)
    {
      throw ogp::ParserException(e, vcov_beg);
    }
  // planner objective parse
  try
    {
      if (plobjective_end > plobjective_beg)
        {
          eqs.parse(stream.substr(plobjective_beg, plobjective_end-plobjective_beg));
          t_plobjective = eqs.pop_last_formula();
        }
    }
  catch (const ogp::ParserException &e)
    {
      throw ogp::ParserException(e, plobjective_beg);
    }
  // planner discount parse
  try
    {
      if (pldiscount_end > pldiscount_beg)
        t_pldiscount = parse_pldiscount(stream.substr(pldiscount_beg, pldiscount_end - pldiscount_beg));
    }
  catch (const ogp::ParserException &e)
    {
      throw ogp::ParserException(e, pldiscount_beg);
    }
  // order parse
  try
    {
      if (order_end > order_beg)
        order = parse_order(stream.substr(order_beg, order_end - order_beg));
    }
  catch (const ogp::ParserException &e)
    {
      throw ogp::ParserException(e, order_beg);
    }

  // check the overridden order
  if (ord != -1)
    order = ord;

  // end parsing job, add planner's FOCs, make substitutions
  DynareModel::final_job();

  // calculate parameters
  calc_params();
  // calculate initial values
  calc_init();

  if (vcov_end > vcov_beg)
    vcov_mat = std::make_unique<ParsedMatrix>(vcov);
  else
    {
      // vcov has not been asserted, set it to unit matrix
      vcov_mat = std::make_unique<TwoDMatrix>(atoms.nexo(), atoms.nexo());
      vcov_mat->unit();
    }

  // check the model
  check_model();

  // differentiate
  if (order >= 1)
    eqs.differentiate(order);
}

DynareParser::DynareParser(const DynareParser &dp)
  : DynareModel(dp),
    pa_atoms(dp.pa_atoms), paramset(dp.paramset, pa_atoms),
    ia_atoms(dp.ia_atoms), initval(dp.initval, ia_atoms), vcov(dp.vcov),
    model_beg(dp.model_beg), model_end(dp.model_end),
    paramset_beg(dp.paramset_beg), paramset_end(dp.paramset_end),
    initval_beg(dp.initval_beg), initval_end(dp.initval_end),
    vcov_beg(dp.vcov_beg), vcov_end(dp.vcov_end),
    order_beg(dp.order_beg), order_end(dp.order_end),
    plobjective_beg(dp.plobjective_beg), plobjective_end(dp.plobjective_end),
    pldiscount_beg(dp.pldiscount_beg), pldiscount_end(dp.pldiscount_end)
{
}

void
DynareParser::add_name(string name, int flag)
{
  DynareModel::add_name(name, flag);
  // register with static atoms used for atom assignements
  if (flag == 1)
    // endogenous
    ia_atoms.register_name(std::move(name));
  else if (flag == 2)
    // exogenous
    ia_atoms.register_name(std::move(name));
  else if (flag == 3)
    {
      // parameter
      pa_atoms.register_name(name);
      ia_atoms.register_name(std::move(name));
    }
  else
    throw DynareException(__FILE__, __LINE__, "Unrecognized flag value.");
}

void
DynareParser::error(string mes)
{
  // throwing zero offset since this exception will be caugth at constructor
  throw ogp::ParserException(std::move(mes), 0);
}

void
DynareParser::print() const
{
  DynareModel::print();
  std::cout << "parameter atoms:\n";
  paramset.print();
  std::cout << "initval atoms:\n";
  initval.print();
  std::cout << "model position: " << model_beg << ' ' << model_end << '\n'
            << "paramset position: " << paramset_beg << ' ' << paramset_end << '\n'
            << "initval position: " << initval_beg << ' ' << initval_end << '\n';
}

/* A global symbol for passing info to the DynareParser from parser. */
DynareParser *dynare_parser;

/* The declarations of functions defined in dynglob_ll.cc and dynglob_tab.cc
   generated from dynglob.lex and dynglob.y */
void *dynglob__scan_string(const char *);
void dynglob__destroy_buffer(void *);
void dynglob_parse();
extern ogp::location_type dynglob_lloc;

void
DynareParser::parse_glob(const string &stream)
{
  void *p = dynglob__scan_string(stream.c_str());
  dynare_parser = this;
  dynglob_parse();
  dynglob__destroy_buffer(p);
}

int
DynareParser::parse_order(const string &str)
{
  return std::stoi(str);
}

int
DynareParser::parse_pldiscount(const string &str)
{
  if (!atoms.is_type(str, DynareDynamicAtoms::atype::param))
    throw ogp::ParserException(std::string{"Name "} + str + " is not a parameter", 0);

  int t = atoms.index(str, 0);
  if (t == -1)
    t = eqs.add_nulary(str);

  return t;
}

void
DynareParser::calc_params()
{
  param_vals = std::make_unique<Vector>(atoms.np());
  ogp::AtomAsgnEvaluator aae(paramset);
  aae.eval();
  for (int i = 0; i < atoms.np(); i++)
    (*param_vals)[i] = aae.get_value(atoms.get_params()[i]);

  for (unsigned int i = 0; i < atoms.get_params().size(); i++)
    if (!std::isfinite((*param_vals)[i]))
      std::cout << "dynare++: warning: value for parameter " << atoms.get_params()[i] << " is not finite\n";
}

void
DynareParser::calc_init()
{
  // update initval atoms assignings according to substitutions
  if (atom_substs)
    initval.apply_subst(atom_substs->get_old2new());

  // calculate the vector of initial values
  init_vals = std::make_unique<Vector>(atoms.ny());
  ogp::AtomAsgnEvaluator aae(initval);
  // set parameters
  for (int ip = 0; ip < atoms.np(); ip++)
    aae.set_user_value(atoms.get_params()[ip], (*param_vals)[ip]);
  // set exogenous to zeros
  for (int ie = 0; ie < atoms.nexo(); ie++)
    aae.set_user_value(atoms.get_exovars()[ie], 0.0);
  // evaluate
  aae.eval();
  // set results to internally ordered vector init_vals
  for (int outer = 0; outer < atoms.ny(); outer++)
    {
      int i = atoms.outer2y_endo()[outer];
      (*init_vals)[i] = aae.get_value(atoms.get_endovars()[outer]);
    }

  /* if the planner's FOCs have been added, then add estimate of Lagrange
     multipliers to the vector */
  if (pbuilder)
    MultInitSS mis(*pbuilder, *param_vals, *init_vals);

  /* if forward substitution builder has been created, we have to its
     substitutions and evaluate them */
  if (fbuilder)
    ogdyn::DynareSteadySubstitutions dss(atoms, eqs.getTree(),
                                         fbuilder->get_aux_map(), *param_vals, *init_vals);

  for (unsigned int i = 0; i < atoms.get_endovars().size(); i++)
    if (!std::isfinite((*init_vals)[i]))
      std::cout << "dynare++: warning: initval for <" << atoms.get_endovars()[atoms.y2outer_endo()[i]] << "> is not finite\n";
}

// this returns false for linear functions
bool
NLSelector::operator()(int t) const
{
  const ogp::Operation &op = model.getParser().getTree().operation(t);
  const DynareDynamicAtoms &atoms = model.getAtoms();
  // if the term is constant, return false
  if (model.is_constant_term(t))
    return false;
  int nary = op.nary();
  if (nary == 0)
    {
      if (atoms.is_type(atoms.name(t), DynareDynamicAtoms::atype::endovar)
          || atoms.is_type(atoms.name(t), DynareDynamicAtoms::atype::exovar))
        return true;
      else
        return false;
    }
  else if (nary == 1)
    {
      if (op.getCode() == ogp::code_t::UMINUS)
        return false;
      else
        return true;
    }
  else
    {
      if (op.getCode() == ogp::code_t::TIMES)
        // if at least one operand is constant, than the TIMES is linear
        if (model.is_constant_term(op.getOp1())
            || model.is_constant_term(op.getOp2()))
          return false;
        else
          return true;
      // both PLUS and MINUS are linear
      if (op.getCode() == ogp::code_t::PLUS
          || op.getCode() == ogp::code_t::MINUS)
        return false;
      // POWER is linear if exponent or base is 0 or one
      if (op.getCode() == ogp::code_t::POWER
          && (op.getOp1() == ogp::OperationTree::zero
              || op.getOp1() == ogp::OperationTree::one
              || op.getOp2() == ogp::OperationTree::zero
              || op.getOp2() == ogp::OperationTree::one))
        return false;
      else
        return true;
      /* DIVIDE is linear if the denominator is constant, or if the nominator
         is zero */
      if (op.getCode() == ogp::code_t::DIVIDE
          && (op.getOp1() == ogp::OperationTree::zero
              || model.is_constant_term(op.getOp2())))
        return false;
      else
        return true;
    }

  throw DynareException(__FILE__, __LINE__,
                        "Wrong operation in operation tree");
  return false;
}

DynareSPModel::DynareSPModel(const std::vector<std::string> &endo,
                             const std::vector<std::string> &exo,
                             const std::vector<std::string> &par,
                             const string &equations,
                             int ord)
  : DynareModel()
{
  // set the order
  order = ord;

  // add names
  for (const auto &it : endo)
    add_name(it, 1);
  for (const auto &it : exo)
    add_name(it, 2);
  for (const auto &it : par)
    add_name(it, 3);

  // parse the equations
  eqs.parse(equations);

  // parsing finished
  atoms.parsing_finished(ogp::VarOrdering::bfspbfpb);

  // create what has to be created from DynareModel
  param_vals = std::make_unique<Vector>(atoms.np());
  init_vals = std::make_unique<Vector>(atoms.ny());
  vcov_mat = std::make_unique<TwoDMatrix>(atoms.nexo(), atoms.nexo());

  // check the model
  check_model();

  // differentiate
  if (order >= 1)
    eqs.differentiate(order);
}

void
ModelSSWriter::write_der0(std::ostream &os)
{
  write_der0_preamble(os);
  write_atom_assignment(os);

  stop_set.clear();
  for (int fi = 0; fi < model.eqs.nformulas(); fi++)
    otree.print_operation_tree(model.eqs.formula(fi), os, *this);

  write_der0_assignment(os);
}

void
ModelSSWriter::write_der1(std::ostream &os)
{
  write_der1_preamble(os);
  write_atom_assignment(os);

  stop_set.clear();

  const vector<int> &variables = model.getAtoms().variables();
  const vector<int> &eam = model.getAtoms().get_endo_atoms_map();
  for (int i = 0; i < model.getParser().nformulas(); i++)
    {
      const ogp::FormulaDerivatives &fder = model.getParser().derivatives(i);
      for (int j : eam)
        {
          int t = fder.derivative(ogp::FoldMultiIndex(variables.size(), 1, j));
          if (t > 0)
            otree.print_operation_tree(t, os, *this);
        }
    }

  write_der1_assignment(os);
}

MatlabSSWriter::MatlabSSWriter(const DynareModel &dm, std::string id_arg)
  : ModelSSWriter(dm), id(std::move(id_arg))
{
}

void
MatlabSSWriter::write_der0_preamble(std::ostream &os) const
{
  os << "% Usage:\n"
     << "%       out = " << id << "_f(params, y)\n"
     << "%   where\n"
     << "%       out    is a (" << model.getAtoms().ny() << ",1) column vector of the residuals\n"
     << "%              of the static system\n";
  write_common1_preamble(os);
  os << "function out = " << id << "_f(params, y)\n";
  write_common2_preamble(os);
}

void
MatlabSSWriter::write_der1_preamble(std::ostream &os) const
{
  os << "% Usage:\n"
     << "%       out = " << id << "_ff(params, y)\n"
     << "%   where\n"
     << "%       out    is a (" << model.getAtoms().ny() << "," << model.getAtoms().ny() << ") matrix of the first order\n"
     << "%              derivatives of the static system residuals\n"
     << "%              columns correspond to endo variables in\n"
     << "%              the ordering as declared\n";
  write_common1_preamble(os);
  os << "function out = " << id << "_ff(params, y)\n";
  write_common2_preamble(os);
}

void
MatlabSSWriter::write_common1_preamble(std::ostream &os) const
{
  os << "%       params is a (" << model.getAtoms().np() << ",1) vector of parameter values\n"
     << "%              in the ordering as declared\n"
     << "%       y      is a (" << model.getAtoms().ny() << ",1) vector of endogenous variables\n"
     << "%              in the ordering as declared\n"
     << "%\n"
     << "% Created by Dynare++ v. " << VERSION << "\n";
  // write ordering of parameters
  os << "\n% params ordering\n% =====================\n";
  for (auto parname : model.getAtoms().get_params())
    os << "% " << parname << "\n";

  // write endogenous variables
  os << "%\n% y ordering\n% =====================\n";
  for (auto endoname : model.getAtoms().get_endovars())
    os << "% " << endoname << "\n";
  os << "\n";
}

void
MatlabSSWriter::write_common2_preamble(std::ostream &os) const
{
  os << "if size(y) ~= [" << model.getAtoms().ny() << ",1]\n"
     << "\terror('Wrong size of y, must be [" << model.getAtoms().ny() << ",1]');\nend\n"
     << "if size(params) ~= [" << model.getAtoms().np() << ",1]\n"
     << "\terror('Wrong size of params, must be [" << model.getAtoms().np() << ",1]');\nend\n\n";
}

void
MatlabSSWriter::write_atom_assignment(std::ostream &os) const
{
  // write OperationTree::num_constants
  os << "% hardwired constants\n";
  ogp::EvalTree etree(model.getParser().getTree(), ogp::OperationTree::num_constants-1);
  for (int i = 0; i < ogp::OperationTree::num_constants; i++)
    {
      format_nulary(i, os);
      double g = etree.eval(i);
      if (std::isnan(g))
        os << " = NaN;\n";
      else
        os << " = " << std::defaultfloat << std::setprecision(8) << etree.eval(i) << ";\n";
    }
  // write numerical constants
  os << "% numerical constants\n";
  const ogp::Constants::Tconstantmap &cmap = model.getAtoms().get_constantmap();
  for (auto it : cmap)
    {
      format_nulary(it.first, os);
      os << " = " << std::defaultfloat << std::setprecision(8) << it.second << ";\n";
    }
  // write parameters
  os << "% parameter values\n";
  for (unsigned int ip = 0; ip < model.getAtoms().get_params().size(); ip++)
    {
      const string &parname = model.getAtoms().get_params()[ip];
      int t = model.getAtoms().index(parname, 0);
      if (t == -1)
        os << "% " << parname << " not used in the model\n";
      else
        {
          format_nulary(t, os);
          os << " = params(" << ip+1 << "); % " << parname << "\n";
        }
    }
  // write exogenous variables
  os << "% exogenous variables to zeros\n";
  for (unsigned int ie = 0; ie < model.getAtoms().get_exovars().size(); ie++)
    {
      const string &exoname = model.getAtoms().get_exovars()[ie];
      try
        {
          const ogp::DynamicAtoms::Tlagmap &lmap = model.getAtoms().lagmap(exoname);
          for (auto it : lmap)
            {
              format_nulary(it.second, os);
              os << " = 0.0; % " << exoname << "\n";
            }
        }
      catch (const ogu::Exception &e)
        {
          // ignore the error of not found variable in the tree
        }
    }
  // write endogenous variables
  os << "% endogenous variables to y\n";
  for (unsigned int ie = 0; ie < model.getAtoms().get_endovars().size(); ie++)
    {
      const string &endoname = model.getAtoms().get_endovars()[ie];
      const ogp::DynamicAtoms::Tlagmap &lmap = model.getAtoms().lagmap(endoname);
      for (auto it : lmap)
        {
          format_nulary(it.second, os);
          os << " = y(" << ie+1 << "); % " << endoname << "\n";
        }
    }
  os << "\n";
}

void
MatlabSSWriter::write_der0_assignment(std::ostream &os) const
{

  // initialize out variable
  os << "% setting the output variable\n"
     << "out = zeros(" << model.getParser().nformulas() << ", 1);\n";

  // fill out with the terms
  for (int i = 0; i < model.getParser().nformulas(); i++)
    {
      os << "out(" << i+1 << ") = ";
      format_term(model.getParser().formula(i), os);
      os << ";\n";
    }
}

void
MatlabSSWriter::write_der1_assignment(std::ostream &os) const
{
  // initialize out variable
  os << "% setting the output variable\n";
  os << "out = zeros(" << model.getParser().nformulas() << ", " << model.getAtoms().ny() << ");\n";

  // fill out with the terms
  const vector<int> &variables = model.getAtoms().variables();
  const vector<int> &eam = model.getAtoms().get_endo_atoms_map();
  for (int i = 0; i < model.getParser().nformulas(); i++)
    {
      const ogp::FormulaDerivatives &fder = model.getParser().derivatives(i);
      for (int j : eam)
        {
          int tvar = variables[j];
          const string &name = model.getAtoms().name(tvar);
          int yi = model.getAtoms().name2outer_endo(name);
          int t = fder.derivative(ogp::FoldMultiIndex(variables.size(), 1, j));
          if (t != ogp::OperationTree::zero)
            {
              os << "out(" << i+1 << "," << yi+1 << ") = out("<< i+1 << "," << yi+1 << ") + ";
              format_term(t, os);
              os <<  "; % " << name << "(" << model.getAtoms().lead(tvar) << ")\n";
            }
        }
    }
}

void
MatlabSSWriter::format_term(int t, std::ostream &os) const
{
  os << 't' << t;
}

void
MatlabSSWriter::format_nulary(int t, std::ostream &os) const
{
  os << 'a' << t;
}

void
DebugOperationFormatter::format_nulary(int t, std::ostream &os) const
{
  const DynareDynamicAtoms &a = model.getAtoms();

  if (t == ogp::OperationTree::zero)
    os << '0';
  else if (t == ogp::OperationTree::one)
    os << '1';
  else if (t == ogp::OperationTree::nan)
    os << "NaN";
  else if (t == ogp::OperationTree::two_over_pi)
    os << "2/sqrt(PI)";
  else if (a.is_constant(t))
    os << a.get_constant_value(t);
  else
    {
      int ll = a.lead(t);
      const std::string &name = a.name(t);
      if (ll == 0)
        os << name;
      else
        os << name << '(' << ll << ')';
    }
}
