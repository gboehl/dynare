/*
 * Copyright © 2004-2011 Ondra Kamenik
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

#include <sstream>
#include <fstream>

#include "dynare3.hh"
#include "dynare_exception.hh"
#include "planner_builder.hh"
#include "forw_subst_builder.hh"

#include "utils/cc/exception.hh"
#include "parser/cc/parser_exception.hh"
#include "parser/cc/atom_substitutions.hh"
#include "../tl/cc/tl_exception.hh"
#include "../kord/kord_exception.hh"

/**************************************************************************************/
/*       DynareNameList class                                                         */
/**************************************************************************************/
std::vector<int>
DynareNameList::selectIndices(const std::vector<std::string> &ns) const
{
  std::vector<int> res;
  for (const auto &n : ns)
    {
      int j = 0;
      while (j < getNum() && n != getName(j))
        j++;
      if (j == getNum())
        throw DynareException(__FILE__, __LINE__,
                              std::string("Couldn't find name for ") + n
                              +" in DynareNameList::selectIndices");
      res.push_back(j);
    }
  return res;
}

/**************************************************************************************/
/*       Dynare class                                                                 */
/**************************************************************************************/

Dynare::Dynare(const std::string &modname, int ord, double sstol, Journal &jr)
  : journal(jr), md(1), ss_tol(sstol)
{
  std::ifstream f{modname};
  if (f.fail())
    throw DynareException(__FILE__, __LINE__, std::string{"Could not open model file "}+modname);

  std::ostringstream buffer;
  buffer << f.rdbuf();
  std::string contents{buffer.str()};

  try
    {
      model = std::make_unique<ogdyn::DynareParser>(contents, ord);
    }
  catch (const ogp::ParserException &pe)
    {
      // Compute line and column, given the offset in the file
      int line = 1;
      int col = 0;
      size_t i = 0;
      while (i < contents.length() && i < static_cast<size_t>(pe.offset()))
        {
          if (contents[i] == '\n')
            {
              line++;
              col = 0;
            }
          i++;
          col++;
        }
      throw DynareException(pe.message(), modname, line, col);
    }
  ysteady = std::make_unique<Vector>(model->getAtoms().ny());
  dnl = std::make_unique<DynareNameList>(*this);
  denl = std::make_unique<DynareExogNameList>(*this);
  dsnl = std::make_unique<DynareStateNameList>(*this, *dnl, *denl);
  fe = std::make_unique<ogp::FormulaEvaluator>(model->getParser());
  fde = std::make_unique<ogp::FormulaDerEvaluator>(model->getParser());
  writeModelInfo(journal);
}

Dynare::Dynare(const std::vector<std::string> &endo,
               const std::vector<std::string> &exo,
               const std::vector<std::string> &par,
               const std::string &equations, int ord,
               double sstol, Journal &jr)
  : journal(jr), md(1), ss_tol(sstol)
{
  try
    {
      model = std::make_unique<ogdyn::DynareSPModel>(endo, exo, par, equations, ord);
    }
  catch (const ogp::ParserException &pe)
    {
      throw DynareException(pe.message(), pe.offset());
    }
  ysteady = std::make_unique<Vector>(model->getAtoms().ny());
  dnl = std::make_unique<DynareNameList>(*this);
  denl = std::make_unique<DynareExogNameList>(*this);
  dsnl = std::make_unique<DynareStateNameList>(*this, *dnl, *denl);
  fe = std::make_unique<ogp::FormulaEvaluator>(model->getParser());
  fde = std::make_unique<ogp::FormulaDerEvaluator>(model->getParser());
  writeModelInfo(journal);
}

Dynare::Dynare(const Dynare &dynare)
  : journal(dynare.journal), md(dynare.md), ss_tol(dynare.ss_tol)
{
  model = dynare.model->clone();
  ysteady = std::make_unique<Vector>(*(dynare.ysteady));
  dnl = std::make_unique<DynareNameList>(*this);
  denl = std::make_unique<DynareExogNameList>(*this);
  dsnl = std::make_unique<DynareStateNameList>(*this, *dnl, *denl);
  fe = std::make_unique<ogp::FormulaEvaluator>(model->getParser());
  fde = std::make_unique<ogp::FormulaDerEvaluator>(model->getParser());
}

void
Dynare::writeMat(mat_t *fd, const std::string &prefix) const
{
  getAllEndoNames().writeMat(fd, prefix + "_vars");
  getAllEndoNames().writeMatIndices(fd, prefix);
  getStateNames().writeMat(fd, prefix + "_state_vars");
  getExogNames().writeMat(fd, prefix + "_shocks");
  getExogNames().writeMatIndices(fd, prefix);
  model->getVcov().writeMat(fd, prefix + "_vcov_exo");
  TwoDMatrix aux(1, 1);
  aux.get(0, 0) = nstat();
  aux.writeMat(fd, prefix + "_nstat");
  aux.get(0, 0) = npred();
  aux.writeMat(fd, prefix + "_npred");
  aux.get(0, 0) = nboth();
  aux.writeMat(fd, prefix + "_nboth");
  aux.get(0, 0) = nforw();
  aux.writeMat(fd, prefix + "_nforw");
}

void
Dynare::writeDump(const std::string &basename) const
{
  std::string fname(basename + ".dump");
  std::ofstream out(fname);
  model->dump_model(out);
  out.close();
}

void
Dynare::solveDeterministicSteady(Vector &steady)
{
  JournalRecordPair pa(journal);
  pa << "Non-linear solver for deterministic steady state" << endrec;
  steady = const_cast<const Vector &>(model->getInit());
  DynareVectorFunction dvf(*this);
  DynareJacobian dj(*this);
  ogu::NLSolver nls(dvf, dj, 500, ss_tol, journal);
  int iter;
  if (!nls.solve(steady, iter))
    throw DynareException(__FILE__, __LINE__,
                          "Could not obtain convergence in non-linear solver");
}

// Evaluate system at given yₜ=yₜ₊₁=yₜ₋₁, and given shocks xₜ
void
Dynare::evaluateSystem(Vector &out, const ConstVector &yy, const Vector &xx)
{
  ConstVector yym(yy, nstat(), nys());
  ConstVector yyp(yy, nstat()+npred(), nyss());
  evaluateSystem(out, yym, yy, yyp, xx);
}

/* Evaluate system at given y*ₜ₋₁, yₜ, y**ₜ₊₁ and at exogenous xₜ, all three
   vectors yym, yy, and yyp have the respective lengths of y*ₜ₋₁, yₜ, y**ₜ₊₁ */
void
Dynare::evaluateSystem(Vector &out, const ConstVector &yym, const ConstVector &yy,
                       const ConstVector &yyp, const Vector &xx)
{
  ogdyn::DynareAtomValues dav(model->getAtoms(), model->getParams(), yym, yy, yyp, xx);
  DynareEvalLoader del(model->getAtoms(), out);
  fe->eval(dav, del);
}

void
Dynare::calcDerivatives(const Vector &yy, const Vector &xx)
{
  ConstVector yym(yy, nstat(), nys());
  ConstVector yyp(yy, nstat()+npred(), nyss());
  ogdyn::DynareAtomValues dav(model->getAtoms(), model->getParams(), yym, yy, yyp, xx);
  DynareDerEvalLoader ddel(model->getAtoms(), md, model->getOrder());
  for (int iord = 1; iord <= model->getOrder(); iord++)
    fde->eval(dav, ddel, iord);
}

void
Dynare::calcDerivativesAtSteady()
{
  Vector xx(nexog());
  xx.zeros();
  calcDerivatives(*ysteady, xx);
}

void
Dynare::writeModelInfo(Journal &jr) const
{
  // write info on variables
  {
    JournalRecordPair rp(journal);
    rp << "Information on variables" << endrec;
    JournalRecord rec1(journal);
    rec1 << "Number of endogenous:            " << ny() << endrec;
    JournalRecord rec2(journal);
    rec2 << "Number of exogenous:             " << nexog() << endrec;
    JournalRecord rec3(journal);
    rec3 << "Number of static:                " << nstat() << endrec;
    JournalRecord rec4(journal);
    rec4 << "Number of predetermined:         " << npred()+nboth() << endrec;
    JournalRecord rec5(journal);
    rec5 << "Number of forward looking:       " << nforw()+nboth() << endrec;
    JournalRecord rec6(journal);
    rec6 << "Number of both:                  " << nboth() << endrec;
  }

  // write info on planner variables
  const ogdyn::PlannerInfo *pinfo = model->get_planner_info();
  if (pinfo)
    {
      JournalRecordPair rp(journal);
      rp << "Information on planner variables" << endrec;
      JournalRecord rec1(journal);
      rec1 << "Number of Lagrange multipliers:  " << pinfo->num_lagrange_mults << endrec;
      JournalRecord rec2(journal);
      rec2 << "Number of auxiliary variables:   " << pinfo->num_aux_variables << endrec;
      JournalRecord rec3(journal);
      rec3 << "Number of new terms in the tree: " << pinfo->num_new_terms << endrec;
    }

  // write info on forward substitutions
  const ogdyn::ForwSubstInfo *finfo = model->get_forw_subst_info();
  if (finfo)
    {
      JournalRecordPair rp(journal);
      rp << "Information on forward substitutions" << endrec;
      JournalRecord rec1(journal);
      rec1 << "Number of affected equations:    " << finfo->num_affected_equations << endrec;
      JournalRecord rec2(journal);
      rec2 << "Number of substituted terms:     " << finfo->num_subst_terms << endrec;
      JournalRecord rec3(journal);
      rec3 << "Number of auxiliary variables:   " << finfo->num_aux_variables << endrec;
      JournalRecord rec4(journal);
      rec4 << "Number of new terms in the tree: " << finfo->num_new_terms << endrec;
    }

  // write info on substitutions
  const ogp::SubstInfo *sinfo = model->get_subst_info();
  if (sinfo)
    {
      JournalRecordPair rp(journal);
      rp << "Information on substitutions" << endrec;
      JournalRecord rec1(journal);
      rec1 << "Number of substitutions:         " << sinfo->num_substs << endrec;
    }
}

DynareNameList::DynareNameList(const Dynare &dynare)
{
  for (int i = 0; i < dynare.ny(); i++)
    {
      int j = dynare.model->getAtoms().y2outer_endo()[i];
      const std::string &name = dynare.model->getAtoms().get_endovars()[j];
      names.push_back(name);
    }
}

DynareStateNameList::DynareStateNameList(const Dynare &dynare, const DynareNameList &dnl,
                                         const DynareExogNameList &denl)
{
  for (int i = 0; i < dynare.nys(); i++)
    names.push_back(dnl.getName(i+dynare.nstat()));
  for (int i = 0; i < dynare.nexog(); i++)
    names.push_back(denl.getName(i));
}

DynareExogNameList::DynareExogNameList(const Dynare &dynare)
{
  for (int i = 0; i < dynare.nexog(); i++)
    {
      int j = dynare.model->getAtoms().y2outer_exo()[i];
      const std::string &name = dynare.model->getAtoms().get_exovars()[j];
      names.push_back(name);
    }
}

DynareEvalLoader::DynareEvalLoader(const ogp::FineAtoms &a, Vector &out)
  : Vector(out)
{
  if (a.ny() != out.length())
    throw DynareException(__FILE__, __LINE__, "Wrong length of out vector in DynareEvalLoader constructor");
}

/* This clears the container of model derivatives and initializes it inserting
   empty sparse tensors up to the given order. */
DynareDerEvalLoader::DynareDerEvalLoader(const ogp::FineAtoms &a,
                                         TensorContainer<FSSparseTensor> &mod_ders,
                                         int order)
  : atoms(a), md(mod_ders)
{
  md.clear();
  for (int iord = 1; iord <= order; iord++)
    {
      auto t = std::make_unique<FSSparseTensor>(iord, atoms.ny()+atoms.nys()+atoms.nyss()+atoms.nexo(), atoms.ny());
      md.insert(std::move(t));
    }
}

void
DynareDerEvalLoader::load(int i, int iord, const int *vars, double res)
{
  FSSparseTensor &t = md.get(Symmetry{iord});
  IntSequence s(iord, 0);
  for (int j = 0; j < iord; j++)
    s[j] = atoms.get_pos_of_all(vars[j]);
  t.insert(s, i, res);
}

DynareJacobian::DynareJacobian(Dynare &dyn)
  : Jacobian(dyn.ny()), d(dyn)
{
  zeros();
}

void
DynareJacobian::eval(const Vector &yy)
{
  ogdyn::DynareSteadyAtomValues
    dav(d.getModel().getAtoms(), d.getModel().getParams(), yy);
  zeros();
  d.fde->eval(dav, *this, 1);
}

void
DynareJacobian::load(int i, int iord, const int *vars, double res)
{
  if (iord != 1)
    throw DynareException(__FILE__, __LINE__,
                          "Derivative order different from order=1 in DynareJacobian::load");

  int t = vars[0];
  int j = d.getModel().getAtoms().get_pos_of_all(t);
  if (j < d.nyss())
    get(i, j+d.nstat()+d.npred()) += res;
  else if (j < d.nyss()+d.ny())
    get(i, j-d.nyss()) += res;
  else if (j < d.nyss()+d.ny()+d.nys())
    get(i, j-d.nyss()-d.ny()+d.nstat()) += res;
}

void
DynareVectorFunction::eval(const ConstVector &in, Vector &out)
{
  check_for_eval(in, out);
  Vector xx(d.nexog());
  xx.zeros();
  d.evaluateSystem(out, in, xx);
}
