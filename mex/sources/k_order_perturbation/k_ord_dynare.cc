/*
 * Copyright (C) 2008-2019 Dynare Team
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

// GP, based on work by O.Kamenik

#include "k_ord_dynare.hh"
#include "dynamic_abstract_class.hh"
#include "dynare_exception.hh"

#include <utility>


KordpDynare::KordpDynare(const std::vector<std::string> &endo,
                         const std::vector<std::string> &exo, int nexog, int npar,
                         Vector &ysteady, TwoDMatrix &vcov, Vector &inParams, int nstat,
                         int npred, int nforw, int nboth, const Vector &nnzd,
                         int nsteps, int norder,
                         Journal &jr, std::unique_ptr<DynamicModelAC> dynamicModelFile_arg,
                         const std::vector<int> &var_order, const TwoDMatrix &llincidence,
                         std::unique_ptr<TwoDMatrix> g1_arg, std::unique_ptr<TwoDMatrix> g2_arg,
                         std::unique_ptr<TwoDMatrix> g3_arg) :
  nStat{nstat}, nBoth{nboth}, nPred{npred}, nForw{nforw}, nExog{nexog}, nPar{npar},
  nYs{npred + nboth}, nYss{nboth + nforw}, nY{nstat + npred + nboth + nforw},
  nJcols{nExog+nY+nYs+nYss}, NNZD{nnzd}, nSteps{nsteps},
  nOrder{norder}, journal{jr}, ySteady{ysteady}, params{inParams}, vCov{vcov},
  md{1}, dnl{*this, endo}, denl{*this, exo}, dsnl{*this, dnl, denl},
  ll_Incidence{llincidence},
  g1p{std::move(g1_arg)}, g2p{std::move(g2_arg)}, g3p{std::move(g3_arg)},
  dynamicModelFile{std::move(dynamicModelFile_arg)}
{
  computeJacobianPermutation(var_order);
}

void
KordpDynare::solveDeterministicSteady()
{
  JournalRecordPair pa(journal);
  pa << "Non-linear solver for deterministic steady state skipped" << endrec;
}

void
KordpDynare::evaluateSystem(Vector &out, const ConstVector &yy, const Vector &xx)
{
  // This method is only called when checking the residuals at steady state (Approximation::check), so return zero residuals
  out.zeros();
}

void
KordpDynare::evaluateSystem(Vector &out, const ConstVector &yym, const ConstVector &yy,
                            const ConstVector &yyp, const Vector &xx)
{
  // This method is only called when checking the residuals at steady state (Approximation::check), so return zero residuals
  out.zeros();
}

void
KordpDynare::calcDerivativesAtSteady()
{
  if (!g1p)
    {
      g1p = std::make_unique<TwoDMatrix>(nY, nJcols);
      g1p->zeros();

      if (nOrder > 1)
        {
          // allocate space for sparse Hessian
          g2p = std::make_unique<TwoDMatrix>(static_cast<int>(NNZD[1]), 3);
          g2p->zeros();
        }

      if (nOrder > 2)
        {
          g3p = std::make_unique<TwoDMatrix>(static_cast<int>(NNZD[2]), 3);
          g3p->zeros();
        }

      Vector xx(nexog());
      xx.zeros();

      Vector out(nY);
      out.zeros();
      Vector llxSteady(nJcols-nExog);
      LLxSteady(ySteady, llxSteady);

      dynamicModelFile->eval(llxSteady, xx, params, ySteady, out, g1p.get(), g2p.get(), g3p.get());
    }

  populateDerivativesContainer(*g1p, 1);

  if (nOrder > 1)
    populateDerivativesContainer(*g2p, 2);

  if (nOrder > 2)
    populateDerivativesContainer(*g3p, 3);
}

void
KordpDynare::populateDerivativesContainer(const TwoDMatrix &g, int ord)
{
  // model derivatives FSSparseTensor instance
  auto mdTi = std::make_unique<FSSparseTensor>(ord, nJcols, nY);

  IntSequence s(ord, 0);

  if (ord == 1)
    for (int i = 0; i < g.ncols(); i++)
      {
        for (int j = 0; j < g.nrows(); j++)
          {
            double x = g.get(j, dynppToDyn[s[0]]);
            if (x != 0.0)
              mdTi->insert(s, j, x);
          }
        s[0]++;
      }
  else if (ord == 2)
    {
      for (int i = 0; i < g.nrows(); i++)
        {
          int j = static_cast<int>(g.get(i, 0))-1; // hessian indices start with 1
          int i1 = static_cast<int>(g.get(i, 1))-1;
          if (j < 0 || i1 < 0)
            continue; // Discard empty entries (see comment in DynamicModelAC::unpackSparseMatrix())
          int s0 = i1 / nJcols;
          int s1 = i1 % nJcols;
          s[0] = dynToDynpp[s0];
          s[1] = dynToDynpp[s1];
          if (s[1] >= s[0])
            {
              double x = g.get(i, 2);
              mdTi->insert(s, j, x);
            }
        }
    }
  else if (ord == 3)
    {
      int nJcols2 = nJcols*nJcols;
      for (int i = 0; i < g.nrows(); i++)
        {
          int j = static_cast<int>(g.get(i, 0))-1;
          int i1 = static_cast<int>(g.get(i, 1))-1;
          if (j < 0 || i1 < 0)
            continue; // Discard empty entries (see comment in DynamicModelAC::unpackSparseMatrix())
          int s0 = i1 / nJcols2;
          int i2 = i1 % nJcols2;
          int s1 = i2 / nJcols;
          int s2 = i2 % nJcols;
          s[0] = dynToDynpp[s0];
          s[1] = dynToDynpp[s1];
          s[2] = dynToDynpp[s2];
          if (s.isSorted())
            {
              double x = g.get(i, 2);
              mdTi->insert(s, j, x);
            }
        }
    }

  md.insert(std::move(mdTi));
}

/* Returns ySteady extended with leads and lags suitable for passing to
   <model>_dynamic */
void
KordpDynare::LLxSteady(const Vector &yS, Vector &llxSteady)
{
  if (yS.length() == nJcols-nExog)
    throw DynareException(__FILE__, __LINE__, "ySteady already of right size");

  /* Create temporary square 2D matrix size nEndo×nEndo (sparse)
     for the lag, current and lead blocks of the jacobian */
  if (llxSteady.length() != nJcols-nExog)
    throw DynareException(__FILE__, __LINE__, "llxSteady has wrong size");

  for (int ll_row = 0; ll_row < ll_Incidence.nrows(); ll_row++)
    // populate (non-sparse) vector with ysteady values
    for (int i = 0; i < nY; i++)
      if (ll_Incidence.get(ll_row, i))
        llxSteady[static_cast<int>(ll_Incidence.get(ll_row, i))-1] = yS[i];
}

/************************************
 * Reorder DynareJacobianIndices of variables in a vector according to
 * given int * varOrder together with lead & lag incidence matrix and
 * any the extra columns for exogenous vars, and then,
 * reorders its blocks given by the varOrder and the Dynare++ expectations:

 * extra	nboth+ npred (t-1) lags
 * varOrder
                static:
    pred
    both
    forward
    * extra both + nforw (t+1) leads, and
    * extra exogen

    * so to match the jacobian organisation expected by the Appoximation class
      both + nforw (t+1) leads
      static
      pred
      both
      forward
      nboth+ npred  (t-1) lags
      exogen
************************************/

void
KordpDynare::computeJacobianPermutation(const std::vector<int> &var_order)
{
  dynppToDyn.resize(nJcols);
  // create temporary square 2D matrix size nEndo×nEndo (sparse)
  // for the lag, current and lead blocks of the jacobian
  std::vector<int> tmp(nY);
  int rjoff = nJcols-nExog-1;

  for (int ll_row = 0; ll_row < ll_Incidence.nrows(); ll_row++)
    {
      // reorder in orde-var order & populate temporary nEndo (sparse) vector with
      // the lag, current and lead blocks of the jacobian respectively
      for (int i = 0; i < nY; i++)
        tmp[i] = static_cast<int>(ll_Incidence.get(ll_row, var_order[i]-1));
      // write the reordered blocks back to the jacobian
      // in reverse order
      for (int j = nY-1; j >= 0; j--)
        if (tmp[j])
          {
            dynppToDyn[rjoff] = tmp[j]-1;
            rjoff--;
            if (rjoff < 0)
              break;
          }
    }

  //add the indices for the nExog exogenous jacobians
  for (int j = nJcols-nExog; j < nJcols; j++)
    dynppToDyn[j] = j;

  // Compute reverse ordering
  dynToDynpp.resize(nJcols);
  for (int i = 0; i < nJcols; i++)
    dynToDynpp[dynppToDyn[i]] = i;
}


DynareNameList::DynareNameList(const KordpDynare &dynare, std::vector<std::string> names_arg)
  : names(std::move(names_arg))
{
}

DynareStateNameList::DynareStateNameList(const KordpDynare &dynare, const DynareNameList &dnl,
                                         const DynareNameList &denl)
{
  for (int i = 0; i < dynare.nYs; i++)
    names.emplace_back(dnl.getName(i+dynare.nstat()));
  for (int i = 0; i < dynare.nexog(); i++)
    names.emplace_back(denl.getName(i));
}
