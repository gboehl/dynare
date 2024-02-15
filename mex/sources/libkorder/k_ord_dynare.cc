/*
 * Copyright © 2008-2024 Dynare Team
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

// GP, based on work by O.Kamenik

#include "k_ord_dynare.hh"
#include "dynamic_abstract_class.hh"
#include "dynare_exception.hh"

#include <cassert>
#include <utility>

KordpDynare::KordpDynare(const std::vector<std::string>& endo, const std::vector<std::string>& exo,
                         int nexog, int npar, Vector& ysteady, TwoDMatrix& vcov, Vector& inParams,
                         int nstat, int npred, int nforw, int nboth, const ConstVector& nnzd,
                         int nsteps, int norder, Journal& jr,
                         std::unique_ptr<DynamicModelAC> dynamicModelFile_arg,
                         const std::vector<int>& dr_order) :
    nStat {nstat},
    nBoth {nboth},
    nPred {npred},
    nForw {nforw},
    nExog {nexog},
    nPar {npar},
    nYs {npred + nboth},
    nYss {nboth + nforw},
    nY {nstat + npred + nboth + nforw},
    nJcols {nExog + nY + nYs + nYss},
    NNZD {nnzd},
    nSteps {nsteps},
    nOrder {norder},
    journal {jr},
    ySteady {ysteady},
    params {inParams},
    vCov {vcov},
    md {1},
    dnl {endo},
    denl {exo},
    dsnl {*this, dnl, denl},
    dynamicModelFile {std::move(dynamicModelFile_arg)}
{
  computeJacobianPermutation(dr_order);
}

void
KordpDynare::solveDeterministicSteady()
{
  JournalRecordPair pa(journal);
  pa << "Non-linear solver for deterministic steady state skipped" << endrec;
}

void
KordpDynare::evaluateSystem(Vector& out, [[maybe_unused]] const ConstVector& yy,
                            [[maybe_unused]] const Vector& xx)
{
  // This method is only called when checking the residuals at steady state (Approximation::check),
  // so return zero residuals
  out.zeros();
}

void
KordpDynare::evaluateSystem(Vector& out, [[maybe_unused]] const ConstVector& yym,
                            [[maybe_unused]] const ConstVector& yy,
                            [[maybe_unused]] const ConstVector& yyp,
                            [[maybe_unused]] const Vector& xx)
{
  // This method is only called when checking the residuals at steady state (Approximation::check),
  // so return zero residuals
  out.zeros();
}

void
KordpDynare::calcDerivativesAtSteady()
{
  assert(md.begin() == md.end());

  Vector xx(nexog());
  xx.zeros();

  Vector out(nY);
  out.zeros();
  Vector llxSteady(3 * nY);
  std::copy_n(ySteady.base(), nY, llxSteady.base());
  std::copy_n(ySteady.base(), nY, llxSteady.base() + nY);
  std::copy_n(ySteady.base(), nY, llxSteady.base() + 2 * nY);

  dynamicModelFile->eval(llxSteady, xx, params, ySteady, out, dynToDynpp, md);
}

/*
   Computes mapping between Dynare and Dynare++ orderings of the (dynamic)
   variable indices in derivatives.

   If one defines:
   – y (resp. x) as the vector of all endogenous (size nY), in DR-order (resp.
     declaration order)
   – y⁻ as the vector of endogenous that appear at previous period (size nYs),
     in DR-order
   – y⁺ as the vector of endogenous that appear at future period (size nYss) in
     DR-order
   – u as the vector of exogenous (size nExog)

   In Dynare++, the vector is of size nY+nYs+nYss+nExog and the ordering is (y⁺, y, y⁻, u).
   In Dynare, the vector is of size 3*nY+nExog and the ordering is (x, x, x, u).
*/
void
KordpDynare::computeJacobianPermutation(const std::vector<int>& dr_order)
{
  // Compute Dynare++ → Dynare ordering
  std::vector<int> dynppToDyn(nJcols);
  int j {0};
  for (; j < nYss; j++)
    dynppToDyn[j] = dr_order[j + nStat + nPred] + 2 * nY; // Forward variables
  for (; j < nYss + nY; j++)
    dynppToDyn[j] = dr_order[j - nYss] + nY; // Variables in current period
  for (; j < nYss + nY + nYs; j++)
    dynppToDyn[j] = dr_order[j - (nY + nYss) + nStat]; // Predetermined variables
  for (; j < nJcols; j++)
    dynppToDyn[j] = j - (nY + nYss + nYs) + 3 * nY; // Exogenous

  // Compute Dynare → Dynare++ ordering
  for (int i {0}; i < nJcols; i++)
    dynToDynpp.emplace(dynppToDyn[i], i);
}

DynareNameList::DynareNameList(std::vector<std::string> names_arg) : names(std::move(names_arg))
{
}

DynareStateNameList::DynareStateNameList(const KordpDynare& dynare, const DynareNameList& dnl,
                                         const DynareNameList& denl)
{
  for (int i = 0; i < dynare.nYs; i++)
    names.emplace_back(dnl.getName(i + dynare.nStat));
  for (int i = 0; i < dynare.nExog; i++)
    names.emplace_back(denl.getName(i));
}
