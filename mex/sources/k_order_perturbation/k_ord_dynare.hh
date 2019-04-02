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

#ifndef K_ORD_DYNARE3_H
#define K_ORD_DYNARE3_H

#include <vector>
#include <string>
#include <memory>

#include "journal.hh"
#include "Vector.hh"
#include "twod_matrix.hh"
#include "t_container.hh"
#include "sparse_tensor.hh"
#include "decision_rule.hh"
#include "dynamic_model.hh"
#include "dynamic_abstract_class.hh"

class KordpDynare;

// instantiations of pure abstract class NameList in dynamic_model.h
class DynareNameList : public NameList
{
  std::vector<std::string> names;
public:
  DynareNameList(const KordpDynare &dynare, std::vector<std::string> names_arg);
  int
  getNum() const override
  {
    return static_cast<int>(names.size());
  }
  const std::string &
  getName(int i) const override
  {
    return names[i];
  }
};

class DynareStateNameList : public NameList
{
  std::vector<std::string> names;
public:
  DynareStateNameList(const KordpDynare &dynare, const DynareNameList &dnl,
                      const DynareNameList &denl);
  int
  getNum() const override
  {
    return static_cast<int>(names.size());
  }
  const std::string &
  getName(int i) const override
  {
    return names[i];
  }
};

// The following only implements DynamicModel with help of ogdyn::DynareModel
// instantiation of pure abstract DynamicModel decl. in dynamic_model.h
class KordpDynare : public DynamicModel
{
public:
  const int nStat;
  const int nBoth;
  const int nPred;
  const int nForw;
  const int nExog;
  const int nPar;
  const int nYs;      // = npred + nboth
  const int nYss;     // = nboth + nforw
  const int nY;       // = nstat + npred + nboth + nforw
  const int nJcols;   // nb of jacobian columns = nExog+nY+nYs+nYss
  const Vector &NNZD; /* the total number of non-zero derivative elements
                         where hessian is 2nd : NNZD(order=2) */
  const int nSteps;
  const int nOrder;
private:
  Journal &journal;
  Vector &ySteady;
  Vector &params;
  TwoDMatrix &vCov;
  TensorContainer<FSSparseTensor> md; // ModelDerivatives
  DynareNameList dnl, denl;
  DynareStateNameList dsnl;
  const TwoDMatrix &ll_Incidence;
  std::vector<int> dynppToDyn; // Maps Dynare++ jacobian variable indices to Dynare ones
  std::vector<int> dynToDynpp; // Maps Dynare jacobian variable indices to Dynare++ ones

  std::unique_ptr<TwoDMatrix> g1p, g2p, g3p;

  std::unique_ptr<DynamicModelAC> dynamicModelFile;
public:
  KordpDynare(const std::vector<std::string> &endo,
              const std::vector<std::string> &exo, int num_exo, int num_par,
              Vector &ySteady, TwoDMatrix &vCov, Vector &params, int nstat, int nPred,
              int nforw, int nboth, const Vector &NNZD,
              int nSteps, int ord,
              Journal &jr, std::unique_ptr<DynamicModelAC> dynamicModelFile_arg,
              const std::vector<int> &varOrder, const TwoDMatrix &ll_Incidence,
              std::unique_ptr<TwoDMatrix> g1_arg,
              std::unique_ptr<TwoDMatrix> g2_arg, std::unique_ptr<TwoDMatrix> g3_arg);

  int
  nstat() const override
  {
    return nStat;
  }
  int
  nboth() const override
  {
    return nBoth;
  }
  int
  npred() const override
  {
    return nPred;
  }
  int
  nforw() const override
  {
    return nForw;
  }
  int
  nexog() const override
  {
    return nExog;
  }
  int
  order() const override
  {
    return nOrder;
  }
  const NameList &
  getAllEndoNames() const override
  {
    return dnl;
  }
  const NameList &
  getStateNames() const override
  {
    return dsnl;
  }
  const NameList &
  getExogNames() const override
  {
    return denl;
  }
  const TwoDMatrix &
  getVcov() const override
  {
    return vCov;
  }

  const TensorContainer<FSSparseTensor> &
  getModelDerivatives() const override
  {
    return md;
  }
  const Vector &
  getSteady() const override
  {
    return ySteady;
  }
  Vector &
  getSteady() override
  {
    return ySteady;
  }

  void solveDeterministicSteady() override;
  void evaluateSystem(Vector &out, const ConstVector &yy, const Vector &xx) override;
  void evaluateSystem(Vector &out, const ConstVector &yym, const ConstVector &yy,
                      const ConstVector &yyp, const Vector &xx) override;
  void calcDerivativesAtSteady() override;
  std::unique_ptr<DynamicModel>
  clone() const override
  {
    std::cerr << "KordpDynare::clone() not implemented" << std::endl;
    exit(EXIT_FAILURE);
  }
private:
  // Given the steady state in yS, returns in llxSteady the steady state extended with leads and lags
  void LLxSteady(const Vector &yS, Vector &llxSteady);
  /* Computes the permutations mapping back and forth between Dynare and
     Dynare++ orderings of variables */
  void computeJacobianPermutation(const std::vector<int> &var_order);
  // Fills tensor container (with reordering) at a given order
  void populateDerivativesContainer(const TwoDMatrix &g, int ord);
};

#endif
