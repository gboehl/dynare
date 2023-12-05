/*
 * Copyright © 2005 Ondra Kamenik
 * Copyright © 2019-2023 Dynare Team
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

// Approximating model solution

/* The class Approximation in this file is a main interface to the
   algorithms calculating approximations to the decision rule about
   deterministic and stochastic steady states.

   The approximation about a deterministic steady state is solved by classes
   FirstOrder and KOrder. The approximation around the stochastic steady state
   is computed by class KOrderStoch together with a method of Approximation
   class.

   The approximation around the stochastic steady state is done with an
   explicit expression of forward derivatives of g**. More formally,
   we have to solve the decision rule g from the implicit system:

    𝔼ₜ(f(g**(g*(y*,uₜ,σ),uₜ₊₁,σ),g(y*,uₜ,σ),yₜ,uₜ)) = 0

   The term within the expectations can be Taylor expanded, and the expectation
   can be driven into the formula. However, when doing this at σ≠0, the term
   g** at σ≠0 is dependent on uₜ₊₁ and thus the integral of its approximation
   includes all derivatives w.r.t. u of g**. Note that for σ=0, the derivatives
   of g** in this context are constant. This is the main difference between the
   approximation at deterministic steady (σ=0), and stochastic steady (σ≠0).
   This means that the k-order derivative of the above equation at σ≠0 depends
   on all derivatives of g** (including those with order greater than k).

   The explicit expression of the forward g** means that the derivatives of g
   are not solved simultaneously, but that the forward derivatives of g** are
   calculated as an extrapolation based on the approximation at lower σ. This
   is exactly what does Approximation::walkStochSteady(). It starts at the
   deterministic steady state, and in a few steps it adds to σ explicitly
   expressing forward g** from a previous step.

   Further details on the both solution methods are given in (todo: put
   references here when they exist).

   Very important note: all classes here used for calculation of decision
   rule approximation are folded. For the time being, it seems that Faà
   Di Bruno formula is quicker for folded tensors, and that is why we
   stick to folded tensors here. However, when the calcs are done, we
   calculate also its unfolded versions, to be available for simulations
   and so on. */

#ifndef APPROXIMATION_H
#define APPROXIMATION_H

#include "decision_rule.hh"
#include "dynamic_model.hh"
#include "journal.hh"
#include "korder.hh"

#include <memory>

/* This class is used to calculate derivatives by Faà Di Bruno of
   f(g**(g*(y*,u,σ),u′,σ),g(y*,u,σ),y*,u) with respect to u′. In order to keep
   it as simple as possible, the class represents an equivalent (with respect
   to u′) container for f(g**(y*,u′,σ),0,0,0). The class is used only for
   evaluation of approximation error in Approximation class, which is
   calculated in Approximation::calcStochShift().

   Since it is a folded version, we inherit from StackContainer<FGSTensor> and
   FoldedStackContainer. To construct it, we need only the g** container and
   size of stacks. */

class ZAuxContainer : public StackContainer<FGSTensor>, public FoldedStackContainer
{
public:
  using _Ctype = StackContainer<FGSTensor>::_Ctype;
  using itype = StackContainer<FGSTensor>::itype;
  ZAuxContainer(const _Ctype* gss, int ngss, int ng, int ny, int nu);
  [[nodiscard]] itype getType(int i, const Symmetry& s) const override;
};

/* This class provides an interface to approximation algorithms. The core
   method is walkStochSteady() which calculates the approximation about
   stochastic steady state in a given number of steps. The number is given as a
   parameter ‘ns’ of the constructor. If the number is equal to zero, the
   resulted approximation is about the deterministic steady state.

   An object is constructed from the DynamicModel, and the number of steps
   ‘ns’. Also, we pass a reference to journal. That's all. The result of the
   core method walkStochSteady() is a decision rule ‘dr’ and a matrix ‘ss’
   whose columns are steady states for increasing σ during the walk. Both can
   be retrived by public methods. The first column of the matrix is the
   deterministic steady state, the last is the stochastic steady state for the
   full size shocks.

   The method walkStochSteady() calls the following methods: approxAtSteady()
   calculates an initial approximation about the deterministic steady,
   saveRuleDerivs() saves derivatives of a rule for the following step in
   ‘rule_ders’ and ‘rule_ders_ss’ (see Approximation::saveRuleDerivs() for
   their description), check() reports an error of the current approximation
   and calcStochShift() (called from check()) calculates a shift of the system
   equations due to uncertainity.

   ‘dr_centralize’ is a new option. Dynare++ was automatically expressing
   results around the fixed point instead of the deterministic steady state.
   ‘dr_centralize’ controls this behavior. */

class Approximation
{
  DynamicModel& model;
  Journal& journal;
  std::unique_ptr<FGSContainer> rule_ders;
  std::unique_ptr<FGSContainer> rule_ders_s;
  std::unique_ptr<FGSContainer> rule_ders_ss;
  std::unique_ptr<FoldDecisionRule> fdr;
  std::unique_ptr<FoldDecisionRule> fdr_pruning;
  std::unique_ptr<UnfoldDecisionRule> udr;
  std::unique_ptr<UnfoldDecisionRule> udr_pruning;
  const PartitionY ypart;
  const FNormalMoments mom;
  IntSequence nvs;
  int steps;
  bool dr_centralize;
  bool pruning;
  double qz_criterium;
  TwoDMatrix ss;

public:
  Approximation(DynamicModel& m, Journal& j, int ns, bool dr_centr, bool pruning, double qz_crit);

  [[nodiscard]] const FoldDecisionRule& getFoldDecisionRule() const;
  [[nodiscard]] const UnfoldDecisionRule& getUnfoldDecisionRulePruning() const;
  [[nodiscard]] const UnfoldDecisionRule& getUnfoldDecisionRule() const;
  [[nodiscard]] const TwoDMatrix&
  getSS() const
  {
    return ss;
  }
  [[nodiscard]] const DynamicModel&
  getModel() const
  {
    return model;
  }

  void walkStochSteady();
  [[nodiscard]] TwoDMatrix calcYCov() const;
  [[nodiscard]] const FGSContainer&
  get_rule_ders() const
  {
    return *rule_ders;
  }
  [[nodiscard]] const FGSContainer&
  get_rule_ders_s() const
  {
    return *rule_ders_s;
  }
  [[nodiscard]] const FGSContainer&
  get_rule_ders_ss() const
  {
    return *rule_ders_ss;
  }

protected:
  void approxAtSteady();
  void calcStochShift(Vector& out, double at_sigma) const;
  void saveRuleDerivs(const FGSContainer& g);
  void check(double at_sigma) const;
};

#endif
