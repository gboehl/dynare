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

#include <utility>

#include "approximation.hh"
#include "first_order.hh"
#include "kord_exception.hh"
#include "korder_stoch.hh"

ZAuxContainer::ZAuxContainer(const _Ctype* gss, int ngss, int ng, int ny, int nu) :
    StackContainer<FGSTensor>(4, 1)
{
  stack_sizes = {ngss, ng, ny, nu};
  conts[0] = gss;
  calculateOffsets();
}

/* The getType() method corresponds to f(g**(y*,u′,σ),0,0,0). For the first
   argument we return ‘matrix’, for other three we return ‘zero’. */

ZAuxContainer::itype
ZAuxContainer::getType(int i, const Symmetry& s) const
{
  if (i == 0)
    {
      if (s[2] > 0)
        return itype::zero;
      else
        return itype::matrix;
    }
  return itype::zero;
}

Approximation::Approximation(DynamicModel& m, Journal& j, int ns, bool dr_centr, bool pruned_dr,
                             double qz_crit) :
    model(m),
    journal(j),
    ypart(model.nstat(), model.npred(), model.nboth(), model.nforw()),
    mom(UNormalMoments(model.order(), model.getVcov())),
    nvs {ypart.nys(), model.nexog(), model.nexog(), 1},
    steps(ns),
    dr_centralize(dr_centr),
    pruning(pruned_dr),
    qz_criterium(qz_crit),
    ss(ypart.ny(), steps + 1)
{
  ss.nans();
}

/* This just returns ‘fdr’ with a check that it is created. */
const FoldDecisionRule&
Approximation::getFoldDecisionRule() const
{
  KORD_RAISE_IF(!fdr,
                "Folded decision rule has not been created in Approximation::getFoldDecisionRule");
  return *fdr;
}

/* This just returns ‘fdr_pruning’ with a check that it is created. */
const UnfoldDecisionRule&
Approximation::getUnfoldDecisionRulePruning() const
{
  KORD_RAISE_IF(
      !udr_pruning,
      "Unfolded decision rule has not been created in Approximation::getUnfoldDecisionRule");
  return *udr_pruning;
}

/* This just returns ‘udr’ with a check that it is created. */
const UnfoldDecisionRule&
Approximation::getUnfoldDecisionRule() const
{
  KORD_RAISE_IF(
      !udr, "Unfolded decision rule has not been created in Approximation::getUnfoldDecisionRule");
  return *udr;
}

/* This methods assumes that the deterministic steady state is
   model.getSteady(). It makes an approximation about it and stores the
   derivatives to ‘rule_ders’ and ‘rule_ders_ss’. Also it runs a check() for
   σ=0. */
void
Approximation::approxAtSteady()
{
  model.calcDerivativesAtSteady();
  FirstOrder fo(model.nstat(), model.npred(), model.nboth(), model.nforw(), model.nexog(),
                model.getModelDerivatives().get(Symmetry {1}), journal, qz_criterium);

  if (model.order() >= 2)
    {
      KOrder korder(model.nstat(), model.npred(), model.nboth(), model.nforw(),
                    model.getModelDerivatives(), fo.getGy(), fo.getGu(), model.getVcov(), journal);
      korder.switchToFolded();
      for (int k = 2; k <= model.order(); k++)
        korder.performStep<Storage::fold>(k);

      saveRuleDerivs(korder.getFoldDers());
    }
  else
    {
      FirstOrderDerivs<Storage::fold> fo_ders(fo);
      saveRuleDerivs(fo_ders);
    }
  check(0.0);
}

/* This is the core routine of Approximation class.

   First we solve for the approximation about the deterministic steady state.
   Then we perform ‘steps’ cycles toward the stochastic steady state. Each
   cycle moves the size of shocks by dsigma=1.0/steps. At the end of a cycle,
   we have ‘rule_ders’ being the derivatives at stochastic steady state for
   σ=sigma_so_far+dsigma and model.getSteady() being the steady state.

   If the number of ‘steps’ is zero, the decision rule ‘dr’ at the bottom is
   created from derivatives about deterministic steady state, with size of σ=1.
   Otherwise, the ‘dr’ is created from the approximation about stochastic
   steady state with σ=0.

   Within each cycle, we first make a backup of the last steady (from
   initialization or from a previous cycle), then we calculate the fix point of
   the last rule with σ=dsigma. This becomes a new steady state at the
   σ=sigma_so_far+dsigma. We calculate expectations of g**(y,σ ηₜ₊₁,σ)
   expressed as a Taylor expansion around the new σ and the new steady state.
   Then we solve for the decision rule with explicit g** at t+1 and save the
   rule.

   After we reached σ=1, the decision rule is formed.

   The biproduct of this method is the matrix ‘ss’, whose columns are steady
   states for subsequent σ’s. The first column is the deterministic steady
   state, the last column is the stochastic steady state for a full size of
   shocks (σ=1). There are ‘steps+1’ columns. */

void
Approximation::walkStochSteady()
{
  // initial approximation at deterministic steady
  /* Here we solve for the deterministic steady state, calculate
     approximation at the deterministic steady and save the steady state
     to ‘ss’. */
  model.solveDeterministicSteady();
  approxAtSteady();
  Vector steady0 {ss.getCol(0)};
  steady0 = model.getSteady();

  double sigma_so_far = 0.0;
  double dsigma = (steps == 0) ? 0.0 : 1.0 / steps;
  for (int i = 1; i <= steps; i++)
    {
      JournalRecordPair pa(journal);
      pa << "Approximation about stochastic steady for sigma=" << sigma_so_far + dsigma << endrec;

      Vector last_steady(const_cast<const Vector&>(model.getSteady()));

      // calculate fix-point of the last rule for ‘dsigma’
      /* We form the DRFixPoint object from the last rule with σ=dsigma. Then
         we save the steady state to ‘ss’. The new steady is also put to
         model.getSteady(). */
      DRFixPoint<Storage::fold> fp(*rule_ders, ypart, model.getSteady(), dsigma);
      bool converged = fp.calcFixPoint(model.getSteady());
      JournalRecord rec(journal);
      rec << "Fix point calcs: iter=" << fp.getNumIter()
          << ", newton_iter=" << fp.getNewtonTotalIter()
          << ", last_newton_iter=" << fp.getNewtonLastIter() << ".";
      if (converged)
        rec << " Converged." << endrec;
      else
        {
          rec << " Not converged!!" << endrec;
          KORD_RAISE_X("Fix point calculation not converged", KORD_FP_NOT_CONV);
        }
      Vector steadyi {ss.getCol(i)};
      steadyi = model.getSteady();

      // calculate ‘hh’ as expectations of the last g**
      /* We form the steady state shift ‘dy’, which is the new steady state
         minus the old steady state. Then we create StochForwardDerivs object,
         which calculates the derivatives of g** expectations at new sigma and
         new steady. */
      Vector dy(const_cast<const Vector&>(model.getSteady()));
      dy.add(-1.0, last_steady);

      StochForwardDerivs<Storage::fold> hh(ypart, model.nexog(), *rule_ders_ss, mom, dy, dsigma,
                                           sigma_so_far);
      JournalRecord rec1(journal);
      rec1 << "Calculation of g** expectations done" << endrec;

      // form KOrderStoch, solve and save
      /* We calculate derivatives of the model at the new steady, form
         KOrderStoch object and solve, and save the rule. */
      model.calcDerivativesAtSteady();
      KOrderStoch korder_stoch(ypart, model.nexog(), model.getModelDerivatives(), hh, journal);
      for (int d = 1; d <= model.order(); d++)
        korder_stoch.performStep<Storage::fold>(d);

      saveRuleDerivs(korder_stoch.getFoldDers());

      check(sigma_so_far + dsigma);
      sigma_so_far += dsigma;
    }

  // construct the resulting decision rules
  udr.reset();
  fdr = std::make_unique<FoldDecisionRule>(*rule_ders, ypart, model.nexog(), model.getSteady(),
                                           1.0 - sigma_so_far);
  if (pruning)
    {
      fdr_pruning = std::make_unique<FoldDecisionRule>(
          *rule_ders, ypart, model.nexog(), model.getSteady(), 1.0 - sigma_so_far, pruning);
      udr_pruning = std::make_unique<UnfoldDecisionRule>(*fdr_pruning);
    }
  if (steps == 0 && dr_centralize)
    {
      // centralize decision rule for zero steps
      DRFixPoint<Storage::fold> fp(*rule_ders, ypart, model.getSteady(), 1.0);
      bool converged = fp.calcFixPoint(model.getSteady());
      JournalRecord rec(journal);
      rec << "Fix point calcs: iter=" << fp.getNumIter()
          << ", newton_iter=" << fp.getNewtonTotalIter()
          << ", last_newton_iter=" << fp.getNewtonLastIter() << ".";
      if (converged)
        rec << " Converged." << endrec;
      else
        {
          rec << " Not converged!!" << endrec;
          KORD_RAISE_X("Fix point calculation not converged", KORD_FP_NOT_CONV);
        }

      {
        JournalRecordPair recp(journal);
        recp << "Centralizing about fix-point." << endrec;
        fdr = std::make_unique<FoldDecisionRule>(*fdr, model.getSteady());
      }
    }
}

/* Here we simply make a new hardcopy of the given rule ‘rule_ders’,
   and make a new container of in-place subtensors of the derivatives
   corresponding to forward looking variables. The given container comes
   from a temporary object and will be destroyed. */

void
Approximation::saveRuleDerivs(const FGSContainer& g)
{
  rule_ders = std::make_unique<FGSContainer>(g);
  rule_ders_s = std::make_unique<FGSContainer>(4);
  rule_ders_ss = std::make_unique<FGSContainer>(4);
  for (auto& run : *rule_ders)
    {
      auto ten_s = std::make_unique<FGSTensor>(ypart.nstat, ypart.nys(), *(run.second));
      rule_ders_s->insert(std::move(ten_s));
      auto ten_ss
          = std::make_unique<FGSTensor>(ypart.nstat + ypart.npred, ypart.nyss(), *(run.second));
      rule_ders_ss->insert(std::move(ten_ss));
    }
}

/* This method calculates a shift of the system equations due to integrating
   shocks at a given σ and current steady state. More precisely, if

    F(y,u,u′,σ)=f(g**(g*(y,u,σ),u′,σ),g(y,u,σ),y,u)

   then the method returns a vector

       σᵈ
    ∑  ── [F_u′ᵈ]_α₁…α_d Σ^α₁…α_d
   d=1 d!

   For a calculation of [F_u′ᵈ] we use ZAuxContainer, so we create its object.
   In each cycle we calculate [F_u′ᵈ], and then multiply with the shocks, and
   add the σᵈ/d! multiple to the result. */

void
Approximation::calcStochShift(Vector& out, double at_sigma) const
{
  KORD_RAISE_IF(out.length() != ypart.ny(),
                "Wrong length of output vector for Approximation::calcStochShift");
  out.zeros();

  ZAuxContainer zaux(rule_ders_ss.get(), ypart.nyss(), ypart.ny(), ypart.nys(), model.nexog());

  int dfac = 1;
  for (int d = 1; d <= rule_ders->getMaxDim(); d++, dfac *= d)
    if (KOrder::is_even(d))
      {
        Symmetry sym {0, d, 0, 0};

        // calculate F_u′ᵈ via ZAuxContainer
        auto ten = std::make_unique<FGSTensor>(ypart.ny(), TensorDimens(sym, nvs));
        ten->zeros();
        for (int l = 1; l <= d; l++)
          {
            const FSSparseTensor& f = model.getModelDerivatives().get(Symmetry {l});
            zaux.multAndAdd(f, *ten);
          }

        // multiply with shocks and add to result
        auto tmp
            = std::make_unique<FGSTensor>(ypart.ny(), TensorDimens(Symmetry {0, 0, 0, 0}, nvs));
        tmp->zeros();
        ten->contractAndAdd(1, *tmp, mom.get(Symmetry {d}));

        out.add(pow(at_sigma, d) / dfac, tmp->getData());
      }
}

/* This method calculates and reports

              σᵈ
    f(ȳ) + ∑  ── [F_u′ᵈ]_α₁…α_d Σ^α₁…α_d
          d=1 d!

   at ȳ, zero shocks and σ. This number should be zero.

   We evaluate the error both at a given σ and σ=1.0. */

void
Approximation::check(double at_sigma) const
{
  Vector stoch_shift(ypart.ny());
  Vector system_resid(ypart.ny());
  Vector xx(model.nexog());
  xx.zeros();
  model.evaluateSystem(system_resid, model.getSteady(), xx);
  calcStochShift(stoch_shift, at_sigma);
  stoch_shift.add(1.0, system_resid);
  JournalRecord rec1(journal);
  rec1 << "Error of current approximation for shocks at sigma " << at_sigma << " is "
       << stoch_shift.getMax() << endrec;
  calcStochShift(stoch_shift, 1.0);
  stoch_shift.add(1.0, system_resid);
  JournalRecord rec2(journal);
  rec2 << "Error of current approximation for full shocks is " << stoch_shift.getMax() << endrec;
}

/* The method returns unconditional variance of endogenous variables
   based on the first order. The first order approximation looks like

    ŷₜ = g_y* ŷ*ₜ₋₁ + gᵤ uₜ

   where ŷ denotes a deviation from the steady state. It can be written as

    ŷₜ = (0 g_y* 0) ŷₜ₋₁ + gᵤ uₜ

   which yields unconditional covariance V for which

    V = GVGᵀ + gᵤ Σ gᵤᵀ

   where G=(0 g_y* 0) and Σ is the covariance of the shocks.

   For solving this Lyapunov equation we use the Sylvester module, which
   solves equation of the type

    AX + BX(C⊗…⊗C) = D

   So we invoke the Sylvester solver for the first dimension with A = I,
   B = −G, C = Gᵀ and D = gᵤ Σ gᵤᵀ. */

TwoDMatrix
Approximation::calcYCov() const
{
  const TwoDMatrix& gy = rule_ders->get(Symmetry {1, 0, 0, 0});
  const TwoDMatrix& gu = rule_ders->get(Symmetry {0, 1, 0, 0});
  TwoDMatrix G(model.numeq(), model.numeq());
  G.zeros();
  G.place(gy, 0, model.nstat());
  TwoDMatrix B(const_cast<const TwoDMatrix&>(G));
  B.mult(-1.0);
  TwoDMatrix C(transpose(G));
  TwoDMatrix A(model.numeq(), model.numeq());
  A.zeros();
  for (int i = 0; i < model.numeq(); i++)
    A.get(i, i) = 1.0;

  TwoDMatrix X((gu * model.getVcov()) * transpose(gu));

  GeneralSylvester gs(1, model.numeq(), model.numeq(), 0, A.getData(), B.getData(), C.getData(),
                      X.getData());
  gs.solve();

  return X;
}
