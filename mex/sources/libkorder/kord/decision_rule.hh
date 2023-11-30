/*
 * Copyright © 2004 Ondra Kamenik
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

// Decision rule and simulation

/* The main purpose of this file is a decision rule representation which can
   run a simulation. So we define an interface for classes providing
   realizations of random shocks, and define the class DecisionRule. The latter
   basically takes tensor container of derivatives of policy rules, and adds
   them up with respect to σ. The class allows to specify the σ different from
   1.

   In addition, we provide classes for running simulations and storing the
   results, calculating some statistics and generating IRF. The class
   DRFixPoint allows for calculation of the fix point of a given decision
   rule. */

#ifndef DECISION_RULE_H
#define DECISION_RULE_H

#include "kord_exception.hh"
#include "korder.hh"

#include <memory>
#include <random>
#include <string>

/* This class is an abstract interface to decision rule. Its main purpose is to
   define a common interface for simulation of a decision rule. We need only a
   simulate, evaluate, centralized clone and output method. */
class DecisionRule
{
public:
  enum class emethod
  {
    horner,
    trad
  };
  virtual ~DecisionRule() = default;

  /* primitive evaluation (it takes a vector of state variables (predetermined,
     both and shocks) and returns the next period variables. Both input and
     output are in deviations from the rule's steady. */
  virtual void eval(emethod em, Vector& out, const ConstVector& v) const = 0;

  /* makes only one step of simulation (in terms of absolute values, not
     deviations) */
  virtual void evaluate(emethod em, Vector& out, const ConstVector& ys, const ConstVector& u) const
      = 0;

  /* returns a new copy of the decision rule, which is centralized about
     provided fix-point */
  virtual std::unique_ptr<DecisionRule> centralizedClone(const Vector& fixpoint) const = 0;

  virtual const Vector& getSteady() const = 0;
  virtual int nexog() const = 0;
  virtual const PartitionY& getYPart() const = 0;
};

/* The main purpose of this class is to implement DecisionRule interface, which
   is a simulation. To be able to do this we have to know the partitioning of
   state vector y since we will need to pick only predetermined part y*. Also,
   we need to know the steady state.

   The decision rule will take the form:

             ₙ                   ᵢ  ⎡y*ₜ₋₁ − ȳ*⎤αₘ
    yₜ − ȳ = ∑  [g_(yu)ⁱ]_α₁…αᵢ  ∏  ⎢          ⎥
            ⁱ⁼⁰                 ᵐ⁼¹ ⎣    uₜ    ⎦

   where the tensors [g_(yu)ⁱ] are tensors of the constructed container, and ȳ
   is the steady state.

   If we know the fix point of the rule (conditional zero shocks) ỹ, the rule
   can be transformed to so called “centralized” form. This is very similar to
   the form above but the zero dimensional tensor is zero:

             ₙ                   ᵢ  ⎡y*ₜ₋₁ − ỹ*⎤αₘ
    yₜ − ỹ = ∑  [g_(yu)ⁱ]_α₁…αᵢ  ∏  ⎢          ⎥
            ⁱ⁼¹                 ᵐ⁼¹ ⎣    uₜ    ⎦

   We provide a method and a constructor to transform a rule to the centralized
   form.

   The class is templated, the template argument is either Storage::fold or
   Storage::unfold. So, there are two implementations of the DecisionRule
   interface. */

template<Storage t>
class DecisionRuleImpl : public ctraits<t>::Tpol, public DecisionRule
{
protected:
  using _Tpol = typename ctraits<t>::Tpol;
  using _Tg = typename ctraits<t>::Tg;
  using _TW = typename ctraits<t>::TW;
  using _Ttensor = typename ctraits<t>::Ttensor;
  using _Ttensym = typename ctraits<t>::Ttensym;
  const Vector ysteady;
  const PartitionY ypart;
  const int nu;

public:
  DecisionRuleImpl(const _Tpol& pol, const PartitionY& yp, int nuu, const ConstVector& ys) :
      ctraits<t>::Tpol(pol), ysteady(ys), ypart(yp), nu(nuu)
  {
  }
  DecisionRuleImpl(_Tpol& pol, const PartitionY& yp, int nuu, const ConstVector& ys) :
      ctraits<t>::Tpol(0, yp.ny(), pol), ysteady(ys), ypart(yp), nu(nuu)
  {
  }
  DecisionRuleImpl(const _Tg& g, const PartitionY& yp, int nuu, const ConstVector& ys,
                   double sigma) :
      ctraits<t>::Tpol(yp.ny(), yp.nys() + nuu), ysteady(ys), ypart(yp), nu(nuu)
  {
    fillTensors(g, sigma);
  }
  DecisionRuleImpl(const _Tg& g, const PartitionY& yp, int nuu, const ConstVector& ys, double sigma,
                   bool pruning) :
      ctraits<t>::Tpol(yp.ny(), yp.nys() + nuu), ysteady(ys), ypart(yp), nu(nuu)
  {
    if (pruning)
      fillTensorsPruning(g);
    else
      fillTensors(g, sigma);
  }

  DecisionRuleImpl(const _TW& W, int nys, int nuu, const ConstVector& ys) :
      ctraits<t>::Tpol(1, nys + nuu), ysteady(ys), nu(nuu)
  {
    fillTensors(W, nys);
  }
  DecisionRuleImpl(const DecisionRuleImpl<t>& dr, const ConstVector& fixpoint) :
      ctraits<t>::Tpol(dr.ypart.ny(), dr.ypart.nys() + dr.nu),
      ysteady(fixpoint),
      ypart(dr.ypart),
      nu(dr.nu)
  {
    centralize(dr);
  }
  const Vector&
  getSteady() const override
  {
    return ysteady;
  }
  void evaluate(emethod em, Vector& out, const ConstVector& ys,
                const ConstVector& u) const override;
  std::unique_ptr<DecisionRule> centralizedClone(const Vector& fixpoint) const override;

  int
  nexog() const override
  {
    return nu;
  }
  const PartitionY&
  getYPart() const override
  {
    return ypart;
  }

protected:
  void fillTensors(const _Tg& g, double sigma);
  void fillTensorsPruning(const _Tg& g);
  void fillTensors(const _TW& W, int nys);
  void centralize(const DecisionRuleImpl& dr);

public:
  void eval(emethod em, Vector& out, const ConstVector& v) const override;
};

/* Here we have to fill the tensor polynomial. This involves two separated
   actions. The first is to evaluate the approximation at a given σ, the second
   is to compile the tensors [g_(yu)ⁱ⁺ʲ] from [g_yⁱuʲ]. The first action is
   done here, the second is done by method addSubTensor() of a full symmetry
   tensor.

   The way how the evaluation is done is described here:

   The q-order approximation to the solution can be written as:

                  ⎡                                                               ⎤
             q  1 ⎢       ⎛  l  ⎞⎡        ⎤            ᵢ ⎡          ⎤αₘ ⱼ ⎡  ⎤βₘ  ⎥
    yₜ − ȳ = ∑  ──⎢   ∑   ⎢     ⎥⎢g_yⁱuʲσᵏ⎥            ∏ ⎢y*ₜ₋₁ − ȳ*⎥   ∏ ⎢uₜ⎥  σᵏ⎥
            ˡ⁼¹ l!⎢ⁱ⁺ʲ⁺ᵏ⁼ˡ⎝i,j,k⎠⎣        ⎦α₁…αⱼβ₁…βⱼ ᵐ⁼¹⎣          ⎦  ᵐ⁼¹⎣  ⎦    ⎥
                  ⎣                                                               ⎦

               ⎡           ⎡                                   ⎤                           ⎤
             q ⎢      ⎛i+j⎞⎢ₗ₋ᵢ₋ⱼ 1  ⎛l⎞ ⎡        ⎤            ⎥  ᵢ ⎡          ⎤αₘ ⱼ ⎡  ⎤βₘ⎥
           = ∑ ⎢  ∑   ⎢   ⎥⎢  ∑   ── ⎢ ⎥ ⎢g_yⁱuʲσᵏ⎥          σᵏ⎥  ∏ ⎢y*ₜ₋₁ − ȳ*⎥   ∏ ⎢uₜ⎥  ⎥
            ˡ⁼¹⎢i+j≤l ⎝ i ⎠⎢ ᵏ⁼⁰  l! ⎝k⎠ ⎣        ⎦α₁…αⱼβ₁…βⱼ  ⎥ ᵐ⁼¹⎣          ⎦  ᵐ⁼¹⎣  ⎦  ⎥
               ⎣           ⎣                                   ⎦                           ⎦

   This means that for each i+j+k=l we have to add

    1  ⎛l⎞                     1
    ── ⎢ ⎥ [g_yⁱuʲσᵏ]·σᵏ = ──────── [g_yⁱuʲσᵏ]·σᵏ
    l! ⎝k⎠                 (i+j)!k!

   to [g_(yu)ⁱ⁺ʲ].
                                         ⎛i+j⎞
   In addition, note that the multiplier ⎝ k ⎠ is applied when the fully symmetric
   tensor [g_(yu)ⁱ⁺ʲ] is evaluated.

   So we go through i+j=d=0…q and in each loop we form the fully symmetric
   tensor [g_(yu)ᵈ] and insert it to the container. */

template<Storage t>
void
DecisionRuleImpl<t>::fillTensors(const _Tg& g, double sigma)
{
  IntSequence tns {ypart.nys(), nu};
  int dfact = 1;
  for (int d = 0; d <= g.getMaxDim(); d++, dfact *= d)
    {
      auto g_yud = std::make_unique<_Ttensym>(ypart.ny(), ypart.nys() + nu, d);
      g_yud->zeros();

      // fill tensor of ‘g_yud’ of dimension ‘d’
      /* Here we have to fill the tensor [g_(yu)ᵈ]. So we go through all pairs
         (i,j) such that i+j=d, and through all k from zero up to maximal
         dimension minus d. In this way we go through all symmetries of
         [g_yⁱuʲσᵏ] which will be added to [g_(yu)ᵈ].

         Note that at the beginning, ‘dfact’ is a factorial of ‘d’. We
         calculate ‘kfact’ is equal to k!. As indicated in
         DecisionRuleImpl::fillTensors(), the added tensor is thus multiplied
         with 1/(d!k!)·σᵏ. */

      for (int i = 0; i <= d; i++)
        {
          int j = d - i;
          int kfact = 1;
          _Ttensor tmp(ypart.ny(), TensorDimens(Symmetry {i, j}, tns));
          tmp.zeros();
          for (int k = 0; k + d <= g.getMaxDim(); k++, kfact *= k)
            if (Symmetry sym {i, j, 0, k}; g.check(sym))
              {
                double mult = pow(sigma, k) / dfact / kfact;
                tmp.add(mult, g.get(sym));
              }
          g_yud->addSubTensor(tmp);
        }

      this->insert(std::move(g_yud));
    }
}

template<Storage t>
void
DecisionRuleImpl<t>::fillTensorsPruning(const _Tg& g)
{
  IntSequence tns {ypart.nys(), nu, 1};
  int dfact = 1;
  for (int d = 0; d <= g.getMaxDim(); d++, dfact *= d)
    {
      auto g_yusd = std::make_unique<_Ttensym>(ypart.ny(), ypart.nys() + nu + 1, d);
      g_yusd->zeros();
      // fill tensor of ‘g_yusd’ of dimension ‘d’
      /*
        Here we have to fill the tensor [g_(yuσ)ᵈ]. So we go through all pairs
        (i,j,k) such that i+j+k=d. We weight it with 1/(i+j+k)! The factorial
        is denoted dfact.
      */
      for (int i = 0; i <= d; i++)
        {
          for (int j = 0; j <= d - i; j++)
            {
              int k = d - i - j;
              _Ttensor tmp(ypart.ny(), TensorDimens(Symmetry {i, j, k}, tns));
              tmp.zeros();
              if (Symmetry sym {i, j, 0, k}; g.check(sym))
                {
                  double mult = 1.0 / dfact;
                  // mexPrintf("Symmetry found: %d %d %d %.2f %d\n", i, j, k, mult, kfact);
                  tmp.add(mult, g.get(sym));
                }
              g_yusd->addSubTensor(tmp);
            }
        }
      this->insert(std::move(g_yusd));
    }
}

template<Storage t>
void
DecisionRuleImpl<t>::fillTensors(const _TW& W, int nys)
{
  IntSequence tns {nys, nu};
  int dfact = 1;
  for (int d = 0; d <= W.getMaxDim(); d++, dfact *= d)
    {
      auto W_yud = std::make_unique<_Ttensym>(1, nys + nu, d);
      W_yud->zeros();

      // fill tensor of ‘g_yud’ of dimension ‘d’
      /* Here we have to fill the tensor [g_(yu)ᵈ]. So we go through all pairs
         (i,j) such that i+j=d, and through all k from zero up to maximal
         dimension minus d. In this way we go through all symmetries of
         [g_yⁱuʲσᵏ] which will be added to [g_(yu)ᵈ].

         Note that at the beginning, ‘dfact’ is a factorial of ‘d’. We
         calculate ‘kfact’ is equal to k!. As indicated in
         DecisionRuleImpl::fillTensors(), the added tensor is thus multiplied
         with 1/(d!k!)·σᵏ. */

      for (int i = 0; i <= d; i++)
        {
          int j = d - i;
          int kfact = 1;
          _Ttensor tmp(1, TensorDimens(Symmetry {i, j}, tns));
          tmp.zeros();
          for (int k = 0; k + d <= W.getMaxDim(); k++, kfact *= k)
            if (Symmetry sym {i, j, 0, k}; W.check(sym))
              {
                double mult = 1.0 / dfact / kfact;
                tmp.add(mult, W.get(sym));
              }
          W_yud->addSubTensor(tmp);
        }

      this->insert(std::move(W_yud));
    }
}

/* The centralization is straightforward. We suppose here that the object’s
   steady state is the fix point ỹ. It is clear that the new derivatives
   [g~_(yu)ⁱ] will be equal to the derivatives of the original decision rule
   ‘dr’ at the new steady state ỹ. So, the new derivatives are obtained by
   derivating the given decision rule ‘dr’ and evaluating its polynomial at:

             ⎡ỹ* − ȳ*⎤
    dstate = ⎢       ⎥,
             ⎣   0   ⎦

   where ȳ is the steady state of the original rule ‘dr’. */

template<Storage t>
void
DecisionRuleImpl<t>::centralize(const DecisionRuleImpl& dr)
{
  Vector dstate(ypart.nys() + nu);
  dstate.zeros();
  Vector dstate_star(dstate, 0, ypart.nys());
  ConstVector newsteady_star(ysteady, ypart.nstat, ypart.nys());
  ConstVector oldsteady_star(dr.ysteady, ypart.nstat, ypart.nys());
  dstate_star.add(1.0, newsteady_star);
  dstate_star.add(-1.0, oldsteady_star);

  _Tpol pol(dr);
  int dfac = 1;
  for (int d = 1; d <= dr.getMaxDim(); d++, dfac *= d)
    {
      pol.derivative(d - 1);
      auto der = pol.evalPartially(d, dstate);
      der->mult(1.0 / dfac);
      this->insert(std::move(der));
    }
}

/* This is one period evaluation of the decision rule. The simulation is a
   sequence of repeated one period evaluations with a difference, that the
   steady state (fix point) is cancelled and added once. Hence we have two
   special methods. */

template<Storage t>
void
DecisionRuleImpl<t>::evaluate(emethod em, Vector& out, const ConstVector& ys,
                              const ConstVector& u) const
{
  KORD_RAISE_IF(ys.length() != ypart.nys() || u.length() != nu,
                "Wrong dimensions of input vectors in DecisionRuleImpl::evaluate");
  KORD_RAISE_IF(out.length() != ypart.ny(),
                "Wrong dimension of output vector in DecisionRuleImpl::evaluate");
  ConstVector ysteady_pred(ysteady, ypart.nstat, ypart.nys());
  Vector ys_u(ypart.nys() + nu);
  Vector ys_u1(ys_u, 0, ypart.nys());
  ys_u1 = ys;
  ys_u1.add(-1.0, ysteady_pred);
  Vector ys_u2(ys_u, ypart.nys(), nu);
  ys_u2 = u;
  eval(em, out, ys_u);
  out.add(1.0, ysteady);
}

/* This is easy. We just return the newly created copy using the centralized
   constructor. */

template<Storage t>
std::unique_ptr<DecisionRule>
DecisionRuleImpl<t>::centralizedClone(const Vector& fixpoint) const
{
  return std::make_unique<DecisionRuleImpl<t>>(*this, fixpoint);
}

/* Here we only encapsulate two implementations to one, deciding according to
   the parameter. */

template<Storage t>
void
DecisionRuleImpl<t>::eval(emethod em, Vector& out, const ConstVector& v) const
{
  if (em == emethod::horner)
    _Tpol::evalHorner(out, v);
  else
    _Tpol::evalTrad(out, v);
}

/* This is exactly the same as DecisionRuleImpl<Storage::fold>. The only
   difference is that we have a conversion from UnfoldDecisionRule, which is
   exactly DecisionRuleImpl<Storage::unfold>. */

class UnfoldDecisionRule;
class FoldDecisionRule : public DecisionRuleImpl<Storage::fold>
{
  friend class UnfoldDecisionRule;

public:
  FoldDecisionRule(const ctraits<Storage::fold>::Tpol& pol, const PartitionY& yp, int nuu,
                   const ConstVector& ys) :
      DecisionRuleImpl<Storage::fold>(pol, yp, nuu, ys)
  {
  }
  FoldDecisionRule(ctraits<Storage::fold>::Tpol& pol, const PartitionY& yp, int nuu,
                   const ConstVector& ys) :
      DecisionRuleImpl<Storage::fold>(pol, yp, nuu, ys)
  {
  }
  FoldDecisionRule(const ctraits<Storage::fold>::Tg& g, const PartitionY& yp, int nuu,
                   const ConstVector& ys, double sigma) :
      DecisionRuleImpl<Storage::fold>(g, yp, nuu, ys, sigma)
  {
  }
  FoldDecisionRule(const ctraits<Storage::fold>::Tg& g, const PartitionY& yp, int nuu,
                   const ConstVector& ys, double sigma, bool pruning) :
      DecisionRuleImpl<Storage::fold>(g, yp, nuu, ys, sigma, pruning)
  {
  }
  FoldDecisionRule(const ctraits<Storage::fold>::TW& W, int nys, int nuu, const ConstVector& ys) :
      DecisionRuleImpl<Storage::fold>(W, nys, nuu, ys)
  {
  }
  FoldDecisionRule(const DecisionRuleImpl<Storage::fold>& dr, const ConstVector& fixpoint) :
      DecisionRuleImpl<Storage::fold>(dr, fixpoint)
  {
  }
  FoldDecisionRule(const UnfoldDecisionRule& udr);
};

/* This is exactly the same as DecisionRuleImpl<Storage::unfold>, but with a
   conversion from FoldDecisionRule, which is exactly
   DecisionRuleImpl<Storage::fold>. */

class UnfoldDecisionRule : public DecisionRuleImpl<Storage::unfold>
{
  friend class FoldDecisionRule;

public:
  UnfoldDecisionRule(const ctraits<Storage::unfold>::Tpol& pol, const PartitionY& yp, int nuu,
                     const ConstVector& ys) :
      DecisionRuleImpl<Storage::unfold>(pol, yp, nuu, ys)
  {
  }
  UnfoldDecisionRule(ctraits<Storage::unfold>::Tpol& pol, const PartitionY& yp, int nuu,
                     const ConstVector& ys) :
      DecisionRuleImpl<Storage::unfold>(pol, yp, nuu, ys)
  {
  }
  UnfoldDecisionRule(const ctraits<Storage::unfold>::Tg& g, const PartitionY& yp, int nuu,
                     const ConstVector& ys, double sigma) :
      DecisionRuleImpl<Storage::unfold>(g, yp, nuu, ys, sigma)
  {
  }
  UnfoldDecisionRule(const DecisionRuleImpl<Storage::unfold>& dr, const ConstVector& fixpoint) :
      DecisionRuleImpl<Storage::unfold>(dr, fixpoint)
  {
  }
  UnfoldDecisionRule(const FoldDecisionRule& udr);
};

/* This class serves for calculation of the fix point of the decision rule
   given that the shocks are zero. The class is very similar to the
   DecisionRuleImpl. Besides the calculation of the fix point, the only
   difference between DRFixPoint and DecisionRuleImpl is that the derivatives
   wrt. shocks are ignored (since shocks are zero during the calculations).
   That is why have a different fillTensor() method.

   The solution algorithm is Newton and is described in
   DRFixPoint::solveNewton(). It solves F(y)=0, where F=g(y,0)−y. The function
   F is given by its derivatives ‘bigf’. The Jacobian of the solved system is
   given by derivatives stored in ‘bigfder’. */

template<Storage t>
class DRFixPoint : public ctraits<t>::Tpol
{
  using _Tpol = typename ctraits<t>::Tpol;
  using _Tg = typename ctraits<t>::Tg;
  using _Ttensor = typename ctraits<t>::Ttensor;
  using _Ttensym = typename ctraits<t>::Ttensym;
  constexpr static int max_iter = 10000;
  constexpr static int max_newton_iter = 50;
  constexpr static int newton_pause = 100;
  constexpr static double tol = 1e-10;
  const Vector ysteady;
  const PartitionY ypart;
  std::unique_ptr<_Tpol> bigf;
  std::unique_ptr<_Tpol> bigfder;

public:
  using emethod = typename DecisionRule::emethod;
  DRFixPoint(const _Tg& g, const PartitionY& yp, const Vector& ys, double sigma);

  bool calcFixPoint(Vector& out);

  int
  getNumIter() const
  {
    return iter;
  }
  int
  getNewtonLastIter() const
  {
    return newton_iter_last;
  }
  int
  getNewtonTotalIter() const
  {
    return newton_iter_total;
  }

protected:
  void fillTensors(const _Tg& g, double sigma);
  bool solveNewton(Vector& y);

private:
  int iter;
  int newton_iter_last;
  int newton_iter_total;
};

/* Here we have to setup the function F=g(y,0)−y and ∂F/∂y. The former is taken
   from the given derivatives of g where a unit matrix is subtracted from the
   first derivative (Symmetry{1}). Then the derivative of the F polynomial is
   calculated. */

template<Storage t>
DRFixPoint<t>::DRFixPoint(const _Tg& g, const PartitionY& yp, const Vector& ys, double sigma) :
    ctraits<t>::Tpol(yp.ny(), yp.nys()), ysteady(ys), ypart(yp)
{
  fillTensors(g, sigma);
  _Tpol yspol(ypart.nstat, ypart.nys(), *this);
  bigf = std::make_unique<_Tpol>(const_cast<const _Tpol&>(yspol));
  _Ttensym& frst = bigf->get(Symmetry {1});
  for (int i = 0; i < ypart.nys(); i++)
    frst.get(i, i) = frst.get(i, i) - 1;
  bigfder = std::make_unique<_Tpol>(*bigf, 0);
}

/* Here we fill the tensors for the DRFixPoint class. We ignore the derivatives
   [g_yⁱuʲσᵏ] for which j>0. So we go through all dimensions ‘d’, and all ‘k’
   such that ‘d+k’ is between the maximum dimension and ‘d’, and add
   σᵏ/(d!k!)[g_yᵈσᵏ] to the tensor [g_yᵈ]. */

template<Storage t>
void
DRFixPoint<t>::fillTensors(const _Tg& g, double sigma)
{
  int dfact = 1;
  for (int d = 0; d <= g.getMaxDim(); d++, dfact *= d)
    {
      auto g_yd = std::make_unique<_Ttensym>(ypart.ny(), ypart.nys(), d);
      g_yd->zeros();
      int kfact = 1;
      for (int k = 0; d + k <= g.getMaxDim(); k++, kfact *= k)
        if (g.check(Symmetry {d, 0, 0, k}))
          {
            const _Ttensor& ten = g.get(Symmetry {d, 0, 0, k});
            double mult = pow(sigma, k) / dfact / kfact;
            g_yd->add(mult, ten);
          }
      this->insert(std::move(g_yd));
    }
}

/* This tries to solve polynomial equation F(y)=0, where F polynomial is ‘bigf’
   and its derivative is in ‘bigfder’. It returns true if the Newton converged.
   The method takes the given vector as initial guess, and rewrites it with a
   solution. The method guarantees to return the vector, which has smaller norm
   of the residual. That is why the input/output vector ‘y’ is always changed.

   The method proceeds with a Newton step, if the Newton step improves the
   residual error. So we track residual errors in ‘flastnorm’ and ‘fnorm’
   (former and current). In addition, at each step we search for an
   underrelaxation parameter ‘urelax’, which improves the residual. If ‘urelax’
   is less that ‘urelax_threshold’, we stop searching and stop the Newton. */

template<Storage t>
bool
DRFixPoint<t>::solveNewton(Vector& y)
{
  const double urelax_threshold = 1.e-5;
  Vector sol(const_cast<const Vector&>(y));
  Vector delta(y.length());
  newton_iter_last = 0;
  bool delta_finite = true;
  double flastnorm = 0.0;
  double fnorm = 0.0;
  bool converged = false;
  double urelax = 1.0;

  do
    {
      auto jacob = bigfder->evalPartially(1, sol);
      bigf->evalHorner(delta, sol);
      if (newton_iter_last == 0)
        flastnorm = delta.getNorm();
      delta_finite = delta.isFinite();
      if (delta_finite)
        {
          ConstTwoDMatrix(*jacob).multInvLeft(delta);

          // find ‘urelax’ improving residual
          /* Here we find the ‘urelax’. We cycle as long as the new residual
             size ‘fnorm’ is greater than last residual size ‘flastnorm’. If
             the urelax is less than ‘urelax_threshold’ we give up. The
             ‘urelax’ is damped by the ratio of ‘flastnorm’ and ‘fnorm’. It the
             ratio is close to one, we damp by one half. */
          bool urelax_found = false;
          urelax = 1.0;
          while (!urelax_found && urelax > urelax_threshold)
            {
              Vector soltmp(const_cast<const Vector&>(sol));
              soltmp.add(-urelax, delta);
              Vector f(sol.length());
              bigf->evalHorner(f, soltmp);
              fnorm = f.getNorm();
              if (fnorm <= flastnorm)
                urelax_found = true;
              else
                urelax *= std::min(0.5, flastnorm / fnorm);
            }

          sol.add(-urelax, delta);
          delta_finite = delta.isFinite();
        }
      newton_iter_last++;
      converged = delta_finite && fnorm < tol;
      flastnorm = fnorm;
    }
  while (!converged && newton_iter_last < max_newton_iter && urelax > urelax_threshold);

  newton_iter_total += newton_iter_last;
  if (!converged)
    newton_iter_last = 0;
  y = const_cast<const Vector&>(sol);
  return converged;
}

/* This method solves the fix point of the no-shocks rule yₜ₊₁=f(yₜ). It
   combines dull steps with Newton attempts. The dull steps correspond to
   evaluations setting yₜ₊₁=f(yₜ). For reasonable models the dull steps
   converge to the fix-point but very slowly. That is why we make Newton
   attempt from time to time. The frequency of the Newton attempts is given by
   ‘newton_pause’. We perform the calculations in deviations from the steady
   state. So, at the end, we have to add the steady state.

   The method also sets the members ‘iter’, ‘newton_iter_last’ and
   ‘newton_iter_total’. These numbers can be examined later.

   The ‘out’ vector is not touched if the algorithm has not convered. */

template<Storage t>
bool
DRFixPoint<t>::calcFixPoint(Vector& out)
{
  KORD_RAISE_IF(out.length() != ypart.ny(), "Wrong length of out in DRFixPoint::calcFixPoint");

  Vector delta(ypart.nys());
  Vector ystar(ypart.nys());
  ystar.zeros();

  iter = 0;
  newton_iter_last = 0;
  newton_iter_total = 0;
  bool converged = false;
  do
    {
      if ((iter / newton_pause) * newton_pause == iter)
        converged = solveNewton(ystar);
      if (!converged)
        {
          bigf->evalHorner(delta, ystar);
          KORD_RAISE_IF_X(!delta.isFinite(), "NaN or Inf asserted in DRFixPoint::calcFixPoint",
                          KORD_FP_NOT_FINITE);
          ystar.add(1.0, delta);
          converged = delta.getNorm() < tol;
        }
      iter++;
    }
  while (iter < max_iter && !converged);

  if (converged)
    {
      _Tpol::evalHorner(out, ystar);
      out.add(1.0, ysteady);
    }

  return converged;
}

#endif
