// Copyright 2007, Ondra Kamenik

// Conjugate family for normal distribution

/* The main purpose here is to implement a class representing conjugate
   distributions for mean and variance of the normal distribution. The class
   has two main methods: the first one is to update itself with respect to one
   observation, the second one is to update itself with respect to anothe
   object of the class. In the both methods, the previous state of the class
   corresponds to the prior distribution, and the final state corresponds to
   the posterior distribution.

   The algebra can be found in Gelman, Carlin, Stern, Rubin (p.87). It goes as
   follows. Prior conjugate distribution takes the following form:

     Σ  ↝ InvWishart_ν₀(Λ₀⁻¹)
    μ|Σ ↝ 𝒩(μ₀,Σ/κ₀)

   If the observations are y₁…yₙ, then the posterior distribution has the same
   form with the following parameters:

          κ₀         n
    μₙ = ──── μ₀ + ──── ȳ
         κ₀+n      κ₀+n

    κₙ = κ₀ + n

    νₙ = ν₀ + n

                  κ₀·n
    Λₙ = Λ₀ + S + ──── (ȳ − μ₀)(ȳ − μ₀)ᵀ
                  κ₀+n

   where

        1  ₙ
    ȳ = ─  ∑ yᵢ
        n ⁱ⁼¹

        ₙ
    S = ∑ (yᵢ − ȳ)(yᵢ − ȳ)ᵀ
       ⁱ⁼¹
*/

#ifndef NORMAL_CONJUGATE_H
#define NORMAL_CONJUGATE_H

#include "twod_matrix.hh"

/* The class is described by the four parameters: μ, κ, ν and Λ. */

class NormalConj
{
protected:
  Vector mu;
  int kappa;
  int nu;
  TwoDMatrix lambda;
public:
  /* We provide the following constructors: The first constructs diffuse
     (Jeffrey’s) prior. It sets κ and Λ to zeros, ν to −1 and also the mean μ
     to zero (it should not be referenced). The second constructs the posterior
     using the diffuse prior and the observed data (columnwise). The third is a
     copy constructor. */
  NormalConj(int d);
  NormalConj(const ConstTwoDMatrix &ydata);
  NormalConj(const NormalConj &) = default;
  NormalConj(NormalConj &&) = default;

  virtual ~NormalConj() = default;
  void update(const ConstVector &y);
  void update(const ConstTwoDMatrix &ydata);
  void update(const NormalConj &nc);
  int
  getDim() const
  {
    return mu.length();
  }
  const Vector &
  getMean() const
  {
    return mu;
  }
  void getVariance(TwoDMatrix &v) const;
};

#endif
