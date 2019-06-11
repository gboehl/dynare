// Copyright 2007, Ondra Kamenik

#include "normal_conjugate.hh"
#include "kord_exception.hh"

// NormalConj diffuse prior constructor
NormalConj::NormalConj(int d)
  : mu(d), kappa(0), nu(-1), lambda(d, d)
{
  mu.zeros();
  lambda.zeros();
}

// NormalConj data update constructor
NormalConj::NormalConj(const ConstTwoDMatrix &ydata)
  : mu(ydata.nrows()), kappa(ydata.ncols()), nu(ydata.ncols()-1),
    lambda(ydata.nrows(), ydata.nrows())
{
  mu.zeros();
  for (int i = 0; i < ydata.ncols(); i++)
    mu.add(1.0/ydata.ncols(), ydata.getCol(i));

  lambda.zeros();
  for (int i = 0; i < ydata.ncols(); i++)
    {
      Vector diff{ydata.getCol(i)};
      diff.add(-1, mu);
      lambda.addOuter(diff);
    }
}

// NormalConj::update() one observation code
/* The method performs the following:

          κ₀         1
    μ₁ = ──── μ₀ + ──── y
         κ₀+1      κ₀+1

    κ₁ = κ₀ + 1

    ν₁ = ν₀ + 1

               κ₀
    Λ₁ = Λ₀ + ──── (y − μ₀)(y − μ₀)ᵀ
              κ₀+1
*/
void
NormalConj::update(const ConstVector &y)
{
  KORD_RAISE_IF(y.length() != mu.length(),
                "Wrong length of a vector in NormalConj::update");

  mu.mult(kappa/(1.0+kappa));
  mu.add(1.0/(1.0+kappa), y);

  Vector diff(y);
  diff.add(-1, mu);
  lambda.addOuter(diff, kappa/(1.0+kappa));

  kappa++;
  nu++;
}

// NormalConj::update() multiple observations code
/* The method evaluates the formula in the header file. */
void
NormalConj::update(const ConstTwoDMatrix &ydata)
{
  NormalConj nc(ydata);
  update(nc);
}

// NormalConj::update() with NormalConj code
void
NormalConj::update(const NormalConj &nc)
{
  double wold = static_cast<double>(kappa)/(kappa+nc.kappa);
  double wnew = 1-wold;

  mu.mult(wold);
  mu.add(wnew, nc.mu);

  Vector diff(nc.mu);
  diff.add(-1, mu);
  lambda.add(1.0, nc.lambda);
  lambda.addOuter(diff);

  kappa = kappa + nc.kappa;
  nu = nu + nc.kappa;
}

/* This returns 1/(ν−d−1)·Λ, which is the mean of the variance in the posterior
   distribution. If the number of degrees of freedom is less than d, then NaNs
   are returned. */
void
NormalConj::getVariance(TwoDMatrix &v) const
{
  if (nu > getDim()+1)
    {
      v = const_cast<const TwoDMatrix &>(lambda);
      v.mult(1.0/(nu-getDim()-1));
    }
  else
    v.nans();
}
