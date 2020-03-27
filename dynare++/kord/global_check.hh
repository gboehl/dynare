/*
 * Copyright © 2005 Ondra Kamenik
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

// Global check

/* The purpose of this file is to provide classes for checking error of
   approximation. If yₜ=g(y*ₜ₋₁,u) is an approximate solution, then we check
   for the error of residuals of the system equations. Let
   F(y*,u,u′)=f(g**(g*(y*,u′),u),g(y*,u),y*,u), then we calculate the integral:

    𝔼ₜ[F(y*,u,u′)]

   which we want to be zero for all y* and u.

   There are a few possibilities how and where the integral is evaluated.
   Currently we offer the following ones:

   — Along shocks. The y* is set to the steady state, and u is set to zero but
     one element is going from minus through plus shocks in few steps. The user
     gives the scaling factor, for instance the interval [−3σ,3σ] (where σ is
     one standard error of the shock), and a number of steps. This is repeated
     for each shock (element of the u vector).

   — Along simulation. Some random simulation is run, and for each realization
     of y* and u along the path we evaluate the residual.

   — On ellipse. Let V=AAᵀ be a covariance matrix of the predetermined
     variables y* based on linear approximation, then we calculate integral for
     points on the ellipse { Ax | ‖x‖₂=1 }. The points are selected by means of
     low discrepancy method and polar transformation. The shock u are zeros.

   — Unconditional distribution.
*/

#ifndef GLOBAL_CHECK_H
#define GLOBAL_CHECK_H

#include <matio.h>

#include <memory>

#include "vector_function.hh"
#include "quadrature.hh"

#include "dynamic_model.hh"
#include "journal.hh"
#include "approximation.hh"

/* This is a class for implementing the VectorFunction interface evaluating the
   residual of equations, this is F(y*,u,u′)=f(g**(g*(y*,u),u′),y*,u) is
   written as a function of u′.

   When the object is constructed, one has to specify (y*,u), this is done by
   the setYU() method. The object has basically two states. One is after
   construction and before the call to setYU(). The second is after the call to
   setYU(). We distinguish between the two states, an object in the second
   state contains ‘yplus’, ‘ystar’, ‘u’, and ‘hss’.

   The vector ‘yplus’ is g*(y*,u). ‘ystar’ is y*, and polynomial ‘hss’ is
   partially evaluated g**(yplus, u).

   The pointer to DynamicModel is important, since the DynamicModel evaluates
   the function f. When copying the object, we have to make also a copy of
   DynamicModel. */

class ResidFunction : public VectorFunction
{
protected:
  const Approximation &approx;
  std::unique_ptr<DynamicModel> model;
  std::unique_ptr<Vector> yplus, ystar, u;
  std::unique_ptr<FTensorPolynomial> hss;
public:
  ResidFunction(const Approximation &app);
  ResidFunction(const ResidFunction &rf);

  std::unique_ptr<VectorFunction>
  clone() const override
  {
    return std::make_unique<ResidFunction>(*this);
  }
  void eval(const Vector &point, const ParameterSignal &sig, Vector &out) override;
  void setYU(const ConstVector &ys, const ConstVector &xx);
};

/* This is a ResidFunction wrapped with GaussConverterFunction. */

class GResidFunction : public GaussConverterFunction
{
public:
  GResidFunction(const Approximation &app)
    : GaussConverterFunction(std::make_unique<ResidFunction>(app), app.getModel().getVcov())
  {
  }
  std::unique_ptr<VectorFunction>
  clone() const override
  {
    return std::make_unique<GResidFunction>(*this);
  }
  void
  setYU(const ConstVector &ys, const ConstVector &xx)
  {
    dynamic_cast<ResidFunction *>(func)->setYU(ys, xx);
  }
};

/* This is a class encapsulating checking algorithms. Its core routine is
   check(), which calculates integral 𝔼[F(y*,u,u′) | y*,u] for given
   realizations of y* and u. The both are given in matrices. The methods
   checking along shocks, on ellipse and anlong a simulation path, just fill
   the matrices and call the core check().

   The method checkUnconditionalAndSave() evaluates unconditional 𝔼[F(y,u,u′)].

   The object also maintains a set of GResidFunction functions ‘vfs’ in order
   to save (possibly expensive) copying of DynamicModel’s. */

class GlobalChecker
{
  const Approximation &approx;
  const DynamicModel &model;
  Journal &journal;
  GResidFunction rf;
  VectorFunctionSet vfs;
public:
  GlobalChecker(const Approximation &app, int n, Journal &jr)
    : approx(app), model(approx.getModel()), journal(jr),
      rf(approx), vfs(rf, n)
  {
  }
  void check(int max_evals, const ConstTwoDMatrix &y,
             const ConstTwoDMatrix &x, TwoDMatrix &out);
  void checkAlongShocksAndSave(mat_t *fd, const std::string &prefix,
                               int m, double mult, int max_evals);
  void checkOnEllipseAndSave(mat_t *fd, const std::string &prefix,
                             int m, double mult, int max_evals);
  void checkAlongSimulationAndSave(mat_t *fd, const std::string &prefix,
                                   int m, int max_evals);
  void checkUnconditionalAndSave(mat_t *fd, const std::string &prefix,
                                 int m, int max_evals);
protected:
  void check(const Quadrature &quad, int level,
             const ConstVector &y, const ConstVector &x, Vector &out);
};

/* Signalled resid function. Not implemented yet. todo: */

class ResidFunctionSig : public ResidFunction
{
public:
  ResidFunctionSig(const Approximation &app, const Vector &ys, const Vector &xx);
};

#endif
