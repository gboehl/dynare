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
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "SymSchurDecomp.hh"

#include "global_check.hh"
#include "seed_generator.hh"

#include "smolyak.hh"
#include "product.hh"
#include "quasi_mcarlo.hh"

#include <utility>
#include <cmath>

/* Here we just set a reference to the approximation, and create a new
   DynamicModel. */

ResidFunction::ResidFunction(const Approximation &app)
  : VectorFunction(app.getModel().nexog(), app.getModel().numeq()), approx(app),
    model(app.getModel().clone())
{
}

ResidFunction::ResidFunction(const ResidFunction &rf)
  : VectorFunction(rf), approx(rf.approx), model(rf.model->clone())
{
  if (rf.yplus)
    yplus = std::make_unique<Vector>(*(rf.yplus));
  if (rf.ystar)
    ystar = std::make_unique<Vector>(*(rf.ystar));
  if (rf.u)
    u = std::make_unique<Vector>(*(rf.u));
  if (rf.hss)
    hss = std::make_unique<FTensorPolynomial>(*(rf.hss));
}

/* This sets y* and u. We have to create ‘ystar’, ‘u’, ‘yplus’ and ‘hss’. */

void
ResidFunction::setYU(const ConstVector &ys, const ConstVector &xx)
{
  ystar = std::make_unique<Vector>(ys);
  u = std::make_unique<Vector>(xx);
  yplus = std::make_unique<Vector>(model->numeq());
  approx.getFoldDecisionRule().evaluate(DecisionRule::emethod::horner,
                                        *yplus, *ystar, *u);

  // make a tensor polynomial of in-place subtensors from decision rule
  /* Note that the non-const polynomial will be used for a construction of
     ‘hss’ and will be used in a const context. So this const cast is safe.

     Note, that there is always a folded decision rule in Approximation. */
  const FoldDecisionRule &dr = approx.getFoldDecisionRule();
  FTensorPolynomial dr_ss(model->nstat()+model->npred(), model->nboth()+model->nforw(),
                          const_cast<FoldDecisionRule &>(dr));

  // make ‘ytmp_star’ be a difference of ‘yplus’ from steady
  Vector ytmp_star(ConstVector(*yplus, model->nstat(), model->npred()+model->nboth()));
  ConstVector ysteady_star(dr.getSteady(), model->nstat(),
                           model->npred()+model->nboth());
  ytmp_star.add(-1.0, ysteady_star);

  // make ‘hss’ and add steady to it
  /* Here is the const context of ‘dr_ss’. */
  hss = std::make_unique<FTensorPolynomial>(dr_ss, ytmp_star);
  ConstVector ysteady_ss(dr.getSteady(), model->nstat()+model->npred(),
                         model->nboth()+model->nforw());
  if (hss->check(Symmetry{0}))
    hss->get(Symmetry{0}).getData().add(1.0, ysteady_ss);
  else
    {
      auto ten = std::make_unique<FFSTensor>(hss->nrows(), hss->nvars(), 0);
      ten->getData() = ysteady_ss;
      hss->insert(std::move(ten));
    }
}

/* Here we evaluate the residual F(y*,u,u′). We have to evaluate ‘hss’ for
   u′=point and then we evaluate the system f. */

void
ResidFunction::eval(const Vector &point, const ParameterSignal &sig, Vector &out)
{
  KORD_RAISE_IF(point.length() != hss->nvars(),
                "Wrong dimension of input vector in ResidFunction::eval");
  KORD_RAISE_IF(out.length() != model->numeq(),
                "Wrong dimension of output vector in ResidFunction::eval");
  Vector yss(hss->nrows());
  hss->evalHorner(yss, point);
  model->evaluateSystem(out, *ystar, *yplus, yss, *u);
}

/* This checks the 𝔼[F(y*,u,u′)] for a given y* and u by integrating with a
   given quadrature. Note that the input ‘ys’ is y* not whole y. */

void
GlobalChecker::check(const Quadrature &quad, int level,
                     const ConstVector &ys, const ConstVector &x, Vector &out)
{
  for (int ifunc = 0; ifunc < vfs.getNum(); ifunc++)
    dynamic_cast<GResidFunction &>(vfs.getFunc(ifunc)).setYU(ys, x);
  quad.integrate(vfs, level, out);
}

/* This method is a bulk version of GlobalChecker::check() vector code. It
   decides between Smolyak and product quadrature according to ‘max_evals’
   constraint.

   Note that ‘y’ can be either full (all endogenous variables including static
   and forward looking), or just y* (state variables). The method is able to
   recognize it. */

void
GlobalChecker::check(int max_evals, const ConstTwoDMatrix &y,
                     const ConstTwoDMatrix &x, TwoDMatrix &out)
{
  JournalRecordPair pa(journal);
  pa << "Checking approximation error for " << y.ncols()
     << " states with at most " << max_evals << " evaluations" << endrec;

  // Decide about which type of quadrature
  GaussHermite gh;

  SmolyakQuadrature dummy_sq(model.nexog(), 1, gh);
  int smol_evals;
  int smol_level;
  dummy_sq.designLevelForEvals(max_evals, smol_level, smol_evals);

  ProductQuadrature dummy_pq(model.nexog(), gh);
  int prod_evals;
  int prod_level;
  dummy_pq.designLevelForEvals(max_evals, prod_level, prod_evals);

  bool take_smolyak = (smol_evals < prod_evals) && (smol_level >= prod_level-1);

  std::unique_ptr<Quadrature> quad;
  int lev;

  // Create the quadrature and report the decision
  if (take_smolyak)
    {
      quad = std::make_unique<SmolyakQuadrature>(model.nexog(), smol_level, gh);
      lev = smol_level;
      JournalRecord rec(journal);
      rec << "Selected Smolyak (level,evals)=(" << smol_level << ","
          << smol_evals << ") over product (" << prod_level << ","
          << prod_evals << ")" << endrec;
    }
  else
    {
      quad = std::make_unique<ProductQuadrature>(model.nexog(), gh);
      lev = prod_level;
      JournalRecord rec(journal);
      rec << "Selected product (level,evals)=(" << prod_level << ","
          << prod_evals << ") over Smolyak (" << smol_level << ","
          << smol_evals << ")" << endrec;
    }

  // Check all columns of ‘y’ and ‘x’
  int first_row = (y.nrows() == model.numeq()) ? model.nstat() : 0;
  ConstTwoDMatrix ysmat(y, first_row, 0, model.npred()+model.nboth(), y.ncols());
  for (int j = 0; j < y.ncols(); j++)
    {
      ConstVector yj{ysmat.getCol(j)};
      ConstVector xj{x.getCol(j)};
      Vector outj{out.getCol(j)};
      check(*quad, lev, yj, xj, outj);
    }
}

/* This method checks an error of the approximation by evaluating residual
   𝔼[F(y*,u,u′) | y*,u] for y* equal to the steady state, and changing u. We go
   through all elements of u and vary them from −mult·σ to mult·σ in ‘m’
   steps. */

void
GlobalChecker::checkAlongShocksAndSave(mat_t *fd, const std::string &prefix,
                                       int m, double mult, int max_evals)
{
  JournalRecordPair pa(journal);
  pa << "Calculating errors along shocks +/- "
     << mult << " std errors, granularity " << m << endrec;

  // Setup ‘y_mat’ of steady states for checking
  TwoDMatrix y_mat(model.numeq(), 2*m*model.nexog()+1);
  for (int j = 0; j < 2*m*model.nexog()+1; j++)
    {
      Vector yj{y_mat.getCol(j)};
      yj = model.getSteady();
    }

  // Setup ‘exo_mat’ for checking
  TwoDMatrix exo_mat(model.nexog(), 2*m*model.nexog()+1);
  exo_mat.zeros();
  for (int ishock = 0; ishock < model.nexog(); ishock++)
    {
      double max_sigma = sqrt(model.getVcov().get(ishock, ishock));
      for (int j = 0; j < 2*m; j++)
        {
          int jmult = (j < m) ? j-m : j-m+1;
          exo_mat.get(ishock, 1+2*m*ishock+j) = mult*jmult*max_sigma/m;
        }
    }

  TwoDMatrix errors(model.numeq(), 2*m*model.nexog()+1);
  check(max_evals, y_mat, exo_mat, errors);

  // Report errors along shock and save them
  TwoDMatrix res(model.nexog(), 2*m+1);
  JournalRecord rec(journal);
  rec << "Shock    value         error" << endrec;
  ConstVector err0{errors.getCol(0)};
  for (int ishock = 0; ishock < model.nexog(); ishock++)
    {
      TwoDMatrix err_out(model.numeq(), 2*m+1);
      for (int j = 0; j < 2*m+1; j++)
        {
          int jj;
          Vector error{err_out.getCol(j)};
          if (j != m)
            {
              if (j < m)
                jj = 1 + 2*m*ishock+j;
              else
                jj = 1 + 2*m*ishock+j-1;
              ConstVector coljj{errors.getCol(jj)};
              error = coljj;
            }
          else
            {
              jj = 0;
              error = err0;
            }
          JournalRecord rec1(journal);
          std::string shockname{model.getExogNames().getName(ishock)};
          shockname.resize(8, ' ');
          rec1 << shockname << ' ' << exo_mat.get(ishock, jj)
               << "\t" << error.getMax() << endrec;
        }
      err_out.writeMat(fd, prefix + "_shock_" + model.getExogNames().getName(ishock) + "_errors");
    }
}

/* This method checks errors on ellipse of endogenous states (predetermined
   variables). The ellipse is shaped according to covariance matrix of
   endogenous variables based on the first order approximation and scaled by
   ‘mult’. The points on the ellipse are chosen as polar images of the low
   discrepancy grid in a cube.

   The method works as follows. First we calculate symmetric Schur factor of
   covariance matrix of the states. Second we generate low discrepancy points
   on the unit sphere. Third we transform the sphere with the
   variance-covariance matrix factor and multiplier ‘mult’ and initialize
   matrix of uₜ to zeros. Fourth we run the check() method and save the
   results. */

void
GlobalChecker::checkOnEllipseAndSave(mat_t *fd, const std::string &prefix,
                                     int m, double mult, int max_evals)
{
  JournalRecordPair pa(journal);
  pa << "Calculating errors at " << m
     << " ellipse points scaled by " << mult << endrec;

  // Make factor of covariance of variables
  /* Here we set ‘ycovfac’ to the symmetric Schur decomposition factor of a
     submatrix of covariances of all endogenous variables. The submatrix
     corresponds to state variables (predetermined plus both). */
  TwoDMatrix ycov{approx.calcYCov()};
  TwoDMatrix ycovpred(const_cast<const TwoDMatrix &>(ycov), model.nstat(), model.nstat(),
                      model.npred()+model.nboth(), model.npred()+model.nboth());
  SymSchurDecomp ssd(ycovpred);
  ssd.correctDefinitness(1.e-05);
  TwoDMatrix ycovfac(ycovpred.nrows(), ycovpred.ncols());
  KORD_RAISE_IF(!ssd.isPositiveSemidefinite(),
                "Covariance matrix of the states not positive \
				  semidefinite in GlobalChecker::checkOnEllipseAndSave");
  ssd.getFactor(ycovfac);

  // Put low discrepancy sphere points to ‘ymat’
  /* Here we first calculate dimension ‘d’ of the sphere, which is a number of
     state variables minus one. We go through the ‘d’-dimensional cube [0,1]ᵈ
     by QMCarloCubeQuadrature and make a polar transformation to the sphere.
     The polar transformation fⁱ can be written recursively w.r.t. the
     dimension i as:

      f⁰() = [1]

                    ⎡cos(2πxᵢ)·fⁱ⁻¹(x₁,…,xᵢ₋₁)⎤
      fⁱ(x₁,…,xᵢ) = ⎣        sin(2πxᵢ)        ⎦
  */
  int d = model.npred()+model.nboth()-1;
  TwoDMatrix ymat(model.npred()+model.nboth(), (d == 0) ? 2 : m);
  if (d == 0)
    {
      ymat.get(0, 0) = 1;
      ymat.get(0, 1) = -1;
    }
  else
    {
      int icol = 0;
      ReversePerScheme ps;
      QMCarloCubeQuadrature qmc(d, m, ps);
      qmcpit beg = qmc.start(m);
      qmcpit end = qmc.end(m);
      for (qmcpit run = beg; run != end; ++run, icol++)
        {
          Vector ycol{ymat.getCol(icol)};
          Vector x(run.point());
          x.mult(2*M_PI);
          ycol[0] = 1;
          for (int i = 0; i < d; i++)
            {
              Vector subsphere(ycol, 0, i+1);
              subsphere.mult(cos(x[i]));
              ycol[i+1] = sin(x[i]);
            }
        }
    }

  // Transform sphere ‘ymat’ and prepare ‘umat’ for checking
  /* Here we multiply the sphere points in ‘ymat’ with the Cholesky factor to
     obtain the ellipse, scale the ellipse by the given ‘mult’, and initialize
     matrix of shocks ‘umat’ to zero. */
  TwoDMatrix umat(model.nexog(), ymat.ncols());
  umat.zeros();
  ymat.mult(mult);
  ymat.multLeft(ycovfac);
  ConstVector ys(model.getSteady(), model.nstat(),
                 model.npred()+model.nboth());
  for (int icol = 0; icol < ymat.ncols(); icol++)
    {
      Vector ycol{ymat.getCol(icol)};
      ycol.add(1.0, ys);
    }

  // Check on ellipse and save
  /* Here we check the points and save the results to MAT-4 file. */
  TwoDMatrix out(model.numeq(), ymat.ncols());
  check(max_evals, ymat, umat, out);

  ymat.writeMat(fd, prefix + "_ellipse_points");
  out.writeMat(fd, prefix + "_ellipse_errors");
}

/* Here we check the errors along a simulation. We simulate, then set ‘x’ to
   zeros, check and save results. */

void
GlobalChecker::checkAlongSimulationAndSave(mat_t *fd, const std::string &prefix,
                                           int m, int max_evals)
{
  JournalRecordPair pa(journal);
  pa << "Calculating errors at " << m
     << " simulated points" << endrec;
  RandomShockRealization sr(model.getVcov(), seed_generator::get_new_seed());
  TwoDMatrix y{approx.getFoldDecisionRule().simulate(DecisionRule::emethod::horner,
                                                     m, model.getSteady(), sr)};
  TwoDMatrix x(model.nexog(), m);
  x.zeros();
  TwoDMatrix out(model.numeq(), m);
  check(max_evals, y, x, out);

  y.writeMat(fd, prefix + "_simul_points");
  out.writeMat(fd, prefix + "_simul_errors");
}
