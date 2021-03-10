/*
 * Copyright © 2019-2021 Dynare Team
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

#include <memory>
#include <string>
#include <utility>

#include <dynmex.h>

#include "tl_static.hh"
#include "decision_rule.hh"
#include "sthread.hh"

/* The class that does the real job. It computes the next iteration for a given
   range of particles. There will be as many instances as there are parallel
   threads.

   Note that we can’t use OpenMP since it is incompatible with MKL, which is
   used internally by Dynare++ routines under MATLAB. Hence we fall back to
   the Dynare++ multithreading abstraction. */
struct ParticleWorker : public sthread::detach_thread
{
  const int npred_both, exo_nbr;
  const std::pair<size_t, size_t> particle_range;
  const ConstGeneralMatrix &yhat, &epsilon;
  const Vector &ys_reordered;
  const UnfoldDecisionRule &dr;
  const ConstVector &restrict_var_list;
  GeneralMatrix &ynext;

  ParticleWorker(int npred_both_arg, int exo_nbr_arg, std::pair<size_t, size_t> particle_range_arg,
                 const ConstGeneralMatrix &yhat_arg, const ConstGeneralMatrix &epsilon_arg,
                 const Vector &ys_reordered_arg, const UnfoldDecisionRule &dr_arg,
                 const ConstVector &restrict_var_list_arg, GeneralMatrix &ynext_arg)
    : npred_both{npred_both_arg}, exo_nbr{exo_nbr_arg}, particle_range{std::move(particle_range_arg)},
      yhat{yhat_arg}, epsilon{epsilon_arg}, ys_reordered{ys_reordered_arg}, dr{dr_arg},
      restrict_var_list{restrict_var_list_arg}, ynext{ynext_arg}
  {
  }
  void
  operator()(std::mutex &mut) override
  {
    Vector dyu(npred_both+exo_nbr);
    Vector dy(dyu, 0, npred_both);
    Vector u(dyu, npred_both, exo_nbr);
    Vector ynext_col_allvars(ys_reordered.length());

    for (size_t i = particle_range.first; i < particle_range.second; i++)
      {
        dy = yhat.getCol(i);
        u = epsilon.getCol(i);

        dr.eval(DecisionRule::emethod::horner, ynext_col_allvars, dyu);

        ynext_col_allvars.add(1.0, ys_reordered);

        /* Select only the variables in restrict_var_list, and copy back to the
           result matrix */
        Vector ynext_col{ynext.getCol(i)};
        for (int j = 0; j < restrict_var_list.length(); j++)
          ynext_col[j] = ynext_col_allvars[restrict_var_list[j]-1];
      }
  }
};

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs != 5 || nlhs != 1)
    mexErrMsgTxt("Must have 5 input arguments and 1 output argument");

  // Give explicit names to input arguments
  const mxArray *yhat_mx = prhs[0];
  const mxArray *epsilon_mx = prhs[1];
  const mxArray *dr_mx = prhs[2];
  const mxArray *M_mx = prhs[3];
  const mxArray *options_mx = prhs[4];

  auto get_int_field = [](const mxArray *struct_mx, const std::string &fieldname)
                       {
                         mxArray *field_mx = mxGetField(struct_mx, 0, fieldname.c_str());
                         if (!(field_mx && mxIsScalar(field_mx) && mxIsNumeric(field_mx)))
                           mexErrMsgTxt(("Field `" + fieldname + "' should be a numeric scalar").c_str());
                         return static_cast<int>(mxGetScalar(field_mx));
                       };

  int nstatic = get_int_field(M_mx, "nstatic");
  int npred = get_int_field(M_mx, "npred");
  int nboth = get_int_field(M_mx, "nboth");
  int nfwrd = get_int_field(M_mx, "nfwrd");
  int endo_nbr = nstatic + npred + nboth + nfwrd;
  int exo_nbr = get_int_field(M_mx, "exo_nbr");
  int order = get_int_field(options_mx, "order");

  const mxArray *order_var_mx = mxGetField(dr_mx, 0, "order_var");
  if (!(order_var_mx && mxIsDouble(order_var_mx) && mxGetNumberOfElements(order_var_mx) == static_cast<size_t>(endo_nbr)))
    mexErrMsgTxt("Field dr.order_var should be a double precision vector with endo_nbr elements");
  const mxArray *ys_mx = mxGetField(dr_mx, 0, "ys");
  if (!ys_mx || !mxIsDouble(ys_mx) || mxGetNumberOfElements(ys_mx) != static_cast<size_t>(endo_nbr))
    mexErrMsgTxt("Field dr.ys should be a double precision vector with endo_nbr elements");
  const mxArray *restrict_var_list_mx = mxGetField(dr_mx, 0, "restrict_var_list");
  if (!(restrict_var_list_mx && mxIsDouble(restrict_var_list_mx)))
    mexErrMsgTxt("Field dr.restrict_var_list should be a double precision vector");

  size_t nparticles = mxGetN(yhat_mx);
  if (mxGetN(epsilon_mx) != nparticles)
    mexErrMsgTxt("epsilon and yhat don't have the same number of columns");
  if (!(mxIsDouble(yhat_mx) && mxGetM(yhat_mx) == static_cast<size_t>(npred + nboth)))
    mexErrMsgTxt("yhat should be a double precision matrix with npred+nboth rows");
  if (!(mxIsDouble(epsilon_mx) && mxGetM(epsilon_mx) == static_cast<size_t>(exo_nbr)))
    mexErrMsgTxt("epsilon should be a double precision matrix with exo_nbr rows");

  const mxArray *threads_mx = mxGetField(options_mx, 0, "threads");
  if (!threads_mx)
    mexErrMsgTxt("Can't find field `threads' in options_");
  int num_threads = get_int_field(threads_mx, "local_state_space_iteration_k");

  ConstGeneralMatrix yhat{yhat_mx};
  ConstGeneralMatrix epsilon{epsilon_mx};
  ConstVector ys{ys_mx};
  const double *order_var = mxGetPr(order_var_mx);
  ConstVector restrict_var_list{restrict_var_list_mx};

  try
    {
      TLStatic::init(order, npred+nboth+exo_nbr);
    }
  catch (TLException &e)
    {
      mexErrMsgTxt(("Dynare++ error: " + e.message).c_str());
    }

  // Form the polynomial (copied from dynare_simul_.cc)
  UTensorPolynomial pol(endo_nbr, npred+nboth+exo_nbr);
  for (int dim = 0; dim <= order; dim++)
    {
      const mxArray *gk_m = mxGetField(dr_mx, 0, ("g_" + std::to_string(dim)).c_str());
      if (!gk_m)
        mexErrMsgTxt(("Can't find field `g_" + std::to_string(dim) + "' in dr structure").c_str());
      ConstTwoDMatrix gk{gk_m};
      FFSTensor ft{endo_nbr, npred+nboth+exo_nbr, dim};
      if (ft.ncols() != gk.ncols())
        mexErrMsgTxt(("Wrong number of columns for folded tensor: got " + std::to_string(gk.ncols()) + " but i want " + std::to_string(ft.ncols()) + '\n').c_str());
      if (ft.nrows() != gk.nrows())
        mexErrMsgTxt(("Wrong number of rows for folded tensor: got " + std::to_string(gk.nrows()) + " but i want " + std::to_string(ft.nrows()) + '\n').c_str());
      ft.zeros();
      ft.add(1.0, gk);
      pol.insert(std::make_unique<UFSTensor>(ft));
    }

  // Construct the reordered steady state (dr.ys(dr.order_var))
  Vector ys_reordered(endo_nbr);
  for (int i = 0; i < endo_nbr; i++)
    ys_reordered[i] = ys[static_cast<int>(order_var[i])-1];

  // Form the decision rule
  UnfoldDecisionRule dr(pol, PartitionY(nstatic, npred, nboth, nfwrd),
                        exo_nbr, ys_reordered);

  // Create the result matrix
  plhs[0] = mxCreateDoubleMatrix(restrict_var_list.length(), nparticles, mxREAL);
  GeneralMatrix ynext{plhs[0]};

  // Run the real job in parallel
  sthread::detach_thread_group::max_parallel_threads = num_threads;
  sthread::detach_thread_group group;
  // The following is equivalent to ceil(nparticles/num_threads), but with integer arithmetic
  int part_by_thread = nparticles / num_threads + (nparticles % num_threads > 0);
  for (size_t i = 0; i < nparticles; i += part_by_thread)
    group.insert(std::make_unique<ParticleWorker>(npred+nboth, exo_nbr,
                                                  std::make_pair(i, std::min(i+part_by_thread, nparticles)),
                                                  yhat, epsilon, ys_reordered, dr, restrict_var_list,
                                                  ynext));
  group.run();
}
