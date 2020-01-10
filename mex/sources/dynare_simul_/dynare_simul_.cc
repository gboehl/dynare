/*
 * Copyright © 2005-2011 Ondra Kamenik
 * Copyright © 2019-2020 Dynare Team
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

// This is the mexFunction providing interface to
// DecisionRule<>::simulate(). It takes the following input
// parameters:
//      order    the order of approximation, needs order+1 derivatives
//      nstat
//      npred
//      nboth
//      nforw
//      nexog
//      ystart   starting value (full vector of endogenous)
//      shocks   matrix of shocks (nexog x number of period)
//      vcov     covariance matrix of shocks (nexog x nexog)
//      seed     integer seed
//      ysteady  full vector of decision rule's steady
//      dr       structure containing matrices of derivatives (g_0, g_1,…)

// output:
//      res      simulated results

#include "dynmex.h"
#include "mex.h"

#include "decision_rule.hh"
#include "fs_tensor.hh"
#include "SylvException.hh"

#include <string>

extern "C" {
  void
  mexFunction(int nlhs, mxArray *plhs[],
              int nhrs, const mxArray *prhs[])
  {
    if (nhrs != 12 || nlhs != 1)
      mexErrMsgTxt("dynare_simul_ must have at exactly 12 input parameters and 1 output argument.");

    int order = static_cast<int>(mxGetScalar(prhs[0]));
    int nstat = static_cast<int>(mxGetScalar(prhs[1]));
    int npred = static_cast<int>(mxGetScalar(prhs[2]));
    int nboth = static_cast<int>(mxGetScalar(prhs[3]));
    int nforw = static_cast<int>(mxGetScalar(prhs[4]));
    int nexog = static_cast<int>(mxGetScalar(prhs[5]));

    const mxArray *const ystart = prhs[6];
    const mxArray *const shocks = prhs[7];
    const mxArray *const vcov = prhs[8];
    int seed = static_cast<int>(mxGetScalar(prhs[9]));
    const mxArray *const ysteady = prhs[10];
    const mxArray *const dr = prhs[11];
    const mwSize *const ystart_dim = mxGetDimensions(ystart);
    const mwSize *const shocks_dim = mxGetDimensions(shocks);
    const mwSize *const vcov_dim = mxGetDimensions(vcov);
    const mwSize *const ysteady_dim = mxGetDimensions(ysteady);

    int ny = nstat + npred + nboth + nforw;
    if (ny != static_cast<int>(ystart_dim[0]))
      mexErrMsgTxt("ystart has wrong number of rows.\n");
    if (1 != ystart_dim[1])
      mexErrMsgTxt("ystart has wrong number of cols.\n");
    int nper = shocks_dim[1];
    if (nexog != static_cast<int>(shocks_dim[0]))
      mexErrMsgTxt("shocks has a wrong number of rows.\n");
    if (nexog != static_cast<int>(vcov_dim[0]))
      mexErrMsgTxt("vcov has a wrong number of rows.\n");
    if (nexog != static_cast<int>(vcov_dim[1]))
      mexErrMsgTxt("vcov has a wrong number of cols.\n");
    if (ny != static_cast<int>(ysteady_dim[0]))
      mexErrMsgTxt("ysteady has wrong number of rows.\n");
    if (1 != ysteady_dim[1])
      mexErrMsgTxt("ysteady has wrong number of cols.\n");

    plhs[0] = mxCreateDoubleMatrix(ny, nper, mxREAL);

    try
      {
        // initialize tensor library
        TLStatic::init(order, npred+nboth+nexog);

        // form the polynomial
        UTensorPolynomial pol(ny, npred+nboth+nexog);
        for (int dim = 0; dim <= order; dim++)
          {
            const mxArray *gk_m = mxGetField(dr, 0, ("g_" + std::to_string(dim)).c_str());
            if (!gk_m)
              mexErrMsgTxt(("Can't find field `g_" + std::to_string(dim) + "' in structured passed as last argument").c_str());
            ConstTwoDMatrix gk(gk_m);
            FFSTensor ft(ny, npred+nboth+nexog, dim);
            if (ft.ncols() != gk.ncols())
              mexErrMsgTxt(("Wrong number of columns for folded tensor: got " + std::to_string(gk.ncols()) + " but i want " + std::to_string(ft.ncols()) + '\n').c_str());
            if (ft.nrows() != gk.nrows())
              mexErrMsgTxt(("Wrong number of rows for folded tensor: got " + std::to_string(gk.nrows()) + " but i want " + std::to_string(ft.nrows()) + '\n').c_str());
            ft.zeros();
            ft.add(1.0, gk);
            pol.insert(std::make_unique<UFSTensor>(ft));
          }
        // form the decision rule
        UnfoldDecisionRule dr(pol, PartitionY(nstat, npred, nboth, nforw),
                              nexog, ConstVector{ysteady});
        // form the shock realization
        ConstTwoDMatrix shocks_mat(nexog, nper, ConstVector{shocks});
        ConstTwoDMatrix vcov_mat(nexog, nexog, ConstVector{vcov});
        GenShockRealization sr(vcov_mat, shocks_mat, seed);
        // simulate and copy the results
        TwoDMatrix res_mat{dr.simulate(DecisionRule::emethod::horner, nper,
                                       ConstVector{ystart}, sr)};
        TwoDMatrix res_tmp_mat{plhs[0]};
        res_tmp_mat = const_cast<const TwoDMatrix &>(res_mat);
      }
    catch (const KordException &e)
      {
        mexErrMsgTxt("Caught Kord exception.");
      }
    catch (const TLException &e)
      {
        mexErrMsgTxt("Caught TL exception.");
      }
    catch (SylvException &e)
      {
        mexErrMsgTxt("Caught Sylv exception.");
      }
  }
};
