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

/* Defines the entry point for the k-order perturbation application DLL.

   See matlab/mex/k_order_perturbation.m for a description of inputs and
   outputs.
*/

#include "dynamic_m.hh"
#include "dynamic_dll.hh"
#include "k_ord_dynare.hh"

#include "approximation.hh"
#include "exception.hh"
#include "dynare_exception.hh"
#include "kord_exception.hh"
#include "tl_exception.hh"
#include "SylvException.hh"

#include <algorithm>
#include <cassert>

#include "dynmex.h"

/* Vector for storing field names like “g_0”, “g_1”, …
   A static structure is needed since MATLAB apparently does not create its own
   copy of the strings (contrary to what is said at:
    https://fr.mathworks.com/matlabcentral/answers/315937-mxcreatestructarray-and-mxcreatestructmatrix-field-name-memory-management
   ) */
std::vector<std::string> g_fieldnames;

/* Convert MATLAB Dynare endo and exo names cell array to a vector<string> array of
   string pointers. */
std::vector<std::string>
DynareMxArrayToString(const mxArray *mxFldp)
{
  assert(mxIsCell(mxFldp));
  std::vector<std::string> r;
  for (size_t i = 0; i < mxGetNumberOfElements(mxFldp); i++)
    r.emplace_back(mxArrayToString(mxGetCell(mxFldp, i)));

  return r;
}

void
copy_derivatives(mxArray *destin, const Symmetry &sym, const FGSContainer &derivs, const char *fieldname)
{
  const FGSTensor &x = derivs.get(sym);
  auto x_unfolded = x.unfold();
  int n = x_unfolded->numRows();
  int m = x_unfolded->numCols();
  mxArray *tmp = mxCreateDoubleMatrix(n, m, mxREAL);
  std::copy_n(x_unfolded->getData().base(), n*m, mxGetPr(tmp));
  mxSetField(destin, 0, fieldname, tmp);
}

extern "C" {

  void
  mexFunction(int nlhs, mxArray *plhs[],
              int nrhs, const mxArray *prhs[])
  {
    if (nrhs < 3 || nlhs < 2)
      DYN_MEX_FUNC_ERR_MSG_TXT("Must have at least 3 input parameters and takes at least 2 output parameters.");

    const mxArray *dr = prhs[0];
    const mxArray *M_ = prhs[1];
    const mxArray *options_ = prhs[2];
    bool use_dll = mxGetScalar(mxGetField(options_, 0, "use_dll")) != 0;

    mxArray *mFname = mxGetField(M_, 0, "fname");
    if (!mxIsChar(mFname))
      DYN_MEX_FUNC_ERR_MSG_TXT("Input must be of type char.");

    std::string fName = mxArrayToString(mFname);

    int kOrder;
    mxArray *mxFldp = mxGetField(options_, 0, "order");
    if (!mxIsNumeric(mxFldp))
      DYN_MEX_FUNC_ERR_MSG_TXT("options_.order must be a numeric value");
    kOrder = static_cast<int>(mxGetScalar(mxFldp));
    if (kOrder < 1 || kOrder > 3)
      DYN_MEX_FUNC_ERR_MSG_TXT("options_.order must be between 1 and 3");

    double qz_criterium = 1+1e-6;
    mxFldp = mxGetField(options_, 0, "qz_criterium");
    if (mxGetNumberOfElements(mxFldp) > 0 && mxIsNumeric(mxFldp))
      qz_criterium = mxGetScalar(mxFldp);

    mxFldp = mxGetField(M_, 0, "params");
    Vector modParams{mxFldp};
    if (!modParams.isFinite())
      DYN_MEX_FUNC_ERR_MSG_TXT("The parameters vector contains NaN or Inf");

    mxFldp = mxGetField(M_, 0, "Sigma_e");
    int dim = static_cast<int>(mxGetN(mxFldp));
    TwoDMatrix vCov(dim, dim, Vector{mxFldp});
    if (!vCov.isFinite())
      DYN_MEX_FUNC_ERR_MSG_TXT("The covariance matrix of shocks contains NaN or Inf");

    mxFldp = mxGetField(dr, 0, "ys");  // and not in order of dr.order_var
    Vector ySteady{mxFldp};
    if (!ySteady.isFinite())
      DYN_MEX_FUNC_ERR_MSG_TXT("The steady state vector contains NaN or Inf");

    mxFldp = mxGetField(M_, 0, "nstatic");
    const int nStat = static_cast<int>(mxGetScalar(mxFldp));
    mxFldp = mxGetField(M_, 0, "npred");
    const int nPred = static_cast<int>(mxGetScalar(mxFldp));
    mxFldp = mxGetField(M_, 0, "nboth");
    const int nBoth = static_cast<int>(mxGetScalar(mxFldp));
    mxFldp = mxGetField(M_, 0, "nfwrd");
    const int nForw = static_cast<int>(mxGetScalar(mxFldp));

    mxFldp = mxGetField(M_, 0, "exo_nbr");
    const int nExog = static_cast<int>(mxGetScalar(mxFldp));
    mxFldp = mxGetField(M_, 0, "endo_nbr");
    const int nEndo = static_cast<int>(mxGetScalar(mxFldp));
    mxFldp = mxGetField(M_, 0, "param_nbr");
    const int nPar = static_cast<int>(mxGetScalar(mxFldp));

    mxFldp = mxGetField(dr, 0, "order_var");
    dim = static_cast<int>(mxGetM(mxFldp));
    if (dim != nEndo)
      DYN_MEX_FUNC_ERR_MSG_TXT("Incorrect size of dr.order_var");
    std::vector<int> dr_order(nEndo);
    std::transform(mxGetPr(mxFldp), mxGetPr(mxFldp)+dim, dr_order.begin(),
                   [](double x) { return static_cast<int>(x)-1; });

    // the lag, current and lead blocks of the jacobian respectively
    TwoDMatrix llincidence(mxGetField(M_, 0, "lead_lag_incidence"));
    if (llincidence.nrows() != 3 || llincidence.ncols() != nEndo)
      DYN_MEX_FUNC_ERR_MSG_TXT("Incorrect size of M_.lead_lag_incidence");

    mxFldp = mxGetField(M_, 0, "NNZDerivatives");
    Vector NNZD{mxFldp};
    if (NNZD[kOrder-1] == -1)
      DYN_MEX_FUNC_ERR_MSG_TXT("The derivatives were not computed for the required order. Make sure that you used the right order option inside the `stoch_simul' command");

    mxFldp = mxGetField(M_, 0, "endo_names");
    std::vector<std::string> endoNames = DynareMxArrayToString(mxFldp);

    mxFldp = mxGetField(M_, 0, "exo_names");
    std::vector<std::string> exoNames = DynareMxArrayToString(mxFldp);

    if (nEndo != static_cast<int>(endoNames.size()) || nExog != static_cast<int>(exoNames.size()))
      DYN_MEX_FUNC_ERR_MSG_TXT("Incorrect size of M_.endo_names or M_.exo_names");

    const int nSteps = 0; // Dynare++ solving steps, for time being default to 0 = deterministic steady state

    try
      {
        Journal journal(fName + ".jnl");

        std::unique_ptr<DynamicModelAC> dynamicModelFile;
        if (use_dll)
          dynamicModelFile = std::make_unique<DynamicModelDLL>(fName);
        else
          dynamicModelFile = std::make_unique<DynamicModelMFile>(fName);

        // intiate tensor library
        TLStatic::init(kOrder, nStat+2*nPred+3*nBoth+2*nForw+nExog);

        // make KordpDynare object
        KordpDynare dynare(endoNames, exoNames, nExog, nPar,
                           ySteady, vCov, modParams, nStat, nPred, nForw, nBoth,
                           NNZD, nSteps, kOrder, journal, std::move(dynamicModelFile),
                           dr_order, llincidence);

        // If model derivatives have been passed as arguments
        if (nrhs > 3)
          {
            dynare.push_back_md(prhs[3]);
            if (nrhs > 4)
              dynare.push_back_md(prhs[4]);
            if (nrhs > 5)
              dynare.push_back_md(prhs[5]);
          }

        // construct main K-order approximation class
        Approximation app(dynare, journal, nSteps, false, qz_criterium);
        // run stochastic steady
        app.walkStochSteady();

        const FoldDecisionRule &fdr = app.getFoldDecisionRule();

        // Add possibly missing field names
        for (int i = static_cast<int>(g_fieldnames.size()); i <= kOrder; i++)
          g_fieldnames.emplace_back("g_" + std::to_string(i));
        // Create structure for storing derivatives in Dynare++ format
        const char *g_fieldnames_c[kOrder+1];
        for (int i = 0; i <= kOrder; i++)
          g_fieldnames_c[i] = g_fieldnames[i].c_str();
        plhs[1] = mxCreateStructMatrix(1, 1, kOrder+1, g_fieldnames_c);

        // Fill that structure
        for (int i = 0; i <= kOrder; i++)
          {
            const FFSTensor &t = fdr.get(Symmetry{i});
            mxArray *tmp = mxCreateDoubleMatrix(t.numRows(), t.numCols(), mxREAL);
            const ConstVector &vec = t.getData();
            assert(vec.skip() == 1);
            std::copy_n(vec.base(), vec.length(), mxGetPr(tmp));
            mxSetField(plhs[1], 0, ("g_" + std::to_string(i)).c_str(), tmp);
          }

        if (nlhs > 2)
          {
            /* Return as 3rd argument a struct containing derivatives in Dynare
               format (unfolded matrices, without Taylor coefficient) up to 3rd
               order */
            const FGSContainer &derivs = app.get_rule_ders();

            size_t nfields = (kOrder == 1 ? 2 : (kOrder == 2 ? 6 : 12));
            const char *c_fieldnames[] = { "gy", "gu", "gyy", "gyu", "guu", "gss",
                                           "gyyy", "gyyu", "gyuu", "guuu", "gyss", "guss" };
            plhs[2] = mxCreateStructMatrix(1, 1, nfields, c_fieldnames);

            copy_derivatives(plhs[2], Symmetry{1, 0, 0, 0}, derivs, "gy");
            copy_derivatives(plhs[2], Symmetry{0, 1, 0, 0}, derivs, "gu");
            if (kOrder >= 2)
              {
                copy_derivatives(plhs[2], Symmetry{2, 0, 0, 0}, derivs, "gyy");
                copy_derivatives(plhs[2], Symmetry{0, 2, 0, 0}, derivs, "guu");
                copy_derivatives(plhs[2], Symmetry{1, 1, 0, 0}, derivs, "gyu");
                copy_derivatives(plhs[2], Symmetry{0, 0, 0, 2}, derivs, "gss");
              }
            if (kOrder >= 3)
              {
                copy_derivatives(plhs[2], Symmetry{3, 0, 0, 0}, derivs, "gyyy");
                copy_derivatives(plhs[2], Symmetry{0, 3, 0, 0}, derivs, "guuu");
                copy_derivatives(plhs[2], Symmetry{2, 1, 0, 0}, derivs, "gyyu");
                copy_derivatives(plhs[2], Symmetry{1, 2, 0, 0}, derivs, "gyuu");
                copy_derivatives(plhs[2], Symmetry{1, 0, 0, 2}, derivs, "gyss");
                copy_derivatives(plhs[2], Symmetry{0, 1, 0, 2}, derivs, "guss");
              }
          }
      }
    catch (const KordException &e)
      {
        e.print();
        DYN_MEX_FUNC_ERR_MSG_TXT(("dynare:k_order_perturbation: Caught Kord exception: " + e.get_message()).c_str());
      }
    catch (const TLException &e)
      {
        e.print();
        DYN_MEX_FUNC_ERR_MSG_TXT("dynare:k_order_perturbation: Caught TL exception");
      }
    catch (SylvException &e)
      {
        e.printMessage();
        DYN_MEX_FUNC_ERR_MSG_TXT("dynare:k_order_perturbation: Caught Sylv exception");
      }
    catch (const DynareException &e)
      {
        DYN_MEX_FUNC_ERR_MSG_TXT(("dynare:k_order_perturbation: Caught KordDynare exception: " + e.message()).c_str());
      }
    catch (const ogu::Exception &e)
      {
        DYN_MEX_FUNC_ERR_MSG_TXT(("dynare:k_order_perturbation: Caught general exception: " + e.message()).c_str());
      }
    plhs[0] = mxCreateDoubleScalar(0);
  } // end of mexFunction()
} // end of extern C

