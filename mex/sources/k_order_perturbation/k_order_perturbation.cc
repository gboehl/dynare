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

/* Convert MATLAB Dynare endo and exo names array to a vector<string> array of
   string pointers. MATLAB “mx” function returns a long string concatenated by
   columns rather than rows hence a rather low level approach is needed. */
void
DynareMxArrayToString(const mxArray *mxFldp, int len, int width, std::vector<std::string> &out)
{
  char *cNamesCharStr = mxArrayToString(mxFldp);

  out.resize(len);

  for (int i = 0; i < width; i++)
    for (int j = 0; j < len; j++)
      // Allow alphanumeric and underscores "_" only:
      if (std::isalnum(cNamesCharStr[j+i*len]) || (cNamesCharStr[j+i*len] == '_'))
        out[j] += cNamesCharStr[j+i*len];
}

void
copy_derivatives(mxArray *destin, const Symmetry &sym, const FGSContainer &derivs, const std::string &fieldname)
{
  const FGSTensor &x = derivs.get(sym);
  auto x_unfolded = x.unfold();
  int n = x_unfolded->numRows();
  int m = x_unfolded->numCols();
  mxArray *tmp = mxCreateDoubleMatrix(n, m, mxREAL);
  std::copy_n(x_unfolded->getData().base(), n*m, mxGetPr(tmp));
  mxSetField(destin, 0, fieldname.c_str(), tmp);
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
    if (mxIsNumeric(mxFldp))
      kOrder = static_cast<int>(mxGetScalar(mxFldp));
    else
      kOrder = 1;

    double qz_criterium = 1+1e-6;
    mxFldp = mxGetField(options_, 0, "qz_criterium");
    if (mxGetNumberOfElements(mxFldp) > 0 && mxIsNumeric(mxFldp))
      qz_criterium = mxGetScalar(mxFldp);

    mxFldp = mxGetField(M_, 0, "params");
    Vector modParams{mxFldp};
    if (!modParams.isFinite())
      DYN_MEX_FUNC_ERR_MSG_TXT("The parameters vector contains NaN or Inf");

    mxFldp = mxGetField(M_, 0, "Sigma_e");
    int npar = static_cast<int>(mxGetN(mxFldp));
    TwoDMatrix vCov(npar, npar, Vector{mxFldp});
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
    auto dparams = mxGetPr(mxFldp);
    npar = static_cast<int>(mxGetM(mxFldp));
    if (npar != nEndo)
      DYN_MEX_FUNC_ERR_MSG_TXT("Incorrect number of input var_order vars.");

    std::vector<int> var_order_vp(nEndo);
    for (int v = 0; v < nEndo; v++)
      var_order_vp[v] = static_cast<int>(*(dparams++));

    // the lag, current and lead blocks of the jacobian respectively
    mxFldp = mxGetField(M_, 0, "lead_lag_incidence");
    npar = static_cast<int>(mxGetN(mxFldp));
    int nrows = static_cast<int>(mxGetM(mxFldp));

    TwoDMatrix llincidence(nrows, npar, Vector{mxFldp});
    if (npar != nEndo)
      DYN_MEX_FUNC_ERR_MSG_TXT(("dynare:k_order_perturbation: Incorrect length of lead lag incidences: ncol="
                                + std::to_string(npar) + " != nEndo=" + std::to_string(nEndo)).c_str());

    mxFldp = mxGetField(M_, 0, "NNZDerivatives");
    Vector NNZD{mxFldp};
    if (NNZD[kOrder-1] == -1)
      DYN_MEX_FUNC_ERR_MSG_TXT("The derivatives were not computed for the required order. Make sure that you used the right order option inside the 'stoch_simul' command");

    mxFldp = mxGetField(M_, 0, "var_order_endo_names");
    const int nendo = static_cast<int>(mxGetM(mxFldp));
    const int widthEndo = static_cast<int>(mxGetN(mxFldp));
    std::vector<std::string> endoNames;
    DynareMxArrayToString(mxFldp, nendo, widthEndo, endoNames);

    mxFldp = mxGetField(M_, 0, "exo_names");
    const int nexo = static_cast<int>(mxGetM(mxFldp));
    const int widthExog = static_cast<int>(mxGetN(mxFldp));
    std::vector<std::string> exoNames;
    DynareMxArrayToString(mxFldp, nexo, widthExog, exoNames);

    if (nEndo != nendo || nExog != nexo)
      DYN_MEX_FUNC_ERR_MSG_TXT("Incorrect number of input parameters.");

    std::unique_ptr<TwoDMatrix> g1m, g2m, g3m;
    if (nrhs > 3)
      {
        // Derivatives have been passed as arguments
        const mxArray *g1 = prhs[3];
        int m = static_cast<int>(mxGetM(g1));
        int n = static_cast<int>(mxGetN(g1));
        g1m = std::make_unique<TwoDMatrix>(m, n, Vector{ConstVector{g1}});
        if (nrhs > 4)
          {
            const mxArray *g2 = prhs[4];
            int m = static_cast<int>(mxGetM(g2));
            int n = static_cast<int>(mxGetN(g2));
            g2m = std::make_unique<TwoDMatrix>(m, n, Vector{ConstVector{g2}});
            if (nrhs > 5)
              {
                const mxArray *g3 = prhs[5];
                int m = static_cast<int>(mxGetM(g3));
                int n = static_cast<int>(mxGetN(g3));
                g3m = std::make_unique<TwoDMatrix>(m, n, Vector{ConstVector{g3}});
              }
          }
      }

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
                           var_order_vp, llincidence,
                           std::move(g1m), std::move(g2m), std::move(g3m));

        // construct main K-order approximation class

        Approximation app(dynare, journal, nSteps, false, qz_criterium);
        // run stochastic steady
        app.walkStochSteady();

        /* Write derivative outputs into memory map */
        std::map<std::string, ConstTwoDMatrix> mm;
        app.getFoldDecisionRule().writeMMap(mm, "");

        // get latest ysteady
        ySteady = dynare.getSteady();

        if (kOrder == 1)
          {
            /* Set the output pointer to the output matrix ysteady. */
            auto cit = mm.begin();
            ++cit;
            plhs[1] = mxCreateDoubleMatrix(cit->second.numRows(), cit->second.numCols(), mxREAL);

            // Copy Dynare++ matrix into MATLAB matrix
            const ConstVector &vec = cit->second.getData();
            assert(vec.skip() == 1);
            std::copy_n(vec.base(), vec.length(), mxGetPr(plhs[1]));
          }
        if (kOrder >= 2)
          {
            int ii = 1;
            for (auto cit = mm.begin(); cit != mm.end() && ii < nlhs; ++cit)
              {
                plhs[ii] = mxCreateDoubleMatrix(cit->second.numRows(), cit->second.numCols(), mxREAL);

                // Copy Dynare++ matrix into MATLAB matrix
                const ConstVector &vec = cit->second.getData();
                assert(vec.skip() == 1);
                std::copy_n(vec.base(), vec.length(), mxGetPr(plhs[ii]));

                ++ii;

              }
            if (kOrder == 3 && nlhs > 5)
              {
                /* Return as 5th argument a struct containing *unfolded* matrices
                   for 3rd order decision rule */
                const FGSContainer &derivs = app.get_rule_ders();
                const std::string fieldnames[] = {"gy", "gu", "gyy", "gyu", "guu", "gss",
                                                  "gyyy", "gyyu", "gyuu", "guuu", "gyss", "guss"};
                // creates the char** expected by mxCreateStructMatrix()
                const char *c_fieldnames[12];
                for (int i = 0; i < 12; ++i)
                  c_fieldnames[i] = fieldnames[i].c_str();
                plhs[ii] = mxCreateStructMatrix(1, 1, 12, c_fieldnames);
                copy_derivatives(plhs[ii], Symmetry{1, 0, 0, 0}, derivs, "gy");
                copy_derivatives(plhs[ii], Symmetry{0, 1, 0, 0}, derivs, "gu");
                copy_derivatives(plhs[ii], Symmetry{2, 0, 0, 0}, derivs, "gyy");
                copy_derivatives(plhs[ii], Symmetry{0, 2, 0, 0}, derivs, "guu");
                copy_derivatives(plhs[ii], Symmetry{1, 1, 0, 0}, derivs, "gyu");
                copy_derivatives(plhs[ii], Symmetry{0, 0, 0, 2}, derivs, "gss");
                copy_derivatives(plhs[ii], Symmetry{3, 0, 0, 0}, derivs, "gyyy");
                copy_derivatives(plhs[ii], Symmetry{0, 3, 0, 0}, derivs, "guuu");
                copy_derivatives(plhs[ii], Symmetry{2, 1, 0, 0}, derivs, "gyyu");
                copy_derivatives(plhs[ii], Symmetry{1, 2, 0, 0}, derivs, "gyuu");
                copy_derivatives(plhs[ii], Symmetry{1, 0, 0, 2}, derivs, "gyss");
                copy_derivatives(plhs[ii], Symmetry{0, 1, 0, 2}, derivs, "guss");
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

