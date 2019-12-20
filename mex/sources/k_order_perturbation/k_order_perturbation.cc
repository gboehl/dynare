/*
 * Copyright © 2008-2019 Dynare Team
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
    {
      const mxArray *cell_mx = mxGetCell(mxFldp, i);
      if (!(cell_mx && mxIsChar(cell_mx)))
        mexErrMsgTxt("Cell is not a character array");
      r.emplace_back(mxArrayToString(cell_mx));
    }

  return r;
}

void
copy_derivatives(mxArray *destin, const Symmetry &sym, const FGSContainer &derivs, const char *fieldname)
{
  const FGSTensor &x = derivs.get(sym);
  auto x_unfolded = x.unfold();
  int n = x_unfolded->nrows();
  int m = x_unfolded->ncols();
  mxArray *tmp = mxCreateDoubleMatrix(n, m, mxREAL);
  std::copy_n(x_unfolded->getData().base(), n*m, mxGetPr(tmp));
  mxSetField(destin, 0, fieldname, tmp);
}

extern "C" {

  void
  mexFunction(int nlhs, mxArray *plhs[],
              int nrhs, const mxArray *prhs[])
  {
    if (nrhs < 3 || nlhs < 2 || nlhs > 3)
      DYN_MEX_FUNC_ERR_MSG_TXT("Must have at least 3 input parameters and takes 2 or 3 output parameters.");

    // Give explicit names to input arguments
    const mxArray *dr_mx = prhs[0];
    const mxArray *M_mx = prhs[1];
    const mxArray *options_mx = prhs[2];

    auto get_int_field = [](const mxArray *struct_mx, const std::string &fieldname)
                         {
                           mxArray *field_mx = mxGetField(struct_mx, 0, fieldname.c_str());
                           if (!(field_mx && mxIsScalar(field_mx) && mxIsNumeric(field_mx)))
                             mexErrMsgTxt(("Field `" + fieldname + "' should be a numeric scalar").c_str());
                           return static_cast<int>(mxGetScalar(field_mx));
                         };

    // Extract various fields from options_
    const int kOrder = get_int_field(options_mx, "order");
    if (kOrder < 1)
      DYN_MEX_FUNC_ERR_MSG_TXT("options_.order must be at least 1");

    const mxArray *use_dll_mx = mxGetField(options_mx, 0, "use_dll");
    if (!(use_dll_mx && mxIsLogicalScalar(use_dll_mx)))
      DYN_MEX_FUNC_ERR_MSG_TXT("options_.use_dll should be a logical scalar");
    bool use_dll = static_cast<bool>(mxGetScalar(use_dll_mx));

    double qz_criterium = 1+1e-6;
    const mxArray *qz_criterium_mx = mxGetField(options_mx, 0, "qz_criterium");
    if (qz_criterium_mx && mxIsScalar(qz_criterium_mx) && mxIsNumeric(qz_criterium_mx))
      qz_criterium = mxGetScalar(qz_criterium_mx);

    const mxArray *threads_mx = mxGetField(options_mx, 0, "threads");
    if (!threads_mx)
      DYN_MEX_FUNC_ERR_MSG_TXT("Can't find field options_.threads");
    const mxArray *num_threads_mx = mxGetField(threads_mx, 0, "k_order_perturbation");
    if (!(num_threads_mx && mxIsScalar(num_threads_mx) && mxIsNumeric(num_threads_mx)))
      DYN_MEX_FUNC_ERR_MSG_TXT("options_.threads.k_order_perturbation be a numeric scalar");
    int num_threads = static_cast<int>(mxGetScalar(num_threads_mx));

    // Extract various fields from M_
    const mxArray *fname_mx = mxGetField(M_mx, 0, "fname");
    if (!(fname_mx && mxIsChar(fname_mx) && mxGetM(fname_mx) == 1))
      DYN_MEX_FUNC_ERR_MSG_TXT("M_.fname should be a character string");
    std::string fname{mxArrayToString(fname_mx)};

    const mxArray *params_mx = mxGetField(M_mx, 0, "params");
    if (!(params_mx && mxIsDouble(params_mx)))
      DYN_MEX_FUNC_ERR_MSG_TXT("M_.params should be a double precision array");
    Vector modParams{ConstVector{params_mx}};
    if (!modParams.isFinite())
      DYN_MEX_FUNC_ERR_MSG_TXT("M_.params contains NaN or Inf");

    const mxArray *sigma_e_mx = mxGetField(M_mx, 0, "Sigma_e");
    if (!(sigma_e_mx && mxIsDouble(sigma_e_mx) && mxGetM(sigma_e_mx) == mxGetN(sigma_e_mx)))
      DYN_MEX_FUNC_ERR_MSG_TXT("M_.Sigma_e should be a double precision square matrix");
    TwoDMatrix vCov{ConstTwoDMatrix{sigma_e_mx}};
    if (!vCov.isFinite())
      DYN_MEX_FUNC_ERR_MSG_TXT("M_.Sigma_e contains NaN or Inf");

    const int nStat = get_int_field(M_mx, "nstatic");
    const int nPred = get_int_field(M_mx, "npred");
    const int nBoth = get_int_field(M_mx, "nboth");
    const int nForw = get_int_field(M_mx, "nfwrd");

    const int nExog = get_int_field(M_mx, "exo_nbr");
    const int nEndo = get_int_field(M_mx, "endo_nbr");
    const int nPar = get_int_field(M_mx, "param_nbr");

    const mxArray *lead_lag_incidence_mx = mxGetField(M_mx, 0, "lead_lag_incidence");
    if (!(lead_lag_incidence_mx && mxIsDouble(lead_lag_incidence_mx) && mxGetM(lead_lag_incidence_mx) == 3
          && mxGetN(lead_lag_incidence_mx) == static_cast<size_t>(nEndo)))
      DYN_MEX_FUNC_ERR_MSG_TXT("M_.lead_lag_incidence should be a double precision matrix with 3 rows and M_.endo_nbr columns");
    ConstTwoDMatrix llincidence{lead_lag_incidence_mx};

    const mxArray *nnzderivatives_mx = mxGetField(M_mx, 0, "NNZDerivatives");
    if (!(nnzderivatives_mx && mxIsDouble(nnzderivatives_mx)))
      DYN_MEX_FUNC_ERR_MSG_TXT("M_.NNZDerivatives should be a double precision array");
    ConstVector NNZD{nnzderivatives_mx};
    if (NNZD.length() < kOrder || NNZD[kOrder-1] == -1)
      DYN_MEX_FUNC_ERR_MSG_TXT("The derivatives were not computed for the required order. Make sure that you used the right order option inside the `stoch_simul' command");

    const mxArray *endo_names_mx = mxGetField(M_mx, 0, "endo_names");
    if (!(endo_names_mx && mxIsCell(endo_names_mx) && mxGetNumberOfElements(endo_names_mx) == static_cast<size_t>(nEndo)))
      DYN_MEX_FUNC_ERR_MSG_TXT("M_.endo_names should be a cell array of M_.endo_nbr elements");
    std::vector<std::string> endoNames = DynareMxArrayToString(endo_names_mx);

    const mxArray *exo_names_mx = mxGetField(M_mx, 0, "exo_names");
    if (!(exo_names_mx && mxIsCell(exo_names_mx) && mxGetNumberOfElements(exo_names_mx) == static_cast<size_t>(nExog)))
      DYN_MEX_FUNC_ERR_MSG_TXT("M_.exo_names should be a cell array of M_.exo_nbr elements");
    std::vector<std::string> exoNames = DynareMxArrayToString(exo_names_mx);

    const mxArray *dynamic_tmp_nbr_mx = mxGetField(M_mx, 0, "dynamic_tmp_nbr");
    if (!(dynamic_tmp_nbr_mx && mxIsDouble(dynamic_tmp_nbr_mx) && mxGetNumberOfElements(dynamic_tmp_nbr_mx) >= static_cast<size_t>(kOrder+1)))
      DYN_MEX_FUNC_ERR_MSG_TXT("M_.dynamic_tmp_nbr should be a double precision array with strictly more elements than the order of derivation");
    int ntt = std::accumulate(mxGetPr(dynamic_tmp_nbr_mx), mxGetPr(dynamic_tmp_nbr_mx)+kOrder+1, 0);

    // Extract various fields from dr
    const mxArray *ys_mx = mxGetField(dr_mx, 0, "ys"); // and not in order of dr.order_var
    if (!(ys_mx && mxIsDouble(ys_mx)))
      DYN_MEX_FUNC_ERR_MSG_TXT("dr.ys should be a double precision array");
    Vector ySteady{ConstVector{ys_mx}};
    if (!ySteady.isFinite())
      DYN_MEX_FUNC_ERR_MSG_TXT("dr.ys contains NaN or Inf");

    const mxArray *order_var_mx = mxGetField(dr_mx, 0, "order_var");
    if (!(order_var_mx && mxIsDouble(order_var_mx) && mxGetNumberOfElements(order_var_mx) == static_cast<size_t>(nEndo)))
      DYN_MEX_FUNC_ERR_MSG_TXT("dr.order_var should be a double precision array of M_.endo_nbr elements");
    std::vector<int> dr_order(nEndo);
    std::transform(mxGetPr(order_var_mx), mxGetPr(order_var_mx)+nEndo, dr_order.begin(),
                   [](double x) { return static_cast<int>(x)-1; });

    const int nSteps = 0; // Dynare++ solving steps, for time being default to 0 = deterministic steady state

    try
      {
        Journal journal(fname + ".jnl");

        std::unique_ptr<DynamicModelAC> dynamicModelFile;
        if (use_dll)
          dynamicModelFile = std::make_unique<DynamicModelDLL>(fname, ntt, kOrder);
        else
          dynamicModelFile = std::make_unique<DynamicModelMFile>(fname, ntt);

        // intiate tensor library
        TLStatic::init(kOrder, nStat+2*nPred+3*nBoth+2*nForw+nExog);

        // Set number of parallel threads
        sthread::detach_thread_group::max_parallel_threads = num_threads;

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
            mxArray *tmp = mxCreateDoubleMatrix(t.nrows(), t.ncols(), mxREAL);
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
