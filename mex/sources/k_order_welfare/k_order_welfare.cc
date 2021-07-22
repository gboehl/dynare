/*
 * Copyright © 2021 Dynare Team
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

#include "dynamic_m.hh"
#include "dynamic_dll.hh"
#include "objective_m.hh"

#include "approximation.hh"
#include "approximation_welfare.hh"
#include "exception.hh"
#include "dynare_exception.hh"
#include "kord_exception.hh"
#include "tl_exception.hh"
#include "SylvException.hh"

#include <algorithm>
#include <cassert>

#include "dynmex.h"

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

/* Vector for storing field names like “W_0”, “W_1”, …
   A static structure is needed since MATLAB apparently does not create its own
   copy of the strings (contrary to what is said at:
    https://fr.mathworks.com/matlabcentral/answers/315937-mxcreatestructarray-and-mxcreatestructmatrix-field-name-memory-management
   )
*/

std::vector<std::string> U_fieldnames;
std::vector<std::string> W_fieldnames;

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

/* The k_order_welfare function computes the conditional welfare starting from the non-stochastic steady state and the derivatives of the felicity and welfare functions U(h(yₜ₋₁, uₜ, σ)) and W(yₜ₋₁, uₜ, σ). The unconditional welfare shall be introduced in a future version with the update of the dynare_simul_.

The routine proceeds in several steps:
1. Computation of the kth-order policy functions: the code is a copy-paste from the k_order_perturbation.cc file
2. Computation of the kth-order derivatives of the felicity and of the welfare functions. The call to the approxAtSteady method in the ApproximationWelfare class carries out the necessary operation 
   - Importing the derivatives of the felicity function with the calcDerivativesAtSteady() method of the KordwDynare class. It relies on the Matlab-generated files, which are handled by the ObjectiveAC and ObjectiveMFile classes
   - Pinpointing the derivatives of the felicity and welfare functions. The performStep method of the KOrderWelfare class carries out the calculations,resorting to the FaaDiBruno class and its methods to get the needed intermediary results.
*/

extern "C" {

  void mexFunction(int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
  {

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
    const mxArray *ramsey_policy_mx = mxGetField(options_mx, 0, "ramsey_policy");
    if (!(ramsey_policy_mx && mxIsLogicalScalar(ramsey_policy_mx)))
      mexErrMsgTxt("options_.ramsey_policy should be a logical scalar");
    bool ramsey_policy = static_cast<bool>(mxGetScalar(ramsey_policy_mx));

    if (ramsey_policy == false)
      mexErrMsgTxt("The considered model must be a Ramsey-typed model !");


    const int kOrder = get_int_field(options_mx, "order");
    if (kOrder < 1)
      mexErrMsgTxt("options_.order must be at least 1");

    const mxArray *use_dll_mx = mxGetField(options_mx, 0, "use_dll");
    if (!(use_dll_mx && mxIsLogicalScalar(use_dll_mx)))
      mexErrMsgTxt("options_.use_dll should be a logical scalar");
    bool use_dll = static_cast<bool>(mxGetScalar(use_dll_mx));

    double qz_criterium = 1+1e-6;
    const mxArray *qz_criterium_mx = mxGetField(options_mx, 0, "qz_criterium");
    if (qz_criterium_mx && mxIsScalar(qz_criterium_mx) && mxIsNumeric(qz_criterium_mx))
      qz_criterium = mxGetScalar(qz_criterium_mx);

    const mxArray *threads_mx = mxGetField(options_mx, 0, "threads");
    if (!threads_mx)
      mexErrMsgTxt("Can't find field options_.threads");
    const mxArray *num_threads_mx = mxGetField(threads_mx, 0, "k_order_perturbation");
    if (!(num_threads_mx && mxIsScalar(num_threads_mx) && mxIsNumeric(num_threads_mx)))
      mexErrMsgTxt("options_.threads.k_order_perturbation be a numeric scalar");
    int num_threads = static_cast<int>(mxGetScalar(num_threads_mx));

    const mxArray *debug_mx = mxGetField(options_mx, 0, "debug");
    if (!(debug_mx && mxIsLogicalScalar(debug_mx)))
      mexErrMsgTxt("options_.debug should be a logical scalar");
    bool debug = static_cast<bool>(mxGetScalar(debug_mx));

    // Extract various fields from M_
    const mxArray *fname_mx = mxGetField(M_mx, 0, "fname");
    if (!(fname_mx && mxIsChar(fname_mx) && mxGetM(fname_mx) == 1))
      mexErrMsgTxt("M_.fname should be a character string");
    std::string fname{mxArrayToString(fname_mx)};

    const mxArray *params_mx = mxGetField(M_mx, 0, "params");
    if (!(params_mx && mxIsDouble(params_mx)))
      mexErrMsgTxt("M_.params should be a double precision array");
    Vector modParams{ConstVector{params_mx}};
    if (!modParams.isFinite())
      mexErrMsgTxt("M_.params contains NaN or Inf");

    const mxArray *param_names_mx = mxGetField(M_mx, 0, "param_names");
    if (!(param_names_mx && mxIsCell(param_names_mx)))
      mexErrMsgTxt("M_.param_names should be a cell array");
    std::vector<std::string> paramNames = DynareMxArrayToString(param_names_mx);
    auto it = std::find(paramNames.begin(), paramNames.end(), "optimal_policy_discount_factor");
    double discount_factor;
    if (it != paramNames.end())
      discount_factor = modParams[std::distance(paramNames.begin(), it)];
    else
      {
        mexErrMsgTxt("M_.param_names does not contain any \"optimal_policy_discount_factor\".");
        return; // To silence a GCC warning about discount_factor unitialized
      }

    const mxArray *sigma_e_mx = mxGetField(M_mx, 0, "Sigma_e");
    if (!(sigma_e_mx && mxIsDouble(sigma_e_mx) && mxGetM(sigma_e_mx) == mxGetN(sigma_e_mx)))
      mexErrMsgTxt("M_.Sigma_e should be a double precision square matrix");
    TwoDMatrix vCov{ConstTwoDMatrix{sigma_e_mx}};
    if (!vCov.isFinite())
      mexErrMsgTxt("M_.Sigma_e contains NaN or Inf");

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
      mexErrMsgTxt("M_.lead_lag_incidence should be a double precision matrix with 3 rows and M_.endo_nbr columns");
    ConstTwoDMatrix llincidence{lead_lag_incidence_mx};

    const mxArray *nnzderivatives_mx = mxGetField(M_mx, 0, "NNZDerivatives");
    if (!(nnzderivatives_mx && mxIsDouble(nnzderivatives_mx)))
      mexErrMsgTxt("M_.NNZDerivatives should be a double precision array");
    ConstVector NNZD{nnzderivatives_mx};
    if (NNZD.length() < kOrder || NNZD[kOrder-1] == -1)
      mexErrMsgTxt("The derivatives were not computed for the required order. Make sure that you used the right order option inside the `stoch_simul' command");

    const mxArray *nnzderivatives_obj_mx = mxGetField(M_mx, 0, "NNZDerivatives_objective");
    if (!(nnzderivatives_obj_mx && mxIsDouble(nnzderivatives_obj_mx)))
      mexErrMsgTxt("M_.NNZDerivatives should be a double precision array");
    ConstVector NNZD_obj{nnzderivatives_obj_mx};
    if (NNZD.length() < kOrder || NNZD_obj[kOrder-1] == -1)
      mexErrMsgTxt("The derivatives were not computed for the required order. Make sure that you used the right order option inside the `stoch_simul' command");

    const mxArray *endo_names_mx = mxGetField(M_mx, 0, "endo_names");
    if (!(endo_names_mx && mxIsCell(endo_names_mx) && mxGetNumberOfElements(endo_names_mx) == static_cast<size_t>(nEndo)))
      mexErrMsgTxt("M_.endo_names should be a cell array of M_.endo_nbr elements");
    std::vector<std::string> endoNames = DynareMxArrayToString(endo_names_mx);

    const mxArray *exo_names_mx = mxGetField(M_mx, 0, "exo_names");
    if (!(exo_names_mx && mxIsCell(exo_names_mx) && mxGetNumberOfElements(exo_names_mx) == static_cast<size_t>(nExog)))
      mexErrMsgTxt("M_.exo_names should be a cell array of M_.exo_nbr elements");
    std::vector<std::string> exoNames = DynareMxArrayToString(exo_names_mx);

    const mxArray *dynamic_tmp_nbr_mx = mxGetField(M_mx, 0, "dynamic_tmp_nbr");
    if (!(dynamic_tmp_nbr_mx && mxIsDouble(dynamic_tmp_nbr_mx) && mxGetNumberOfElements(dynamic_tmp_nbr_mx) >= static_cast<size_t>(kOrder+1)))
      mexErrMsgTxt("M_.dynamic_tmp_nbr should be a double precision array with strictly more elements than the order of derivation");
    int ntt = std::accumulate(mxGetPr(dynamic_tmp_nbr_mx), mxGetPr(dynamic_tmp_nbr_mx)+kOrder+1, 0);

    // Extract various fields from dr
    const mxArray *ys_mx = mxGetField(dr_mx, 0, "ys"); // and not in order of dr.order_var
    if (!(ys_mx && mxIsDouble(ys_mx)))
      mexErrMsgTxt("dr.ys should be a double precision array");
    Vector ySteady{ConstVector{ys_mx}};
    if (!ySteady.isFinite())
      mexErrMsgTxt("dr.ys contains NaN or Inf");

    const mxArray *order_var_mx = mxGetField(dr_mx, 0, "order_var");
    if (!(order_var_mx && mxIsDouble(order_var_mx) && mxGetNumberOfElements(order_var_mx) == static_cast<size_t>(nEndo)))
      mexErrMsgTxt("dr.order_var should be a double precision array of M_.endo_nbr elements");
    std::vector<int> dr_order(nEndo);
    std::transform(mxGetPr(order_var_mx), mxGetPr(order_var_mx)+nEndo, dr_order.begin(),
                   [](double x) { return static_cast<int>(x)-1; });

    const int nSteps = 0; // Dynare++ solving steps, for time being default to 0 = deterministic steady state

    // Journal is not written on-disk, unless options_.debug = true (see #1735)
    Journal journal;
    if (debug)
      journal = Journal{fname + ".jnl"};

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


    // construct main K-order approximation class
    Approximation app(dynare, journal, nSteps, false, qz_criterium);
    // run stochastic steady
    app.walkStochSteady();

    const mxArray *objective_tmp_nbr_mx = mxGetField(M_mx, 0, "objective_tmp_nbr");
    if (!(objective_tmp_nbr_mx && mxIsDouble(objective_tmp_nbr_mx) && mxGetNumberOfElements(objective_tmp_nbr_mx) >= static_cast<size_t>(kOrder+1)))
      mexErrMsgTxt("M_.objective_tmp_nbr should be a double precision array with strictly more elements than the order of derivation");
    int ntt_objective = std::accumulate(mxGetPr(objective_tmp_nbr_mx), mxGetPr(objective_tmp_nbr_mx)+kOrder+1, 0);

    //Getting derivatives of the planner's objective function
    std::unique_ptr<ObjectiveAC> objectiveFile;
    objectiveFile = std::make_unique<ObjectiveMFile>(fname, ntt_objective);

    // make KordwDynare object
    KordwDynare welfare(dynare, NNZD_obj, journal, modParams, std::move(objectiveFile), dr_order);

    // construct main K-order approximation class of welfare
    ApproximationWelfare appwel(welfare, discount_factor, app.get_rule_ders(), app.get_rule_ders_s(), journal);
    appwel.approxAtSteady();

    // Conditional welfare approximation at non-stochastic steady state when stochastic shocks are enabled
    const Vector &cond_approx = appwel.getCond();
    plhs[0] = mxCreateDoubleMatrix(cond_approx.length(), 1, mxREAL);
    std::copy_n(cond_approx.base(), cond_approx.length(), mxGetPr(plhs[0]));

    if (nlhs > 1)
      {
        const FoldDecisionRule &unc_fdr = appwel.getFoldUncWel();
        // Add possibly missing field names
        for (int i = static_cast<int>(U_fieldnames.size()); i <= kOrder; i++)
          U_fieldnames.emplace_back("U_" + std::to_string(i));
        // Create structure for storing derivatives in Dynare++ format
        const char *U_fieldnames_c[kOrder+1];
        for (int i = 0; i <= kOrder; i++)
          U_fieldnames_c[i] = U_fieldnames[i].c_str();
        plhs[1] = mxCreateStructMatrix(1, 1, kOrder+1, U_fieldnames_c);

        // Fill that structure
        for (int i = 0; i <= kOrder; i++)
          {
            const FFSTensor &t = unc_fdr.get(Symmetry{i});
            mxArray *tmp = mxCreateDoubleMatrix(t.nrows(), t.ncols(), mxREAL);
            const ConstVector &vec = t.getData();
            assert(vec.skip() == 1);
            std::copy_n(vec.base(), vec.length(), mxGetPr(tmp));
            mxSetField(plhs[1], 0, ("U_" + std::to_string(i)).c_str(), tmp);
          }

        if (nlhs > 2)
          {
            const FoldDecisionRule &cond_fdr = appwel.getFoldCondWel();
            // Add possibly missing field names
            for (int i = static_cast<int>(W_fieldnames.size()); i <= kOrder; i++)
              W_fieldnames.emplace_back("W_" + std::to_string(i));
            // Create structure for storing derivatives in Dynare++ format
            const char *W_fieldnames_c[kOrder+1];
            for (int i = 0; i <= kOrder; i++)
              W_fieldnames_c[i] = W_fieldnames[i].c_str();
            plhs[2] = mxCreateStructMatrix(1, 1, kOrder+1, W_fieldnames_c);
      
            // Fill that structure
            for (int i = 0; i <= kOrder; i++)
              {
                const FFSTensor &t = cond_fdr.get(Symmetry{i});
                mxArray *tmp = mxCreateDoubleMatrix(t.nrows(), t.ncols(), mxREAL);
                const ConstVector &vec = t.getData();
                assert(vec.skip() == 1);
                std::copy_n(vec.base(), vec.length(), mxGetPr(tmp));
                mxSetField(plhs[2], 0, ("W_" + std::to_string(i)).c_str(), tmp);
              }
            if (nlhs > 3)
              {

                const FGSContainer &U = appwel.get_unc_ders();

                size_t nfields = (kOrder == 1 ? 2 : (kOrder == 2 ? 6 : 12));
                const char *c_fieldnames[] = { "Uy", "Uu", "Uyy", "Uyu", "Uuu", "Uss", "Uyyy", "Uyyu", "Uyuu", "Uuuu", "Uyss", "Uuss" };
                plhs[3] = mxCreateStructMatrix(1, 1, nfields, c_fieldnames);

                copy_derivatives(plhs[3], Symmetry{1, 0, 0, 0}, U, "Uy");
                copy_derivatives(plhs[3], Symmetry{0, 1, 0, 0}, U, "Uu");
                if (kOrder >= 2)
                  {
                    copy_derivatives(plhs[3], Symmetry{2, 0, 0, 0}, U, "Uyy");
                    copy_derivatives(plhs[3], Symmetry{0, 2, 0, 0}, U, "Uuu");
                    copy_derivatives(plhs[3], Symmetry{1, 1, 0, 0}, U, "Uyu");
                    copy_derivatives(plhs[3], Symmetry{0, 0, 0, 2}, U, "Uss");
                  }
                if (kOrder >= 3)
                  {
                    copy_derivatives(plhs[3], Symmetry{3, 0, 0, 0}, U, "Uyyy");
                    copy_derivatives(plhs[3], Symmetry{0, 3, 0, 0}, U, "Uuuu");
                    copy_derivatives(plhs[3], Symmetry{2, 1, 0, 0}, U, "Uyyu");
                    copy_derivatives(plhs[3], Symmetry{1, 2, 0, 0}, U, "Uyuu");
                    copy_derivatives(plhs[3], Symmetry{1, 0, 0, 2}, U, "Uyss");
                    copy_derivatives(plhs[3], Symmetry{0, 1, 0, 2}, U, "Uuss");
                  }

                if (nlhs > 4)
                  {
                    const FGSContainer &W = appwel.get_cond_ders();

                    size_t nfields = (kOrder == 1 ? 2 : (kOrder == 2 ? 6 : 12));
                    const char *c_fieldnames[] = { "Wy", "Wu", "Wyy", "Wyu", "Wuu", "Wss", "Wyyy", "Wyyu", "Wyuu", "Wuuu", "Wyss", "Wuss" };
                    plhs[4] = mxCreateStructMatrix(1, 1, nfields, c_fieldnames);

                    copy_derivatives(plhs[4], Symmetry{1, 0, 0, 0}, W, "Wy");
                    copy_derivatives(plhs[4], Symmetry{0, 1, 0, 0}, W, "Wu");
                    if (kOrder >= 2)
                      {
                        copy_derivatives(plhs[4], Symmetry{2, 0, 0, 0}, W, "Wyy");
                        copy_derivatives(plhs[4], Symmetry{0, 2, 0, 0}, W, "Wuu");
                        copy_derivatives(plhs[4], Symmetry{1, 1, 0, 0}, W, "Wyu");
                        copy_derivatives(plhs[4], Symmetry{0, 0, 0, 2}, W, "Wss");
                      }
                    if (kOrder >= 3)
                      {
                        copy_derivatives(plhs[4], Symmetry{3, 0, 0, 0}, W, "Wyyy");
                        copy_derivatives(plhs[4], Symmetry{0, 3, 0, 0}, W, "Wuuu");
                        copy_derivatives(plhs[4], Symmetry{2, 1, 0, 0}, W, "Wyyu");
                        copy_derivatives(plhs[4], Symmetry{1, 2, 0, 0}, W, "Wyuu");
                        copy_derivatives(plhs[4], Symmetry{1, 0, 0, 2}, W, "Wyss");
                        copy_derivatives(plhs[4], Symmetry{0, 1, 0, 2}, W, "Wuss");
                      }
                  }
              }
          }
      }
  } // end of mexFunction()
} // end of extern C
