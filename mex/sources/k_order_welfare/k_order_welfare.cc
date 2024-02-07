/*
 * Copyright © 2021-2024 Dynare Team
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

#include "dynamic_dll.hh"
#include "dynamic_m.hh"
#include "objective_m.hh"

#include "SylvException.hh"
#include "approximation.hh"
#include "approximation_welfare.hh"
#include "dynare_exception.hh"
#include "exception.hh"
#include "kord_exception.hh"
#include "tl_exception.hh"

#include <algorithm>
#include <cassert>
#include <string>

#include "dynmex.h"

/* Convert MATLAB Dynare endo and exo names cell array to a vector<string> array of
   string pointers. */
std::vector<std::string>
DynareMxArrayToString(const mxArray* mxFldp)
{
  assert(mxIsCell(mxFldp));
  std::vector<std::string> r;
  for (size_t i = 0; i < mxGetNumberOfElements(mxFldp); i++)
    {
      const mxArray* cell_mx = mxGetCell(mxFldp, i);
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

std::vector<std::string> W_fieldnames;

void
copy_derivatives(mxArray* destin, const Symmetry& sym, const FGSContainer& derivs,
                 const char* fieldname)
{
  const FGSTensor& x = derivs.get(sym);
  auto x_unfolded = x.unfold();
  int n = x_unfolded->nrows();
  int m = x_unfolded->ncols();
  mxArray* tmp = mxCreateDoubleMatrix(n, m, mxREAL);
  std::copy_n(x_unfolded->getData().base(), n * m, mxGetPr(tmp));
  mxSetField(destin, 0, fieldname, tmp);
}

/*
The k_order_welfare function computes the conditional welfare starting from the non-stochastic
steady state and the derivatives of the felicity and welfare functions U(h(yₜ₋₁, uₜ, σ)) and W(yₜ₋₁,
uₜ, σ). The unconditional welfare shall be introduced in a future version with the update of the
dynare_simul_.

The routine proceeds in several steps:
1. Computation of the kth-order policy functions: the code is a copy-paste from the
   k_order_perturbation.cc file
2. Computation of the kth-order derivatives of the felicity and of the welfare functions. The call
   to the approxAtSteady method in the ApproximationWelfare class carries out the necessary
   operation
   - Importing the derivatives of the felicity function with the calcDerivativesAtSteady() method of
     the KordwDynare class. It relies on the MATLAB-generated files, which are handled by the
     ObjectiveMFile class
   - Pinpointing the derivatives of the felicity and welfare functions. The performStep method of
     the KOrderWelfare class carries out the calculations,resorting to the FaaDiBruno class and its
     methods to get the needed intermediary results.
*/

extern "C"
{

  void
  mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
  {
    if (nlhs != 1 || nrhs != 3)
      mexErrMsgTxt("Must have exactly 3 input arguments and 1 output argument");

    const mxArray* dr_mx = prhs[0];
    const mxArray* M_mx = prhs[1];
    const mxArray* options_mx = prhs[2];

    auto get_int_field = [](const mxArray* struct_mx, const std::string& fieldname) {
      mxArray* field_mx = mxGetField(struct_mx, 0, fieldname.c_str());
      if (!(field_mx && mxIsScalar(field_mx) && mxIsNumeric(field_mx)))
        mexErrMsgTxt(("Field `" + fieldname + "' should be a numeric scalar").c_str());
      return static_cast<int>(mxGetScalar(field_mx));
    };

    // Extract various fields from options_
    const mxArray* ramsey_policy_mx = mxGetField(options_mx, 0, "ramsey_policy");
    if (!(ramsey_policy_mx && mxIsLogicalScalar(ramsey_policy_mx)))
      mexErrMsgTxt("options_.ramsey_policy should be a logical scalar");
    bool ramsey_policy = static_cast<bool>(mxGetScalar(ramsey_policy_mx));

    if (ramsey_policy == false)
      mexErrMsgTxt("The considered model must be a Ramsey-typed model !");

    const int kOrder = get_int_field(options_mx, "order");
    if (kOrder < 1)
      mexErrMsgTxt("options_.order must be at least 1");

    const mxArray* use_dll_mx = mxGetField(options_mx, 0, "use_dll");
    if (!(use_dll_mx && mxIsLogicalScalar(use_dll_mx)))
      mexErrMsgTxt("options_.use_dll should be a logical scalar");
    bool use_dll = static_cast<bool>(mxGetScalar(use_dll_mx));

    double qz_criterium = 1 + 1e-6;
    const mxArray* qz_criterium_mx = mxGetField(options_mx, 0, "qz_criterium");
    if (qz_criterium_mx && mxIsScalar(qz_criterium_mx) && mxIsNumeric(qz_criterium_mx))
      qz_criterium = mxGetScalar(qz_criterium_mx);

    const mxArray* threads_mx = mxGetField(options_mx, 0, "threads");
    if (!threads_mx)
      mexErrMsgTxt("Can't find field options_.threads");
    const mxArray* num_threads_mx = mxGetField(threads_mx, 0, "k_order_perturbation");
    if (!(num_threads_mx && mxIsScalar(num_threads_mx) && mxIsNumeric(num_threads_mx)))
      mexErrMsgTxt("options_.threads.k_order_perturbation be a numeric scalar");
    int num_threads = static_cast<int>(mxGetScalar(num_threads_mx));

    const mxArray* debug_mx = mxGetField(options_mx, 0, "debug");
    if (!(debug_mx && mxIsLogicalScalar(debug_mx)))
      mexErrMsgTxt("options_.debug should be a logical scalar");
    bool debug = static_cast<bool>(mxGetScalar(debug_mx));

    const mxArray* pruning_mx = mxGetField(options_mx, 0, "pruning");
    if (!(pruning_mx && mxIsLogicalScalar(pruning_mx)))
      mexErrMsgTxt("options_.pruning should be a logical scalar");
    bool pruning = static_cast<bool>(mxGetScalar(pruning_mx));

    // Extract various fields from M_
    const mxArray* fname_mx = mxGetField(M_mx, 0, "fname");
    if (!(fname_mx && mxIsChar(fname_mx) && mxGetM(fname_mx) == 1))
      mexErrMsgTxt("M_.fname should be a character string");
    std::string fname {mxArrayToString(fname_mx)};

    const mxArray* params_mx = mxGetField(M_mx, 0, "params");
    if (!(params_mx && mxIsDouble(params_mx)))
      mexErrMsgTxt("M_.params should be a double precision array");
    Vector modParams {ConstVector {params_mx}};
    if (!modParams.isFinite())
      mexErrMsgTxt("M_.params contains NaN or Inf");

    const mxArray* param_names_mx = mxGetField(M_mx, 0, "param_names");
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

    const mxArray* sigma_e_mx = mxGetField(M_mx, 0, "Sigma_e");
    if (!(sigma_e_mx && mxIsDouble(sigma_e_mx) && mxGetM(sigma_e_mx) == mxGetN(sigma_e_mx)))
      mexErrMsgTxt("M_.Sigma_e should be a double precision square matrix");
    TwoDMatrix vCov {ConstTwoDMatrix {sigma_e_mx}};
    if (!vCov.isFinite())
      mexErrMsgTxt("M_.Sigma_e contains NaN or Inf");

    const int nStat = get_int_field(M_mx, "nstatic");
    const int nPred = get_int_field(M_mx, "npred");
    const int nBoth = get_int_field(M_mx, "nboth");
    const int nForw = get_int_field(M_mx, "nfwrd");

    const int nExog = get_int_field(M_mx, "exo_nbr");
    const int nEndo = get_int_field(M_mx, "endo_nbr");
    const int nPar = get_int_field(M_mx, "param_nbr");

    const mxArray* lead_lag_incidence_mx = mxGetField(M_mx, 0, "lead_lag_incidence");
    if (!(lead_lag_incidence_mx && mxIsDouble(lead_lag_incidence_mx)
          && mxGetM(lead_lag_incidence_mx) == 3
          && mxGetN(lead_lag_incidence_mx) == static_cast<size_t>(nEndo)))
      mexErrMsgTxt("M_.lead_lag_incidence should be a double precision matrix with 3 rows and "
                   "M_.endo_nbr columns");
    ConstTwoDMatrix llincidence {lead_lag_incidence_mx};

    const mxArray* nnzderivatives_mx = mxGetField(M_mx, 0, "NNZDerivatives");
    if (!(nnzderivatives_mx && mxIsDouble(nnzderivatives_mx)))
      mexErrMsgTxt("M_.NNZDerivatives should be a double precision array");
    ConstVector NNZD {nnzderivatives_mx};
    if (NNZD.length() < kOrder || NNZD[kOrder - 1] == -1)
      mexErrMsgTxt("The derivatives were not computed for the required order. Make sure that you "
                   "used the right order option inside the `stoch_simul' command");

    const mxArray* endo_names_mx = mxGetField(M_mx, 0, "endo_names");
    if (!(endo_names_mx && mxIsCell(endo_names_mx)
          && mxGetNumberOfElements(endo_names_mx) == static_cast<size_t>(nEndo)))
      mexErrMsgTxt("M_.endo_names should be a cell array of M_.endo_nbr elements");
    std::vector<std::string> endoNames = DynareMxArrayToString(endo_names_mx);

    const mxArray* exo_names_mx = mxGetField(M_mx, 0, "exo_names");
    if (!(exo_names_mx && mxIsCell(exo_names_mx)
          && mxGetNumberOfElements(exo_names_mx) == static_cast<size_t>(nExog)))
      mexErrMsgTxt("M_.exo_names should be a cell array of M_.exo_nbr elements");
    std::vector<std::string> exoNames = DynareMxArrayToString(exo_names_mx);

    const mxArray* dynamic_tmp_nbr_mx = mxGetField(M_mx, 0, "dynamic_tmp_nbr");
    if (!(dynamic_tmp_nbr_mx && mxIsDouble(dynamic_tmp_nbr_mx) && !mxIsComplex(dynamic_tmp_nbr_mx)
          && !mxIsSparse(dynamic_tmp_nbr_mx)
          && mxGetNumberOfElements(dynamic_tmp_nbr_mx) >= static_cast<size_t>(kOrder + 1)))
      mexErrMsgTxt("M_.dynamic_tmp_nbr should be a double precision array with strictly more "
                   "elements than the order of derivation");
    int ntt
        = std::accumulate(mxGetPr(dynamic_tmp_nbr_mx), mxGetPr(dynamic_tmp_nbr_mx) + kOrder + 1, 0);

    // Extract various fields from dr
    const mxArray* ys_mx = mxGetField(dr_mx, 0, "ys"); // and not in order of dr.order_var
    if (!(ys_mx && mxIsDouble(ys_mx) && !mxIsComplex(ys_mx) && !mxIsSparse(ys_mx)))
      mexErrMsgTxt("dr.ys should be a real dense array");
    Vector ySteady {ConstVector {ys_mx}};
    if (!ySteady.isFinite())
      mexErrMsgTxt("dr.ys contains NaN or Inf");

    const mxArray* order_var_mx = mxGetField(dr_mx, 0, "order_var");
    if (!(order_var_mx && mxIsDouble(order_var_mx) && !mxIsComplex(order_var_mx)
          && !mxIsSparse(order_var_mx)
          && mxGetNumberOfElements(order_var_mx) == static_cast<size_t>(nEndo)))
      mexErrMsgTxt("dr.order_var should be a real dense array of M_.endo_nbr elements");
    std::vector<int> dr_order(nEndo);
    std::transform(mxGetPr(order_var_mx), mxGetPr(order_var_mx) + nEndo, dr_order.begin(),
                   [](double x) { return static_cast<int>(x) - 1; });

    const mxArray* objective_g1_sparse_rowval_mx
        = mxGetField(M_mx, 0, "objective_g1_sparse_rowval");
    if (!(objective_g1_sparse_rowval_mx && mxIsInt32(objective_g1_sparse_rowval_mx)))
      mexErrMsgTxt("M_.objective_g1_sparse_rowval should be an int32 array");

    const mxArray* objective_g1_sparse_colval_mx
        = mxGetField(M_mx, 0, "objective_g1_sparse_colval");
    if (!(objective_g1_sparse_colval_mx && mxIsInt32(objective_g1_sparse_colval_mx)))
      mexErrMsgTxt("M_.objective_g1_sparse_colval should be an int32 array");

    const mxArray* objective_g1_sparse_colptr_mx
        = mxGetField(M_mx, 0, "objective_g1_sparse_colptr");
    if (!(objective_g1_sparse_colptr_mx && mxIsInt32(objective_g1_sparse_colptr_mx)))
      mexErrMsgTxt("M_.objective_g1_sparse_colptr should be an int32 array");

    std::vector<const mxArray*> objective_gN_sparse_indices;
    for (int o {2}; o <= kOrder; o++)
      {
        using namespace std::string_literals;
        auto fieldname {"objective_g"s + std::to_string(o) + "_sparse_indices"};
        const mxArray* indices = mxGetField(M_mx, 0, fieldname.c_str());
        if (!(indices && mxIsInt32(indices)))
          mexErrMsgTxt(("M_."s + fieldname + " should be an int32 array").c_str());
        objective_gN_sparse_indices.push_back(indices);
      }

    const int nSteps
        = 0; // Dynare++ solving steps, for time being default to 0 = deterministic steady state

    // Journal is not written on-disk, unless options_.debug = true (see #1735)
    Journal journal;
    if (debug)
      journal = Journal {fname + ".jnl"};

    std::unique_ptr<DynamicModelAC> dynamicModelFile;
    if (use_dll)
      dynamicModelFile = std::make_unique<DynamicModelDLL>(fname, ntt, kOrder);
    else
      dynamicModelFile = std::make_unique<DynamicModelMFile>(fname, ntt);

    // intiate tensor library
    TLStatic::init(kOrder, nStat + 2 * nPred + 3 * nBoth + 2 * nForw + nExog);

    // Set number of parallel threads
    sthread::detach_thread_group::max_parallel_threads = num_threads;

    // make KordpDynare object
    KordpDynare dynare(endoNames, exoNames, nExog, nPar, ySteady, vCov, modParams, nStat, nPred,
                       nForw, nBoth, NNZD, nSteps, kOrder, journal, std::move(dynamicModelFile),
                       dr_order, llincidence);

    // construct main K-order approximation class
    Approximation app(dynare, journal, nSteps, false, pruning, qz_criterium);
    // run stochastic steady
    app.walkStochSteady();

    // Getting derivatives of the planner's objective function
    std::unique_ptr<ObjectiveMFile> objectiveFile;
    objectiveFile = std::make_unique<ObjectiveMFile>(
        fname, kOrder, objective_g1_sparse_rowval_mx, objective_g1_sparse_colval_mx,
        objective_g1_sparse_colptr_mx, objective_gN_sparse_indices);

    // make KordwDynare object
    KordwDynare welfare(dynare, journal, modParams, std::move(objectiveFile), dr_order);

    // construct main K-order approximation class of welfare
    ApproximationWelfare appwel(welfare, discount_factor, app.get_rule_ders(),
                                app.get_rule_ders_s(), journal);
    appwel.approxAtSteady();

    const FoldDecisionRule& cond_fdr = appwel.getFoldCondWel();
    // Add possibly missing field names
    for (int i = static_cast<int>(W_fieldnames.size()); i <= kOrder; i++)
      W_fieldnames.emplace_back("W_" + std::to_string(i));
    // Create structure for storing derivatives in Dynare++ format
    const char* W_fieldnames_c[kOrder + 1];
    for (int i = 0; i <= kOrder; i++)
      W_fieldnames_c[i] = W_fieldnames[i].c_str();
    plhs[0] = mxCreateStructMatrix(1, 1, kOrder + 1, W_fieldnames_c);

    // Fill that structure
    for (int i = 0; i <= kOrder; i++)
      {
        const FFSTensor& t = cond_fdr.get(Symmetry {i});
        mxArray* tmp = mxCreateDoubleMatrix(t.nrows(), t.ncols(), mxREAL);
        const ConstVector& vec = t.getData();
        assert(vec.skip() == 1);
        std::copy_n(vec.base(), vec.length(), mxGetPr(tmp));
        mxSetField(plhs[0], 0, ("W_" + std::to_string(i)).c_str(), tmp);
      }
  } // end of mexFunction()
} // end of extern C
