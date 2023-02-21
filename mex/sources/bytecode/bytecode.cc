/*
 * Copyright © 2007-2023 Dynare Team
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

#include <ctime>
#include <cmath>
#include <cstring>
#include <type_traits>

#include "Interpreter.hh"
#include "ErrorHandling.hh"

string
Get_Argument(const mxArray *prhs)
{
  const mxArray *mxa = prhs;
  auto buflen = mwSize(mxGetM(mxa) * mxGetN(mxa) + 1);
  char *first_argument;
  first_argument = static_cast<char *>(mxCalloc(buflen, sizeof(char)));
  size_t status = mxGetString(mxa, first_argument, buflen);
  if (status != 0)
    mexWarnMsgTxt("Not enough space. The first argument is truncated.");
  string f(first_argument);
  mxFree(first_argument);
  return f;
}

string
deblank(string x)
{
  for (int i = 0; i < static_cast<int>(x.length()); i++)
    if (x[i] == ' ')
      x.erase(i--, 1);
  return x;
}

void
Get_Arguments_and_global_variables(int nrhs,
                                   const mxArray *prhs[],
                                   int &count_array_argument,
                                   double *yd[], size_t &row_y, size_t &col_y,
                                   double *xd[], size_t &row_x, size_t &col_x,
                                   double *params[],
                                   double *steady_yd[], size_t &steady_row_y, size_t &steady_col_y,
                                   unsigned int &periods,
                                   mxArray *block_structur[],
                                   bool &steady_state, bool &block_decomposed,
                                   bool &evaluate, int &block,
                                   mxArray *M_[], mxArray *oo_[], mxArray *options_[], bool &global_temporary_terms,
                                   bool &print,
                                   bool &print_error,
                                   mxArray *GlobalTemporaryTerms[],
                                   string *plan_struct_name, string *pfplan_struct_name, bool *extended_path, mxArray *ep_struct[])
{
  size_t pos;
  *extended_path = false;
  for (int i = 0; i < nrhs; i++)
    {
#ifdef DEBUG
      if (mxIsChar(prhs[i]))
        mexPrintf("Arg %d: %s\n", i, Get_Argument(prhs[i]).c_str());
#endif
      if (!mxIsChar(prhs[i]))
        {
          switch (count_array_argument)
            {
            case 0:
              *yd = mxGetPr(prhs[i]);
              row_y = mxGetM(prhs[i]);
              col_y = mxGetN(prhs[i]);
              break;
            case 1:
              *xd = mxGetPr(prhs[i]);
              row_x = mxGetM(prhs[i]);
              col_x = mxGetN(prhs[i]);
              break;
            case 2:
              *params = mxGetPr(prhs[i]);
              break;
            case 3:
              *steady_yd = mxGetPr(prhs[i]);
              steady_row_y = mxGetM(prhs[i]);
              steady_col_y = mxGetN(prhs[i]);
              break;
            case 4:
              periods = static_cast<int>(mxGetScalar(prhs[i]));
              break;
            case 5:
              *block_structur = mxDuplicateArray(prhs[i]);
              break;
            case 6:
              global_temporary_terms = true;
              *GlobalTemporaryTerms = mxDuplicateArray(prhs[i]);
              break;
            default:
              mexPrintf("Unknown argument count_array_argument=%d\n", count_array_argument);
              break;
            }
          count_array_argument++;
        }
      else
        if (Get_Argument(prhs[i]) == "static")
          steady_state = true;
        else if (Get_Argument(prhs[i]) == "dynamic")
          steady_state = false;
        else if (Get_Argument(prhs[i]) == "block_decomposed")
          block_decomposed = true;
        else if (Get_Argument(prhs[i]) == "evaluate")
          evaluate = true;
        else if (Get_Argument(prhs[i]) == "global_temporary_terms")
          global_temporary_terms = true;
        else if (Get_Argument(prhs[i]) == "print")
          print = true;
        else if (Get_Argument(prhs[i]) == "no_print_error")
          print_error = false;
        else
          {
            pos = 0;
            if (Get_Argument(prhs[i]).substr(0, 5) == "block")
              {
                size_t pos1 = Get_Argument(prhs[i]).find("=", pos + 5);
                if (pos1 != string::npos)
                  pos = pos1 + 1;
                else
                  pos += 5;
                block = atoi(Get_Argument(prhs[i]).substr(pos, string::npos).c_str())-1;
              }
            else if (Get_Argument(prhs[i]).substr(0, 13) == "extended_path")
              {
                *extended_path = true;
                if ((i+1) >= nrhs)
                  *ep_struct = nullptr;
                else
                  {
                    *ep_struct = mxDuplicateArray(prhs[i + 1]);
                    i++;
                  }
              }
            else if (Get_Argument(prhs[i]).substr(0, 6) == "pfplan")
              {
                size_t pos1 = Get_Argument(prhs[i]).find("=", pos + 6);
                if (pos1 != string::npos)
                  pos = pos1 + 1;
                else
                  pos += 6;
                *pfplan_struct_name = deblank(Get_Argument(prhs[i]).substr(pos, string::npos));
              }
            else if (Get_Argument(prhs[i]).substr(0, 4) == "plan")
              {
                size_t pos1 = Get_Argument(prhs[i]).find("=", pos + 4);
                if (pos1 != string::npos)
                  pos = pos1 + 1;
                else
                  pos += 4;
                *plan_struct_name = deblank(Get_Argument(prhs[i]).substr(pos, string::npos));
              }
            else
              throw FatalException{"In main, unknown argument : " + Get_Argument(prhs[i])};
          }
    }
  if (count_array_argument > 0 && count_array_argument < 5)
    {
      if (count_array_argument == 3 && steady_state)
        periods = 1;
      else
        throw FatalException{"In main, missing arguments. All the following arguments have to be indicated y, x, params, it_, ys"};
    }
  *M_ = mexGetVariable("global", "M_");
  if (!*M_)
    throw FatalException{"In main, global variable not found: M_"};

  /* Gets variables and parameters from global workspace of Matlab */
  *oo_ = mexGetVariable("global", "oo_");
  if (!*oo_)
    throw FatalException{"In main, global variable not found: oo_"};

  *options_ = mexGetVariable("global", "options_");
  if (!*options_)
    throw FatalException{"In main, global variable not found: options_"};
}

/* The gateway routine */
void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *M_, *oo_, *options_;
  mxArray *GlobalTemporaryTerms;
  mxArray *block_structur = nullptr;
  mxArray *pfplan_struct = nullptr;
  size_t i, row_y = 0, col_y = 0, row_x = 0, col_x = 0;
  size_t steady_row_y, steady_col_y;
  int y_kmin = 0, y_kmax = 0, y_decal = 0;
  unsigned int periods = 1;
  double *direction;
  bool steady_state = false;
  bool block_decomposed {false};
  bool evaluate = false;
  int block = -1;
  double *params = nullptr;
  double *yd = nullptr, *xd = nullptr;
  int count_array_argument = 0;
  bool global_temporary_terms = false;
  bool print = false, print_error = true, print_it = false;
  double *steady_yd = nullptr;
  string plan, pfplan;
  bool extended_path;
  mxArray *extended_path_struct;

  table_conditional_local_type conditional_local;
  vector<s_plan> splan, spfplan, sextended_path, sconditional_extended_path;
  vector_table_conditional_local_type vector_conditional_local;
  table_conditional_global_type table_conditional_global;

  int max_periods = 0;

#ifdef DEBUG
  mexPrintf("**************************************\n");
  mexPrintf("ENTERING BYTECODE: nargin=%d, nargout=%d\n", nrhs, nlhs);
#endif

  try
    {
      Get_Arguments_and_global_variables(nrhs, prhs, count_array_argument,
                                         &yd, row_y, col_y,
                                         &xd, row_x, col_x,
                                         &params,
                                         &steady_yd, steady_row_y, steady_col_y,
                                         periods,
                                         &block_structur,
                                         steady_state, block_decomposed, evaluate, block,
                                         &M_, &oo_, &options_, global_temporary_terms,
                                         print, print_error, &GlobalTemporaryTerms,
                                         &plan, &pfplan, &extended_path, &extended_path_struct);
    }
  catch (GeneralException &feh)
    {
      mexErrMsgTxt(feh.message.c_str());
    }
#ifdef DEBUG
  mexPrintf("**************************************\n");
#endif
  if (!count_array_argument)
    {
      int field = mxGetFieldNumber(M_, "params");
      if (field < 0)
        mexErrMsgTxt("params is not a field of M_");
      params = mxGetPr(mxGetFieldByNumber(M_, 0, field));
    }

  BasicSymbolTable symbol_table;
  vector<string> dates;

  if (extended_path)
    {
      if (!extended_path_struct)
        mexErrMsgTxt("The 'extended_path' option must be followed by the extended_path descriptor");
      mxArray *date_str = mxGetField(extended_path_struct, 0, "date_str");
      if (!date_str)
        mexErrMsgTxt("The extended_path description structure does not contain the member: date_str");
      int nb_periods = mxGetM(date_str) * mxGetN(date_str);

      mxArray *constrained_vars_ = mxGetField(extended_path_struct, 0, "constrained_vars_");
      if (!constrained_vars_)
        mexErrMsgTxt("The extended_path description structure does not contain the member: constrained_vars_");
      mxArray *constrained_paths_ = mxGetField(extended_path_struct, 0, "constrained_paths_");
      if (!constrained_paths_)
        mexErrMsgTxt("The extended_path description structure does not contain the member: constrained_paths_");
      mxArray *constrained_int_date_ = mxGetField(extended_path_struct, 0, "constrained_int_date_");
      if (!constrained_int_date_)
        mexErrMsgTxt("The extended_path description structure does not contain the member: constrained_int_date_");
      mxArray *constrained_perfect_foresight_ = mxGetField(extended_path_struct, 0, "constrained_perfect_foresight_");
      if (!constrained_perfect_foresight_)
        mexErrMsgTxt("The extended_path description structure does not contain the member: constrained_perfect_foresight_");
      mxArray *shock_var_ = mxGetField(extended_path_struct, 0, "shock_vars_");
      if (!shock_var_)
        mexErrMsgTxt("The extended_path description structure does not contain the member: shock_vars_");
      mxArray *shock_paths_ = mxGetField(extended_path_struct, 0, "shock_paths_");
      if (!shock_paths_)
        mexErrMsgTxt("The extended_path description structure does not contain the member: shock_paths_");
      mxArray *shock_int_date_ = mxGetField(extended_path_struct, 0, "shock_int_date_");
      if (!shock_int_date_)
        mexErrMsgTxt("The extended_path description structure does not contain the member: shock_int_date_");
      mxArray *shock_str_date_ = mxGetField(extended_path_struct, 0, "shock_str_date_");
      if (!shock_str_date_)
        mexErrMsgTxt("The extended_path description structure does not contain the member: shock_str_date_");
      int nb_constrained = mxGetM(constrained_vars_) * mxGetN(constrained_vars_);
      int nb_controlled = 0;
      mxArray *options_cond_fcst_ = mxGetField(extended_path_struct, 0, "options_cond_fcst_");
      mxArray *controlled_varexo = nullptr;
      if (options_cond_fcst_)
        {
          controlled_varexo = mxGetField(options_cond_fcst_, 0, "controlled_varexo");
          nb_controlled = mxGetM(controlled_varexo) * mxGetN(controlled_varexo);
          if (nb_controlled != nb_constrained)
            mexErrMsgTxt("The number of exogenized variables and the number of exogenous controlled variables should be equal.");
        }
      double *controlled_varexo_value = nullptr;
      if (controlled_varexo)
        controlled_varexo_value = mxGetPr(controlled_varexo);
      double *constrained_var_value = mxGetPr(constrained_vars_);
      sconditional_extended_path.resize(nb_constrained);
      max_periods = 0;
      if (nb_constrained)
        {
          conditional_local.is_cond = false;
          conditional_local.var_exo = 0;
          conditional_local.var_endo = 0;
          conditional_local.constrained_value = 0;
          for (int i = 0; i < nb_periods; i++)
            {
              vector_conditional_local.clear();
              for (unsigned int j = 0; j < row_y; j++)
                {
                  conditional_local.var_endo = j;
                  vector_conditional_local.push_back(conditional_local);
                }
              table_conditional_global[i] = vector_conditional_local;
            }
        }

      vector_table_conditional_local_type vv3 = table_conditional_global[0];
      for (int i = 0; i < nb_constrained; i++)
        {
          sconditional_extended_path[i].exo_num = ceil(constrained_var_value[i]);
          sconditional_extended_path[i].var_num = ceil(controlled_varexo_value[i]);
          mxArray *Array_constrained_paths_ = mxGetCell(constrained_paths_, i);
          double *specific_constrained_paths_ = mxGetPr(Array_constrained_paths_);
          double *specific_constrained_int_date_ = mxGetPr(mxGetCell(constrained_int_date_, i));
          int nb_local_periods = mxGetM(Array_constrained_paths_) * mxGetN(Array_constrained_paths_);
          int *constrained_int_date = static_cast<int *>(mxMalloc(nb_local_periods * sizeof(int)));
          test_mxMalloc(constrained_int_date, __LINE__, __FILE__, __func__, nb_local_periods * sizeof(int));
          if (nb_periods < nb_local_periods)
            mexErrMsgTxt(("The total number of simulation periods (" + to_string(nb_periods)
                          + ") is lesser than the number of periods in the shock definitions ("
                          + to_string(nb_local_periods)).c_str());

          sconditional_extended_path[i].per_value.resize(nb_local_periods);
          sconditional_extended_path[i].value.resize(nb_periods);
          for (int j = 0; j < nb_periods; j++)
            sconditional_extended_path[i].value[j] = 0;
          for (int j = 0; j < nb_local_periods; j++)
            {
              constrained_int_date[j] = static_cast<int>(specific_constrained_int_date_[j]) - 1;
              conditional_local.is_cond = true;
              conditional_local.var_exo = sconditional_extended_path[i].var_num - 1;
              conditional_local.var_endo = sconditional_extended_path[i].exo_num - 1;
              conditional_local.constrained_value = specific_constrained_paths_[j];
              table_conditional_global[constrained_int_date[j]][sconditional_extended_path[i].exo_num - 1] = conditional_local;
              sconditional_extended_path[i].per_value[j] = { constrained_int_date[j], specific_constrained_paths_[j] };
              sconditional_extended_path[i].value[constrained_int_date[j]] = specific_constrained_paths_[j];
              max_periods = max(max_periods, constrained_int_date[j] + 1);
            }
          mxFree(constrained_int_date);
        }
      vector_table_conditional_local_type vv = table_conditional_global[0];
      double *shock_var_value = mxGetPr(shock_var_);
      int nb_shocks = mxGetM(shock_var_) * mxGetN(shock_var_);
      sextended_path.resize(nb_shocks);
      for (int i = 0; i < nb_shocks; i++)
        {
          sextended_path[i].exo_num = ceil(shock_var_value[i]);
          mxArray *Array_shock_paths_ = mxGetCell(shock_paths_, i);
          double *specific_shock_paths_ = mxGetPr(Array_shock_paths_);
          double *specific_shock_int_date_ = mxGetPr(mxGetCell(shock_int_date_, i));
          int nb_local_periods = mxGetM(Array_shock_paths_) * mxGetN(Array_shock_paths_);
          if (nb_periods < nb_local_periods)
            mexErrMsgTxt(("The total number of simulation periods (" + to_string(nb_periods)
                          + ") is lesser than the number of periods in the shock definitions ("
                          + to_string(nb_local_periods)).c_str());
          sextended_path[i].per_value.resize(nb_local_periods);
          sextended_path[i].value.resize(nb_periods);
          for (int j = 0; j < nb_periods; j++)
            sextended_path[i].value[j] = 0;
          for (int j = 0; j < nb_local_periods; j++)
            {
              sextended_path[i].per_value[j] = { static_cast<int>(specific_shock_int_date_[j]), specific_shock_paths_[j] };
              sextended_path[i].value[static_cast<int>(specific_shock_int_date_[j]-1)] = specific_shock_paths_[j];
              max_periods = max(max_periods, static_cast<int>(specific_shock_int_date_[j]));
            }
        }
      for (int i = 0; i < nb_periods; i++)
        {
          int buflen = mxGetNumberOfElements(mxGetCell(date_str, i)) + 1;
          char *buf = static_cast<char *>(mxCalloc(buflen, sizeof(char)));
          int info = mxGetString(mxGetCell(date_str, i), buf, buflen);
          if (info)
            mexErrMsgTxt("Can not allocated memory to store the date_str in the extended path descriptor");
          dates.emplace_back(buf); //string(Dates[i]);
          mxFree(buf);
        }
    }
  if (plan.length() > 0)
    {
      mxArray *plan_struct = mexGetVariable("base", plan.c_str());
      if (!plan_struct)
        mexErrMsgTxt(("Can't find the plan: " + plan).c_str());
      size_t n_plan = mxGetN(plan_struct);
      splan.resize(n_plan);
      for (int i = 0; i < static_cast<int>(n_plan); i++)
        {
          splan[i].var = "";
          splan[i].exo = "";
          mxArray *tmp = mxGetField(plan_struct, i, "exo");
          if (tmp)
            {
              char name[100];
              mxGetString(tmp, name, 100);
              splan[i].var = name;
              auto [variable_type, exo_num] = symbol_table.getIDAndType(name);
              if (variable_type == SymbolType::exogenous)
                splan[i].var_num = exo_num;
              else
                mexErrMsgTxt(("The variable '"s + name + "'  defined as var in plan is not an exogenous").c_str());
            }
          tmp = mxGetField(plan_struct, i, "var");
          if (tmp)
            {
              char name[100];
              mxGetString(tmp, name, 100);
              splan[i].exo = name;
              auto [variable_type, exo_num] = symbol_table.getIDAndType(name);
              if (variable_type == SymbolType::endogenous)
                splan[i].exo_num = exo_num;
              else
                mexErrMsgTxt(("The variable '"s + name + "'  defined as exo in plan is not an endogenous variable").c_str());
            }
          tmp = mxGetField(plan_struct, i, "per_value");
          if (tmp)
            {
              size_t num_shocks = mxGetM(tmp);
              splan[i].per_value.resize(num_shocks);
              double *per_value = mxGetPr(tmp);
              for (int j = 0; j < static_cast<int>(num_shocks); j++)
                splan[i].per_value[j] = { ceil(per_value[j]), per_value[j + num_shocks] };
            }
        }
      int i = 0;
      for (auto & it : splan)
        {
          mexPrintf("----------------------------------------------------------------------------------------------------\n");
          mexPrintf("surprise #%d\n", i+1);
          if (it.exo.length())
            mexPrintf(" plan fliping var=%s (%d) exo=%s (%d) for the following periods and with the following values:\n", it.var.c_str(), it.var_num, it.exo.c_str(), it.exo_num);
          else
            mexPrintf(" plan shocks on var=%s for the following periods and with the following values:\n", it.var.c_str());
          for (auto &[period, value]: it.per_value)
            mexPrintf("  %3d %10.5f\n", period, value);
          i++;
        }
    }

  if (pfplan.length() > 0)
    {
      pfplan_struct = mexGetVariable("base", pfplan.c_str());
      if (!pfplan_struct)
        mexErrMsgTxt(("Can't find the pfplan: " + pfplan).c_str());
      size_t n_plan = mxGetN(pfplan_struct);
      spfplan.resize(n_plan);
      for (int i = 0; i < static_cast<int>(n_plan); i++)
        {
          spfplan[i].var = "";
          spfplan[i].exo = "";
          mxArray *tmp = mxGetField(pfplan_struct, i, "var");
          if (tmp)
            {
              char name[100];
              mxGetString(tmp, name, 100);
              spfplan[i].var = name;
              auto [variable_type, exo_num] = symbol_table.getIDAndType(name);
              if (variable_type == SymbolType::exogenous)
                splan[i].var_num = exo_num;
              else
                mexErrMsgTxt(("The variable '"s + name + "' defined as var in pfplan is not an exogenous").c_str());
            }
          tmp = mxGetField(pfplan_struct, i, "exo");
          if (tmp)
            {
              char name[100];
              mxGetString(tmp, name, 100);
              spfplan[i].exo = name;
              auto [variable_type, exo_num] = symbol_table.getIDAndType(name);
              if (variable_type == SymbolType::endogenous)
                spfplan[i].exo_num = exo_num;
              else
                mexErrMsgTxt(("The variable '"s + name + "' defined as exo in pfplan  is not an endogenous variable").c_str());
            }
          tmp = mxGetField(pfplan_struct, i, "per_value");
          if (tmp)
            {
              size_t num_shocks = mxGetM(tmp);
              double *per_value = mxGetPr(tmp);
              spfplan[i].per_value.resize(num_shocks);
              for (int j = 0; j < static_cast<int>(num_shocks); j++)
                spfplan[i].per_value[j] = { ceil(per_value[j]), per_value[j+ num_shocks] };
            }
        }
      int i = 0;
      for (auto & it : spfplan)
        {
          mexPrintf("----------------------------------------------------------------------------------------------------\n");
          mexPrintf("perfect foresight #%d\n", i+1);
          if (it.exo.length())
            mexPrintf(" plan flipping var=%s (%d) exo=%s (%d) for the following periods and with the following values:\n", it.var.c_str(), it.var_num, it.exo.c_str(), it.exo_num);
          else
            mexPrintf(" plan shocks on var=%s (%d) for the following periods and with the following values:\n", it.var.c_str(), it.var_num);
          for (auto &[period, value] : it.per_value)
            mexPrintf("  %3d %10.5f\n", period, value);
          i++;
        }
    }

  int field_steady_state = mxGetFieldNumber(oo_, "steady_state");
  if (field_steady_state < 0)
    mexErrMsgTxt("steady_state is not a field of oo_");
  int field_exo_steady_state = mxGetFieldNumber(oo_, "exo_steady_state");
  if (field_exo_steady_state < 0)
    mexErrMsgTxt("exo_steady_state is not a field of oo_");

  if (!steady_state)
    {
      int field_endo_simul = mxGetFieldNumber(oo_, "endo_simul");
      if (field_endo_simul < 0)
        mexErrMsgTxt("endo_simul is not a field of oo_");

      int field_exo_simul = mxGetFieldNumber(oo_, "exo_simul");
      if (field_exo_simul < 0)
        mexErrMsgTxt("exo_simul is not a field of oo_");

      if (!count_array_argument)
        {
          mxArray *endo_sim_arr = mxGetFieldByNumber(oo_, 0, field_endo_simul);
          yd = mxGetPr(endo_sim_arr);
          row_y = mxGetM(endo_sim_arr);
          col_y = mxGetN(endo_sim_arr);
          mxArray *exo_sim_arr = mxGetFieldByNumber(oo_, 0, field_exo_simul);
          xd = mxGetPr(exo_sim_arr);
          row_x = mxGetM(exo_sim_arr);
          col_x = mxGetN(exo_sim_arr);
        }
      int field = mxGetFieldNumber(M_, "maximum_lag");
      if (field >= 0)
        y_kmin = static_cast<int>(floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, field)))));
      else
        mexErrMsgTxt("maximum_lag is not a field of M_");
      field = mxGetFieldNumber(M_, "maximum_lead");
      if (field >= 0)
        y_kmax = static_cast<int>(floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, field)))));
      else
        mexErrMsgTxt("maximum_lead is not a field of M_");
      field = mxGetFieldNumber(M_, "maximum_endo_lag");
      if (field >= 0)
        y_decal = max(0, y_kmin-static_cast<int>(floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, field))))));
      else
        mexErrMsgTxt("maximum_endo_lag is not a field of M_");

      if (!count_array_argument)
        {
          int field = mxGetFieldNumber(options_, "periods");
          if (field >= 0)
            periods = static_cast<int>(floor(*(mxGetPr(mxGetFieldByNumber(options_, 0, field)))));
          else
            mexErrMsgTxt("options_ is not a field of options_");
        }

      if (!steady_yd)
        {
          mxArray *steady_state_arr = mxGetFieldByNumber(oo_, 0, field_steady_state);
          steady_yd = mxGetPr(steady_state_arr);
          steady_row_y = mxGetM(steady_state_arr);
          steady_col_y = mxGetN(steady_state_arr);
        }
    }
  else
    {
      if (!count_array_argument)
        {
          mxArray *steady_state_arr = mxGetFieldByNumber(oo_, 0, field_steady_state);
          yd = mxGetPr(steady_state_arr);
          row_y = mxGetM(steady_state_arr);
          col_y = mxGetN(steady_state_arr);

          mxArray *exo_steady_state_arr = mxGetFieldByNumber(oo_, 0, field_exo_steady_state);
          xd = mxGetPr(exo_steady_state_arr);
          row_x = mxGetM(exo_steady_state_arr);
          col_x = mxGetN(exo_steady_state_arr);
        }
    }
  int field = mxGetFieldNumber(options_, "verbosity");
  int verbose = 0;
  if (field >= 0)
    verbose = static_cast<int>(*mxGetPr((mxGetFieldByNumber(options_, 0, field))));
  else
    mexErrMsgTxt("verbosity is not a field of options_");
  if (verbose)
    print_it = true;
  if (!steady_state)
    field = mxGetFieldNumber(options_, "simul");
  else
    field = mxGetFieldNumber(options_, "steady");
  mxArray *temporaryfield;
  if (field >= 0)
    temporaryfield = mxGetFieldByNumber(options_, 0, field);
  else
    {
      if (!steady_state)
        mexErrMsgTxt("simul is not a field of options_");
      else
        mexErrMsgTxt("steady is not a field of options_");
    }
  field = mxGetFieldNumber(temporaryfield, "maxit");
  if (field < 0)
    {
      if (!steady_state)
        mexErrMsgTxt("maxit is not a field of options_.simul");
      else
        mexErrMsgTxt("maxit is not a field of options_.steady");
    }
  int maxit_ = static_cast<int>(floor(*mxGetPr(mxGetFieldByNumber(temporaryfield, 0, field))));
  field = mxGetFieldNumber(options_, "markowitz");
  if (field < 0)
    mexErrMsgTxt("markowitz is not a field of options_");
  auto markowitz_c = static_cast<double>(*mxGetPr(mxGetFieldByNumber(options_, 0, field)));
  field = mxGetFieldNumber(options_, "minimal_solving_periods");
  if (field < 0)
    mexErrMsgTxt("minimal_solving_periods is not a field of options_");
  int minimal_solving_periods = static_cast<int>(*mxGetPr(mxGetFieldByNumber(options_, 0, field)));
  field = mxGetFieldNumber(options_, "stack_solve_algo");
  if (field < 0)
    mexErrMsgTxt("stack_solve_algo is not a field of options_");
  int stack_solve_algo = static_cast<int>(*mxGetPr(mxGetFieldByNumber(options_, 0, field)));
  int solve_algo;
  double solve_tolf;

  if (steady_state)
    {
      int field = mxGetFieldNumber(options_, "solve_algo");
      if (field >= 0)
        solve_algo = static_cast<int>(*mxGetPr(mxGetFieldByNumber(options_, 0, field)));
      else
        mexErrMsgTxt("solve_algo is not a field of options_");
      field = mxGetFieldNumber(options_, "solve_tolf");
      if (field >= 0)
        solve_tolf = *mxGetPr(mxGetFieldByNumber(options_, 0, field));
      else
        mexErrMsgTxt("solve_tolf is not a field of options_");
    }
  else
    {
      solve_algo = stack_solve_algo;
      int field = mxGetFieldNumber(options_, "dynatol");
      mxArray *dynatol;
      if (field >= 0)
        dynatol = mxGetFieldByNumber(options_, 0, field);
      else
        mexErrMsgTxt("dynatol is not a field of options_");
      field = mxGetFieldNumber(dynatol, "f");
      if (field >= 0)
        solve_tolf = *mxGetPr(mxGetFieldByNumber(dynatol, 0, field));
      else
        mexErrMsgTxt("f is not a field of options_.dynatol");
    }
  field = mxGetFieldNumber(M_, "fname");
  mxArray *mxa;
  if (field >= 0)
    mxa = mxGetFieldByNumber(M_, 0, field);
  else
    mexErrMsgTxt("fname is not a field of M_");
  size_t buflen = mxGetM(mxa) * mxGetN(mxa) + 1;
  char *fname = static_cast<char *>(mxCalloc(buflen+1, sizeof(char)));
  size_t status = mxGetString(mxa, fname, static_cast<int>(buflen));
  fname[buflen] = ' ';
  if (status != 0)
    mexWarnMsgTxt("Not enough space. Filename is truncated.");
  string file_name = fname;
  mxFree(fname);

  if (stack_solve_algo == 7 && !steady_state)
    mexErrMsgTxt("Bytecode: Can't use option stack_solve_algo=7");

  if (steady_state && !evaluate && (solve_algo < 5 || solve_algo > 8))
    mexErrMsgTxt("Bytecode: solve_algo must be between 5 and 8 when using the internal steady state solver");

  size_t size_of_direction = col_y*row_y*sizeof(double);
  auto *y = static_cast<double *>(mxMalloc(size_of_direction));
  test_mxMalloc(y, __LINE__, __FILE__, __func__, size_of_direction);
  auto *ya = static_cast<double *>(mxMalloc(size_of_direction));
  test_mxMalloc(ya, __LINE__, __FILE__, __func__, size_of_direction);
  direction = static_cast<double *>(mxMalloc(size_of_direction));
  test_mxMalloc(direction, __LINE__, __FILE__, __func__, size_of_direction);
  memset(direction, 0, size_of_direction);
  auto *x = static_cast<double *>(mxMalloc(col_x*row_x*sizeof(double)));
  test_mxMalloc(x, __LINE__, __FILE__, __func__, col_x*row_x*sizeof(double));
  for (i = 0; i < row_x*col_x; i++)
    x[i] = static_cast<double>(xd[i]);
  for (i = 0; i < row_y*col_y; i++)
    {
      y[i] = static_cast<double>(yd[i]);
      ya[i] = static_cast<double>(yd[i]);
    }
  size_t y_size = row_y;
  size_t nb_row_x = row_x;

  clock_t t0 = clock();
  Interpreter interprete {params, y, ya, x, steady_yd, direction, y_size, nb_row_x,
                          periods, y_kmin, y_kmax, maxit_, solve_tolf, size_of_direction, y_decal,
                          markowitz_c, file_name, minimal_solving_periods, stack_solve_algo,
                          solve_algo, global_temporary_terms, print, print_error, GlobalTemporaryTerms,
                          steady_state, block_decomposed, print_it, col_x, col_y, symbol_table};
  int nb_blocks = 0;
  double *pind;

  if (extended_path)
    {
      try
        {
          interprete.extended_path(file_name, evaluate, block, nb_blocks, max_periods, sextended_path, sconditional_extended_path, dates, table_conditional_global);
        }
      catch (GeneralException &feh)
        {
          // Release the lock on dynamic.bin for MATLAB+Windows, see #1815
          interprete.Close_SaveCode();
          mexErrMsgTxt(feh.message.c_str());
        }
    }
  else
    {
      try
        {
          interprete.compute_blocks(file_name, evaluate, block, nb_blocks);
        }
      catch (GeneralException &feh)
        {
          // Release the lock on dynamic.bin for MATLAB+Windows, see #1815
          interprete.Close_SaveCode();
          mexErrMsgTxt(feh.message.c_str());
        }
    }

  clock_t t1 = clock();
  if (!steady_state && !evaluate && print)
    mexPrintf("Simulation Time=%f milliseconds\n",
              1000.0*(static_cast<double>(t1)-static_cast<double>(t0))/static_cast<double>(CLOCKS_PER_SEC));
  bool dont_store_a_structure = false;
  if (nlhs > 0)
    {
      if (block >= 0)
        {
          if (evaluate)
            {
              vector<double> residual = interprete.get_residual();
              plhs[0] = mxCreateDoubleMatrix(static_cast<int>(residual.size()/static_cast<double>(col_y)),
                                             static_cast<int>(col_y), mxREAL);
              pind = mxGetPr(plhs[0]);
              for (i = 0; i < residual.size(); i++)
                pind[i] = residual[i];
            }
          else
            {
              int out_periods = extended_path ? max_periods + y_kmin : row_y;
              plhs[0] = mxCreateDoubleMatrix(out_periods, static_cast<int>(col_y), mxREAL);
              pind = mxGetPr(plhs[0]);
              for (i = 0; i < out_periods*col_y; i++)
                pind[i] = y[i];
            }
        }
      else
        {
          int out_periods = extended_path ? max_periods + y_kmin : col_y;
          plhs[0] = mxCreateDoubleMatrix(static_cast<int>(row_y), out_periods, mxREAL);
          pind = mxGetPr(plhs[0]);
          if (evaluate)
            {
              vector<double> residual = interprete.get_residual();
              for (i = 0; i < residual.size(); i++)
                pind[i] = residual[i];
            }
          else
            for (i = 0; i < row_y*out_periods; i++)
              pind[i] = y[i];
        }
      if (nlhs > 1)
        {
          if (evaluate)
            {
              int jacob_field_number = 0, jacob_exo_field_number = 0,
                jacob_exo_det_field_number = 0;
              if (!block_structur)
                {
                  const char *field_names[] = {"g1", "g1_x", "g1_xd"};
                  jacob_field_number = 0;
                  jacob_exo_field_number = 1;
                  jacob_exo_det_field_number = 2;
                  mwSize dims[1] = { static_cast<mwSize>(nb_blocks) };
                  plhs[1] = mxCreateStructArray(1, dims, std::extent_v<decltype(field_names)>, field_names);
                }
              else if (!mxIsStruct(block_structur))
                {
                  plhs[1] = interprete.get_jacob(0);
                  dont_store_a_structure = true;
                }
              else
                {
                  plhs[1] = block_structur;
                  jacob_field_number = mxAddField(plhs[1], "g1");
                  if (jacob_field_number == -1)
                    mexErrMsgTxt("Fatal error in bytecode: in main, cannot add extra field jacob to the structArray");
                  jacob_exo_field_number = mxAddField(plhs[1], "g1_x");
                  if (jacob_exo_field_number == -1)
                    mexErrMsgTxt("Fatal error in bytecode: in main, cannot add extra field jacob_exo to the structArray");
                  jacob_exo_det_field_number = mxAddField(plhs[1], "g1_xd");
                  if (jacob_exo_det_field_number == -1)
                    mexErrMsgTxt("Fatal error in bytecode: in main, cannot add extra field jacob_exo_det to the structArray");
                }
              if (!dont_store_a_structure)
                for (int i = 0; i < nb_blocks; i++)
                  {
                    mxSetFieldByNumber(plhs[1], i, jacob_field_number, interprete.get_jacob(i));
                    if (!steady_state)
                      {
                        mxSetFieldByNumber(plhs[1], i, jacob_exo_field_number, interprete.get_jacob_exo(i));
                        mxSetFieldByNumber(plhs[1], i, jacob_exo_det_field_number, interprete.get_jacob_exo_det(i));
                      }
                  }
            }
          else
            {
              plhs[1] = mxCreateDoubleMatrix(static_cast<int>(row_x), static_cast<int>(col_x), mxREAL);
              pind = mxGetPr(plhs[1]);
              for (i = 0; i < row_x*col_x; i++)
                pind[i] = x[i];
            }
          if (nlhs > 2)
            {
              plhs[2] = mxCreateDoubleMatrix(static_cast<int>(row_y), static_cast<int>(col_y), mxREAL);
              pind = mxGetPr(plhs[2]);
              for (i = 0; i < row_y*col_y; i++)
                pind[i] = y[i];
              if (nlhs > 3)
                {
                  mxArray *GlobalTemporaryTerms = interprete.get_Temporary_Terms();
                  size_t nb_temp_terms = mxGetM(GlobalTemporaryTerms);
                  plhs[3] = mxCreateDoubleMatrix(static_cast<int>(nb_temp_terms), 1, mxREAL);
                  pind = mxGetPr(plhs[3]);
                  double *tt = mxGetPr(GlobalTemporaryTerms);
                  for (i = 0; i < nb_temp_terms; i++)
                    pind[i] = tt[i];
                }
            }

        }
    }
  if (x)
    mxFree(x);
  if (y)
    mxFree(y);
  if (ya)
    mxFree(ya);
  if (direction)
    mxFree(direction);
}
