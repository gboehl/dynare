function plan = init_plan(date)
% Creates and initializes a new forecast scenario
%
% INPUTS
%  o date                 [dates]           The period of the forecast
%
%
% OUTPUTS
%  plan                   [structure]       Returns a structure containing a new forecast scenario
%
%
% Copyright (C) 2013-2018 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
global M_
plan = struct();
plan.date = date;
plan.date_str = strings(date);
plan.endo_names = M_.endo_names(1:M_.orig_endo_nbr);
plan.exo_names = M_.exo_names(1:M_.exo_nbr);
plan.constrained_vars_ = [];
plan.constrained_paths_ = [];
plan.constrained_date_ = [];
plan.constrained_int_date_ = [];
plan.constrained_str_date_ = [];
plan.constrained_perfect_foresight_ = [];
plan.shock_vars_ = [];
plan.shock_paths_ = [];
plan.shock_date_ = [];
plan.shock_int_date_ = [];
plan.shock_str_date_ = [];
plan.shock_perfect_foresight_ = [];
plan.options_cond_fcst_ = struct();
plan.options_cond_fcst_.parameter_set = 'calibration';
plan.options_cond_fcst_.simulation_type = 'deterministic';
plan.options_cond_fcst_.controlled_varexo = [];
