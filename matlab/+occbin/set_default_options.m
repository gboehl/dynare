function options_occbin_ = set_default_options(options_occbin_,M_,flag)
% function options_occbin_ = set_default_options(options_occbin_,M_,flag)
% Sets default options for Occbin
%
% INPUTS
% - options_occbin_ [structure]     Matlab's structure describing the current options
% - M_              [structure]     Matlab's structure describing the model
% - flag            [cell]          govern what/how much to initialize
%
% OUTPUTS
% - options_occbin_ [structure]     Matlab's structure describing the current options

% Copyright © 2021 Dynare Team
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

if nargin<3
    flag='all';
end

if ismember(flag,{'all'})
    options_occbin_.solver.solve_algo=3; %solver for match_function: csolve
    options_occbin_.solver.solve_tolx=1e-10;
    options_occbin_.solver.solve_tolf=1e-5;
    options_occbin_.solver.maxit=10;
    options_occbin_.write_regimes.periods=[];
    options_occbin_.write_regimes.type='simul';
    options_occbin_.write_regimes.filename=[M_.fname '_occbin_regimes'];
end

if ismember(flag,{'filter','all'})
    options_occbin_.filter.state_covariance = false;
    options_occbin_.filter.guess_regime = false;
    options_occbin_.filter.periodic_solution = true;
    options_occbin_.filter.use_relaxation = false;
end

if ismember(flag,{'forecast','all'})
    options_occbin_.forecast.check_ahead_periods=30;
    options_occbin_.forecast.debug_flag=false;
    options_occbin_.forecast.frcst_regimes=[];
    options_occbin_.forecast.maxit=30;
    options_occbin_.forecast.qmc=0;
    options_occbin_.forecast.replic=0;
    options_occbin_.forecast.SHOCKS0=[];
end

if ismember(flag,{'irf','all'})    
    options_occbin_.irf.check_ahead_periods=30;
    options_occbin_.irf.exo_names=M_.exo_names;
    options_occbin_.irf.init_regime=[];
    options_occbin_.irf.maxit=30;
%     options_occbin_.irf.periods=options_.irf;
    options_occbin_.irf.shocksize=[];
    options_occbin_.irf.shocksigns = {'pos','neg'}; 
    options_occbin_.irf.t0=0;
end

if ismember(flag,{'likelihood','all'})
    options_occbin_.likelihood.brute_force_regime_guess = true;
    options_occbin_.likelihood.curb_retrench = false;
    options_occbin_.likelihood.first_period_occbin_update = 1;
    options_occbin_.likelihood.full_output = false;
    options_occbin_.likelihood.IF_likelihood = false;
    options_occbin_.likelihood.init_regime_history = [];
    options_occbin_.likelihood.init_binding_indicator = false(0);
    options_occbin_.likelihood.inversion_filter = false;
    options_occbin_.likelihood.IVF_shock_observable_mapping = [];
    options_occbin_.likelihood.loss_function_regime_guess = false;
    options_occbin_.likelihood.maxit = 30; % this is for occbin solver algo
    options_occbin_.likelihood.max_number_of_iterations = 10; % this is for occbin_kalman_update loop
    options_occbin_.likelihood.max_check_ahead_periods=inf;
    options_occbin_.likelihood.number_of_initial_periods_with_extra_regime_guess=0;
    options_occbin_.likelihood.periods = 3;
    options_occbin_.likelihood.check_ahead_periods=200;
    options_occbin_.likelihood.periodic_solution=false;
    options_occbin_.likelihood.piecewise_only = true;
    options_occbin_.likelihood.restrict_state_space = true;
    options_occbin_.likelihood.status=true; %initialized to false in default_option_values
    options_occbin_.likelihood.use_updated_regime = true;
    options_occbin_.likelihood.waitbar=false;
end

if ismember(flag,{'plot_irf','all'})
    options_occbin_.plot_irf.add_steadystate = 0;
    options_occbin_.plot_irf.endo_scaling_factor = [];
    options_occbin_.plot_irf.grid            = true;
    options_occbin_.plot_irf.ncols            = 3;
    options_occbin_.plot_irf.nrows            = 3;
    options_occbin_.plot_irf.shocksigns = {'pos','neg'}; 
    options_occbin_.plot_irf.simulname='';
end

if ismember(flag,{'plot_shock_decomp','all'})
    options_occbin_.plot_shock_decomp.add_steadystate = false;
    options_occbin_.plot_shock_decomp.add_zero_line = false;
    options_occbin_.plot_shock_decomp.decomp_type='qoq';
    options_occbin_.plot_shock_decomp.figure_size = [200 100 650 850];
    options_occbin_.plot_shock_decomp.grid = false;
    options_occbin_.plot_shock_decomp.graph_line=false;
    options_occbin_.plot_shock_decomp.graph_regime=false;
    options_occbin_.plot_shock_decomp.graph_simul=char('total','linear');
    options_occbin_.plot_shock_decomp.graph_simul_zoom=char('total','piecewise','linear');
    options_occbin_.plot_shock_decomp.graph_zoom=false;
    options_occbin_.plot_shock_decomp.init_names_=[];
    options_occbin_.plot_shock_decomp.lw=2;
    options_occbin_.plot_shock_decomp.mystyles = {':','--',':','-.'};
    options_occbin_.plot_shock_decomp.ncol = 3;
    options_occbin_.plot_shock_decomp.no_legend = false;
    options_occbin_.plot_shock_decomp.no_others = false;
    options_occbin_.plot_shock_decomp.T0 = []; % initial date in plot (must be >= TINIT)
    options_occbin_.plot_shock_decomp.TINIT = dates();  % date of initialized states for shock decomp
%     options_occbin_.plot_shock_decomp.use_shock_groups=options_.plot_shock_decomp.use_shock_groups;
end

if ismember(flag,{'plot_simul','all'})
    options_occbin_.plot_simul.add_steadystate  = false;
    options_occbin_.plot_simul.add_vertical_line = false;
    options_occbin_.plot_simul.cutoff           = false;
    options_occbin_.plot_simul.figure_size      = [50 50 850 850];
    options_occbin_.plot_simul.labels           = {'Linear','Piecewise','Simulation 3','Simulation 4'};
    options_occbin_.plot_simul.legend           = true;
    options_occbin_.plot_simul.length_simul     = [];
    options_occbin_.plot_simul.linewidth        = 2;
    options_occbin_.plot_simul.log_normalize_graph = false;
    options_occbin_.plot_simul.marg_h(1)        = 0.08;
    options_occbin_.plot_simul.marg_h(2)        = 0.055;
    options_occbin_.plot_simul.mycolors         = get(groot,'DefaultAxesColorOrder');     
    options_occbin_.plot_simul.my_dir           = 'OccBinSimul';
    options_occbin_.plot_simul.mystyles         = {'-','--',':','-.'};
    options_occbin_.plot_simul.mystst_simul_pos = false; 
    options_occbin_.plot_simul.ncols            = 3;
    options_occbin_.plot_simul.normalization_point = 1;
    options_occbin_.plot_simul.nrows            = 3;
    options_occbin_.plot_simul.print_emf        = false;
    options_occbin_.plot_simul.scale_y          = false;
    options_occbin_.plot_simul.scale_y1         = true;
    options_occbin_.plot_simul.scale_y2         = true;
    options_occbin_.plot_simul.simulname        = 'occbin_simul';
    options_occbin_.plot_simul.subplot_gap      = 0.07;
    options_occbin_.plot_simul.threshold        = 10^-6;
    options_occbin_.plot_simul.timeaxis         = [];
    options_occbin_.plot_simul.use_grid         = true;
    
end

if ismember(flag,{'shock_decomp','all'})
    options_occbin_.shock_decomp.additive=false;
    options_occbin_.shock_decomp.curb_retrench=false;
    options_occbin_.shock_decomp.debug=false;
    options_occbin_.shock_decomp.init_in_others=false;
    options_occbin_.shock_decomp.init_names_=[];
    options_occbin_.shock_decomp.init_total=false;
    options_occbin_.shock_decomp.init2shocks= false;
    options_occbin_.shock_decomp.main_effect=false;
    options_occbin_.shock_decomp.main_effect_init=false;
    options_occbin_.shock_decomp.maxit = 100;
    options_occbin_.shock_decomp.nfrcst=0;
    options_occbin_.shock_decomp.periods = 60;
    options_occbin_.shock_decomp.check_ahead_periods=200;
    options_occbin_.shock_decomp.shocks_only=false;
    options_occbin_.shock_decomp.total_effect=false;
    options_occbin_.shock_decomp.conditional_only=true;
    options_occbin_.shock_decomp.TINIT = dates(); % date to initialize states for shock decomp
%     options_occbin_.shock_decomp.use_shock_groups=options_.plot_shock_decomp.use_shock_groups;
end

if ismember(flag,{'simul','all'})
    options_occbin_.simul.algo_truncation = 1;
    options_occbin_.simul.debug = false;
    options_occbin_.simul.curb_retrench=false;
    options_occbin_.simul.endo_init=zeros(M_.endo_nbr,1);
    options_occbin_.simul.full_output=true;
    options_occbin_.simul.init_regime=[];
    options_occbin_.simul.init_binding_indicator=false(0);
    options_occbin_.simul.exo_pos=1:M_.exo_nbr;
    options_occbin_.simul.maxit=30;
    options_occbin_.simul.max_check_ahead_periods=inf;
    options_occbin_.simul.periods=100;
    options_occbin_.simul.check_ahead_periods=200;
    options_occbin_.simul.periodic_solution=false;
    options_occbin_.simul.periodic_solution_threshold=1;
    options_occbin_.simul.periodic_solution_strict=true;
    options_occbin_.simul.piecewise_only = false;
    options_occbin_.simul.reset_check_ahead_periods_in_new_period = false;
    options_occbin_.simul.reset_regime_in_new_period = false;
    options_occbin_.simul.restrict_state_space=false;
    options_occbin_.simul.SHOCKS=zeros(1,M_.exo_nbr);
    options_occbin_.simul.waitbar=true;
end

if ismember(flag,{'smoother','all'})
    options_occbin_.smoother.curb_retrench = false;
    options_occbin_.smoother.debug = false;
    options_occbin_.smoother.fast = false;
    options_occbin_.smoother.first_period_occbin_update = 1;
    options_occbin_.smoother.full_output = false;
%     options.occbin.smoother.init_mode = 1; % 0 = standard;  1 = unconditional frcsts zero shocks+smoothed states in each period
    options_occbin_.smoother.init_regime_history = [];
    options_occbin_.smoother.init_binding_indicator = false(0);
    options_occbin_.smoother.inversion_filter = false;
    options_occbin_.smoother.linear_smoother = true;
    options_occbin_.smoother.maxit = 30; % this is for occbin solver algo
    options_occbin_.smoother.max_check_ahead_periods=inf;
    options_occbin_.smoother.max_number_of_iterations = 10; % this is for smoother loop
    options_occbin_.smoother.periods = 3;
    options_occbin_.smoother.check_ahead_periods=200;
    options_occbin_.smoother.periodic_solution=false;
    options_occbin_.smoother.piecewise_only = true;
    options_occbin_.smoother.plot = true;
    options_occbin_.smoother.status=true;
    options_occbin_.smoother.waitbar=true;
%     options.occbin.smoother.restrict_state_space = 1;
end

if ismember(flag,{'graph','all'})
    options_occbin_.graph.steady_state=true;
end