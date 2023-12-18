function oo_ = initial_condition_decomposition(M_,oo_,options_,varlist,bayestopt_,estim_params_)
% function oo_ = initial_condition_decomposition(M_,oo_,options_,varlist,bayestopt_,estim_params_)
% Computes initial condition contribution to a simulated trajectory. The field set is
% oo_.initval_decomposition. It is a n_var by n_var+2 by nperiods array. The
% first n_var columns store the respective endogenous initval contribution, column n+1
% stores the role of the shocks, while column n+2 stores the
% value of the smoothed variables.  Variables are stored
% in the order of declaration, i.e. M_.endo_names.
%
% INPUTS
%    M_:            [structure]                Definition of the model
%    oo_:           [structure]                Storage of results
%    options_:      [structure]                Options
%    varlist:       [cell of char array]       List of variables
%    bayestopt_:    [structure]                Description of the priors
%    estim_params_: [structure]                Estimated parameters
%
% OUTPUTS
%    oo_:           [structure]                Storage of results
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2017-2018 Dynare Team
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

options_.plot_shock_decomp.colormap = options_.initial_condition_decomp.colormap;
options_.plot_shock_decomp.nodisplay = options_.initial_condition_decomp.nodisplay;
options_.plot_shock_decomp.graph_format = options_.initial_condition_decomp.graph_format;
options_.plot_shock_decomp.fig_name = options_.initial_condition_decomp.fig_name;
options_.plot_shock_decomp.detail_plot = options_.initial_condition_decomp.detail_plot;
options_.plot_shock_decomp.init2shocks = options_.initial_condition_decomp.init2shocks;
options_.plot_shock_decomp.steadystate = options_.initial_condition_decomp.steadystate;
options_.plot_shock_decomp.write_xls = options_.initial_condition_decomp.write_xls;
options_.plot_shock_decomp.type = options_.initial_condition_decomp.type;
options_.plot_shock_decomp.plot_init_date = options_.initial_condition_decomp.plot_init_date;
options_.plot_shock_decomp.plot_end_date = options_.initial_condition_decomp.plot_end_date;
options_.plot_shock_decomp.diff = options_.initial_condition_decomp.diff;
options_.plot_shock_decomp.flip = options_.initial_condition_decomp.flip;
options_.plot_shock_decomp.max_nrows = options_.initial_condition_decomp.max_nrows;

if isfield(options_.initial_condition_decomp,'init2shocks') % private trap for uimenu calls
    init2shocks=options_.initial_condition_decomp.init2shocks;
else
    init2shocks=[];
end
% indices of endogenous variables
if isempty(varlist)
    varlist = M_.endo_names(1:M_.orig_endo_nbr);
end

if ~isequal(varlist,0)
    [~, ~, index_uniques] = varlist_indices(varlist, M_.endo_names);
    varlist = varlist(index_uniques);
end

% number of variables
endo_nbr = M_.endo_nbr;

% parameter set
parameter_set = options_.parameter_set;
if isempty(parameter_set)
    if isfield(oo_,'posterior_mean')
        parameter_set = 'posterior_mean';
    elseif isfield(oo_,'mle_mode')
        parameter_set = 'mle_mode';
    elseif isfield(oo_,'posterior')
        parameter_set = 'posterior_mode';
    else
        error(['shock_decomposition: option parameter_set is not specified ' ...
               'and posterior mode is not available'])
    end
end

if ~isfield(oo_,'initval_decomposition') || isequal(varlist,0)
    if isfield(oo_,'shock_decomposition_info') && isfield(oo_.shock_decomposition_info,'i_var')
        if isfield (oo_,'realtime_conditional_shock_decomposition') ...
                || isfield (oo_,'realtime_forecast_shock_decomposition') ...
                || isfield (oo_,'realtime_shock_decomposition') ...
                || isfield (oo_,'conditional_shock_decomposition') ...
                || isfield (oo_,'shock_decomposition')
            error('initval_decomposition::squeezed shock decompositions are already stored in oo_')
        end
    end
    with_epilogue = options_.initial_condition_decomp.with_epilogue;
    options_.selected_variables_only = 0; %make sure all variables are stored
    options_.plot_priors=0;
    [oo_local,~,~,~,Smoothed_Variables_deviation_from_mean] = evaluate_smoother(parameter_set,varlist,M_,oo_,options_,bayestopt_,estim_params_);

    % reduced form
    dr = oo_local.dr;

    % data reordering
    order_var = dr.order_var;
    inv_order_var = dr.inv_order_var;


    % coefficients
    A = dr.ghx;
    B = dr.ghu;

    % initialization
    gend = length(oo_local.SmoothedShocks.(M_.exo_names{1})); %+options_.forecast;
    z = zeros(endo_nbr,endo_nbr+2,gend);
    z(:,end,:) = Smoothed_Variables_deviation_from_mean;

    for i=1:endo_nbr
        z(i,i,1) = Smoothed_Variables_deviation_from_mean(i,1);
    end

    maximum_lag = M_.maximum_lag;

    k2 = dr.kstate(find(dr.kstate(:,2) <= maximum_lag+1),[1 2]);
    i_state = order_var(k2(:,1))+(min(i,maximum_lag)+1-k2(:,2))*M_.endo_nbr;
    for i=1:gend
        if i > 1 && i <= maximum_lag+1
            lags = min(i-1,maximum_lag):-1:1;
        end

        if i > 1
            tempx = permute(z(:,1:endo_nbr,lags),[1 3 2]);
            m = min(i-1,maximum_lag);
            tempx = [reshape(tempx,endo_nbr*m,endo_nbr); zeros(endo_nbr*(maximum_lag-i+1),endo_nbr)];
            z(:,1:endo_nbr,i) = A(inv_order_var,:)*tempx(i_state,:);
            lags = lags+1;
        end
        z(:,endo_nbr+1,i) = z(:,endo_nbr+2,i) - sum(z(:,1:endo_nbr,i),2);

    end

    if with_epilogue
        [z, epilogue_steady_state] = epilogue_shock_decomposition(z, M_, oo_);
        if ~isfield(oo_,'shock_decomposition_info') || ~isfield(oo_.shock_decomposition_info,'epilogue_steady_state')
            oo_.shock_decomposition_info.epilogue_steady_state = epilogue_steady_state;
        end
    end
    oo_.initval_decomposition = z;
end

% when varlist==0, we only store results in oo_ and do not make any plot
if ~isequal(varlist,0)

    % if ~options_.no_graph.shock_decomposition
    oo_local=oo_;
    oo_local.shock_decomposition = oo_.initval_decomposition;
    if ~isempty(init2shocks)
        init2shocks = M_.init2shocks.(init2shocks);
        n=size(init2shocks,1);
        for i=1:n
            j=strmatch(init2shocks{i}{1},M_.endo_names,'exact');
            oo_local.shock_decomposition(:,end-1,:)=oo_local.shock_decomposition(:,j,:)+oo_local.shock_decomposition(:,end-1,:);
            oo_local.shock_decomposition(:,j,:)=0;
        end
    end
    M_.exo_names = M_.endo_names;
    M_.exo_nbr = M_.endo_nbr;
    options_.plot_shock_decomp.realtime=0;
    options_.plot_shock_decomp.screen_shocks=1;
    options_.plot_shock_decomp.use_shock_groups = '';
    options_.plot_shock_decomp.init_cond_decomp = 1; % private flag to plotting utilities

    plot_shock_decomposition(M_,oo_local,options_,varlist);
end
