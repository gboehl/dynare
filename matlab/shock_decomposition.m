function [oo_,M_] = shock_decomposition(M_,oo_,options_,varlist,bayestopt_,estim_params_)
% function z = shock_decomposition(M_,oo_,options_,varlist)
% Computes shocks contribution to a simulated trajectory. The field set is
% oo_.shock_decomposition. It is a n_var by nshock+2 by nperiods array. The
% first nshock columns store the respective shock contributions, column n+1
% stores the role of the initial conditions, while column n+2 stores the
% value of the smoothed variables.  Both the variables and shocks are stored
% in the order of declaration, i.e. M_.endo_names and M_.exo_names, respectively.
%
% INPUTS
%    M_:          [structure]  Definition of the model
%    oo_:         [structure]  Storage of results
%    options_:    [structure]  Options
%    varlist:     [char]       List of variables
%    bayestopt_:  [structure]  describing the priors
%    estim_params_: [structure] characterizing parameters to be estimated
%
% OUTPUTS
%    oo_:         [structure]  Storage of results
%    M_:          [structure]  Definition of the model; makes sure that
%                               M_.params is correctly updated
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2009-2020 Dynare Team
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

% indices of endogenous variables

if isfield(oo_,'shock_decomposition_info') && isfield(oo_.shock_decomposition_info,'i_var')
    if isfield (oo_,'realtime_conditional_shock_decomposition') ...
            || isfield (oo_,'realtime_forecast_shock_decomposition') ...
            || isfield (oo_,'realtime_shock_decomposition') ...
            || isfield (oo_,'conditional_shock_decomposition') ...
            || isfield (oo_,'initval_decomposition')
        error('shock_decomposition::squeezed shock decompositions are already stored in oo_')
    end
end
with_epilogue = options_.shock_decomp.with_epilogue;

if isempty(varlist)
    varlist = M_.endo_names(1:M_.orig_endo_nbr);
end

[~, ~,index_uniques] = varlist_indices(varlist, M_.endo_names);
varlist = varlist(index_uniques);

% number of variables
endo_nbr = M_.endo_nbr;

% number of shocks
nshocks = M_.exo_nbr;

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


options_.selected_variables_only = 0; %make sure all variables are stored
options_.plot_priors=0;
[oo_, M_, ~, ~, Smoothed_Variables_deviation_from_mean, initial_date] = evaluate_smoother(parameter_set, varlist, M_, oo_, options_, bayestopt_, estim_params_);

options_.initial_date=initial_date;
% reduced form
dr = oo_.dr;

% data reordering
order_var = dr.order_var;
inv_order_var = dr.inv_order_var;


% coefficients
A = dr.ghx;
B = dr.ghu;

% initialization
gend = size(oo_.SmoothedShocks.(M_.exo_names{1}),1);
epsilon=NaN(nshocks,gend);
for i=1:nshocks
    epsilon(i,:) = oo_.SmoothedShocks.(M_.exo_names{i});
end

z = zeros(endo_nbr,nshocks+2,gend);

z(:,end,:) = Smoothed_Variables_deviation_from_mean;

maximum_lag = M_.maximum_lag;

k2 = dr.kstate(find(dr.kstate(:,2) <= maximum_lag+1),[1 2]);
i_state = order_var(k2(:,1))+(min(i,maximum_lag)+1-k2(:,2))*M_.endo_nbr;
for i=1:gend
    if i > 1 && i <= maximum_lag+1
        lags = min(i-1,maximum_lag):-1:1;
    end

    if i > 1
        tempx = permute(z(:,1:nshocks,lags),[1 3 2]);
        m = min(i-1,maximum_lag);
        tempx = [reshape(tempx,endo_nbr*m,nshocks); zeros(endo_nbr*(maximum_lag-i+1),nshocks)];
        z(:,1:nshocks,i) = A(inv_order_var,:)*tempx(i_state,:);
        lags = lags+1;
    end

    if i > options_.shock_decomp.init_state
        z(:,1:nshocks,i) = z(:,1:nshocks,i) + B(inv_order_var,:).*repmat(epsilon(:,i)',endo_nbr,1);
    end
    z(:,nshocks+1,i) = z(:,nshocks+2,i) - sum(z(:,1:nshocks,i),2);
end

if with_epilogue
    [z, oo_.shock_decomposition_info.epilogue_steady_state] = epilogue_shock_decomposition(z, M_, oo_);
end

oo_.shock_decomposition = z;

if ~options_.no_graph.shock_decomposition
    oo_ = plot_shock_decomposition(M_,oo_,options_,varlist);
end

oo_.gui.ran_shock_decomposition = true;
