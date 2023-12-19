function [initial_conditions, innovations, pfm, ep, verbosity, options_, oo_] = extended_path_initialization(initial_conditions, sample_size, exogenousvariables, options_, M_, oo_)
% [initial_conditions, innovations, pfm, ep, verbosity, options_, oo_] = extended_path_initialization(initial_conditions, sample_size, exogenousvariables, options_, M_, oo_)
% Initialization of the extended path routines.
%
% INPUTS
%  o initial_conditions     [double]    m*1 array, where m is the number of endogenous variables in the model.
%  o sample_size            [integer]   scalar, size of the sample to be simulated.
%  o exogenousvariables     [double]    T*n array, values for the structural innovations.
%  o options_               [struct]    Dynare's options structure
%  o M_                     [struct]    Dynare's model structure
%  o oo_                    [struct]    Dynare's result structure
%
% OUTPUTS
%
% ALGORITHM
%
% SPECIAL REQUIREMENTS

% Copyright © 2016-2023 Dynare Team
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

ep  = options_.ep;

% Set verbosity levels.
options_.verbosity = ep.verbosity;
verbosity = ep.verbosity+ep.debug;

% Set maximum number of iterations for the deterministic solver.
options_.simul.maxit = ep.maxit;

% Prepare a structure needed by the matlab implementation of the perfect foresight model solver
pfm = setup_stochastic_perfect_foresight_model_solver(M_, options_, oo_);

% Check that the user did not use varexo_det
if M_.exo_det_nbr~=0
    error('Extended path does not support varexo_det.')
end

% Set default initial conditions.
if isempty(initial_conditions)
    if isempty(M_.endo_histval)
        initial_conditions = oo_.steady_state;
    else
        initial_conditions = M_.endo_histval;
    end
end

% Set the number of periods for the (stochastic) perfect foresight model
pfm.periods = ep.periods;

pfm.i_upd = pfm.ny+(1:pfm.periods*pfm.ny);

pfm.block = options_.block;

% Set the algorithm for the perfect foresight solver
options_.stack_solve_algo = ep.stack_solve_algo;

% Compute the first order reduced form if needed.
dr = struct();
if ep.init
    options_.order = 1;
    oo_.dr=set_state_space(dr,M_);
    [oo_.dr,info,M_.params] = resol(0,M_,options_,oo_.dr,oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    if info(1)
        print_info(info,options_.noprint,options_);
    end
end

% Do not use a minimal number of perdiods for the perfect foresight solver (with bytecode and blocks)
options_.minimal_solving_period = options_.ep.periods;

% Set the covariance matrix of the structural innovations.
if isempty(exogenousvariables)
    innovations = struct();
    innovations.positive_var_indx = find(diag(M_.Sigma_e)>0);
    innovations.effective_number_of_shocks = length(innovations.positive_var_indx);
    innovations.covariance_matrix = M_.Sigma_e(innovations.positive_var_indx,innovations.positive_var_indx);
    innovations.covariance_matrix_upper_cholesky = chol(innovations.covariance_matrix);
else
    innovations = struct();
end

% Set seed.
if ep.set_dynare_seed_to_default
    options_=set_dynare_seed_local_options(options_,'default');
end

% hybrid correction
pfm.hybrid_order = ep.stochastic.hybrid_order;
if pfm.hybrid_order
    oo_.dr = set_state_space(oo_.dr, M_);
    options = options_;
    options.order = pfm.hybrid_order;
    [pfm.dr, M_.params] = resol(0, M_, options, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
else
    pfm.dr = [];
end

% number of nonzero derivatives
pfm.nnzA = M_.NNZDerivatives(1);

% setting up integration nodes if order > 0
if ep.stochastic.order > 0
    [nodes,weights,nnodes] = setup_integration_nodes(options_.ep,pfm);
    pfm.nodes = nodes;
    pfm.weights = weights;
    pfm.nnodes = nnodes;
    % compute number of blocks
    [block_nbr,pfm.world_nbr] = get_block_world_nbr(ep.stochastic.algo,nnodes,ep.stochastic.order,ep.periods);
else
    block_nbr = ep.periods;
end

% set boundaries if mcp
[lb,ub,pfm.eq_index] = get_complementarity_conditions(M_, options_.ramsey_policy);
if options_.ep.solve_algo == 10
    options_.lmmcp.lb = repmat(lb,block_nbr,1);
    options_.lmmcp.ub = repmat(ub,block_nbr,1);
elseif options_.ep.solve_algo == 11
    options_.mcppath.lb = repmat(lb,block_nbr,1);
    options_.mcppath.ub = repmat(ub,block_nbr,1);
end
pfm.block_nbr = block_nbr;
