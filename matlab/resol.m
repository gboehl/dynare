function [dr, info, M, oo] = resol(check_flag, M, options, oo)

% Computes the perturbation based decision rules of the DSGE model (orders 1 to 3)
%
% INPUTS
% - check_flag    [integer]       scalar, equal to 0 if all the approximation is required, equal to 1 if only the eigenvalues are to be computed.
% - M             [structure]     Matlab's structure describing the model (M_).
% - options       [structure]     Matlab's structure describing the current options (options_).
% - oo            [structure]     Matlab's structure containing the results (oo_).
%
% OUTPUTS
% - dr            [structure]     Reduced form model.
% - info          [integer]       scalar or vector, error code.
% - M             [structure]     Matlab's structure describing the model (M_).
% - oo            [structure]     Matlab's structure containing the results (oo_).
%
% REMARKS
% Possible values for the error codes are:
%
%   info(1)=0     ->    No error.
%   info(1)=1     ->    The model doesn't determine the current variables uniquely.
%   info(1)=2     ->    MJDGGES returned an error code.
%   info(1)=3     ->    Blanchard & Kahn conditions are not satisfied: no stable equilibrium.
%   info(1)=4     ->    Blanchard & Kahn conditions are not satisfied: indeterminacy.
%   info(1)=5     ->    Blanchard & Kahn conditions are not satisfied: indeterminacy due to rank failure.
%   info(1)=6     ->    The jacobian evaluated at the deterministic steady state is complex.
%   info(1)=19    ->    The steadystate routine has thrown an exception (inconsistent deep parameters).
%   info(1)=20    ->    Cannot find the steady state, info(2) contains the sum of square residuals (of the static equations).
%   info(1)=21    ->    The steady state is complex, info(2) contains the sum of square of imaginary parts of the steady state.
%   info(1)=22    ->    The steady has NaNs.
%   info(1)=23    ->    M_.params has been updated in the steadystate routine and has complex valued scalars.
%   info(1)=24    ->    M_.params has been updated in the steadystate routine and has some NaNs.
%   info(1)=30    ->    Ergodic variance can't be computed.

% Copyright (C) 2001-2018 Dynare Team
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

if isfield(oo,'dr')
    if isfield(oo.dr,'kstate')
        dr.kstate = oo.dr.kstate;
    end
    if isfield(oo.dr,'inv_order_var')
        dr.inv_order_var = oo.dr.inv_order_var;
    end
    if isfield(oo.dr,'order_var')
        dr.order_var = oo.dr.order_var;
    end
    if isfield(oo.dr,'restrict_var_list')
        dr.restrict_var_list = oo.dr.restrict_var_list;
    end
    if isfield(oo.dr,'restrict_columns')
        dr.restrict_columns = oo.dr.restrict_columns;
    end
    if isfield(oo.dr,'obs_var')
        dr.obs_var = oo.dr.obs_var;
    end
end

if M.exo_nbr == 0
    oo.exo_steady_state = [] ;
end

[dr.ys,M.params,info] = evaluate_steady_state(oo.steady_state,M,options,oo,~options.steadystate.nocheck);

if info(1)
    oo.dr = dr;
    return
end

if options.loglinear
    threshold = 1e-16;
    % Find variables with non positive steady state. Skip auxiliary
    % variables for lagges/leaded exogenous variables
    idx = find(dr.ys(get_all_variables_but_lagged_leaded_exogenous(M))<threshold);
    if length(idx)
        if options.debug
            variables_with_non_positive_steady_state = M.endo_names{idx};
            skipline()
            fprintf('You are attempting to simulate/estimate a loglinear approximation of a model, but\n')
            fprintf('the steady state level of the following variables is not strictly positive:\n')
            for var_iter=1:length(idx)
                fprintf(' - %s (%s)\n',deblank(variables_with_non_positive_steady_state(var_iter,:)), num2str(dr.ys(idx(var_iter))))
            end
            if isinestimationobjective()
                fprintf('You should check that the priors and/or bounds over the deep parameters are such\n')
                fprintf('that the steady state levels of all the variables are strictly positive, or consider\n')
                fprintf('a linearization of the model instead of a log linearization.\n')
            else
                fprintf('You should check that the calibration of the deep parameters is such that the\n')
                fprintf('steady state levels of all the variables are strictly positive, or consider\n')
                fprintf('a linearization of the model instead of a log linearization.\n')
            end
        end
        info(1)=26;
        info(2)=sum(dr.ys(dr.ys<threshold).^2);
        return
    end
end

if options.block
    [dr,info,M,oo] = dr_block(dr,check_flag,M,options,oo);
    dr.edim = nnz(abs(dr.eigval) > options.qz_criterium);
    dr.sdim = dr.nd-dr.edim;
else
    [dr,info] = stochastic_solvers(dr,check_flag,M,options,oo);
end
oo.dr = dr;
