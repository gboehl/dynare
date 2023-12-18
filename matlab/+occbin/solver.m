function [dr, out, ss] = solver(M_, options_, dr ,steady_state, exo_steady_state, exo_det_steady_state)
% [dr, out, ss] = solver(M_,oo_,options_, dr ,steady_state, exo_steady_state, exo_det_steady_state
% Solves the model with an OBC and produces simulations/IRFs
%
% INPUT: 
% - M_                      [structure]     Matlab's structure describing the model
% - options_                [structure]     Matlab's structure containing the options
% - dr                      [structure]     model information structure
% - endo_steady_state       [vector]        steady state value for endogenous variables
% - exo_steady_state        [vector]        steady state value for exogenous variables
% - exo_det_steady_state    [vector]        steady state value for exogenous deterministic variables                                    
%
% OUTPUT: 
% - dr                      [structure]     decision rules
% - out                     [structure]     simulation result containing fields:
%                                               - linear: paths for endogenous variables ignoring OBC (linear solution)
%                                               - piecewise: paths for endogenous variables satisfying the OBC (occbin/piecewise solution)
%                                               - ys: vector of steady state values
%                                               - regime_history: information on number and time of regime transitions
% - ss                      [structure]     State space solution
%                                               - T: [n_vars by n_vars by n_shock_period] array of transition matrices
%                                               - R: [n_vars by n_exo by n_shock_period] array of shock response matrices
%                                               - C: [n_vars by n_shock_period] array of constants

% Copyright © 2021-2023 Dynare Team
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

persistent sto_M sto_dr

% check dr
solve_dr=0;
if isempty(sto_M) || isempty(sto_dr)
    solve_dr=1;
else
    inan = find(~isnan(M_.params));
    inan0 = find(~isnan(sto_M.params));
    if ~isequal(inan,inan0) || ~isequal(sto_M.params(inan),M_.params(inan))
        solve_dr=1;
    end
end

ss=struct();
if solve_dr
    if isempty(options_.qz_criterium)
        options_.qz_criterium = 1+1e-6;
    end
    if options_.order>1
        options_.order = 1;
    end

    [dr,error_flag,M_.params] = resol(0,M_,options_,dr,steady_state,exo_steady_state,exo_det_steady_state);
    out.error_flag=error_flag;
    if error_flag
        print_info(error_flag, options_.noprint, options_)
        return;
    end
    sto_dr=dr;
    sto_M=M_;
else
    dr=sto_dr;
end

if options_.occbin.simul.check_ahead_periods>options_.occbin.simul.max_check_ahead_periods
    options_.occbin.simul.check_ahead_periods=options_.occbin.simul.max_check_ahead_periods;
    disp(['occbin::options::' simul '_check_ahead_periods cannot exceed ' simul '_max_check_ahead_periods'])
    disp(['occbin::options::' simul '_check_ahead_periods is re-set to be equal to ' simul '_max_check_ahead_periods'])
end

if M_.occbin.constraint_nbr==1
    [out, ss, error_flag  ] = occbin.solve_one_constraint(M_,dr,options_.occbin.simul,solve_dr);
elseif M_.occbin.constraint_nbr==2
    [out, ss, error_flag  ] = occbin.solve_two_constraints(M_,dr,options_.occbin.simul,solve_dr);
end

out.error_flag=error_flag;
if error_flag
    print_info(error_flag, options_.noprint, options_)
    return;
end

% add back steady state
if ~options_.occbin.simul.piecewise_only
    out.linear = out.linear + out.ys';
end
out.piecewise = out.piecewise + out.ys';
out.exo_pos = options_.occbin.simul.exo_pos;
