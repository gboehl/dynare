function [oo_, out, ss] = solver(M_,oo_,options_)
% function [oo_, out, ss] = solver(M_,oo_,options_)
% Solves the model with an OBC and produces simulations/IRFs
%
% INPUT: 
% - M_                  [structure]     Matlab's structure describing the model
% - oo_                 [structure]     Matlab's structure containing the results
% - options_            [structure]     Matlab's structure containing the options
%
% OUTPUT: 
% - oo_                 [structure]     Matlab's structure containing the results
% - out                 [structure]     simulation result containing fields:
%                                           - linear: paths for endogenous variables ignoring OBC (linear solution)
%                                           - piecewise: paths for endogenous variables satisfying the OBC (occbin/piecewise solution)
%                                           - ys: vector of steady state values
%                                           - regime_history: information on number and time of regime transitions
% - ss                  [structure]     State space solution
%                                           - T: [n_vars by n_vars by n_shock_period] array of transition matrices
%                                           - R: [n_vars by n_exo by n_shock_period] array of shock response matrices
%                                           - C: [n_vars by n_shock_period] array of constants

% Copyright (C) 2021 Dynare Team
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

if solve_dr
    if isempty(options_.qz_criterium)
        options_.qz_criterium = 1+1e-6;
    end

    [dr,error_flag,M_,oo_] = resol(0,M_,options_,oo_);
    out.error_flag=error_flag;
    if error_flag
        print_info(error_flag, options_.noprint, options_)
        return;
    end
    oo_.dr = dr;
    sto_dr=dr;
    sto_M=M_;
else
    oo_.dr=sto_dr;
end

if M_.occbin.constraint_nbr==1
    [out, ss, error_flag  ] = occbin.solve_one_constraint(M_,oo_.dr,options_.occbin.simul,solve_dr);
elseif M_.occbin.constraint_nbr==2
    [out, ss, error_flag  ] = occbin.solve_two_constraints(M_,oo_.dr,options_.occbin.simul,solve_dr);
end

out.error_flag=error_flag;
if error_flag
    print_info(error_flag, options_.noprint, options_)
    return;
end

% add back steady state
if ~options_.occbin.simul.piecewise_only
    if ~isoctave && matlab_ver_less_than('9.1') % Automatic broadcasting was introduced in MATLAB R2016b
        out.linear = bsxfun(@plus, out.linear, out.ys');
    else
        out.linear = out.linear + out.ys';
    end
end
if ~isoctave && matlab_ver_less_than('9.1') % Automatic broadcasting was introduced in MATLAB R2016b
    out.piecewise = bsxfun(@plus, out.piecewise, out.ys');
else
    out.piecewise = out.piecewise + out.ys';
end
out.exo_pos = options_.occbin.simul.exo_pos;

oo_.occbin.simul=out;
