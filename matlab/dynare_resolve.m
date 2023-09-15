function [A,B,ys,info,dr,params] = dynare_resolve(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,mode)
% function [A,B,ys,info,dr,params] = dynare_resolve(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,mode)
% Computes the linear approximation and the matrices A and B of the transition equation.
%
% Inputs:
% - M_                  [structure]     Matlab's structure describing the model
% - options_            [structure]     Matlab's structure containing the options
% - dr                  [structure]     Reduced form model.
% - endo_steady_state   [vector]        steady state value for endogenous variables
% - exo_steady_state    [vector]        steady state value for exogenous variables
% - exo_det_steady_state    [vector]    steady state value for exogenous deterministic variables
% - mode                [string]        if provided, use restricted state space
%
% Outputs:
% - A                   [double]        State transition matrix (potentially for restricted state space)
% - B                   [double]        shock impact matrix (potentially for restricted state space)
% - ys                  [double]        vector of steady state values
% - info                [double]        4 by 1 vector with exit flag and information
% - dr                  [structure]     Reduced form model.
% - params              [double]        vector of potentially updated parameters

% Copyright © 2001-2023 Dynare Team
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

[dr,info,params] =compute_decision_rules(M_,options_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);

if info(1) > 0
    A = [];
    if nargout>1
        B = [];
        if nargout>2
            ys = [];
        end
    end
    return
end

switch nargin
  case 6
    endo_nbr = M_.endo_nbr;
    nstatic = M_.nstatic;
    nspred = M_.nspred;
    iv = (1:endo_nbr)';
    ic = [ nstatic+(1:nspred) endo_nbr+(1:size(dr.ghx,2)-nspred) ]';
  case 7
    iv = dr.restrict_var_list;
    ic = dr.restrict_columns;
  otherwise
    error('dynare_resolve:: Error in the calling sequence!')
end

if nargout==1
    A = kalman_transition_matrix(dr,iv,ic);
    return
end

[A,B] = kalman_transition_matrix(dr,iv,ic);
ys = dr.ys;
