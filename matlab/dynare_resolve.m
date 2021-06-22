function [A,B,ys,info,M_,oo_] = dynare_resolve(M_,options_,oo_,mode)
% function [A,B,ys,info,M_,options_,oo_] = dynare_resolve(M_,options_,oo_,mode)
% Computes the linear approximation and the matrices A and B of the transition equation.
%
% Inputs:
% - M_                  [structure]     Matlab's structure describing the model
% - options_            [structure]     Matlab's structure containing the options
% - oo_                 [structure]     Matlab's structure containing the results
% - mode                [string]        if provided, use restricted state space
%
% Outputs:
% - A                   [double]        State transition matrix (potentially for restricted state space)
% - B                   [double]        shock impact matrix (potentially for restricted state space)
% - ys                  [double]        vector of steady state values
% - info                [double]        4 by 1 vector with exit flag and information
% - M_                  [structure]     Matlab's structure describing the model
% - oo_                 [structure]     Matlab's structure containing the results

% Copyright (C) 2001-2021 Dynare Team
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

[dr,info,M_,oo_] =compute_decision_rules(M_,options_,oo_);

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
  case 3
    endo_nbr = M_.endo_nbr;
    nstatic = M_.nstatic;
    nspred = M_.nspred;
    iv = (1:endo_nbr)';
    if ~options_.block
        ic = [ nstatic+(1:nspred) endo_nbr+(1:size(oo_.dr.ghx,2)-nspred) ]';
    else
        ic = oo_.dr.restrict_columns;
    end
  case 4
    iv = oo_.dr.restrict_var_list;
    ic = oo_.dr.restrict_columns;
  otherwise
    error('dynare_resolve:: Error in the calling sequence!')
end

if nargout==1
    A = kalman_transition_matrix(oo_.dr,iv,ic,M_.exo_nbr);
    return
end

[A,B] = kalman_transition_matrix(oo_.dr,iv,ic,M_.exo_nbr);
ys = oo_.dr.ys;
