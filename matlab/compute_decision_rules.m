function [dr,info,params] =compute_decision_rules(M_,options_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% function [dr,info,params] =compute_decision_rules(M_,options_,oo_)
% INPUTS
% - M_            [structure]     Matlab's structure describing the model (M_).
% - options_      [structure]     Matlab's structure describing the current options (options_).
% - dr            [structure]     Reduced form model.
% - endo_steady_state       [vector]     steady state value for endogenous variables
% - exo_steady_state        [vector]     steady state value for exogenous variables
% - exo_det_steady_state    [vector]     steady state value for exogenous deterministic variables                                    
%
% OUTPUTS
% - dr            [structure]     Reduced form model.
% - info          [integer]       scalar or vector, error code.
% - params        [double]        vector of potentially updated parameters

% Copyright © 2020-2023 Dynare Team
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

if options_.discretionary_policy
    [dr,info,params] = discretionary_policy_1(M_,options_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
else
    [dr,info,params] = resol(0,M_,options_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
end
