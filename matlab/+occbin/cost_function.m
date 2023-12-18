function [cost, out] = cost_function(err_0, current_obs, weights, opts_simul,...
                                     M_, dr,endo_steady_state,exo_steady_state,exo_det_steady_state, options_)
% [cost, out] = cost_function(err_0, current_obs, opts_simul,...
%                             M_, dr,endo_steady_state,exo_steady_state,exo_det_steady_state, options_)
% Outputs:
%  - cost               [double]        penalty 
%  - out                [structure]     Occbin's results structure
%
% Inputs
% - err_0               [double]        value of shocks 
% - current_obs         [double]        [1 by n_obs] current value of observables
% - weights             [double]        [1 by n_obs] variance of observables,
% - opts_simul          [structure]     Structure with simulation options
%                                       used in cost function
% - M_                  [structure]     Matlab's structure describing the model (M_).
% - dr_                 [structure]     model information structure
% - endo_steady_state   [vector]        steady state value for endogenous variables
% - exo_steady_state    [vector]        steady state value for exogenous variables
% - exo_det_steady_state    [vector]        steady state value for exogenous deterministic variables                                    
% - options_            [structure]     Matlab's structure describing the current options (options_).

% Copyright © 2023 Dynare Team
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


opts_simul.SHOCKS = err_0';
options_.occbin.simul=opts_simul;
options_.occbin.simul.full_output=1;
options_.noprint = 1;
[~, out] = occbin.solver(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state);

cost = 0;
if ~out.error_flag
    cost = mean((out.piecewise(1,opts_simul.varobs_id)'-current_obs').^2./weights);
else
    cost = cost+1.e10;
end