function [y, out, cost] = findmin(d_index, a0, P1, Qt, Y, ZZ, opts_simul,M_, dr,endo_steady_state,exo_steady_state,exo_det_steady_state, options_)
% [y, out, cost] = findmin(d_index, a0, P1, Qt, Y, ZZ, opts_simul,M_, dr,endo_steady_state,exo_steady_state,exo_det_steady_state, options_)
% Outputs:
%  - cost               [double]        penalty 
%  - out                [structure]     Occbin's results structure
%
% Inputs
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

current_obs = Y(d_index,2)'+dr.ys(options_.varobs_id(d_index))';
err_index = find(diag(Qt(:,:,2))~=0);
F = ZZ(d_index,:)*P1(:,:,2)*ZZ(d_index,:)' ;

weights=diag(F);

filtered_errs_init = zeros(1,length(err_index));
opts_simul.varobs_id=options_.varobs_id(d_index)';
opts_simul.exo_pos=err_index; %err_index is predefined mapping from observables to shocks
opts_simul.SHOCKS = filtered_errs_init;
if opts_simul.restrict_state_space
    tmp=zeros(M_.endo_nbr,1);
    tmp(dr.restrict_var_list,1)=a0(:,1);  %updated state
    opts_simul.endo_init = tmp(dr.inv_order_var,1);
else
    opts_simul.endo_init = a0(dr.inv_order_var,1);
end


[y] = fminsearch(@cost_function,filtered_errs_init');

[cost, out] = occbin.cost_function(y, current_obs, weights, opts_simul,...
            M_, dr,endo_steady_state,exo_steady_state,exo_det_steady_state, options_);

    function cost = cost_function(x)
        cost = occbin.cost_function(x, current_obs, weights, opts_simul,...
            M_, dr,endo_steady_state,exo_steady_state,exo_det_steady_state, options_);
    end


end