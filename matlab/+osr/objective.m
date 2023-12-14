function [loss,info,exit_flag,df,vx]=objective(x,M_, oo_, options_,i_params,i_var,weights)
% [loss,info,exit_flag,df,vx]=objective(x,M_, oo_, options_,i_params,i_var,weights)
% Objective function for optimal simple rules (OSR)
% INPUTS
%   x                         vector           values of the parameters
%                                              over which to optimize
%   M_                        [structure]      Dynare's model structure
%   oo_                       [structure]      Dynare's results structure
%   options_                  [structure]      Dynare's options structure
%   i_params                  vector           index of optimizing parameters in M_.params
%   i_var                     vector           variables indices
%   weights                   vector           weights in the OSRs
%
% OUTPUTS
%   loss                      scalar           loss function returned to solver
%   info                      vector           info vector returned by resol
%   exit_flag                 scalar           exit flag returned to solver
%   df                        vectcor          Analytic Jacobian
%   vx                        vector           variances of the endogenous variables
%
% SPECIAL REQUIREMENTS
%   none
% Copyright © 2005-2023 Dynare Team
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

exit_flag = 1;
vx = [];
df=NaN(length(i_params),1);
% set parameters of the policy rule
M_.params(i_params) = x;

[oo_.dr,info,M_.params] = resol(0,M_,options_,oo_.dr ,oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);

if info(1)
    if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 ||...
            info(1) == 20 || info(1) == 21 || info(1) == 23 || info(1) == 26 || ...
            info(1) == 81 || info(1) == 84 ||  info(1) == 85
        loss = 1e8;
        info(4)=info(2);
        return
    else
        loss = 1e8;
        info(4)=0.1;
        return
    end
end

if ~options_.analytic_derivation
    vx = osr.get_variance_of_endogenous_variables(M_,options_,oo_.dr,i_var);
    loss = full(weights(:)'*vx(:));
else
    totparam_nbr=length(i_params);
    oo_.dr.derivs = identification.get_perturbation_params_derivs(M_, options_, [], oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state, i_params, [], [], 0); %analytic derivatives of perturbation matrices

    pruned_state_space = pruned_SS.pruned_state_space_system(M_, options_, oo_.dr, i_var, 0, 0, 1);
    vx = pruned_state_space.Var_y + pruned_state_space.E_y*pruned_state_space.E_y';
    dE_yy = pruned_state_space.dVar_y;
    for jp=1:length(i_params)
        dE_yy(:,:,jp) = dE_yy(:,:,jp) + pruned_state_space.dE_y(:,jp)*pruned_state_space.E_y' + pruned_state_space.E_y*pruned_state_space.dE_y(:,jp)';
    end

    model_moments_params_derivs = reshape(dE_yy,length(i_var)^2,totparam_nbr);

    df = NaN(totparam_nbr,1);
    loss = full(weights(:)'*vx(:));

    for jp=1:length(i_params)
        df(jp,1) = sum(weights(:).*model_moments_params_derivs(:,jp));
    end
end