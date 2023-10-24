function optimize_prior(options_, M_, oo_, bayestopt_, estim_params_)
% optimize_prior(options_, M_, oo_, bayestopt_, estim_params_)
% This routine computes the mode of the prior density using an optimization algorithm.
%
% INPUTS
%   options_            [structure] describing the options
%   M_                  [structure] describing the model
%   oo_                 [structure] storing the results
%   bayestopt_          [structure] describing the priors
%   estim_params_       [structure] characterizing parameters to be estimated

% Copyright © 2015-2023 Dynare Team
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

% Initialize to the prior mean
oo_.dr = set_state_space(oo_.dr,M_);
xparam1 = bayestopt_.p1;

% Pertubation of the initial condition.
look_for_admissible_initial_condition = 1; scale = 1.0; iter  = 0;
while look_for_admissible_initial_condition
    xinit = xparam1+scale*randn(size(xparam1));
    if all(xinit(:)>bayestopt_.p3) && all(xinit(:)<bayestopt_.p4)
        M_ = set_all_parameters(xinit,estim_params_,M_);
        [oo_.dr,INFO,M_.params] = resol(0,M_,options_,oo_.dr,oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
        if ~INFO(1)
            look_for_admissible_initial_condition = 0;
        end
    else
        if iter == 2000
            scale = scale/1.1;
            iter  = 0;
        else
            iter = iter+1;
        end
    end
end

% Maximization of the prior density
[xparams, lpd, hessian_mat] = ...
    maximize_prior_density(xinit, bayestopt_.pshape, ...
                           bayestopt_.p6, ...
                           bayestopt_.p7, ...
                           bayestopt_.p3, ...
                           bayestopt_.p4,options_,M_,bayestopt_,estim_params_,oo_);

% Display the results.
skipline(2)
disp('------------------')
disp('PRIOR OPTIMIZATION')
disp('------------------')
skipline()
for i = 1:length(xparams)
    disp(['deep parameter ' int2str(i) ': ' get_the_name(i,0,M_,estim_params_,options_) '.'])
    disp(['  Initial condition ....... ' num2str(xinit(i)) '.'])
    disp(['  Prior mode .............. ' num2str(bayestopt_.p5(i)) '.'])
    disp(['  Optimized prior mode .... ' num2str(xparams(i)) '.'])
    skipline()
end
skipline()