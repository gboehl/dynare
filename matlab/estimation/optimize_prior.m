function optimize_prior(options_, M_, oo_, Prior, estim_params_, pnames)

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

oo_.dr = set_state_space(oo_.dr, M_, options_);

% Initialize to the prior mean
xparam1 = Prior.p1;

% Pertubation of the initial condition.
look_for_admissible_initial_condition = true; scale = 1.0; iter  = 0;
while look_for_admissible_initial_condition
    xinit = xparam1+scale*randn(size(xparam1));
    if all(xinit>Prior.p3) && all(xinit<Prior.p4)
        M_ = set_all_parameters(xinit, estim_params_, M_);
        [~, INFO, M_, oo_] = resol(0, M_, options_, oo_);
        if ~INFO(1)
            look_for_admissible_initial_condition = false;
        end
    else
        if iter==2000
            scale = scale/1.1;
            iter = 0;
        else
            iter = iter+1;
        end
    end
end

% Maximization of the prior density
xparams = maximize_prior_density(xinit, pnames, options_, M_, Prior, estim_params_, oo_);

% Display results.
skipline(2)
disp('------------------')
disp('PRIOR OPTIMIZATION')
disp('------------------')
skipline()
for i = 1:length(xparams)
    dprintf('deep parameter %u: %s.', i, get_the_name(i, 0, M_, estim_params_, options_.varobs))
    dprintf('  Initial condition ........ %s.', num2str(xinit(i)))
    dprintf('  Prior mode ............... %s.', num2str(Prior.p5(i)))
    dprintf('  Optimized prior mode ..... %s.', num2str(xparams(i)))
    skipline()
end
skipline()
