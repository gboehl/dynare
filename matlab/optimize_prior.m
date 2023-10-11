function optimize_prior(DynareOptions, ModelInfo, DynareResults, BayesInfo, EstimationInfo)

% This routine computes the mode of the prior density using an optimization algorithm.

% Copyright © 2015-2017 Dynare Team
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
DynareResults.dr = set_state_space(DynareResults.dr,ModelInfo);
xparam1 = BayesInfo.p1;

% Pertubation of the initial condition.
look_for_admissible_initial_condition = 1; scale = 1.0; iter  = 0;
while look_for_admissible_initial_condition
    xinit = xparam1+scale*randn(size(xparam1));
    if all(xinit(:)>BayesInfo.p3) && all(xinit(:)<BayesInfo.p4)
        ModelInfo = set_all_parameters(xinit,EstimationInfo,ModelInfo);
        [DynareResults.dr,INFO,ModelInfo.params] = resol(0,ModelInfo,DynareOptions,DynareResults.dr,DynareResults.steady_state, DynareResults.exo_steady_state, DynareResults.exo_det_steady_state);
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
    maximize_prior_density(xinit, BayesInfo.pshape, ...
                           BayesInfo.p6, ...
                           BayesInfo.p7, ...
                           BayesInfo.p3, ...
                           BayesInfo.p4,DynareOptions,ModelInfo,BayesInfo,EstimationInfo,DynareResults);

% Display the results.
skipline(2)
disp('------------------')
disp('PRIOR OPTIMIZATION')
disp('------------------')
skipline()
for i = 1:length(xparams)
    disp(['deep parameter ' int2str(i) ': ' get_the_name(i,0,ModelInfo,EstimationInfo,DynareOptions) '.'])
    disp(['  Initial condition ....... ' num2str(xinit(i)) '.'])
    disp(['  Prior mode .............. ' num2str(BayesInfo.p5(i)) '.'])
    disp(['  Optimized prior mode .... ' num2str(xparams(i)) '.'])
    skipline()
end
skipline()