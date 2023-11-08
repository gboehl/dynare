function [endogenousvariables, success, err, iter] = sim1_lbj(endogenousvariables, exogenousvariables, steadystate, M_, options_)

% Performs deterministic simulations with lead or lag on one period using the historical LBJ algorithm
%
% INPUTS
%   ...
%
% OUTPUTS
%   endogenousvariables [matrix]        All endogenous variables of the model
%   success             [logical]       Whether a solution was found
%   err                 [double]        ∞-norm of Δendogenousvariables
%   iter                [integer]       Number of iterations
%
% ALGORITHM
%   Laffargue, Boucekkine, Juillard (LBJ)
%   see Juillard (1996) Dynare: A program for the resolution and
%   simulation of dynamic models with forward variables through the use
%   of a relaxation algorithm. CEPREMAP. Couverture Orange. 9602.
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright © 1996-2023 Dynare Team
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

lead_lag_incidence = M_.lead_lag_incidence;

ny = size(endogenousvariables,1);
nyp = nnz(lead_lag_incidence(1,:));
nyf = nnz(lead_lag_incidence(3,:));
nrs = ny+nyp+nyf+1;
nrc = nyf+1;
iyf = find(lead_lag_incidence(3,:)>0);
iyp = find(lead_lag_incidence(1,:)>0);
isp = 1:nyp;
is = nyp+1:ny+nyp;
isf = iyf+nyp;
isf1 = nyp+ny+1:nyf+nyp+ny+1;
success = false;
iz = 1:ny+nyp+nyf;

function [r, g1] = bytecode_wrapper(y, xpath, params, ys, it_)
    ypath = NaN(ny, 3);
    ypath(find(lead_lag_incidence')) = y;
    [r, s] = bytecode('dynamic', 'evaluate', M_, options_, ypath, xpath(it_+(-1:1), :), params, ys, 1);
    g1 = s.g1;
end

if options_.bytecode
    dynamicmodel = @bytecode_wrapper;
else
    dynamicmodel = str2func(sprintf('%s.dynamic', M_.fname));
end

verbose = options_.verbosity;

if verbose
    printline(56)
    fprintf('MODEL SIMULATION :\n')
end

h1 = clock;

for iter = 1:options_.simul.maxit
    h2 = clock;
    c = zeros(ny*options_.periods, nrc);
    it_ = M_.maximum_lag+1;
    z = [endogenousvariables(iyp,it_-1) ; endogenousvariables(:,it_) ; endogenousvariables(iyf,it_+1)];
    [d1, jacobian] = dynamicmodel(z, exogenousvariables, M_.params, steadystate, it_);
    jacobian = [jacobian(:,iz), -d1];
    ic = 1:ny;
    icp = iyp;
    c (ic,:) = jacobian(:,is)\jacobian(:,isf1);
    for it_ = M_.maximum_lag+(2:options_.periods)
        z = [endogenousvariables(iyp,it_-1) ; endogenousvariables(:,it_) ; endogenousvariables(iyf,it_+1)];
        [d1, jacobian] = dynamicmodel(z, exogenousvariables, M_.params, steadystate, it_);
        jacobian = [jacobian(:,iz), -d1];
        jacobian(:,[isf nrs]) = jacobian(:,[isf nrs])-jacobian(:,isp)*c(icp,:);
        ic = ic + ny;
        icp = icp + ny;
        c (ic,:) = jacobian(:,is)\jacobian(:,isf1);
    end
    c = back_subst_lbj(c, ny, iyf, options_.periods);
    endogenousvariables(:,M_.maximum_lag+(1:options_.periods)) = endogenousvariables(:,M_.maximum_lag+(1:options_.periods))+c;
    err = max(max(abs(c)));
    if verbose
        fprintf('Iter: %s,\t err. = %s, \t time = %s\n', num2str(iter), num2str(err), num2str(etime(clock, h2)));
    end
    if err < options_.dynatol.x
        success = true; % Convergency obtained.
        break
    end
end

if verbose
    fprintf('\nTotal time of simulation: %s\n', num2str(etime(clock,h1)))
    if success
        printline(56)
    else
        disp('Maximum number of iterations is reached (modify option maxit).')
        printline(62)
    end
    skipline()
end

end
