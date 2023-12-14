function [vdec, cc, ac] = monte_carlo_moments(mm, ss, dr, M_, options_, estim_params_)
% [vdec, cc, ac] = monte_carlo_moments(mm, ss, dr, M_, options_,estim_params_)
% Conduct Monte Carlo simulation of second moments for GSA
% Inputs:
%  - dr                 [structure]     decision rules
%  - M_                 [structure]     model structure
%  - options_           [structure]     Matlab's structure describing the current options
%  - estim_params_      [structure]     characterizing parameters to be estimated
%
% Outputs:
% - vdec                [double]        variance decomposition matrix
% - cc                  [double]        vector of unique elements of cross correlation matrix
% - ac                  [cell]          autocorrelation matrix


% Copyright © 2012-2023 Dynare Team
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

[~, nc1, nsam] = size(mm);
nobs=length(options_.varobs);
disp('monte_carlo_moments: Computing theoretical moments ...')
h = dyn_waitbar(0,'Theoretical moments ...');
vdec = zeros(nobs,M_.exo_nbr,nsam);
cc = zeros(nobs,nobs,nsam);
ac = zeros(nobs,nobs*options_.ar,nsam);

for j=1:nsam
    dr.ghx = mm(:, 1:(nc1-M_.exo_nbr),j);
    dr.ghu = mm(:, (nc1-M_.exo_nbr+1):end, j);
    if ~isempty(ss)
        M_=gsa.set_shocks_param(M_,estim_params_,ss(j,:));
    end
    [vdec(:,:,j), corr, autocorr] = gsa.th_moments(dr,options_,M_);
    cc(:,:,j)=triu(corr);
    dum=NaN(nobs,nobs*options_.ar);
    for i=1:options_.ar
        dum(:,(i-1)*nobs+1:i*nobs)=autocorr{i};
    end
    ac(:,:,j)=dum;
    if mod(j,3)==0
        dyn_waitbar(j/nsam,h)
    end
end
dyn_waitbar_close(h)
skipline()
disp('... done !')
