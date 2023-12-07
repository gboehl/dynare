function [vdec, corr, autocorr, z, zz] = th_moments(dr,options_,M_)
% [vdec, corr, autocorr, z, zz] = th_moments(dr,options_,M_)
% Computes theoretical moments for GSA
%
% INPUTS
% - dr            [structure]     model information structure
% - options_      [structure]     Matlab's structure describing the current options
% - M_            [structure]     Matlab's structure describing the model
%
% OUTPUTS
% - vdec          [double]        variance decomposition matrix
% - corr          [double]        correlation matrix
% - autocorr      [cell]          contains autocorrelation or
%                                 auto- and cross-covariance matrices
% - z             [double]        matrix containing mean, standard
%                                 deviation, and variance vector
% - zz            [double]        autocorrelation matrix

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

nvar = length(options_.varobs);
if nvar == 0
    nvar = length(dr.order_var);
    ivar = [1:nvar]';
else
    ivar=zeros(nvar,1);
    for i=1:nvar
        i_tmp = strmatch(options_.varobs{i}, M_.endo_names, 'exact');
        if isempty(i_tmp)
            error('th_moments: One of the variables specified does not exist');
        else
            ivar(i) = i_tmp;
        end
    end
end

[gamma_y,stationary_vars] = th_autocovariances(dr,ivar,M_, options_);
m = dr.ys(ivar(stationary_vars));

i1 = 1:length(ivar);
s2 = diag(gamma_y{1});
sd = sqrt(s2);

z = [ m sd s2 ];

%'VARIANCE DECOMPOSITION (in percent)';
if M_.exo_nbr>1
    vdec = 100*gamma_y{options_.ar+2}(i1,:);
else
    vdec = 100*ones(size(gamma_y{1}(i1,1)));
end
%'MATRIX OF CORRELATIONS';
if options_.opt_gsa.useautocorr
    corr = gamma_y{1}(i1,i1)./(sd(i1)*sd(i1)');
    corr = corr-diag(diag(corr))+diag(diag(gamma_y{1}(i1,i1)));
else
    corr = gamma_y{1}(i1,i1);
end
if options_.ar > 0
    zz=NaN(length(ivar),options_.ar);
    %'COEFFICIENTS OF AUTOCORRELATION';
    for i=1:options_.ar
        if options_.opt_gsa.useautocorr
            autocorr{i} = gamma_y{i+1}(i1,i1);
        else
            autocorr{i} = gamma_y{i+1}(i1,i1).*(sd(i1)*sd(i1)');
        end
        zz(:,i) = diag(gamma_y{i+1}(i1,i1));
    end
end
