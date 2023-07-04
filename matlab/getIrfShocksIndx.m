function irf_shocks_indx=getIrfShocksIndx(M_, options_)
% irf_shocks_indx=getIrfShocksIndx(M_, options_)
% returns the unique indices of the exogenous shocks specified for IRF
% generation using the irf_shocks-command
%
% Inputs:
% - M_            [structure]     Matlab's structure describing the model (M_).
% - options_      [structure]     Matlab's structure describing the current options (options_).
% Outputs:
% - irf_shocks_indx: [1 by n_irf_shocks] vector storing the indices
%
% Copyright © 2011-2022 Dynare Team
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

if (isfield(options_,'irf_shocks')==0)
    irf_shocks_indx = 1:M_.exo_nbr;
else
    irf_shocks_indx = zeros(1,size(options_.irf_shocks,1));
    for i=1:size(options_.irf_shocks,1)
        irf_shocks_indx(i) = find(strcmp(deblank(options_.irf_shocks(i,:)), M_.exo_names));
    end
    irf_shocks_indx_unique=unique(irf_shocks_indx);
    if options_.debug && (length(irf_shocks_indx_unique) ~= length(irf_shocks_indx))
        fprintf('\nSTOCH_SIMUL: Warning: The IRFs for some shocks have been requested twice.\n')
        fprintf('STOCH_SIMUL: The redundant entries will be ignored.\n')
    end
    irf_shocks_indx=irf_shocks_indx_unique;
end
