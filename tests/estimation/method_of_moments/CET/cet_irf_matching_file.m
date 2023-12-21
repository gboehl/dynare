function [modelIrf, error_indicator] = cet_irf_matching_file(modelIrf, M_, options_mom_, ys_)
% [modelIrf, error_indicator] = cet_irf_matching_file(modelIrf, M_, options_mom_, ys_)
% -------------------------------------------------------------------------
% This file manipulates model IRFs to be consistent with empirical IRFS
% -------------------------------------------------------------------------
% INPUTS
% - modelIrf:        [options_mom_.irf by M_.endo_nbr by M_.exo_nbr]
%                                 array of IRFs for all model variables and all shocks
% - M_:              [structure]  Dynare model structure
% - options_mom_:    [structure]  Dynare options structure
% - ys_:             [double]     steady state values of all endogenous variables
% -------------------------------------------------------------------------
% OUTPUTS
% - modelIrf:        [options_mom_.irf by M_.endo_nbr by M_.exo_nbr]
%                                 modified array of IRFs for all model variables and all shocks
% - error_indicator: [boolean]    indicator of success (0) or failure (1)
% -------------------------------------------------------------------------
% This function is called by
% - mom.run
% -------------------------------------------------------------------------

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
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

% initialize error indicator
error_indicator = 0;

% get indices of variables
idmuF        = ismember(M_.endo_names,'muF');
idmupsiF     = ismember(M_.endo_names,'mupsiF');
idGDPAGG     = ismember(M_.endo_names,'GDPAGG');
idpiAGG      = ismember(M_.endo_names,'piAGG');
idRAGG       = ismember(M_.endo_names,'RAGG');
idukAGG      = ismember(M_.endo_names,'ukAGG');
idlAGG       = ismember(M_.endo_names,'lAGG');
idwAGG       = ismember(M_.endo_names,'wAGG');
idcAGG       = ismember(M_.endo_names,'cAGG');
idiAGG       = ismember(M_.endo_names,'iAGG');
idpinvestAGG = ismember(M_.endo_names,'pinvestAGG');
idpunempAGG  = ismember(M_.endo_names,'unempAGG');
idvTotAGG    = ismember(M_.endo_names,'vTotAGG');
idfAGG       = ismember(M_.endo_names,'fAGG');

modelIrf = 100.*modelIrf; % convert to percent deviations
for jexo=1:M_.exo_nbr
    if jexo==1
        mulev_vec = 0;
        mupsilev_vec = 0;
    else
        mulev_vec = cumsum(modelIrf(:,idmuF,jexo));
        mupsilev_vec = cumsum(modelIrf(:,idmupsiF,jexo));
    end
    modelIrf(:,idGDPAGG,jexo) = modelIrf(:,idGDPAGG,jexo) + mulev_vec; % gdp = GDPAGG + cumsum(muF)
    modelIrf(:,idpiAGG,jexo) = modelIrf(:,idpiAGG,jexo); % inflation = piAGG
    modelIrf(:,idRAGG,jexo) = modelIrf(:,idRAGG,jexo); % Nominal interest rate = RAGG
    modelIrf(:,idukAGG,jexo) = modelIrf(:,idukAGG,jexo); % capital utilization = ukAGG
    modelIrf(:,idlAGG,jexo) = modelIrf(:,idlAGG,jexo); % labor = lAGG
    modelIrf(:,idwAGG,jexo) = modelIrf(:,idwAGG,jexo) + mulev_vec; % real wage = wAGG + cumsum(muF)
    modelIrf(:,idcAGG,jexo) = modelIrf(:,idcAGG,jexo) + mulev_vec; % consumption = cAGG + cumsum(muF)
    modelIrf(:,idiAGG,jexo) = modelIrf(:,idiAGG,jexo) + mulev_vec + mupsilev_vec; % investment = iAGG + cumsum(muF) + cumsum(mupsiF)
    modelIrf(:,idpinvestAGG,jexo) = cumsum(modelIrf(:,idpinvestAGG,jexo)); % rel. price investment = cumsum(pinvestAGG)
    modelIrf(:,idpunempAGG,jexo) = modelIrf(:,idpunempAGG,jexo)*M_.params(ismember(M_.param_names,'u')); % vacancies = vTotAGG*u
    modelIrf(:,idvTotAGG,jexo) = modelIrf(:,idvTotAGG,jexo); % aggregated total vacancies = vTotAGG
    modelIrf(:,idfAGG,jexo) = modelIrf(:,idfAGG,jexo)*M_.params(ismember(M_.param_names,'f')); % job finding rate = fAGG*f
end