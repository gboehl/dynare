function [modelIrf, check] = cet_irf_matching_file(modelIrf, M_, options_mom_, ys_)
% Based on replication codes for Christiano, Eichenbaum, Trabandt (2016, Econometrica) - Unemployment and the Business Cycle
% This showcases how to manipulate model IRFs in a irf_matching_file
% =========================================================================
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
% =========================================================================
% - column  1: gdp                    corresponds to  GDPAGG + cumsum(muF)
% - column  6: real wage              corresponds to  wAGG + cumsum(muF)
% - column  7: consumption            corresponds to  cAGG + cumsum(muF)
% - column  7: investment             corresponds to  iAGG + cumsum(muF) + cumsum(mupsiF)
% - column  8: rel. price investment  corresponds to  cumsum(pinvestAGG)
% - column 11: vacancies              corresponds to  vTotAGG*u
% - column 12: labor force
% - column 13: separation rate
% - column 14: job finding rate       corresponds to  fAGG*f

% initialize indicator
check = 0;

modelIrf = 100.*modelIrf;

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

for jexo=1:M_.exo_nbr
    if jexo==1
        mulev_vec = 0;
        mupsilev_vec = 0;
    else
        mulev_vec = cumsum(modelIrf(:,idmuF,jexo));
        mupsilev_vec = cumsum(modelIrf(:,idmupsiF,jexo));
    end    
    modelIrf(:,idGDPAGG,jexo) = modelIrf(:,idGDPAGG,jexo) + mulev_vec;
    modelIrf(:,idpiAGG,jexo) = modelIrf(:,idpiAGG,jexo);
    modelIrf(:,idRAGG,jexo) = modelIrf(:,idRAGG,jexo);
    modelIrf(:,idukAGG,jexo) = modelIrf(:,idukAGG,jexo);
    modelIrf(:,idlAGG,jexo) = modelIrf(:,idlAGG,jexo);
    modelIrf(:,idwAGG,jexo) = modelIrf(:,idwAGG,jexo) + mulev_vec;
    modelIrf(:,idcAGG,jexo) = modelIrf(:,idcAGG,jexo) + mulev_vec;
    modelIrf(:,idiAGG,jexo) = modelIrf(:,idiAGG,jexo) + mulev_vec + mupsilev_vec;
    modelIrf(:,idpinvestAGG,jexo) = cumsum(modelIrf(:,idpinvestAGG,jexo));
    modelIrf(:,idpunempAGG,jexo) = modelIrf(:,idpunempAGG,jexo)*M_.params(ismember(M_.param_names,'u'));
    modelIrf(:,idvTotAGG,jexo) = modelIrf(:,idvTotAGG,jexo);
    modelIrf(:,idfAGG,jexo) = modelIrf(:,idfAGG,jexo)*M_.params(ismember(M_.param_names,'f'));
end

