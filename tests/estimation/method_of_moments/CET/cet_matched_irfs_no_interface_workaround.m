function [matched_irfs, matched_irfs_weights] = cet_matched_irfs_no_interface_workaround(endo_names,exo_names)
% Based on replication codes for Christiano, Eichenbaum, Trabandt (2016, Econometrica) - Unemployment and the Business Cycle
% This currently replaces the interface for the IRF Matching capabilities of the method_of_moments toolbox.
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

%% settings
irf_horizon = 15;        % horizon of impulse responses to match
do_monetary_shock_only = 0; % if = 0 all shocks used in estimation

%% load VAR impulse responses from Christiano, Trabandt and Walentin (2010)- Handbook of Monetary Economics Chapter( see main text for reference).
% IRFFF: impulse responses with respect to monetary policy shock
% IRFFz: impulse responses with respect to neutral tech shock
% IRFFu: impulse responses with respect to invest tech shock
% IRFFFSE, IRFzSE, IRFuSE contain the corresponding estimated standard errors
% dimensions: rows are periods, columns are variables:
% - column  1: gdp                    corresponds to  GDPAGG + cumsum(muF)
% - column  2: inflation              corresponds to  piAGG
% - column  3: federal funds rate     corresponds to  RAGG
% - column  4: capacity utilization   corresponds to  ukAGG
% - column  5: total hours            corresponds to  lAGG
% - column  6: real wage              corresponds to  wAGG + cumsum(muF)
% - column  7: consumption            corresponds to  cAGG + cumsum(muF)
% - column  7: investment             corresponds to  iAGG + cumsum(muF) + cumsum(mupsiF)
% - column  8: rel. price investment  corresponds to  cumsum(pinvestAGG)
% - column 10: unemployment rate      corresponds to  unempAGG
% - column 11: vacancies              corresponds to  vTotAGG*u
% - column 12: labor force
% - column 13: separation rate
% - column 14: job finding rate       corresponds to  fAGG*f
load('cet_data','IRFz','IRFzSE','IRFFF','IRFFFSE','IRFu','IRFuSE');

%% map empirical irf data to a model variable
% note that any further required transformations or manipulations to model variables (such as cumsum, adding muF and mupsiF)
% as well as selection of which periods to match occurs in an extra function cet_irf_matching.m
% if no such function is given then the mapping is exact and the whole horizon will be considered

% irfs with respect to monetary shock
if do_monetary_shock_only
    SHOCKNAMES = {'epsR_eps'};
else
    SHOCKNAMES = {'epsR_eps', 'muz_eps', 'mupsi_eps'};
end
VARNAMES = {'GDPAGG' 'piAGG' 'RAGG' 'ukAGG' 'lAGG' 'wAGG' 'cAGG' 'iAGG' 'pinvestAGG' 'unempAGG' 'vTotAGG' 'fAGG'};
RESCALE  = [1        400     400    1       1      1      1      1      1            100        1         100   ];
idx=1;
for jexo = 1:length(SHOCKNAMES)
    id_shock = strmatch(SHOCKNAMES{jexo},exo_names);
    if SHOCKNAMES{jexo}=="epsR_eps"
        IRF = -1*IRFFF(:,[1:11 14]); IRFSE = IRFFFSE(:,[1:11 14]);
    elseif SHOCKNAMES{jexo}=="muz_eps"
        IRF = IRFz(:,[1:11 14]); IRFSE = IRFzSE(:,[1:11 14]);
    elseif SHOCKNAMES{jexo}=="mupsi_eps"
        IRF = IRFu(:,[1:11 14]); IRFSE = IRFuSE(:,[1:11 14]);
    end
    for jvar = 1:length(VARNAMES)
        id_var = strmatch(VARNAMES{jvar},endo_names);
        for jirf=1:irf_horizon
            if IRF(jirf,jvar) ~= 0
                matched_irfs(idx,:) = {[id_var, id_shock, jirf], IRF(jirf,jvar)/RESCALE(jvar)};
                matched_irfs_weights(idx,:) = {[id_var, id_shock, jirf], [id_var, id_shock, jirf], 1/(IRFSE(jirf,jvar)/RESCALE(jvar))^2 };
                idx = idx+1;
            end            
        end
    end
end

end
