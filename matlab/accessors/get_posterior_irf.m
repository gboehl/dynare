function r = get_posterior_irf(endo, exo)

% Copyright © 2020 Dynare Team
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

global M_ options_

endo_id = find(strcmp(endo, M_.endo_names));
exo_id = find(strcmp(exo, M_.exo_names));
if isempty(exo_id)
    error('Unknown exogenous: %s', exo)
end

ivar = get_variables_list(options_,M_);

if ~ismember(endo_id, ivar)
    error('Posterior IRFs have not been computed for %s', endo)
end

r = struct('draw', [], 'irf', []);


% Fetch parameter draws

NumberDrawsFiles = length(dir([M_.dname '/metropolis/' M_.fname '_param_irf*' ]));
if NumberDrawsFiles == 0
    error('Can''t find posterior draws file(s)')
end

idx = 1;
for file = 1:NumberDrawsFiles
    load([M_.dname '/metropolis/' M_.fname '_param_irf' int2str(file) ],'stock');
    for i = 1:size(stock, 1)
        r(idx).draw = stock(i, :);
        idx = idx+1;
    end
end


% Fetch IRFs

% We use the raw files, i.e. those that are *not* the output of
% ReshapeMatFiles.m. Their filename contains 'irf_dsge' in *lowercase*.
filesUpperAndLower = dir([M_.dname '/metropolis/' M_.fname '_irf_dsge*']);
filesLower = cellfun(@any, regexp({filesUpperAndLower.name}, [M_.fname '_irf_dsge\d\.mat']));
NumberIRFsFiles = sum(filesLower);
if NumberIRFsFiles == 0
    error('Can''t find posterior IRFs file(s)')
end

endo_id_varlist = find(ivar == endo_id);

idx = 1;
for file = 1:NumberIRFsFiles
    load([M_.dname '/metropolis/' M_.fname '_irf_dsge' int2str(file) ],'stock_irf_dsge');
    for i = 1:size(stock_irf_dsge, 4)
        r(idx).irf = stock_irf_dsge(:, endo_id_varlist, exo_id, i);
        idx = idx+1;
    end
end
