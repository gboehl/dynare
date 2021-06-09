function r = get_posterior_covariance(endo1, endo2)

% Copyright (C) 2020 Dynare Team
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

id1 = find(strcmp(endo1, M_.endo_names));
if nargin < 2
    id2 = id1;
else
    id2 = find(strcmp(endo2, M_.endo_names));
end

ivar = get_variables_list(options_,M_);

if ~ismember(id1, ivar)
    error('Posterior moments have not been computed for %s', endo1)
end
if ~ismember(id2, ivar)
    error('Posterior moments have not been computed for %s', endo2)
end

r = struct('draw', [], 'moment', []);


% Fetch parameter draws

NumberDrawsFiles = length(dir([M_.dname '/metropolis/' M_.fname '_posterior_draws*' ]));
if NumberDrawsFiles == 0
    error('Can''t find posterior draws file(s)')
end

idx = 1;
for file = 1:NumberDrawsFiles
    load([M_.dname '/metropolis/' M_.fname '_posterior_draws' int2str(file) ],'pdraws');
    for i = 1:size(pdraws, 1)
        r(idx).draw = pdraws{i, 1};
        idx = idx+1;
    end
end


% Fetch moments

symidx = symmetric_matrix_index(find(ivar==id1), find(ivar==id2), length(ivar));

NumberMomentsFiles = length(dir([M_.dname '/metropolis/' M_.fname '_Posterior2ndOrderMoments*']));
if NumberMomentsFiles == 0
    error('Can''t find posterior 2nd order moments file(s)')
end

idx = 1;
for file = 1:NumberMomentsFiles
    load([M_.dname '/metropolis/' M_.fname '_Posterior2ndOrderMoments' int2str(file) ],'Covariance_matrix');
    for i = 1:size(Covariance_matrix, 1)
        r(idx).moment = Covariance_matrix(i, symidx);
        idx = idx+1;
    end
end
