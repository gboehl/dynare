function read(json)
% function varargout = read(json)
% Read JSON and run perfect foresight solver. Potentially return output as
% JSON
%
% INPUTS
%   json         [string]   JSON string representing options to run perfect
%                           foresight solver
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2019-2020 Dynare Team
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

global M_ options_ oo_

%loading JSON
jm = loadjson_(json, 'SimplifyCell', 1);
data2json=struct();

M_.exo_det_length = 0;
for nshocks = 1:length(jm.stochasticshocksdescription)
	covartype=jm.stochasticshocksdescription{nshocks}.shockattributevalue;
	thisshock=(jm.stochasticshocksdescription{nshocks}.shockindex)+1;
	assoshock=(jm.stochasticshocksdescription{nshocks}.assoshockindex)+1;
	switch covartype
		case 1
			M_.Sigma_e(thisshock, thisshock) = (jm.stochasticshocksdescription{nshocks}.shockvalue)^2;
		case 2
			M_.Sigma_e(thisshock, thisshock) = jm.stochasticshocksdescription{nshocks}.shockvalue;
		case 3
			M_.Sigma_e(thisshock, assoshock) = jm.stochasticshocksdescription{nshocks}.shockvalue;
			M_.Sigma_e(assoshock, thisshock) = M_.Sigma_e(thisshock, assoshock);
			M_.sigma_e_is_diagonal = 0;
		case 4
			M_.Sigma_e(thisshock, assoshock) = 2*sqrt(M_.Sigma_e(thisshock, thisshock)*M_.Sigma_e(assoshock, assoshock));
			M_.Sigma_e(assoshock, thisshock) = M_.Sigma_e(thisshock, assoshock);
			M_.Correlation_matrix(thisshock, assoshock) = jm.stochasticshocksdescription{nshocks}.shockvalue;
			M_.Correlation_matrix(assoshock, thisshock) = M_.Correlation_matrix(thisshock, assoshock);
			M_.sigma_e_is_diagonal = 0;
		end
end

options_.irf = jm.irfperiods;
options_.nograph = 1;
options_.order = jm.taylororder;
% if jm.taylororder==3
% 	options_.k_order_solver = 3;
% end
var_list_ = char();
[~, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_);

irfnames=fieldnames(oo_.irfs);
for jj = 1:numel(fieldnames(oo_.irfs))
	data2json.irfs.(strtrim(char(irfnames(jj))))=oo_.irfs.(irfnames{jj});
end

savejson('',data2json,'stochsimout.JSON');

return;
