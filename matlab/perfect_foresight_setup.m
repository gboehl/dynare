function perfect_foresight_setup()
% Prepares a deterministic simulation, by filling oo_.exo_simul and oo_.endo_simul
%  
% INPUTS
%   None
%  
% OUTPUTS
%   none
%    
% ALGORITHM
%   
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 1996-2014 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global M_ options_ oo_

test_for_deep_parameters_calibration(M_);

if size(M_.lead_lag_incidence,2)-nnz(M_.lead_lag_incidence(M_.maximum_endo_lag+1,:)) > 0
    mess = ['PERFECT_FORESIGHT_SETUP: error in model specification : variable ' M_.endo_names(find(M_.lead_lag_incidence(M_.maximum_lag+1,:)==0),:)];
    mess = [mess ' doesn''t appear as current variable.'];
    error(mess)
end

if options_.periods == 0
    error('PERFECT_FORESIGHT_SETUP: number of periods for the simulation isn''t specified')
end

if ~options_.initval_file
    if isempty(options_.datafile)
        make_ex_;
        make_y_;
    else
        read_data_;
    end
end
