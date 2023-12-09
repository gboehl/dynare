function oo_=perfect_foresight_setup(M_, options_, oo_)
% Prepares a deterministic simulation, by filling oo_.exo_simul and oo_.endo_simul
%
% INPUTS
%   M_                  [structure] describing the model
%   options_            [structure] describing the options
%   oo_                 [structure] storing the results
%
% OUTPUTS
%   oo_                 [structure] storing the results
%
% ALGORITHM
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 1996-2023 Dynare Team
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

test_for_deep_parameters_calibration(M_);

if size(M_.lead_lag_incidence,2)-nnz(M_.lead_lag_incidence(M_.maximum_endo_lag+1,:)) > 0
    mess = 'PERFECT_FORESIGHT_SETUP: error in model specification : the variable(s) ';
    var_list = M_.endo_names(M_.lead_lag_incidence(M_.maximum_endo_lag+1,:)==0);
    for i=1:length(var_list)
        if i<length(var_list)
            mess = [mess, var_list{i} ', '];
        else
            mess = [mess, var_list{i} ];
        end
    end
    mess = [mess ' don''t appear as current period variables.'];
    error(mess)
end

if options_.periods == 0
    error('PERFECT_FORESIGHT_SETUP: number of periods for the simulation isn''t specified')
end

if ~isempty(M_.det_shocks) && options_.periods<max([M_.det_shocks.periods])
    % Some expected shocks happen after the terminal period.
    mess = sprintf('\nPERFECT_FORESIGHT_SETUP: Problem with the declaration of the expected shocks:\n');
    for i=1:length(M_.det_shocks)
        if any(M_.det_shocks(i).periods>options_.periods)
            mess = sprintf('%s\n   At least one expected value for %s has been declared after the terminal period.', mess, M_.exo_names{M_.det_shocks(i).exo_id});
        end
    end
    disp(mess)
    skipline()
    error('PERFECT_FORESIGHT_SETUP: Please check the declaration of the shocks or increase the value of the periods option.')
end

if options_.simul.endval_steady && M_.maximum_lead == 0
    error('PERFECT_FORESIGHT_SETUP: Option endval_steady cannot be used on a purely backward or static model.')
end

oo_ = make_ex_(M_,options_,oo_);
oo_ = make_y_(M_,options_,oo_);
