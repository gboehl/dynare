function M_=discretionary_policy_initialization(M_,options_)
% function M_=discretionary_policy_initialization(M_,options_)
% INPUTS
% - M_            [structure]     Matlab's structure describing the model (M_).
% - options_      [structure]     Matlab's structure describing the current options (options_).
%
% OUTPUTS
% - M_            [structure]     Matlab's structure describing the model (M_).

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


if options_.loglinear
    % Ensure it's ok to ignore options_ returned from stoch_simul. #1197
    error('discretionary_policy is not compatible with `loglinear` option set to 1')
end

% safeguard against issues like running ramsey policy first and then running discretion
if isfield(M_,'orig_model')
    M_.endo_nbr = M_.orig_model.endo_nbr;
    M_.endo_names = M_.orig_model.endo_names;
    M_.lead_lag_incidence = M_.orig_model.lead_lag_incidence;
    M_.maximum_lead = M_.orig_model.maximum_lead;
    M_.maximum_endo_lead = M_.orig_model.maximum_endo_lead;
    M_.maximum_lag = M_.orig_model.maximum_lag;
    M_.maximum_endo_lag = M_.orig_model.maximum_endo_lag;
end

instr_nbr=M_.endo_nbr-M_.eq_nbr;

if instr_nbr==0
    error('discretionary_policy:: There are no available instruments, because the model has as many equations as variables.')
end
if size(options_.instruments,1)< instr_nbr
    error('discretionary_policy:: There are fewer declared instruments than omitted equations.')
elseif size(options_.instruments,1)> instr_nbr
    error('discretionary_policy:: There are more declared instruments than omitted equations.')
end

instr_id=NaN(size(options_.instruments,1),1);
for j=1:size(options_.instruments,1)
    vj=deblank(options_.instruments{j});
    vj_id=strmatch(vj, M_.endo_names, 'exact');
    if ~isempty(vj_id)
        instr_id(j)=vj_id;
    else
        error([mfilename,':: instrument ',vj,' not found'])
    end
end
M_.instr_id=instr_id;
