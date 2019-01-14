function [fp, lp] = get_ols_start_end_dates(Y, lhssub, X, jsonmodel)
% function [fp, lp] = get_ols_start_end_dates(Y, lhssub, X, jsonmodel)
% Find the common first observed date and last observed date for X, Y, and
% lhssub. Impose sample tag if passed to equation
%
% INPUTS
%   Y              [cell array]  dependent variables
%   lhssub         [cell array]  RHS to subtract from Y
%   X              [cell array]  regressors
%   jsonmodel      [cell array]  JSON representation of model block
%
% OUTPUTS
%   fp             [date]        first observed period
%   lp             [date]        last observed period
%
% SPECIAL REQUIREMENTS
%   none
% Copyright (C) 2019 Dynare Team
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

fp = max(Y.firstobservedperiod, X.firstobservedperiod);
lp = min(Y.lastobservedperiod, X.lastobservedperiod);
if ~isempty(lhssub)
    fp = max(fp, lhssub.firstobservedperiod);
    lp = min(lp, lhssub.lastobservedperiod);
end
if isfield(jsonmodel, 'tags') ...
        && isfield(jsonmodel.tags, 'sample') ...
        && ~isempty(jsonmodel.tags.sample)
    colon_idx = strfind(jsonmodel.tags.sample, ':');
    fsd = dates(jsonmodel.tags.sample(1:colon_idx-1));
    lsd = dates(jsonmodel.tags.sample(colon_idx+1:end));
    if fp > fsd
        warning(['The sample over which you want to estimate contains NaNs. '...
            'Adjusting estimation range to begin on: ' fp.char])
    else
        fp = fsd;
    end
    if lp < lsd
        warning(['The sample over which you want to estimate contains NaNs. '...
            'Adjusting estimation range to end on: ' lp.char])
    else
        lp = lsd;
    end
end
end
