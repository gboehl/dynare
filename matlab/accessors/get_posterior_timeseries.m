function r = get_posterior_timeseries(type, endo)
% Returns the posterior distribution of either
% smoothed/updated/filtered/forecast timeseries for a given endogenous variable.
% "type" must be either 'smoothed', 'updated', 'filtered' or 'forecast'.
%
% For filtered variables, returns a matrix with step-ahead in lines and periods
%  in columns.
% For forecasts, returns the "point" forecast (i.e. with uncertainty about both
%  parameters and shocks).

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

global M_

endo_id = find(strcmp(endo, M_.endo_names));

if isempty(endo_id)
    error('Unknown endogenous variable %s', endo)
end

r = struct('draw', [], type, []);


% Fetch parameter draws

DrawsFiles = dir([M_.dname '/metropolis/' M_.fname '_param*' ]);
NumberDrawsFiles = 0;
for i=1:length(DrawsFiles)
    % We need to filter out the _param_irf* files
    if ~isempty(regexp(DrawsFiles(i).name, '.*_param[0-9]+\.mat'))
        NumberDrawsFiles = NumberDrawsFiles + 1;
    end
end
if NumberDrawsFiles == 0
    error('Can''t find posterior draws file(s)')
end

idx = 1;
for file = 1:NumberDrawsFiles
    load([M_.dname '/metropolis/' M_.fname '_param' int2str(file) ],'stock');
    for i = 1:size(stock, 1)
        r(idx).draw = stock(i, :);
        idx = idx+1;
    end
end


% Fetch timeseries

switch type
    case 'smoothed'
        basename = 'smooth';
    case 'updated'
        basename = 'update';
    case 'filtered'
        basename = 'filter_step_ahead';
    case 'forecast'
        basename = 'forc_point';
    otherwise
        error('Unknown type requested. Should be one of: smoothed, updated, filtered, forecast')
end


if strcmp(type, 'smoothed')
    SmoothedFiles = dir([M_.dname '/metropolis/' M_.fname '_smooth*']);
    NumberTimeseriesFiles = 0;
    for i=1:length(SmoothedFiles)
        %% We need to filter out the _smoothed_{constant,trend}* files
        if ~isempty(regexp(SmoothedFiles(i).name, '.*_smooth[0-9]+\.mat'))
            NumberTimeseriesFiles = NumberTimeseriesFiles + 1;
        end
    end
else
    NumberTimeseriesFiles = length(dir([M_.dname '/metropolis/' M_.fname '_' ...
                                          basename '*']));
end
if NumberTimeseriesFiles == 0
    error('Can''t find file(s) with posterior timeseries of requested type')
end

idx = 1;
for file = 1:NumberTimeseriesFiles
    load([M_.dname '/metropolis/' M_.fname '_' basename int2str(file) ],'stock');
    if strcmp(type, 'filtered')
        for i = 1:size(stock, 4)
            r(idx).(type) = reshape(stock(:, endo_id, :, i), size(stock,1), size(stock,3));
            idx = idx+1;
        end
    else
        for i = 1:size(stock, 3)
            r(idx).(type) = stock(endo_id, :, i)';
            idx = idx+1;
        end
    end
end
