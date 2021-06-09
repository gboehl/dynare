function create_sur_report()
%function create_sur_report()
% Creates report for all SUR models estimated in oo_
%
% INPUTS
% none
%
% OUTPUTS
% none
%
% SPECIAL REQUIREMENTS
%   dynare must have been run with the option: json=compute

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
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

global oo_

if ~isfield(oo_, 'sur')
    disp(['create_sur_report: to use this function you must '...
        'have already estimated a SUR model']);
    return
end

%% Begin Report
rep = report();

%% loop through SUR estimationsd
fields = fieldnames(oo_.sur);
for i = 1:length(fields)
    rep = rep.addPage(...
        'title', ['SUR model ' regexprep(fields{i}, '_', '\\_')], ...
        'titleFormat', '\large\bfseries');

    rep = rep.addSection('cols', 1);
    preamble = sprintf('No. Equations: %d\\newline Observations: %d\\newline', ...
        oo_.sur.(fields{i}).neqs, oo_.sur.(fields{i}).dof);
    rep = rep.addParagraph(...
        'indent', false, ...
        'text', preamble);

    rep = rep.addSection('cols', 1);
    column_names = {'', 'Estimates','t-statistic','Std. Error'};
    rep = rep.addTable(...
        'precision', 5, ...
        'column_names', column_names);
    rep = rep.addData('data', {...
        oo_.sur.(fields{i}).pname, ...
        oo_.sur.(fields{i}).beta ...
        oo_.sur.(fields{i}).tstat ...
        oo_.sur.(fields{i}).stderr});
    
    rep = rep.addSection('cols', 1);
    afterward = [sprintf('$R^2$: %f\\newline ', oo_.sur.(fields{i}).R2), ...
        sprintf('$R^2$ Adjusted: %f\\newline ', oo_.sur.(fields{i}).adjR2), ...
        sprintf('$s^2$: %f\\newline ', oo_.sur.(fields{i}).s2), ...
        sprintf('Durbin-Watson: %f\\newline ', oo_.sur.(fields{i}).dw)];
    rep = rep.addParagraph(...
        'indent', false, ...
        'text', afterward);
end

%% Write & Compile Report
rep.write();
rep.compile();
end
