function collapse_figures_in_tabgroup

% Copyright © 2023 Eduard Benet Cerda
% Copyright © 2024 Dynare Team
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

% Create a new figure with results
fig = uifigure(Name = 'Dynare Results');

% Add a grid layout to make sure it spans the entire width
g = uigridlayout(fig, [1,1], Padding = 0);

% Add a tabgroup
tg = uitabgroup(g);

% Find all figures with Dynare Tag 
f = findobj('-regexp','tag','dynare-figure');

% Loop over all figures and reparent them to a tab. Avoid legends, they are
% automatically tied.
for j = 1 : numel(f)
    t = uitab(tg);
    types = arrayfun(@class, f(j).Children, 'UniformOutput', false);
    idx = ismember(types, 'matlab.graphics.illustration.Legend'); % no need to reparent legends
    set(f(j).Children(~idx),'Parent',t)
    t.Title = f(j).Name;
    delete(f(j))
end