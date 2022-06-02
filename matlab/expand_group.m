function expand_group(use_shock_groups,var_list_, ic, no_graph)
% function expand_group(use_shock_groups,var_list_, ic)
% Expands shocks contributions out of a group of shocks
%
% INPUTS
%    use_shock_groups: [char]  name of the group
%    var_list_:        [char]  list of variables
%    ic:               [int]   group # to expand
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2016-2018 Dynare Team
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

if nargin<4
    no_graph=0;
end
filename = get(gcf,'filename');
[filepath, name, ext]=fileparts(filename);
M = evalin('base','M_');
oo = evalin('base','oo_');
options = evalin('base','options_');
mydata=get(findobj(gcf,'tag',['group' int2str(ic)]),'userdata');
if isfield(mydata,'shock_decomp')
    options.shock_decomp=mydata.shock_decomp;
end
options.plot_shock_decomp=mydata.plot_shock_decomp;
options.first_obs=mydata.first_obs;
options.nobs=mydata.nobs;
% define expanded group
label = mydata.shock_group.label;
label = strrep(label,' ','_');
label = strrep(label,'-','_');
label = strrep(label,'(','');
label = strrep(label,')','');
label = strrep(label,'.','');
shocks = mydata.shock_group.shocks;
options.plot_shock_decomp.fig_name = [mydata.fig_name '. Expand'];
options.plot_shock_decomp.use_shock_groups = label; %[use_shock_groups_old int2str(ic)];
options.plot_shock_decomp.filepath = filepath;
for j=1:length(shocks)
    M.shock_groups.(options.plot_shock_decomp.use_shock_groups).(['group' int2str(j)]).label=shocks{j};
    M.shock_groups.(options.plot_shock_decomp.use_shock_groups).(['group' int2str(j)]).shocks=shocks(j);
end
M.exo_names = mydata.exo_names;

options.plot_shock_decomp.interactive=0;
options.plot_shock_decomp.expand=1;
options.plot_shock_decomp.nodisplay=0;
if no_graph
    options.no_graph.plot_shock_decomposition=1;
    options.plot_shock_decomp.write_xls=1;
else
    options.plot_shock_decomp.write_xls=0;
end
%% set optimal colormap
func = @(x) colorspace('RGB->Lab',x);
MAP = distinguishable_colors(length(shocks)+1,'w',func);
MAP(end,:) = [0.7 0.7 0.7];

options.plot_shock_decomp.colormap = MAP;

options.plot_shock_decomp.endo_names = mydata.endo_names;
options.plot_shock_decomp.endo_names_tex = mydata.endo_names_tex;

plot_shock_decomposition(M,oo,options,{var_list_});
