function options = set_default_initial_condition_decomposition_options(options)
%function options = set_default_initial_condition_decomposition_options(options)
% sets the default options for prior_shock_decomposition
%
% INPUTS
%    options
%
% OUTPUTS
%    options
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2017-2019 Dynare Team
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

options.initial_condition_decomp.colormap = '';
options.initial_condition_decomp.nodisplay = false;
options.initial_condition_decomp.graph_format = 'eps';
options.initial_condition_decomp.fig_name = '';
options.initial_condition_decomp.detail_plot = false;
options.initial_condition_decomp.init2shocks = [];
options.initial_condition_decomp.steadystate = false;
options.initial_condition_decomp.with_epilogue = options.shock_decomp.with_epilogue;
options.initial_condition_decomp.write_xls = false;
options.initial_condition_decomp.type = '';
options.initial_condition_decomp.plot_init_date = [];
options.initial_condition_decomp.plot_end_date = [];
options.initial_condition_decomp.diff = false;
options.initial_condition_decomp.flip = false;
options.initial_condition_decomp.max_nrows = 6;
end
