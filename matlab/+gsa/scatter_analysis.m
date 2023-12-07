function scatter_analysis(lpmat, xdata, options_scatter, options_)
% scatter_analysis(lpmat, xdata, options_scatter, options_)
% Plot scatter plot analysis
%
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu
%

% Copyright © 2017 European Commission
% Copyright © 2017-2023 Dynare Team
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

param_names = options_scatter.param_names;

if options_.TeX
    if ~isfield(options_scatter,'param_names_tex')
        param_names_tex = options_scatter.param_names;
    else
        param_names_tex = options_scatter.param_names_tex;
    end
end
amcf_name = options_scatter.amcf_name;
amcf_title = options_scatter.amcf_title;
fname_ = options_scatter.fname_;
xparam1=[];
if isfield(options_scatter,'xparam1')
    xparam1=options_scatter.xparam1;
end
OutputDirectoryName = options_scatter.OutputDirectoryName;

if ~options_.nograph
    skipline()
    xx=[];
    if ~isempty(xparam1)
        xx=xparam1;
    end
    if options_.TeX
        gsa.scatter_plots(lpmat, xdata, param_names_tex, '.', [fname_, '_', amcf_name], OutputDirectoryName, amcf_title, xx, options_)
    else 
        gsa.scatter_plots(lpmat, xdata, param_names, '.', [fname_, '_', amcf_name], OutputDirectoryName, amcf_title, xx, options_)
    end
end
