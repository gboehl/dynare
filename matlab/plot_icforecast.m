function plot_icforecast(Variables,periods,options_,oo_)
% Build plots for the conditional forecasts.
%
% INPUTS
%  o Variables     [cell]        Names of the endogenous variables to be plotted.
%  o periods       [int]         Number of periods to be plotted.
%  o options_      [structure]   Options.
%  o oo_           [structure]   Storage of results.
%
% OUTPUTS
%  None.
%
% SPECIAL REQUIREMENTS
%  This routine has to be called after imcforecast.m.

% Copyright © 2006-2019 Dynare Team
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

if ~isfield(oo_, 'conditional_forecast')
    error('Can''t find conditional forecasts');
else
    forecasts = oo_.conditional_forecast;
end

forecast_periods = length(forecasts.cond.Mean.(Variables{1}));
if nargin==1 || isempty(periods) % Set default number of periods.
    periods=forecast_periods;
else
    periods=periods+1; %take care of initial period that is not forecasted
end

if periods>forecast_periods
    fprintf('\nplot_icforecast:: Number of periods for plotting exceeds forecast horizon. Setting periods to %d.\n',forecast_periods-1)
    periods=forecast_periods;
end

forecasts.graph.OutputDirectoryName = CheckPath('graphs',forecasts.graph.fname);

for i=1:length(Variables)
    ci1 = forecasts.cond.ci.(Variables{i});
    m1 = forecasts.cond.Mean.(Variables{i});
    ci2 = forecasts.uncond.ci.(Variables{i});
    m2 = forecasts.uncond.Mean.(Variables{i});
    build_figure(Variables{i},ci1(:,1:periods),ci2(:,1:periods),m1(1:periods),m2(1:periods),options_,forecasts.graph);
end

function build_figure(name,cci1,cci2,mm1,mm2,options_,graphoptions)
hh = dyn_figure(options_.nodisplay,'Name',['Conditional forecast (' graphoptions.title ,'): ' name '.']);
H = length(mm1);
h1 = area(1:H,cci1(2,1:H),'BaseValue',min([min(cci1(1,:)),min(cci2(1,:))]),'FaceColor',[.9 .9 .9]);
hold on
h2 = area(1:H,cci1(1,1:H),'BaseValue',min([min(cci1(1,:)),min(cci2(1,:))]),'FaceColor',[1 1 1]);
plot(1:H,mm1,'-k','linewidth',3)
plot(1:H,mm2,'--k','linewidth',3)
plot(1:H,cci2(1,:),'--k','linewidth',1)
plot(1:H,cci2(2,:),'--k','linewidth',1)
axis tight
hold off
dyn_saveas(hh,[graphoptions.OutputDirectoryName '/Conditional_forecast_',strrep(deblank(graphoptions.title),' ','_'),'_',name],options_.nodisplay,options_.graph_format)