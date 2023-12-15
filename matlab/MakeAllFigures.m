function MakeAllFigures(NumberOfPlots,Caption,FigureProperties,Info)

% Copyright © 2005-2023 Dynare Team
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

global M_ options_

FigHandle = figure('Name',FigureProperties.Name);

NAMES = cell(NumberOfPlots,1);
if options_.TeX
    TeXNAMES = cell(NumberOfPlots,1);
end

if NumberOfPlots == 9
    nr = 3;
    nc = 3;
elseif NumberOfPlots == 8
    nr = 3;
    nc = 3;
elseif NumberOfPlots == 7
    nr = 3;
    nc = 3;
elseif NumberOfPlots == 6
    nr = 2;
    nc = 3;
elseif NumberOfPlots == 5
    nr = 3;
    nc = 2;
elseif NumberOfPlots == 4
    nr = 2;
    nc = 2;
elseif NumberOfPlots == 3
    nr = 2;
    nc = 2;
elseif NumberOfPlots == 2
    nr = 1;
    nc = 2;
elseif NumberOfPlots == 1
    nr = 1;
    nc = 1;
end

for plt = 1:NumberOfPlots
    field1 = sprintf('Box%u', plt);
    NumberOfCurves = Info.(field1).Number;
    NumberOfObservations = zeros(2,1);
    x = cell(NumberOfCurves,1);
    y = cell(NumberOfCurves,1);
    PltType = cell(NumberofCurves,1);
    top = NaN(NumberOfCurves,1);
    bottom = NaN(NumberOfCurves,1);
    binf = NaN(NumberOfCurves,1);
    bsup = NaN(NumberOfCurves,1);
    for curve = 1:NumberOfCurves
        field2 = sprintf('Curve%u', curve);
        x{curve} = Info.(field1).(field2).xdata;
        y{curve} = Info.(field1).(field2).ydata;
        name = Info.(field1).(field2).variablename;
        PltType{curve} = Info.(field1).(field2).type;
        if length(x{curve})-length(y{curve})
            disp('MakeFigure :: The number of observations in x doesn''t match with ')
            disp(['the number of observation in y for ' name ])
            return
        end
        if Info.PlotProperties.CutTop
            top(curve) = max(y{curve});
        else Info.PlotProperties.CutBottom
            bottom(curve) = min(y{curve});
        end
        binf(curve) = min(x{curve});
        bsup(curve) = max(x{curve});
    end
    ymax = max(top);
    ymin = min(bottom);
    xmin = min(binf);
    xmax = max(bsup);
    if isnan(ymin(plt))
        ymin = 0;
    end
    NAMES{plt} = Info.(field1).name;
    if options_.TeX
        TeXNAMES{plt} = Info.(field1).texname;
    end
    subplot(nr,nc,plt)
    hold on
    for curve = 1:NumberOfCurves
        hh_fig = plot(x{curve},y{curve});
        if strcmpi(PltType{curve},'PriorDensity')
            set(hh_fig,'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',2)
            %
            %
        elseif strcmpi(PltType{curve},'DensityEstimate')
            set(hh_fig,'Color','k','LineStyle','-','LineWidth',2)
            %
            %
        elseif strcmpi(PltType{curve},'ModeEstimate')
            set(hh_fig,'Color','g','LineStyle','--','LineWidth',2)
            %
            %
        elseif strcmpi(PltType{curve},'SmoothVariable')
            set(hh_fig,'Color','k','LineStyle','-','LineWidth',2)
            %
            %
        elseif strcmpi(PltType{curve},'Deciles')
            set(hh_fig,'Color','g','LineStyle','-','LineWidth',1)
            %
            %
        elseif strcmpi(PltType{curve},'Forecasts')
            set(hh_fig,'Color','','LineStyle','-','LineWidth',2)
            %
            %
        elseif strcmpi(PltType{curve},'ForecastsHPD')
            set(hh_fig,'Color','k','LineStyle','-','LineWidth',1)
            %
            %
        elseif strcmpi(PltType{curve},'ForecastsDeciles')
            set(hh_fig,'Color','g','LineStyle','-','LineWidth',1)
            %
            %
        elseif strcmpi(PltType{curve},'DiagnosticWithin')
            set(hh_fig,'Color','b','LineStyle','-','LineWidth',2)
            %
            %
        elseif strcmpi(PltType{curve},'DiagnosticPooled')
            set(hh_fig,'Color','r','LineStyle','-','LineWidth',2)
            %
            %
        end
    end
    axis([xmin xmax ymin ymax])
    title(NAMES{plt})
    drawnow
    hold off
end

if Info.SaveFormat.Eps
    if isempty(Info.SaveFormat.Name)
        print(sprintf('%s%s%u', M_.fname, Info.SaveFormat.GenericName, Info.SaveFormat.Number), '-depsc2')
    else
        print(sprintf('%s%s%s', M_.fname, Info.SaveFormat.GenericName, Info.SaveFormat.Name), '-depsc2')
    end
end

if Info.SaveFormat.Pdf && ~isoctave
    if isempty(Info.SaveFormat.Name)
        print(sprintf('%s%s%u', M_.fname, Info.SaveFormat.GenericName, Info.SaveFormat.Number), '-dpdf')
    else
        print(sprintf('%s%s%s', M_.fname, Info.SaveFormat.GenericName, Info.SaveFormat.Name), '-dpdf')
    end
end

if Info.SaveFormat.Fig && ~isoctave
    if isempty(Info.SaveFormat.Name)
        saveas(FigHandle, sprintf('%s%s%u.fig', M_.fname, Info.SaveFormat.GenericName, Info.SaveFormat.Number));
    else
        saveas(FigHandle, sprintf('%s%s%s.fig', M_.fname, Info.SaveFormat.GenericName, Info.SaveFormat.Name));
    end
end
