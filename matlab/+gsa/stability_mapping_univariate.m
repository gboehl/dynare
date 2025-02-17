function [proba, dproba] = stability_mapping_univariate(lpmat, ibehaviour, inonbehaviour, aname, fname_, options_, parnames, estim_params_, iplot, ipar, dirname, pcrit, atitle)
% [proba, dproba] = stability_mapping_univariate(lpmat, ibehaviour, inonbehaviour, aname, fname_, options_, parnames, estim_params_, iplot, ipar, dirname, pcrit, atitle)
% Inputs:
% - lpmat               [double]        Monte Carlo matrix
% - ibehaviour          [integer]       index of behavioural runs
% - inonbehaviour       [integer]       index of non-behavioural runs
% - aname               [string]        label of the analysis
% - fname_              [string]        file name
% - options_            [structure]     options structure
% - parnames            [char]          parameter name vector
% - estim_params_       [structure]     characterizing parameters to be estimated
% - iplot               [boolean]       1 plot cumulative distributions (default)
%                                       0 no plots
% - ipar                [integer]       index array of parameters to plot
% - dirname             [string]        (OPTIONAL) path of the directory where to save
%                                       (default: current directory)
% - pcrit               [double]        (OPTIONAL) critical value of the pvalue below which show the plots
%
% Plots: dotted lines for BEHAVIOURAL
%        solid lines for NON BEHAVIOURAL
% USES gsa.smirnov_test.m
%
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu

% Copyright © 2012-2016 European Commission
% Copyright © 2012-2023 Dynare Team
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

if nargin<9
    iplot=1;
end
if nargin<11
    dirname='';
end
if nargin<13
    atitle=aname;
end

nshock = estim_params_.nvx;
nshock = nshock + estim_params_.nvn;
nshock = nshock + estim_params_.ncx;
nshock = nshock + estim_params_.ncn;

npar=size(lpmat,2);
ishock= npar>estim_params_.np;

if nargin<10
    ipar=[];
end
if nargin<12 || isempty(pcrit)
    pcrit=1;
end

% Smirnov test for Blanchard
proba=NaN(npar,1);
dproba=NaN(npar,1);
for j=1:npar
    [~,P,KSSTAT] = gsa.smirnov_test(lpmat(ibehaviour,j),lpmat(inonbehaviour,j));
    proba(j)=P;
    dproba(j)=KSSTAT;
end
if isempty(ipar)
    ipar=find(proba<pcrit);
end
nparplot=length(ipar);
if iplot && ~options_.nograph
    lpmat=lpmat(:,ipar);
    ftit=parnames(ipar+nshock*(1-ishock));

    for i=1:ceil(nparplot/12)
        hh_fig=dyn_figure(options_.nodisplay,'name',atitle);
        for j=1+12*(i-1):min(nparplot,12*i)
            subplot(3,4,j-12*(i-1))
            if ~isempty(ibehaviour)
                h=gsa.cumplot(lpmat(ibehaviour,j));
                set(h,'color',[0 0 1], 'linestyle',':','LineWidth',1.5)
            end
            hold on
            if ~isempty(inonbehaviour)
                h=gsa.cumplot(lpmat(inonbehaviour,j));
                set(h,'color',[0 0 0],'LineWidth',1.5)
            end
            title([ftit{j},'. p-value ', num2str(proba(ipar(j)),2)],'interpreter','none')
        end
        if nparplot>12
            dyn_saveas(hh_fig,[dirname,filesep,fname_,'_',aname,'_SA_',int2str(i)],options_.nodisplay,options_.graph_format);
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                fidTeX = fopen([dirname,filesep,fname_,'_',aname,'_SA_',int2str(i) '.tex'],'w');
                fprintf(fidTeX,'%% TeX eps-loader file generated by gsa.stability_mapping_univariate.m (Dynare).\n');
                fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s}\n',strrep([dirname,filesep,fname_,'_',aname,'_SA_',int2str(i)],'\','/'));
                fprintf(fidTeX,'\\caption{%s.}',atitle);
                fprintf(fidTeX,'\\label{Fig:%s:%u}\n',atitle,i);
                fprintf(fidTeX,'\\end{figure}\n\n');
                fprintf(fidTeX,'%% End Of TeX file. \n');
                fclose(fidTeX);
            end
        else
            dyn_saveas(hh_fig,[dirname,filesep,fname_,'_',aname,'_SA'],options_.nodisplay,options_.graph_format);
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                fidTeX = fopen([dirname,filesep,fname_,'_',aname,'_SA.tex'],'w');
                fprintf(fidTeX,'%% TeX eps-loader file generated by gsa.stability_mapping_univariate.m (Dynare).\n');
                fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s}\n',options_.figures.textwidth*min((j-12*(i-1))/3,1),strrep([dirname,filesep,fname_,'_',aname,'_SA'],'\','/'));
                fprintf(fidTeX,'\\caption{%s.}',atitle);
                fprintf(fidTeX,'\\label{Fig:%s}\n',atitle);
                fprintf(fidTeX,'\\end{figure}\n\n');
                fprintf(fidTeX,'%% End Of TeX file. \n');
                fclose(fidTeX);
            end
        end
    end
end
