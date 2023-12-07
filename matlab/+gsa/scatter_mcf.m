function  scatter_mcf(X,Y,vnames,plotsymbol, fnam, dirname, figtitle, xparam1, options_, beha_name, non_beha_name, beha_name_latex, non_beha_name_latex)
% scatter_mcf(X,Y,vnames,plotsymbol, fnam, dirname, figtitle, xparam1, options_, beha_name, non_beha_name, beha_name_latex, non_beha_name_latex)
%
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu

% Copyright © 2014-2016 European Commission
% Copyright © 2014-2023 Dynare Team
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

% PURPOSE: Pairwise scatter plots of the columns of x and y after
% Monte Carlo filtering
%---------------------------------------------------
% USAGE:    scatter_mcf(x,y,vnames,pltsym,diagon)
%        or scatter_mcf(x,y) which relies on defaults
% where:
%        x = an nxk matrix with columns containing behavioural sample
%        y = an mxk matrix with columns containing non-behavioural sample
%   vnames = a vector of variable names
%            (default = numeric labels 1,2,3 etc.)
%   pltsym = a plt symbol
%            (default = '.' for npts > 100, 'o' for npts < 100


Z=[X;Y];
[n,p] = size(Z);
clear Z;

nflag = 0;
if nargin >=3
    nflag = 1;
end

if nargin<4 || isempty(plotsymbol)
    if n*p<100
        plotsymbol = 'o';
    else 
        plotsymbol = '.';
    end
end

if nargin<5
    fnam='';
end
if nargin<6
    dirname='';
    nograph=1;
else
    nograph=0;
end
if nargin<7
    figtitle=fnam;
end
if nargin<8
    xparam1=[];
end
if nargin<10
    beha_name = 'BEHAVIOUR';
    non_beha_name = 'NON-BEHAVIOUR';
end
if nargin==10
    non_beha_name = ['NON-' beha_name];
end

figtitle_tex=strrep(figtitle,'_','\_');

fig_nam_=fnam;
if ~nograph
    hh_fig=dyn_figure(options_.nodisplay,'name',figtitle);
end

bf = 0.1;
ffs = 0.05/(p-1);
ffl = (1-2*bf-0.05)/p;
if p>1
    fL = linspace(bf,1-bf+ffs,p+1);
else
    fL = bf;
end
for i = 1:p
    for j = 1:p
        h = axes('position',[fL(i),fL(p+1-j),ffl,ffl]);
        if i==j
            h1=gsa.cumplot(X(:,j));
            set(h1,'color',[0 0 1],'LineWidth',1.5)
            hold on,
            h2=gsa.cumplot(Y(:,j));
            set(h2,'color',[1 0 0],'LineWidth',1.5)
            if ~isempty(xparam1)
                hold on, plot(xparam1([j j]),[0 1],'k--')
            end
            if j<p
                set(gca,'XTickLabel',[],'XTick',[]);
            else
                grid off
            end
            set(gca,'YTickLabel',[],'YTick',[]);
        else
            if j>i
                plot(X(:,i),X(:,j),[plotsymbol,'b'])
                hold on,
                plot(Y(:,i),Y(:,j),[plotsymbol,'r'])
            else
                plot(Y(:,i),Y(:,j),[plotsymbol,'r'])
                hold on,
                plot(X(:,i),X(:,j),[plotsymbol,'b'])
            end
            if ~isempty(xparam1)
                hold on 
                plot(xparam1(i),xparam1(j),'s','MarkerFaceColor',[0 0.75 0],'MarkerEdgeColor',[0 0.75 0])
            end
            hold off;
            if i<p
                set(gca,'YTickLabel',[],'YTick',[]);
            else
                set(gca,'yaxislocation','right');
            end
            if j<p
                set(gca,'XTickLabel',[],'XTick',[]);
            end
        end
        if nflag == 1
            set(gca,'fontsize',9);
        end
        if i==1
            if nflag == 1
                if options_.TeX
                    ylabel(vnames(j,:),'Rotation',45, ...
                        'HorizontalAlignment','right','VerticalAlignment','middle','Interpreter','latex');
                else
                    ylabel(vnames(j,:),'Rotation',45, ...
                        'HorizontalAlignment','right','VerticalAlignment','middle','Interpreter','none');
                end
            else
                ylabel([num2str(j),' '],'Rotation',90)
            end
        end
        if j==1
            if nflag == 1
                if options_.TeX
                    title(vnames(i,:),'Rotation',45, ...
                      'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex')
                else
                    title(vnames(i,:),'Rotation',45, ...
                      'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','none')
                end
            else
                title(num2str(i))
            end
        end
        drawnow
    end
end
if ~isoctave
    if options_.TeX
        annotation('textbox', [0.1,0,0.35,0.05],'String', beha_name_latex,'Color','Blue','horizontalalignment','center','interpreter','latex');
        annotation('textbox', [0.55,0,0.35,0.05],'String', non_beha_name_latex,'Color','Red','horizontalalignment','center','interpreter','latex');
    else
        annotation('textbox', [0.1,0,0.35,0.05],'String', beha_name,'Color','Blue','horizontalalignment','center','interpreter','none');
        annotation('textbox', [0.55,0,0.35,0.05],'String', non_beha_name,'Color','Red','horizontalalignment','center','interpreter','none');
    end
end

if ~nograph
    dyn_saveas(hh_fig,[dirname,filesep,fig_nam_],options_.nodisplay,options_.graph_format);
    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fidTeX = fopen([dirname,'/',fig_nam_ '.tex'],'w');
        fprintf(fidTeX,'%% TeX eps-loader file generated by scatter_mcf.m (Dynare).\n');
        fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s}\n',strrep([dirname,'/',fig_nam_],'\','/'));
        fprintf(fidTeX,'\\caption{%s.}',figtitle_tex);
        fprintf(fidTeX,'\\label{Fig:%s}\n',fig_nam_);
        fprintf(fidTeX,'\\end{figure}\n\n');
        fprintf(fidTeX,'%% End Of TeX file. \n');
        fclose(fidTeX);
    end
end