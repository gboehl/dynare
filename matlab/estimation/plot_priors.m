function plot_priors(bayestopt_,M_,estim_params_,options_,optional_title)
% plot_priors(bayestopt_,M_,estim_params_,options_,optional_title)
% plots prior density
%
% INPUTS
%    o bayestopt_       [structure]
%    o M_               [structure]
%    o estim_params_    [structure]
%    o options_         [structure]
%    o optional_title   [string]

% OUTPUTS
%    None
%
% SPECIAL REQUIREMENTS
%    None

% Copyright © 2004-2023 Dynare Team
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
latexDirectoryName = CheckPath('latex',M_.dname);
graphDirectoryName = CheckPath('graphs',M_.dname);

TeX = options_.TeX;
if nargin<5
    figurename = 'Priors';
else
    figurename = optional_title;
end
npar = length(bayestopt_.p1);
[nbplt,nr,nc,~,~,nstar] = pltorg(npar);

if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fidTeX = fopen([latexDirectoryName filesep M_.fname '_Priors.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by plot_priors.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
end
for plt = 1:nbplt
    hh_fig = dyn_figure(options_.nodisplay,'Name',figurename);
    if TeX
        TeXNAMES = [];
        NAMES    = [];
    end
    nstar0 = min(nstar,npar-(plt-1)*nstar);
    for index=1:nstar0
        names = [];
        i = (plt-1)*nstar + index;
        [x,f] = draw_prior_density(i,bayestopt_);
        [nam,texnam] = get_the_name(i,TeX,M_,estim_params_,options_.varobs);
        subplot(nr,nc,index)
        hh_plt = plot(x,f,'-k','linewidth',2);
        set(hh_plt,'color',[0.7 0.7 0.7]);
        box on
        if TeX
            title(texnam,'Interpreter','latex')
        else
            title(nam,'Interpreter','none')
        end
        drawnow
    end
    dyn_saveas(hh_fig,[graphDirectoryName filesep M_.fname '_Priors' int2str(plt)],options_.nodisplay,options_.graph_format);
    if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        fprintf(fidTeX,'\\centering\n');
        fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_Priors%s}\n',options_.figures.textwidth*min(index/nc,1),[graphDirectoryName '/' M_.fname],int2str(plt));% don't use filesep as it will create issues with LaTeX on Windows
        fprintf(fidTeX,'\\caption{Priors.}');
        fprintf(fidTeX,'\\label{Fig:Priors:%s}\n',int2str(plt));
        fprintf(fidTeX,'\\end{figure}\n');
    end
end
if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,'%% End of TeX file.\n');
    fclose(fidTeX);
end