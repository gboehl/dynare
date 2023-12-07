function indcorr = stab_map_2(x,alpha2, pvalue_crit, M_,options_,bayestopt_,estim_params_, fnam, dirname,xparam1,figtitle)
% indcorr = stab_map_2(x,alpha2, pvalue_crit, M_,options_,bayestopt_,estim_params_, fnam, dirname,xparam1,figtitle)
% Inputs:
%  - x
%  - alpha2
%  - pvalue_crit
%  - M_                     [structure]     Matlab's structure describing the model
%  - options_               [structure]     Matlab's structure describing the current options
%  - bayestopt_             [structure]     describing the priors
%  - estim_params_          [structure]     characterizing parameters to be estimated
%  - fnam                   [string]        file name
%  - dirnam                 [string]        directory name
%  - xparam1                [double]        parameter vector
%  - figtitle               [string]        figure title
%
% Output:
%  - indcorr
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu

% Copyright © 2011-2016 European Commission
% Copyright © 2011-2023 Dynare Team
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

npar=size(x,2);
ishock= npar>estim_params_.np;
nograph = options_.nograph;
if nargin<8
    fnam='';
end
if nargin<9
    dirname='';
    nograph=1;
end
if nargin<10
    xparam1=[];
end
if nargin<11
    figtitle=fnam;
end

fname_ = M_.fname;
nshock = estim_params_.nvx;
nshock = nshock + estim_params_.nvn;
nshock = nshock + estim_params_.ncx;
nshock = nshock + estim_params_.ncn;

[c0, pvalue] = corrcoef(x);
c00=tril(c0,-1);
fig_nam_=[fname_,'_',fnam,'_corr_'];
fig_nam_tex_table=strrep([fnam,'_corr'],' ','_');
fig_nam_=strrep(fig_nam_,' ','_');

ifig=0;
j2=0;
if ishock==0
    npar=estim_params_.np;
    if ~isempty(xparam1)
        xparam1=xparam1(nshock+1:end);
    end
else
    npar=estim_params_.np+nshock;
end
title_string=['Correlation analysis for ',fnam];
title_string_tex=['Correlation analysis for ',strrep(fnam,'_','\\_')];

indcorr = [];
entry_iter=1;
for j=1:npar
    i2=find(abs(c00(:,j))>alpha2);
    if ~isempty(i2)
        for jx=1:length(i2)
            if pvalue(j,i2(jx))<pvalue_crit
                indcorr = [indcorr; [j i2(jx)]];
                j2=j2+1;
                if ishock
                    if options_.TeX
                        [param_name_temp1, param_name_tex_temp1]= get_the_name(j,options_.TeX,M_,estim_params_,options_);
                        [param_name_temp2, param_name_tex_temp2]= get_the_name(i2(jx),options_.TeX,M_,estim_params_,options_);
                        tmp_name=(['[',param_name_temp1,',',param_name_temp2,']']);
                        tmp_name_tex=(['[',param_name_tex_temp1,',',param_name_tex_temp2,']']);
                        name{entry_iter,1}=tmp_name;
                        name_tex{entry_iter,1}=tmp_name_tex;
                    else
                        [param_name_temp1]= get_the_name(j,options_.TeX,M_,estim_params_,options_);
                        [param_name_temp2]= get_the_name(i2(jx),options_.TeX,M_,estim_params_,options_);
                        tmp_name=(['[',param_name_temp1,',',param_name_temp2,']']);
                        name{entry_iter,1}=tmp_name;
                    end
                else
                    if options_.TeX
                        [param_name_temp1, param_name_tex_temp1]= get_the_name(j+nshock,options_.TeX,M_,estim_params_,options_);
                        [param_name_temp2, param_name_tex_temp2]= get_the_name(i2(jx)+nshock,options_.TeX,M_,estim_params_,options_);
                        tmp_name=(['[',param_name_temp1,',',param_name_temp2,']']);
                        tmp_name_tex=(['[',param_name_tex_temp1,',',param_name_tex_temp2,']']);
                        name{entry_iter,1}=tmp_name;
                        name_tex{entry_iter,1}=tmp_name_tex;
                    else
                        [param_name_temp1]= get_the_name(j+nshock,options_.TeX,M_,estim_params_,options_);
                        [param_name_temp2]= get_the_name(i2(jx)+nshock,options_.TeX,M_,estim_params_,options_);
                        tmp_name=(['[',param_name_temp1,',',param_name_temp2,']']);
                        name{entry_iter,1}=tmp_name;
                    end
                end
                data_mat(entry_iter,1)=c0(i2(jx),j);
                entry_iter=entry_iter+1;

                if ~nograph
                    if mod(j2,12)==1
                        ifig=ifig+1;
                        hh_fig=dyn_figure(options_.nodisplay,'name',[figtitle,' sample bivariate projection ', num2str(ifig)]);
                    end
                    subplot(3,4,j2-(ifig-1)*12)
                    plot(x(:,j),x(:,i2(jx)),'.')
                    if ~isempty(xparam1)
                        hold on, plot(xparam1(j),xparam1(i2(jx)),'ro')
                    end
                    if ishock
                        xlabel(bayestopt_.name{j},'interpreter','none'),
                        ylabel(bayestopt_.name{i2(jx)},'interpreter','none'),
                    else
                        xlabel(bayestopt_.name{j+nshock},'interpreter','none'),
                        ylabel(bayestopt_.name{i2(jx)+nshock},'interpreter','none'),
                    end
                    title(['cc = ',num2str(c0(i2(jx),j))])
                    if (mod(j2,12)==0) && j2>0
                        dyn_saveas(hh_fig,[dirname,filesep,fig_nam_,int2str(ifig)],options_.nodisplay,options_.graph_format);
                        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                            fidTeX = fopen([dirname,filesep,fig_nam_,int2str(ifig),'.tex'],'w');
                            fprintf(fidTeX,'%% TeX eps-loader file generated by stab_map_2.m (Dynare).\n');
                            fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
                            fprintf(fidTeX,'\\begin{figure}[H]\n');
                            fprintf(fidTeX,'\\centering \n');
                            fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s}\n',strrep([dirname,'/',fig_nam_,int2str(ifig)],'\','/'));
                            fprintf(fidTeX,'\\caption{%s.}',[figtitle,' sample bivariate projection ', num2str(ifig)]);
                            fprintf(fidTeX,'\\label{Fig:%s:%u}\n',fig_nam_,ifig);
                            fprintf(fidTeX,'\\end{figure}\n\n');
                            fprintf(fidTeX,'%% End Of TeX file. \n');
                            fclose(fidTeX);
                        end
                    end
                end
            end

        end
    end
    if ~nograph && (j==(npar)) && j2>0 && (mod(j2,12)~=0)
        dyn_saveas(hh_fig,[dirname,filesep,fig_nam_,int2str(ifig)],options_.nodisplay,options_.graph_format);
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fidTeX = fopen([dirname,filesep,fig_nam_,int2str(ifig),'.tex'],'w');
            fprintf(fidTeX,'%% TeX eps-loader file generated by stab_map_2.m (Dynare).\n');
            fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s}\n',options_.figures.textwidth*min((j2-(ifig-1)*12)/3,1),strrep([dirname,'/',fig_nam_,int2str(ifig)],'\','/'));
            fprintf(fidTeX,'\\caption{%s.}',[figtitle,' sample bivariate projection ', num2str(ifig)]);
            fprintf(fidTeX,'\\label{Fig:%s:%u}\n',fig_nam_,ifig);
            fprintf(fidTeX,'\\end{figure}\n\n');
            fprintf(fidTeX,'%% End Of TeX file. \n');
            fclose(fidTeX);
        end
    end
end

if j2==0
    skipline();
    disp(['No correlation term with pvalue <', num2str(pvalue_crit),' and |corr. coef.| >',num2str(alpha2),' found for ',fnam])
else
    headers={'Parameters'; 'corrcoef'};
    if ~options_.noprint
        skipline();
    end
    dyntable(options_,title_string,headers, name, data_mat, 0, 7, 3);
    if options_.TeX
        dyn_latex_table(M_, options_, title_string_tex, fig_nam_tex_table, headers, name_tex, data_mat, 0, 7, 3);
    end
end