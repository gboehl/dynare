function plot(M_, params, idemoments, idehess, idemodel, idelre, advanced, tittxt, name, IdentifDirectoryName, fname, options_, estim_params_, bayestopt_, tit_TeX, name_tex)
% plot(M_, params,idemoments,idehess,idemodel, idelre, advanced, tittxt, name, IdentifDirectoryName, fname, options_, estim_params_, bayestopt_, tit_TeX, name_tex)
%
% INPUTS
%    o M_                   [structure] model
%    o params               [array] parameter values for identification checks
%    o idemoments           [structure] identification results for the moments
%    o idehess              [structure] identification results for the Hessian
%    o idemodel             [structure] identification results for the reduced form solution
%    o idelre               [structure] identification results for the LRE model
%    o advanced             [integer] flag for advanced identification checks
%    o tittxt               [char] name of the results to plot
%    o name                 [char] list of parameter names
%    o IdentifDirectoryName [char] directory name
%    o fname                [char] file name
%    o options_             [structure] structure describing the current options
%    o estim_params_        [structure] characterizing parameters to be estimated
%    o bayestopt_           [structure] describing the priors
%    o tittxt               [char] TeX-name of the results to plot
%    o name_tex             [char] TeX-names of the parameters
%
% OUTPUTS
%    None
%
% SPECIAL REQUIREMENTS
%    None

% Copyright © 2008-2023 Dynare Team
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

if nargin <14 || isempty(tit_TeX)
    tit_TeX=tittxt;
end

if nargin <15
    name_tex=name;
end

[SampleSize, nparam]=size(params);
si_dMOMENTSnorm = idemoments.si_dMOMENTSnorm;
si_dTAUnorm = idemodel.si_dREDUCEDFORMnorm;
si_dLREnorm = idelre.si_dDYNAMICnorm;

tittxt1=regexprep(tittxt, ' ', '_');
tittxt1=strrep(tittxt1, '.', '');
if SampleSize == 1
    hh_fig = dyn_figure(options_.nodisplay,'Name',[tittxt, ' - Identification using info from observables']);
    subplot(211)
    mmm = (idehess.ide_strength_dMOMENTS);
    [~, is] = sort(mmm);
    if ~all(isnan(idehess.ide_strength_dMOMENTS_prior)) ...
       && ~(nparam == 1 && ~isoctave && matlab_ver_less_than('9.7')) % MATLAB < R2019b does not accept bar(1, [2 3])
        bar(1:nparam,log([idehess.ide_strength_dMOMENTS(:,is)' idehess.ide_strength_dMOMENTS_prior(:,is)']))
    else
        bar(1:nparam,log(idehess.ide_strength_dMOMENTS(:,is)'))
    end
    hold on
    plot((1:length(idehess.ide_strength_dMOMENTS(:,is)))-0.15,log(idehess.ide_strength_dMOMENTS(:,is)'),'o','MarkerSize',7,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none')
    plot((1:length(idehess.ide_strength_dMOMENTS_prior(:,is)))+0.15,log(idehess.ide_strength_dMOMENTS_prior(:,is)'),'o','MarkerSize',7,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none')
    if any(isinf(log(idehess.ide_strength_dMOMENTS(idehess.identified_parameter_indices))))
        %-Inf, i.e. 0 strength
        inf_indices=find(isinf(log(idehess.ide_strength_dMOMENTS(idehess.identified_parameter_indices))) & log(idehess.ide_strength_dMOMENTS(idehess.identified_parameter_indices))<0);
        inf_pos=ismember(is,idehess.identified_parameter_indices(inf_indices));
        plot(find(inf_pos)-0.15,zeros(sum(inf_pos),1),'o','MarkerSize',7,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0])
        %+Inf, i.e. Inf strength
        inf_indices=find(isinf(log(idehess.ide_strength_dMOMENTS(idehess.identified_parameter_indices))) & log(idehess.ide_strength_dMOMENTS(idehess.identified_parameter_indices))>0);
        inf_pos=ismember(is,idehess.identified_parameter_indices(inf_indices));
        plot(find(inf_pos)-0.15,zeros(sum(inf_pos),1),'o','MarkerSize',7,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
    end
    if any(isinf(log(idehess.ide_strength_dMOMENTS_prior(idehess.identified_parameter_indices))))
        %-Inf, i.e. 0 strength
        inf_indices=find(isinf(log(idehess.ide_strength_dMOMENTS_prior(idehess.identified_parameter_indices))) & log(idehess.ide_strength_dMOMENTS_prior(idehess.identified_parameter_indices))<0);
        inf_pos=ismember(is,idehess.identified_parameter_indices(inf_indices));
        plot(find(inf_pos)+0.15,zeros(sum(inf_pos),1),'o','MarkerSize',7,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0])
        %+Inf, i.e. 0 strength
        inf_indices=find(isinf(log(idehess.ide_strength_dMOMENTS_prior(idehess.identified_parameter_indices))) & log(idehess.ide_strength_dMOMENTS_prior(idehess.identified_parameter_indices))>0);
        inf_pos=ismember(is,idehess.identified_parameter_indices(inf_indices));
        plot(find(inf_pos)+0.15,zeros(sum(inf_pos),1),'o','MarkerSize',7,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
    end
    set(gca,'xlim',[0 nparam+1])
    set(gca,'xticklabel','')
    dy = get(gca,'ylim');
    for ip=1:nparam
        if options_.TeX
            text(ip,dy(1),name_tex{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','latex')
        else
            text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
        end
    end
    if ~all(isnan(idehess.ide_strength_dMOMENTS_prior))
        legend('relative to param value','relative to prior std','Location','Best')
    else
        legend('relative to param value','Location','Best')
    end
    if  idehess.flag_score
        title('Identification strength with asymptotic Information matrix (log-scale)')
    else
        title('Identification strength with moments Information matrix (log-scale)')
    end

    subplot(212)
    if ~all(isnan(idehess.deltaM_prior)) ...
       && ~(nparam == 1 && ~isoctave && matlab_ver_less_than('9.7')) % MATLAB < R2019b does not accept bar(1, [2 3])
        bar(1:nparam, log([idehess.deltaM(is) idehess.deltaM_prior(is)]))
    else
        bar(1:nparam, log([idehess.deltaM(is)]))
    end
    hold on
    plot((1:length(idehess.deltaM(is)))-0.15,log(idehess.deltaM(is)'),'o','MarkerSize',7,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none')
    plot((1:length(idehess.deltaM_prior(is)))+0.15,log(idehess.deltaM_prior(is)'),'o','MarkerSize',7,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none')
    inf_pos=find(isinf(log(idehess.deltaM)));
    if ~isempty(inf_pos)
        inf_indices=~ismember(inf_pos,idehess.sensitivity_zero_pos);
        inf_pos=ismember(is,inf_pos(inf_indices));
        plot(find(inf_pos)-0.15,zeros(sum(inf_pos),1),'o','MarkerSize',7,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
    end
    inf_pos=find(isinf(log(idehess.deltaM_prior)));
    if ~isempty(inf_pos)
        inf_indices=~ismember(inf_pos,idehess.sensitivity_zero_pos);
        inf_pos=ismember(is,inf_pos(inf_indices));
        plot(find(inf_pos)+0.15,zeros(sum(inf_pos),1),'o','MarkerSize',7,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
    end
    set(gca,'xlim',[0 nparam+1])
    set(gca,'xticklabel','')
    dy = get(gca,'ylim');
    for ip=1:nparam
        if options_.TeX
            text(ip,dy(1),name_tex{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','latex')
        else
            text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
        end
    end
    if ~all(isnan(idehess.deltaM_prior))
        legend('relative to param value','relative to prior std','Location','Best')
    else
        legend('relative to param value','Location','Best')
    end
    if  idehess.flag_score
        title('Sensitivity component with asymptotic Information matrix (log-scale)')
    else
        title('Sensitivity component with moments Information matrix (log-scale)')
    end
    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fidTeX = fopen([IdentifDirectoryName '/' fname '_ident_strength_' tittxt1,'.tex'],'w');
        fprintf(fidTeX,'%% TeX eps-loader file generated by identification.plot.m (Dynare).\n');
        fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s_ident_strength_%s}\n',[IdentifDirectoryName '/' fname],tittxt1);
        fprintf(fidTeX,'\\caption{%s  - Identification using info from observables.}',tit_TeX);
        fprintf(fidTeX,'\\label{Fig:ident:%s}\n',deblank(tittxt));
        fprintf(fidTeX,'\\end{figure}\n\n');
        fprintf(fidTeX,'%% End Of TeX file. \n');
        fclose(fidTeX);
    end
    dyn_saveas(hh_fig,[IdentifDirectoryName '/' fname '_ident_strength_' tittxt1],options_.nodisplay,options_.graph_format);

    if advanced
        if ~options_.nodisplay
            skipline()
            disp('Plotting advanced diagnostics')
        end
        if all(isnan([si_dMOMENTSnorm';si_dTAUnorm';si_dLREnorm']))
            fprintf('\nIDENTIFICATION: Skipping sensitivity plot, because standard deviation of parameters is NaN, possibly due to the use of ML.\n')
        else
            hh_fig = dyn_figure(options_.nodisplay,'Name',[tittxt, ' - Sensitivity plot']);
            subplot(211)
            mmm = (si_dMOMENTSnorm)'./max(si_dMOMENTSnorm);
            mmm1 = (si_dTAUnorm)'./max(si_dTAUnorm);
            mmm=[mmm mmm1];
            mmm1 = (si_dLREnorm)'./max(si_dLREnorm);
            offset=length(si_dTAUnorm)-length(si_dLREnorm);
            mmm1 = [NaN(offset,1); mmm1];
            mmm=[mmm mmm1];

            bar(log(mmm(is,:).*100))
            set(gca,'xlim',[0 nparam+1])
            set(gca,'xticklabel','')
            dy = get(gca,'ylim');
            for ip=1:nparam
                if options_.TeX
                    text(ip,dy(1),name_tex{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','latex')
                else
                    text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
                end
            end
            legend('Moments','Model','LRE model','Location','Best')
            title('Sensitivity bars using derivatives (log-scale)')
            dyn_saveas(hh_fig,[IdentifDirectoryName '/' fname '_sensitivity_' tittxt1 ],options_.nodisplay,options_.graph_format);
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                fidTeX = fopen([IdentifDirectoryName '/' fname '_sensitivity_' tittxt1,'.tex'],'w');
                fprintf(fidTeX,'%% TeX eps-loader file generated by identification.plot.m (Dynare).\n');
                fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s_sensitivity_%s}\n',[IdentifDirectoryName '/' fname],tittxt1);
                fprintf(fidTeX,'\\caption{%s  - Sensitivity plot.}',tit_TeX);
                fprintf(fidTeX,'\\label{Fig:sensitivity:%s}\n',deblank(tittxt));
                fprintf(fidTeX,'\\end{figure}\n\n');
                fprintf(fidTeX,'%% End Of TeX file. \n');
                fclose(fidTeX);
            end
        end
        % identificaton patterns
        for  j=1:size(idemoments.cosndMOMENTS,2)
            pax=NaN(nparam,nparam);
            for i=1:nparam
                namx='';
                for in=1:j
                    dumpindx = idemoments.pars{i,j}(in);
                    if isnan(dumpindx)
                        namx=[namx ' ' sprintf('%-15s','--')];
                    else
                        if options_.TeX
                            namx=[namx ' ' sprintf('%-15s',name_tex{dumpindx})];
                        else
                            namx=[namx ' ' sprintf('%-15s',name{dumpindx})];
                        end
                        pax(i,dumpindx)=idemoments.cosndMOMENTS(i,j);
                    end
                end
            end
            hh_fig = dyn_figure(options_.nodisplay,'Name',[tittxt,' - Collinearity patterns with ', int2str(j) ,' parameter(s)']);
            imagesc(pax,[0 1]);
            set(gca,'xticklabel','')
            set(gca,'yticklabel','')
            for ip=1:nparam
                if options_.TeX
                    text(ip,(0.5),name_tex{ip},'rotation',90,'HorizontalAlignment','left','interpreter','latex')
                    text(0.5,ip,name_tex{ip},'rotation',0,'HorizontalAlignment','right','interpreter','latex')
                else
                    text(ip,(0.5),name{ip},'rotation',90,'HorizontalAlignment','left','interpreter','none')
                    text(0.5,ip,name{ip},'rotation',0,'HorizontalAlignment','right','interpreter','none')
                end
            end
            colorbar;
            colormap('jet');
            ax=colormap;
            ax(1,:)=[0.9 0.9 0.9];
            colormap(ax);
            if nparam>10
                set(gca,'xtick',(5:5:nparam))
                set(gca,'ytick',(5:5:nparam))
            end
            set(gca,'xgrid','on')
            set(gca,'ygrid','on')
            xlabel([tittxt,' - Collinearity patterns with ', int2str(j) ,' parameter(s)'],'interpreter','none')
            dyn_saveas(hh_fig,[ IdentifDirectoryName '/' fname '_ident_collinearity_' tittxt1 '_' int2str(j) ],options_.nodisplay,options_.graph_format);
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                fidTeX = fopen([ IdentifDirectoryName '/' fname '_ident_collinearity_' tittxt1 '_' int2str(j),'.tex'],'w');
                fprintf(fidTeX,'%% TeX eps-loader file generated by identification.plot.m (Dynare).\n');
                fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s_ident_collinearity_%s_%u}\n',[IdentifDirectoryName '/' fname],tittxt1,j);
                fprintf(fidTeX,'\\caption{%s  - Collinearity patterns with %u parameter(s).}',tit_TeX,j);
                fprintf(fidTeX,'\\label{Fig:collinearity:%s:%u_pars}\n',deblank(tittxt),j);
                fprintf(fidTeX,'\\end{figure}\n\n');
                fprintf(fidTeX,'%% End Of TeX file. \n');
                fclose(fidTeX);
            end
        end
        skipline()
        [~,S,V]=svd(idehess.AHess,0);
        S=diag(S);
        if idehess.flag_score
            if nparam<5
                f1 = dyn_figure(options_.nodisplay,'Name',[tittxt,' - Identification patterns (Information matrix)']);
                tex_tit_1=[tittxt,' - Identification patterns (Information matrix)'];
            else
                f1 = dyn_figure(options_.nodisplay,'Name',[tittxt,' - Identification patterns (Information matrix): SMALLEST SV']);
                tex_tit_1=[tittxt,' - Identification patterns (Information matrix): SMALLEST SV'];
                f2 = dyn_figure(options_.nodisplay,'Name',[tittxt,' - Identification patterns (Information matrix): HIGHEST SV']);
                tex_tit_2=[tittxt,' - Identification patterns (Information matrix): HIGHEST SV'];
            end
        else
            if nparam<5
                f1 = dyn_figure(options_.nodisplay,'Name',[tittxt,' - Identification patterns (moments Information matrix)']);
                tex_tit_1=[tittxt,' - Identification patterns (moments Information matrix)'];
            else
                f1 = dyn_figure(options_.nodisplay,'Name',[tittxt,' - Identification patterns (moments Information matrix): SMALLEST SV']);
                tex_tit_1=[tittxt,' - Identification patterns (moments Information matrix): SMALLEST SV'];
                f2 = dyn_figure(options_.nodisplay,'Name',[tittxt,' - Identification patterns (moments Information matrix): HIGHEST SV']);
                tex_tit_2=[tittxt,' - Identification patterns (moments Information matrix): HIGHEST SV'];
            end
        end
        for j=1:min(nparam,8)
            if j<5
                set(0,'CurrentFigure',f1),
                jj=j;
            else
                set(0,'CurrentFigure',f2),
                jj=j-4;
            end
            subplot(4,1,jj)
            if j<5
                bar(abs(V(:,end-j+1)))
                Stit = S(end-j+1);
            else
                bar(abs(V(:,jj))),
                Stit = S(jj);
            end
            set(gca,'xticklabel','')
            if j==4 || j==nparam || j==8
                for ip=1:nparam
                    if options_.TeX
                        text(ip,-0.02,name_tex{ip},'rotation',90,'HorizontalAlignment','right','interpreter','latex')
                    else
                        text(ip,-0.02,name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
                    end
                end
            end
            title(['Singular value ',num2str(Stit)])
        end
        dyn_saveas(f1,[  IdentifDirectoryName '/' fname '_ident_pattern_' tittxt1 '_1' ],options_.nodisplay,options_.graph_format);
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fidTeX = fopen([  IdentifDirectoryName '/' fname '_ident_pattern_' tittxt1 '_1','.tex'],'w');
            fprintf(fidTeX,'%% TeX eps-loader file generated by identification.plot.m (Dynare).\n');
            fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s_ident_pattern_%s_1}\n',[IdentifDirectoryName '/' fname],tittxt1);
            fprintf(fidTeX,'\\caption{%s.}',tex_tit_1);
            fprintf(fidTeX,'\\label{Fig:ident_pattern:%s:1}\n',tittxt1);
            fprintf(fidTeX,'\\end{figure}\n\n');
            fprintf(fidTeX,'%% End Of TeX file. \n');
            fclose(fidTeX);
        end
        if nparam>4
            dyn_saveas(f2,[  IdentifDirectoryName '/' fname '_ident_pattern_' tittxt1 '_2' ],options_.nodisplay,options_.graph_format);
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                fidTeX = fopen([  IdentifDirectoryName '/' fname '_ident_pattern_' tittxt1 '_2.tex'],'w');
                fprintf(fidTeX,'%% TeX eps-loader file generated by identification.plot.m (Dynare).\n');
                fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s_ident_pattern_%s_2}\n',[IdentifDirectoryName '/' fname],tittxt1);
                fprintf(fidTeX,'\\caption{%s.}',tex_tit_2);
                fprintf(fidTeX,'\\label{Fig:ident_pattern:%s:2}\n',tittxt1);
                fprintf(fidTeX,'\\end{figure}\n\n');
                fprintf(fidTeX,'%% End Of TeX file. \n');
                fclose(fidTeX);
            end
        end
    end

else
    hh_fig = dyn_figure(options_.nodisplay,'Name','MC sensitivities');
    subplot(211)
    mmm = (idehess.ide_strength_dMOMENTS);
    [~, is] = sort(mmm);
    mmm = mean(si_dMOMENTSnorm)';
    mmm = mmm./max(mmm);
    if advanced
        mmm1 = mean(si_dTAUnorm)';
        mmm=[mmm mmm1./max(mmm1)];
        mmm1 = mean(si_dLREnorm)';
        offset=size(si_dTAUnorm,2)-size(si_dLREnorm,2);
        mmm1 = [NaN(offset,1); mmm1./max(mmm1)];
        mmm=[mmm mmm1];
    end

    bar(mmm(is,:))
    set(gca,'xlim',[0 nparam+1])
    set(gca,'xticklabel','')
    dy = get(gca,'ylim');
    for ip=1:nparam
        if options_.TeX
            text(ip,dy(1),name_tex{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','latex')
        else
            text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
        end
    end
    if advanced
        legend('Moments','Model','LRE model','Location','Best')
    end
    title('MC mean of sensitivity measures')
    dyn_saveas(hh_fig,[ IdentifDirectoryName '/' fname '_MC_sensitivity' ],options_.nodisplay,options_.graph_format);
    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fidTeX = fopen([ IdentifDirectoryName '/' fname '_MC_sensitivity.tex'],'w');
        fprintf(fidTeX,'%% TeX eps-loader file generated by identification.plot.m (Dynare).\n');
        fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s_MC_sensitivity}\n',[IdentifDirectoryName '/' fname]);
        fprintf(fidTeX,'\\caption{MC mean of sensitivity measures}');
        fprintf(fidTeX,'\\label{Fig:_MC_sensitivity}\n');
        fprintf(fidTeX,'\\end{figure}\n\n');
        fprintf(fidTeX,'%% End Of TeX file. \n');
        fclose(fidTeX);
    end

    if advanced
        if ~options_.nodisplay
            skipline()
            disp('Displaying advanced diagnostics')
        end
        hh_fig = dyn_figure(options_.nodisplay,'Name','MC Condition Number');
        subplot(221)
        if isoctave
            hist(log10(idemodel.cond))
        else
            histogram(log10(idemodel.cond))
        end
        title('log10 of Condition number in the model')
        subplot(222)
        if isoctave
            hist(log10(idemoments.cond))
        else
            histogram(log10(idemoments.cond))
        end
        title('log10 of Condition number in the moments')
        subplot(223)
        if isoctave
            hist(log10(idelre.cond))
        else
            histogram(log10(idelre.cond))
        end
        title('log10 of Condition number in the LRE model')
        dyn_saveas(hh_fig,[IdentifDirectoryName '/' fname '_ident_COND' ],options_.nodisplay,options_.graph_format);
        options_mcf.pvalue_ks = 0.1;
        options_mcf.pvalue_corr = 0.001;
        options_mcf.alpha2 = 0;
        options_mcf.param_names = name;
        options_mcf.param_names_tex = name_tex;
        options_mcf.fname_ = fname;
        options_mcf.OutputDirectoryName = IdentifDirectoryName;
        options_mcf.beha_title = 'LOW condition nbr';
        options_mcf.nobeha_title = 'HIGH condition nbr';
        if options_.TeX
            options_mcf.beha_title_latex = 'LOW condition nbr';
            options_mcf.nobeha_title_latex = 'HIGH condition nbr';
        end
        options_mcf.amcf_name = 'MC_HighestCondNumberLRE';
        options_mcf.amcf_title = 'MC Highest Condition Number LRE Model';
        options_mcf.title = 'MC Highest Condition Number LRE Model';
        ncut=floor(SampleSize/10*9);
        [~,is]=sort(idelre.cond);
        gsa.monte_carlo_filtering_analysis(params, is(1:ncut), is(ncut+1:end), options_mcf, M_, options_, bayestopt_, estim_params_);
        options_mcf.amcf_name = 'MC_HighestCondNumberModel';
        options_mcf.amcf_title = 'MC Highest Condition Number Model Solution';
        options_mcf.title = 'MC Highest Condition Number Model Solution';
        [~,is]=sort(idemodel.cond);
        gsa.monte_carlo_filtering_analysis(params, is(1:ncut), is(ncut+1:end), options_mcf, M_, options_, bayestopt_, estim_params_);
        options_mcf.amcf_name = 'MC_HighestCondNumberMoments';
        options_mcf.amcf_title = 'MC Highest Condition Number Model Moments';
        options_mcf.title = 'MC Highest Condition Number Model Moments';
        [~,is]=sort(idemoments.cond);
        gsa.monte_carlo_filtering_analysis(params, is(1:ncut), is(ncut+1:end), options_mcf, M_, options_, bayestopt_, estim_params_);

        if nparam<5
            f1 = dyn_figure(options_.nodisplay,'Name',[tittxt,' - MC Identification patterns (moments): HIGHEST SV']);
            tex_tit_1=[tittxt,' - MC Identification patterns (moments): HIGHEST SV'];
        else
            f1 = dyn_figure(options_.nodisplay,'Name',[tittxt,' - MC Identification patterns (moments): SMALLEST SV']);
            tex_tit_1=[tittxt,' - MC Identification patterns (moments): SMALLEST SV'];
            f2 = dyn_figure(options_.nodisplay,'Name',[tittxt,' - MC Identification patterns (moments): HIGHEST SV']);
            tex_tit_2=[tittxt,' - MC Identification patterns (moments): HIGHEST SV'];
        end
        nplots=min(nparam,8);
        if nplots>4
            nsubplo=ceil(nplots/2);
        else
            nsubplo=nplots;
        end
        for j=1:nplots
            if (nparam>4 && j<=ceil(nplots/2)) || nparam<5
                set(0,'CurrentFigure',f1),
                jj=j;
                VVV=squeeze(abs(idemoments.V(:,:,end-j+1)));
                SSS = idemoments.S(:,end-j+1);
            else
                set(0,'CurrentFigure',f2),
                jj=j-ceil(nplots/2);
                VVV=squeeze(abs(idemoments.V(:,:,jj)));
                SSS = idemoments.S(:,jj);
            end
            subplot(nsubplo,1,jj)
            post_median=NaN(1,nparam);
            hpd_interval=NaN(nparam,2);
            for i=1:nparam
                [~, post_median(:,i), ~, hpd_interval(i,:)] = posterior_moments(VVV(:,i),0.9);
            end
            bar(post_median)
            hold on, plot(hpd_interval,'--*r'),
            Stit=mean(SSS);

            set(gca,'xticklabel','')
            if j==4 || j==nparam || j==8
                for ip=1:nparam
                    if options_.TeX
                        text(ip,-0.02,name_tex{ip},'rotation',90,'HorizontalAlignment','right','interpreter','latex')
                    else
                        text(ip,-0.02,name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
                    end
                end
            end
            title(['MEAN Singular value ',num2str(Stit)])
        end
        dyn_saveas(f1,[IdentifDirectoryName '/' fname '_MC_ident_pattern_1' ],options_.nodisplay,options_.graph_format);
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fidTeX = fopen([IdentifDirectoryName '/' fname '_MC_ident_pattern_1.tex'],'w');
            fprintf(fidTeX,'%% TeX eps-loader file generated by identification.plot.m (Dynare).\n');
            fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s_MC_ident_pattern_1}\n',[IdentifDirectoryName '/' fname]);
            fprintf(fidTeX,'\\caption{%s.}',tex_tit_1);
            fprintf(fidTeX,'\\label{Fig:MC_ident_pattern:1}\n');
            fprintf(fidTeX,'\\end{figure}\n\n');
            fprintf(fidTeX,'%% End Of TeX file. \n');
            fclose(fidTeX);
        end
        if nparam>4
            dyn_saveas(f2,[  IdentifDirectoryName '/' fname '_MC_ident_pattern_2' ],options_.nodisplay,options_.graph_format);
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                fidTeX = fopen([  IdentifDirectoryName '/' fname '_MC_ident_pattern_2.tex'],'w');
                fprintf(fidTeX,'%% TeX eps-loader file generated by identification.plot.m (Dynare).\n');
                fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s_MC_ident_pattern_2}\n',[IdentifDirectoryName '/' fname]);
                fprintf(fidTeX,'\\caption{%s.}',tex_tit_2);
                fprintf(fidTeX,'\\label{Fig:MC_ident_pattern:2}\n');
                fprintf(fidTeX,'\\end{figure}\n\n');
                fprintf(fidTeX,'%% End Of TeX file. \n');
                fclose(fidTeX);
            end
        end
    end
end
