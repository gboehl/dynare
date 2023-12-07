function reduced_form_screening(dirname, options_gsa_, estim_params_, M_, dr, options_, bayestopt_)
% reduced_form_screening(dirname, options_gsa_, estim_params_, M_, dr, options_, bayestopt_)
% Conduct reduced form screening
% Inputs:
%  - dirname             [string]    name of the output directory
%  - options_gsa_        [structure] GSA options_
%  - estim_params        [structure] describing the estimated parameters
%  - M_                  [structure] describing the model
%  - dr                  [structure] decision rules
%  - options_            [structure] describing the options
%  - bayestopt_          [structure] describing the priors
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

anamendo = options_gsa_.namendo;
anamlagendo = options_gsa_.namlagendo;
anamexo = options_gsa_.namexo;
anamendo_tex = options_gsa_.namendo_tex;
anamlagendo_tex = options_gsa_.namlagendo_tex;
anamexo_tex = options_gsa_.namexo_tex;

nliv = options_gsa_.morris_nliv;

pnames = M_.param_names(estim_params_.param_vals(:,1));
if options_.TeX
    for par_iter=1:size(estim_params_.param_vals(:,1),1)
        [~,tex_names{par_iter,1}]=get_the_name(estim_params_.param_vals(par_iter,1),options_.TeX, M_, estim_params_, options_.varobs);
    end
end
if nargin==0
    dirname='';
end

load([dirname,'/',M_.fname,'_prior'],'lpmat','lpmat0','istable','T');

nspred=M_.nspred;

[kn, np]=size(lpmat);
nshock = length(bayestopt_.pshape)-np;

nsok = length(find(M_.lead_lag_incidence(M_.maximum_lag,:)));

js=0;
for j=1:size(anamendo,1)
    namendo = anamendo{j,:};
    namendo_tex = anamendo_tex{j,:};
    iendo = strmatch(namendo, M_.endo_names(dr.order_var), 'exact');
    iplo=0;
    ifig=0;
    for jx=1:size(anamexo,1)
        namexo = anamexo{jx};
        namexo_tex = anamexo_tex{jx};
        iexo = strmatch(namexo, M_.exo_names, 'exact');
        if ~isempty(iexo)
            y0=gsa.teff(T(iendo,iexo+nspred,:), kn, istable);
            if ~isempty(y0)
                if mod(iplo,9)==0
                    ifig = ifig+1;
                    hh_fig = dyn_figure(options_.nodisplay, 'name', [namendo,' vs. shocks ', int2str(ifig)]);
                    iplo = 0;
                end
                iplo = iplo+1;
                js = js+1;
                subplot(3, 3, iplo)
                [~, SAMorris] = gsa.Morris_Measure_Groups(np+nshock, [lpmat0 lpmat], y0, nliv);
                SAM = squeeze(SAMorris(nshock+1:end,1));
                SA(:,js) = SAM./(max(SAM)+eps);
                [~, iso] = sort(-SA(:,js));
                bar(SA(iso(1:min(np,10)),js))
                set(gca,'xticklabel',' ','fontsize',10)
                set(gca,'xlim',[0.5 10.5])
                for ip=1:min(np,10)
                    if options_.TeX
                        text(ip,-0.02,tex_names(iso(ip)),'rotation',90,'HorizontalAlignment','right','interpreter','latex')
                    else
                        text(ip,-0.02,pnames(iso(ip)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
                    end
                end
                if options_.TeX
                    title([namendo_tex,' vs. ',namexo_tex],'interpreter','latex')
                else
                    title([namendo,' vs. ',namexo],'interpreter','none')
                end
                if iplo==9
                    dyn_saveas(hh_fig,[dirname,'/',M_.fname,'_', namendo,'_vs_shock_',num2str(ifig)],options_.nodisplay,options_.graph_format);
                    create_TeX_loader(options_,[dirname,'/',M_.fname,'_', namendo,'_vs_shock_',num2str(ifig)],ifig,[namendo,' vs. shocks ',int2str(ifig)],[namendo,'_vs_shock'],1)
                end

            end
        end
    end
    if iplo<9 && iplo>0 && ifig
        dyn_saveas(hh_fig,[dirname,'/',M_.fname,'_', namendo,'_vs_shocks_',num2str(ifig)],options_.nodisplay,options_.graph_format);
        create_TeX_loader(options_,[dirname,'/',M_.fname,'_', namendo,'_vs_shock_',num2str(ifig)],ifig,[namendo,' vs. shocks ',int2str(ifig)],[namendo,'_vs_shock'],options_.figures.textwidth*min(iplo/3))
    end

    iplo=0;
    ifig=0;
    for je=1:size(anamlagendo,1)
        namlagendo=anamlagendo{je};
        namlagendo_tex=anamlagendo_tex{je};
        ilagendo=strmatch(namlagendo, M_.endo_names(dr.order_var(M_.nstatic+1:M_.nstatic+nsok)), 'exact');

        if ~isempty(ilagendo)
            y0=gsa.teff(T(iendo,ilagendo,:),kn,istable);
            if ~isempty(y0)
                if mod(iplo,9)==0
                    ifig=ifig+1;
                    hh_fig=dyn_figure(options_.nodisplay,'name',[namendo,' vs. lagged endogenous ',int2str(ifig)]);
                    iplo=0;
                end
                iplo=iplo+1;
                js=js+1;
                subplot(3,3,iplo),
                [~, SAMorris] = gsa.Morris_Measure_Groups(np+nshock, [lpmat0 lpmat], y0,nliv);
                SAM = squeeze(SAMorris(nshock+1:end,1));
                SA(:,js)=SAM./(max(SAM)+eps);
                [~, iso] = sort(-SA(:,js));
                bar(SA(iso(1:min(np,10)),js))
                set(gca,'xticklabel',' ','fontsize',10)
                set(gca,'xlim',[0.5 10.5])
                for ip=1:min(np,10)
                    if options_.TeX
                        text(ip,-0.02,tex_names(iso(ip)),'rotation',90,'HorizontalAlignment','right','interpreter','latex')
                    else
                        text(ip,-0.02,pnames{iso(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
                    end
                end

                if options_.TeX
                    title([namendo_tex,' vs. ',namlagendo_tex,'(-1)'],'interpreter','latex')
                else
                    title([namendo,' vs. ',namlagendo,'(-1)'],'interpreter','none')
                end
                if iplo==9
                    dyn_saveas(hh_fig,[dirname,'/',M_.fname,'_', namendo,'_vs_lags_',num2str(ifig)],options_.nodisplay,options_.graph_format);
                    create_TeX_loader(options_,[dirname,'/',M_.fname,'_', namendo,'_vs_lags_',num2str(ifig)],ifig,[namendo,' vs. lagged endogenous ',int2str(ifig)],[namendo,'_vs_lags'],1)
                end
            end
        end
    end
    if iplo<9 && iplo>0 && ifig
        dyn_saveas(hh_fig,[dirname,'/',M_.fname,'_', namendo,'_vs_lags_',num2str(ifig)],options_.nodisplay,options_.graph_format);
        create_TeX_loader(options_,[dirname,'/',M_.fname,'_', namendo,'_vs_lags_',num2str(ifig)],ifig,[namendo,' vs. lagged endogenous ',int2str(ifig)],[namendo,'_vs_lags'],options_.figures.textwidth*min(iplo/3))
    end
end

hh_fig=dyn_figure(options_.nodisplay,'Name','Reduced form screening');
gsa.boxplot(SA',[],'.',[],10)
set(gca,'xticklabel',' ','fontsize',10,'xtick',1:np)
set(gca,'xlim',[0.5 np+0.5])
set(gca,'ylim',[0 1])
set(gca,'position',[0.13 0.2 0.775 0.7])
for ip=1:np
    if options_.TeX
        text(ip,-0.02,tex_names(ip),'rotation',90,'HorizontalAlignment','right','interpreter','latex')
    else
        text(ip,-0.02,pnames{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
end
xlabel(' ')
ylabel('Elementary Effects')
title('Reduced form screening')
dyn_saveas(hh_fig,[dirname,'/',M_.fname,'_redform_screen'],options_.nodisplay,options_.graph_format);
create_TeX_loader(options_,[dirname,'/',M_.fname,'_redform_screen'],1,'Reduced form screening','redform_screen',1)


function []=create_TeX_loader(options_,figpath,label_number,caption,label_name,scale_factor)
if nargin<6
    scale_factor=1;
end
if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fidTeX = fopen([figpath '.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by reduced_form_screening.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s}\n',scale_factor,strrep(figpath,'\','/'));
    fprintf(fidTeX,'\\caption{%s.}',caption);
    fprintf(fidTeX,'\\label{Fig:%s:%u}\n',label_name,label_number);
    fprintf(fidTeX,'\\end{figure}\n\n');
    fprintf(fidTeX,'%% End Of TeX file. \n');
    fclose(fidTeX);
end
