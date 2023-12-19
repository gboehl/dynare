function map_calibration(OutputDirectoryName, M_, options_, oo_, estim_params_, bayestopt_)
% map_calibration(OutputDirectoryName, M_, options_, oo_, estim_params_, bayestopt_)
% Inputs:
%  - OutputDirectoryName [string]    name of the output directory
%  - M_                  [structure] describing the model
%  - options_            [structure] describing the options
%  - oo_                 [structure] storing the results
%  - estim_params_       [structure] characterizing parameters to be estimated
%  - bayestopt_          [structure] describing the priors
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

fname_ = M_.fname;

np = estim_params_.np;
nshock = estim_params_.nvx + estim_params_.nvn + estim_params_.ncx + estim_params_.ncn;
pnames=cell(np,1);
pnames_tex=cell(np,1);
for jj=1:np
    if options_.TeX
        [param_name_temp, param_name_tex_temp]= get_the_name(nshock+jj, options_.TeX, M_, estim_params_, options_.varobs);
        pnames_tex{jj,1} = param_name_tex_temp;
        pnames{jj,1} = param_name_temp;
    else
        param_name_temp = get_the_name(nshock+jj, options_.TeX, M_, estim_params_, options_.varobs);
        pnames{jj,1} = param_name_temp;
    end
end

indx_irf = [];
indx_moment = [];
init = ~options_.opt_gsa.load_stab;

options_mcf.pvalue_ks = options_.opt_gsa.pvalue_ks;
options_mcf.pvalue_corr = options_.opt_gsa.pvalue_corr;
options_mcf.alpha2 = options_.opt_gsa.alpha2_stab;
options_mcf.param_names = pnames;
if options_.TeX
    options_mcf.param_names_tex = pnames_tex;
end
options_mcf.fname_ = fname_;
options_mcf.OutputDirectoryName = OutputDirectoryName;

skipline()
disp('Sensitivity analysis for calibration criteria')

if options_.opt_gsa.ppost
    filetoload=dir([M_.dname filesep 'metropolis' filesep fname_ '_param_irf*.mat']);
    lpmat=[];
    for j=1:length(filetoload)
        load([M_.dname filesep 'metropolis' filesep fname_ '_param_irf',int2str(j),'.mat'],'stock')
        lpmat = [lpmat; stock];
        clear stock
    end
    type = 'post';
else
    if options_.opt_gsa.pprior
        filetoload=[OutputDirectoryName '/' fname_ '_prior'];
        load(filetoload,'lpmat','lpmat0')
        lpmat = [lpmat0 lpmat];
        type = 'prior';
    else
        filetoload=[OutputDirectoryName '/' fname_ '_mc'];
        load(filetoload,'lpmat','lpmat0')
        lpmat = [lpmat0 lpmat];
        type = 'mc';
    end
end
[Nsam, np] = size(lpmat);
npar = size(pnames,1);
nshock = np - npar;

nbr_irf_restrictions = size(options_.endogenous_prior_restrictions.irf,1);
nbr_moment_restrictions = size(options_.endogenous_prior_restrictions.moment,1);

if init
    mat_irf=cell(nbr_irf_restrictions,1);
    for ij=1:nbr_irf_restrictions
        mat_irf{ij}=NaN(Nsam,length(options_.endogenous_prior_restrictions.irf{ij,3}));
    end

    mat_moment=cell(nbr_moment_restrictions,1);
    for ij=1:nbr_moment_restrictions
        mat_moment{ij}=NaN(Nsam,length(options_.endogenous_prior_restrictions.moment{ij,3}));
    end

    irestrictions = 1:Nsam;
    h = dyn_waitbar(0,'Please wait...');
    for j=1:Nsam
        M_ = set_all_parameters(lpmat(j,:)',estim_params_,M_);
        if nbr_moment_restrictions
            [Tt,Rr,~,info,oo_.dr, M_.params] = dynare_resolve(M_,options_,oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
        else
            [Tt,Rr,~,info,oo_.dr, M_.params] = dynare_resolve(M_,options_,oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state,'restrict');
        end
        if info(1)==0
            [~, info_irf, info_moment, data_irf, data_moment]=endogenous_prior_restrictions(Tt,Rr,M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
            if ~isempty(info_irf)
                for ij=1:nbr_irf_restrictions
                    mat_irf{ij}(j,:)=data_irf{ij}(:,2)';
                end
                indx_irf(j,:)=info_irf(:,1);
            end
            if ~isempty(info_moment)
                for ij=1:nbr_moment_restrictions
                    mat_moment{ij}(j,:)=data_moment{ij}(:,2)';
                end
                indx_moment(j,:)=info_moment(:,1);
            end
        else
            irestrictions(j)=0;
        end
        if mod(j,3)==0
            dyn_waitbar(j/Nsam,h,['MC iteration ',int2str(j),'/',int2str(Nsam)])
        end
    end
    dyn_waitbar_close(h);

    irestrictions=irestrictions(find(irestrictions));
    xmat=lpmat(irestrictions,:);
    skipline()
    endo_prior_restrictions=options_.endogenous_prior_restrictions;
    save([OutputDirectoryName,filesep,fname_,'_',type,'_restrictions'],'xmat','mat_irf','mat_moment','irestrictions','indx_irf','indx_moment','endo_prior_restrictions');
else
    load([OutputDirectoryName,filesep,fname_,'_',type,'_restrictions'],'xmat','mat_irf','mat_moment','irestrictions','indx_irf','indx_moment','endo_prior_restrictions');
end
if ~isempty(indx_irf)
    skipline()
    disp('Deleting old IRF calibration plots ...')
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_irf_calib*.eps']);
    for j=1:length(a)
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_irf_calib*.fig']);
    for j=1:length(a)
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_irf_calib*.pdf']);
    for j=1:length(a)
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_irf_restrictions.eps']);
    for j=1:length(a)
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_irf_restrictions.fig']);
    for j=1:length(a)
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_irf_restrictions.pdf']);
    for j=1:length(a)
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    disp('done !')
    skipline()

    % For single legend search which has maximum nbr of restrictions
    all_irf_couples = cellstr([char(endo_prior_restrictions.irf(:,1)) char(endo_prior_restrictions.irf(:,2))]);
    irf_couples = unique(all_irf_couples);
    nbr_irf_couples = size(irf_couples,1);
    plot_indx = NaN(nbr_irf_couples,1);
    time_matrix=cell(nbr_irf_couples,1);
    indx_irf_matrix=zeros(length(irestrictions),nbr_irf_couples);
    irf_matrix=cell(nbr_irf_couples,1);
    irf_mean=cell(nbr_irf_couples,1);
    irf_median=cell(nbr_irf_couples,1);
    irf_var=cell(nbr_irf_couples,1);
    irf_HPD=cell(nbr_irf_couples,1);
    irf_distrib=cell(nbr_irf_couples,1);
    maxijv=0;
    for ij=1:nbr_irf_restrictions
        if length(endo_prior_restrictions.irf{ij,3})>maxijv
            maxijv=length(endo_prior_restrictions.irf{ij,3});
        end
        plot_indx(ij) = find(strcmp(irf_couples,all_irf_couples(ij,:)));
        time_matrix{plot_indx(ij)} = [time_matrix{plot_indx(ij)} endo_prior_restrictions.irf{ij,3}];
    end
    iplot_indx = ones(size(plot_indx));

    indx_irf = indx_irf(irestrictions,:);
    if ~options_.nograph
        h1=dyn_figure(options_.nodisplay,'name',[type ' evaluation of irf restrictions']);
        nrow=ceil(sqrt(nbr_irf_couples));
        ncol=nrow;
        if nrow*(nrow-1)>nbr_irf_couples
            ncol=nrow-1;
        end
    end
    for ij=1:nbr_irf_restrictions
        mat_irf{ij}=mat_irf{ij}(irestrictions,:);
        irf_matrix{plot_indx(ij)} = [irf_matrix{plot_indx(ij)} mat_irf{ij}];
        indx_irf_matrix(:,plot_indx(ij)) = indx_irf_matrix(:,plot_indx(ij)) + indx_irf(:,ij);
        for ik=1:size(mat_irf{ij},2)
            [Mean,Median,Var,HPD,Distrib] = ...
                posterior_moments(mat_irf{ij}(:,ik),options_.mh_conf_sig);
            irf_mean{plot_indx(ij)} = [irf_mean{plot_indx(ij)}; Mean];
            irf_median{plot_indx(ij)} = [irf_median{plot_indx(ij)}; Median];
            irf_var{plot_indx(ij)} = [irf_var{plot_indx(ij)}; Var];
            irf_HPD{plot_indx(ij)} = [irf_HPD{plot_indx(ij)}; HPD];
            irf_distrib{plot_indx(ij)} = [irf_distrib{plot_indx(ij)}; Distrib'];
        end
        leg = num2str(endo_prior_restrictions.irf{ij,3}(1));
        aleg = num2str(endo_prior_restrictions.irf{ij,3}(1));
        if size(mat_irf{ij},2)>1
            leg = [leg,':' ,num2str(endo_prior_restrictions.irf{ij,3}(end))];
            aleg = [aleg,'-' ,num2str(endo_prior_restrictions.irf{ij,3}(end))];
            iplot_indx(ij)=0;
        end
        if ~options_.nograph && length(time_matrix{plot_indx(ij)})==1
            set(0,'currentfigure',h1),
            subplot(nrow,ncol, plot_indx(ij)),
            hc = gsa.cumplot(mat_irf{ij}(:,ik));
            a=axis;
            delete(hc);
            x1val=max(endo_prior_restrictions.irf{ij,4}(1),a(1));
            x2val=min(endo_prior_restrictions.irf{ij,4}(2),a(2));
            hp = patch([x1val x2val x2val x1val],a([3 3 4 4]),'b');
            hold all,
            set(hp,'FaceColor', [0.7 0.8 1])
            hc = gsa.cumplot(mat_irf{ij}(:,ik));
            set(hc,'color','k','linewidth',2)
            hold off,
            %         hold off,
            title([endo_prior_restrictions.irf{ij,1},' vs ',endo_prior_restrictions.irf{ij,2}, '(', leg,')'],'interpreter','none'),
        end
        indx1 = find(indx_irf(:,ij)==0);
        indx2 = find(indx_irf(:,ij)~=0);
        atitle0=[endo_prior_restrictions.irf{ij,1},' vs ',endo_prior_restrictions.irf{ij,2}, '(', leg,')'];
        fprintf(['%4.1f%% of the ',type,' support matches IRF ',atitle0,' inside [%4.1f, %4.1f]\n'],length(indx1)/length(irestrictions)*100,endo_prior_restrictions.irf{ij,4})
        aname=[type '_irf_calib_',endo_prior_restrictions.irf{ij,1},'_vs_',endo_prior_restrictions.irf{ij,2},'_',aleg];
        atitle=[type ' IRF Calib: Parameter(s) driving ',endo_prior_restrictions.irf{ij,1},' vs ',endo_prior_restrictions.irf{ij,2}, '(', leg,')'];
        options_mcf.amcf_name = aname;
        options_mcf.amcf_title = atitle;
        options_mcf.beha_title = 'IRF restriction';
        options_mcf.nobeha_title = 'NO IRF restriction';
        if options_.TeX
            options_mcf.beha_title_latex = 'IRF restriction';
            options_mcf.nobeha_title_latex = 'NO IRF restriction';
        end
        options_mcf.title = atitle0;
        if ~isempty(indx1) && ~isempty(indx2)
            gsa.monte_carlo_filtering_analysis(xmat(:,nshock+1:end), indx1, indx2, options_mcf, M_, options_, bayestopt_, estim_params_);
        end
    end
    for ij=1:nbr_irf_couples
        if length(time_matrix{ij})>1
            if ~options_.nograph
                set(0,'currentfigure',h1);
                subplot(nrow,ncol, ij)
                itmp = (find(plot_indx==ij));
                htmp = plot(time_matrix{ij},[max(irf_matrix{ij})' min(irf_matrix{ij})'],'k--','linewidth',2);
                a=axis;
                delete(htmp);
                tmp=[];
                for ir=1:length(itmp)
                    for it=1:length(endo_prior_restrictions.irf{itmp(ir),3})
                        temp_index = find(time_matrix{ij}==endo_prior_restrictions.irf{itmp(ir),3}(it));
                        tmp(temp_index,:) = endo_prior_restrictions.irf{itmp(ir),4};
                    end
                end
                tmp(isinf(tmp(:,1)),1)=a(3);
                tmp(isinf(tmp(:,2)),2)=a(4);
                hp = patch([time_matrix{ij} time_matrix{ij}(end:-1:1)],[tmp(:,1); tmp(end:-1:1,2)],'c');
                set(hp,'FaceColor',[0.7 0.8 1])
                hold on,
                plot(time_matrix{ij},[max(irf_matrix{ij})' min(irf_matrix{ij})'],'k--','linewidth',2)
                plot(time_matrix{ij},irf_median{ij},'k','linewidth',2)
                plot(time_matrix{ij},[irf_distrib{ij}],'k-')
                plot(a(1:2),[0 0],'r')
                hold off
                axis([max(1,a(1)) a(2:4)])
                box on
                itmp = min(itmp);
                title([endo_prior_restrictions.irf{itmp,1},' vs ',endo_prior_restrictions.irf{itmp,2}],'interpreter','none'),
            end
            if any(iplot_indx.*plot_indx==ij)
                % MCF of the couples with logical AND
                itmp = min(find(plot_indx==ij));
                indx1 = find(indx_irf_matrix(:,ij)==0);
                indx2 = find(indx_irf_matrix(:,ij)~=0);
                leg = num2str(time_matrix{ij}(1));
                leg = [leg '...' num2str(time_matrix{ij}(end))];
                aleg = 'ALL';
                atitle0=[endo_prior_restrictions.irf{itmp,1},' vs ',endo_prior_restrictions.irf{itmp,2}, '(', leg,')'];
                fprintf(['%4.1f%% of the ',type,' support matches IRF restrictions ',atitle0,'\n'],length(indx1)/length(irestrictions)*100)
                aname=[type '_irf_calib_',endo_prior_restrictions.irf{itmp,1},'_vs_',endo_prior_restrictions.irf{itmp,2},'_',aleg];
                atitle=[type ' IRF Calib: Parameter(s) driving ',endo_prior_restrictions.irf{itmp,1},' vs ',endo_prior_restrictions.irf{itmp,2}, '(', leg,')'];
                options_mcf.amcf_name = aname;
                options_mcf.amcf_title = atitle;
                options_mcf.beha_title = 'IRF restriction';
                options_mcf.nobeha_title = 'NO IRF restriction';
                if options_.TeX
                    options_mcf.beha_title_latex = 'IRF restriction';
                    options_mcf.nobeha_title_latex = 'NO IRF restriction';
                end

                options_mcf.title = atitle0;
                if ~isempty(indx1) && ~isempty(indx2)
                    gsa.monte_carlo_filtering_analysis(xmat(:,nshock+1:end), indx1, indx2, options_mcf, M_, options_, bayestopt_, estim_params_);
                end
            end
        end
    end
    if ~options_.nograph
        dyn_saveas(h1,[OutputDirectoryName,filesep,fname_,'_',type,'_irf_restrictions'],options_.nodisplay,options_.graph_format);
        create_TeX_loader(options_,[OutputDirectoryName,filesep,fname_,'_',type,'_irf_restrictions'],[type ' evaluation of irf restrictions'],'irf_restrictions',type,options_.figures.textwidth*min(ij/ncol,1))
    end
    skipline()
end

if ~isempty(indx_moment)
    skipline()
    disp('Deleting old MOMENT calibration plots ...')
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_moment_calib*.eps']);
    for j=1:length(a)
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_moment_calib*.fig']);
    for j=1:length(a)
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_moment_calib*.pdf']);
    for j=1:length(a)
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_moment_restrictions.eps']);
    for j=1:length(a)
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_moment_restrictions.fig']);
    for j=1:length(a)
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_moment_restrictions.pdf']);
    for j=1:length(a)
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    disp('done !')
    skipline()

    %get parameter names including standard deviations
    np=size(bayestopt_.name,1);
    name=cell(np,1);
    name_tex=cell(np,1);
    for jj=1:np
        if options_.TeX
            [param_name_temp, param_name_tex_temp]= get_the_name(jj,options_.TeX,M_,estim_params_,options_.varobs);
            name_tex{jj,1} = param_name_tex_temp;
            name{jj,1} = param_name_temp;
        else
            param_name_temp = get_the_name(jj,options_.TeX,M_,estim_params_,options_.varobs);
            name{jj,1} = param_name_temp;
        end
    end
    options_mcf.param_names = name;
    if options_.TeX
        options_mcf.param_names_tex = name_tex;
    end
    options_mcf.param_names = bayestopt_.name;
    all_moment_couples = cellstr([char(endo_prior_restrictions.moment(:,1)) char(endo_prior_restrictions.moment(:,2))]);
    moment_couples = unique(all_moment_couples);
    nbr_moment_couples = size(moment_couples,1);
    plot_indx = NaN(nbr_moment_couples,1);
    time_matrix=cell(nbr_moment_couples,1);
    indx_moment_matrix=zeros(length(irestrictions),nbr_moment_couples);
    moment_matrix=cell(nbr_moment_couples,1);
    moment_mean=cell(nbr_moment_couples,1);
    moment_median=cell(nbr_moment_couples,1);
    moment_var=cell(nbr_moment_couples,1);
    moment_HPD=cell(nbr_moment_couples,1);
    moment_distrib=cell(nbr_moment_couples,1);
    % For single legend search which has maximum nbr of restrictions
    maxijv=0;
    for ij=1:nbr_moment_restrictions
        endo_prior_restrictions.moment{ij,3} = sort(endo_prior_restrictions.moment{ij,3});
        if length(endo_prior_restrictions.moment{ij,3})>maxijv
            maxijv=length(endo_prior_restrictions.moment{ij,3});
        end
        plot_indx(ij) = find(strcmp(moment_couples,all_moment_couples(ij,:)));
        time_matrix{plot_indx(ij)} = [time_matrix{plot_indx(ij)} endo_prior_restrictions.moment{ij,3}];
    end
    iplot_indx = ones(size(plot_indx));

    indx_moment = indx_moment(irestrictions,:);
    if ~options_.nograph
        h2=dyn_figure(options_.nodisplay,'name',[type ' evaluation of moment restrictions']);
        nrow=ceil(sqrt(nbr_moment_couples));
        ncol=nrow;
        if nrow*(nrow-1)>nbr_moment_couples
            ncol=nrow-1;
        end
    end

    for ij=1:nbr_moment_restrictions
        mat_moment{ij}=mat_moment{ij}(irestrictions,:);
        moment_matrix{plot_indx(ij)} = [moment_matrix{plot_indx(ij)} mat_moment{ij}];
        indx_moment_matrix(:,plot_indx(ij)) = indx_moment_matrix(:,plot_indx(ij)) + indx_moment(:,ij);
        for ik=1:size(mat_moment{ij},2)
            [Mean,Median,Var,HPD,Distrib] = ...
                posterior_moments(mat_moment{ij}(:,ik),options_.mh_conf_sig);
            moment_mean{plot_indx(ij)} = [moment_mean{plot_indx(ij)}; Mean];
            moment_median{plot_indx(ij)} = [moment_median{plot_indx(ij)}; Median];
            moment_var{plot_indx(ij)} = [moment_var{plot_indx(ij)}; Var];
            moment_HPD{plot_indx(ij)} = [moment_HPD{plot_indx(ij)}; HPD];
            moment_distrib{plot_indx(ij)} = [moment_distrib{plot_indx(ij)}; Distrib'];
        end
        leg = num2str(endo_prior_restrictions.moment{ij,3}(1));
        aleg = num2str(endo_prior_restrictions.moment{ij,3}(1));
        if size(mat_moment{ij},2)>1
            leg = [leg,':' ,num2str(endo_prior_restrictions.moment{ij,3}(end))];
            aleg = [aleg,'_' ,num2str(endo_prior_restrictions.moment{ij,3}(end))];
            iplot_indx(ij)=0;
        end
        if ~options_.nograph && length(time_matrix{plot_indx(ij)})==1
            set(0,'currentfigure',h2);
            subplot(nrow,ncol,plot_indx(ij)),
            hc = gsa.cumplot(mat_moment{ij}(:,ik));
            a=axis; delete(hc),
            %     hist(mat_moment{ij}),
            x1val=max(endo_prior_restrictions.moment{ij,4}(1),a(1));
            x2val=min(endo_prior_restrictions.moment{ij,4}(2),a(2));
            hp = patch([x1val x2val x2val x1val],a([3 3 4 4]),'b');
            set(hp,'FaceColor', [0.7 0.8 1])
            hold all
            hc = gsa.cumplot(mat_moment{ij}(:,ik));
            set(hc,'color','k','linewidth',2)
            hold off
            title([endo_prior_restrictions.moment{ij,1},' vs ',endo_prior_restrictions.moment{ij,2},'(',leg,')'],'interpreter','none'),
        end
        indx1 = find(indx_moment(:,ij)==0);
        indx2 = find(indx_moment(:,ij)~=0);
        atitle0=[endo_prior_restrictions.moment{ij,1},' vs ',endo_prior_restrictions.moment{ij,2}, '(', leg,')'];
        fprintf(['%4.1f%% of the ',type,' support matches MOMENT ',atitle0,' inside [%4.1f, %4.1f]\n'],length(indx1)/length(irestrictions)*100,endo_prior_restrictions.moment{ij,4})
        aname=[type '_moment_calib_',endo_prior_restrictions.moment{ij,1},'_vs_',endo_prior_restrictions.moment{ij,2},'_',aleg];
        atitle=[type ' MOMENT Calib: Parameter(s) driving ',endo_prior_restrictions.moment{ij,1},' vs ',endo_prior_restrictions.moment{ij,2}, '(', leg,')'];
        options_mcf.amcf_name = aname;
        options_mcf.amcf_title = atitle;
        options_mcf.beha_title = 'moment restriction';
        options_mcf.nobeha_title = 'NO moment restriction';
        if options_.TeX
            options_mcf.beha_title_latex = 'moment restriction';
            options_mcf.nobeha_title_latex = 'NO moment restriction';
        end
        options_mcf.title = atitle0;
        if ~isempty(indx1) && ~isempty(indx2)
            gsa.monte_carlo_filtering_analysis(xmat, indx1, indx2, options_mcf, M_, options_, bayestopt_, estim_params_);
        end
    end
    for ij=1:nbr_moment_couples
        time_matrix{ij} = sort(time_matrix{ij});
        if length(time_matrix{ij})>1
            if ~options_.nograph
                itmp = (find(plot_indx==ij));
                set(0,'currentfigure',h2);
                subplot(nrow,ncol, ij)
                htmp = plot(time_matrix{ij},[max(moment_matrix{ij})' min(moment_matrix{ij})'],'k--','linewidth',2);
                a=axis;
                delete(htmp);
                tmp=[];
                for ir=1:length(itmp)
                    for it=1:length(endo_prior_restrictions.moment{itmp(ir),3})
                        temp_index = find(time_matrix{ij}==endo_prior_restrictions.moment{itmp(ir),3}(it));
                        tmp(temp_index,:) = endo_prior_restrictions.moment{itmp(ir),4};
                    end
                end
                tmp(isinf(tmp(:,1)),1)=a(3);
                tmp(isinf(tmp(:,2)),2)=a(4);
                hp = patch([time_matrix{ij} time_matrix{ij}(end:-1:1)],[tmp(:,1); tmp(end:-1:1,2)],'b');
                set(hp,'FaceColor',[0.7 0.8 1])
                hold on
                plot(time_matrix{ij},[max(moment_matrix{ij})' min(moment_matrix{ij})'],'k--','linewidth',2)
                plot(time_matrix{ij},moment_median{ij},'k','linewidth',2)
                plot(time_matrix{ij},[moment_distrib{ij}],'k-')
                plot(a(1:2),[0 0],'r')
                hold off
                axis(a)
                box on
                itmp = min(itmp);
                title([endo_prior_restrictions.moment{itmp,1},' vs ',endo_prior_restrictions.moment{itmp,2}],'interpreter','none'),
            end
            if any(iplot_indx.*plot_indx==ij)
                % MCF of the couples with logical AND
                itmp = min(find(plot_indx==ij));
                indx1 = find(indx_moment_matrix(:,ij)==0);
                indx2 = find(indx_moment_matrix(:,ij)~=0);
                leg = num2str(time_matrix{ij}(1));
                leg = [leg '...' num2str(time_matrix{ij}(end))];
                aleg = 'ALL';
                atitle0=[endo_prior_restrictions.moment{itmp,1},' vs ',endo_prior_restrictions.moment{itmp,2}, '(', leg,')'];
                fprintf(['%4.1f%% of the ',type,' support matches MOMENT restrictions ',atitle0,'\n'],length(indx1)/length(irestrictions)*100)
                aname=[type '_moment_calib_',endo_prior_restrictions.moment{itmp,1},'_vs_',endo_prior_restrictions.moment{itmp,2},'_',aleg];
                atitle=[type ' MOMENT Calib: Parameter(s) driving ',endo_prior_restrictions.moment{itmp,1},' vs ',endo_prior_restrictions.moment{itmp,2}, '(', leg,')'];
                options_mcf.amcf_name = aname;
                options_mcf.amcf_title = atitle;
                options_mcf.beha_title = 'moment restriction';
                options_mcf.nobeha_title = 'NO moment restriction';
                if options_.TeX
                    options_mcf.beha_title_latex = 'moment restriction';
                    options_mcf.nobeha_title_latex = 'NO moment restriction';
                end                
                options_mcf.title = atitle0;
                if ~isempty(indx1) && ~isempty(indx2)
                    gsa.monte_carlo_filtering_analysis(xmat, indx1, indx2, options_mcf, M_, options_, bayestopt_, estim_params_);
                end
            end
        end
    end
    if ~options_.nograph
        dyn_saveas(h2,[OutputDirectoryName,filesep,fname_,'_',type,'_moment_restrictions'],options_.nodisplay,options_.graph_format);
        create_TeX_loader(options_,[OutputDirectoryName,filesep,fname_,'_',type,'_moment_restrictions'],[type ' evaluation of moment restrictions'],'moment_restrictions',type,options_.figures.textwidth*min(ij/ncol,1))
    end

    skipline()
end
return

function []=create_TeX_loader(options,figpath,caption,label_name,label_type,scale_factor)
if options.TeX && any(strcmp('eps',cellstr(options.graph_format)))
    fidTeX = fopen([figpath '.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by map_calibration.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s}\n',scale_factor,strrep(figpath,'\','/'));
    fprintf(fidTeX,'\\caption{%s.}',caption);
    fprintf(fidTeX,'\\label{Fig:%s:%s}\n',label_name,label_type);
    fprintf(fidTeX,'\\end{figure}\n\n');
    fprintf(fidTeX,'%% End Of TeX file. \n');
    fclose(fidTeX);
end
