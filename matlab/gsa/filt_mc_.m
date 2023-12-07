function [rmse_MC, ixx] = filt_mc_(OutDir,options_gsa_,dataset_,dataset_info,M_,oo_,options_,bayestopt_,estim_params_)
% [rmse_MC, ixx] = filt_mc_(OutDir,options_gsa_,dataset_,dataset_info,M_,oo_,options_,bayestopt_,estim_params_
% Inputs:
%  - OutputDirectoryName [string]       name of the output directory
%  - options_gsa_        [structure]    GSA options
%  - dataset_            [dseries]      object storing the dataset
%  - dataset_info        [structure]    storing informations about the sample.
%  - M_                  [structure]    Matlab's structure describing the model
%  - oo_                 [structure]    storing the results
%  - options_            [structure]    Matlab's structure describing the current options
%  - bayestopt_          [structure]    describing the priors
%  - estim_params_       [structure]    characterizing parameters to be estimated
%
% Outputs:
%  - rmse_MC             [double]       RMSE by nvar matrix of the RMSEs
%  - ixx                 [double]       RMSE by nvar matrix of sorting
%                                       indices (descending order of RMSEs)

% inputs (from opt_gsa structure)
% vvarvecm = options_gsa_.var_rmse;
% loadSA   = options_gsa_.load_rmse;
% pfilt    = options_gsa_.pfilt_rmse;
% alpha    = options_gsa_.alpha_rmse;
% alpha2   = options_gsa_.alpha2_rmse;
% istart   = options_gsa_.istart_rmse;
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

vvarvecm = options_gsa_.var_rmse;
if options_.TeX
    vvarvecm_tex = options_gsa_.var_rmse_tex;
else
    vvarvecm_tex = {};
end
loadSA   = options_gsa_.load_rmse;
pfilt    = options_gsa_.pfilt_rmse;
alpha    = options_gsa_.alpha_rmse;
alpha2 = 0;
pvalue   = options_gsa_.alpha2_rmse;
istart   = max(2,options_gsa_.istart_rmse);

fname_ = M_.fname;

skipline(2)
disp('Starting sensitivity analysis')
disp('for the fit of EACH observed series ...')
skipline()
if ~options_.nograph
    disp('Deleting old SA figures...')
    a=dir([OutDir,filesep,'*.*']);
    tmp1='0';
    if options_.opt_gsa.ppost
        tmp='_rmse_post';
    else
        if options_.opt_gsa.pprior
            tmp='_rmse_prior';
        else
            tmp='_rmse_mc';
        end
        if options_gsa_.lik_only
            tmp1 = [tmp,'_post_SA'];
            tmp = [tmp,'_lik_SA'];
        end
    end
    for j=1:length(a)
        if strmatch([fname_,tmp],a(j).name)
            if options_.debug
                disp(a(j).name)
            end
            delete([OutDir,filesep,a(j).name])
        end
        if strmatch([fname_,tmp1],a(j).name)
            if options_.debug
                disp(a(j).name)
            end
            delete([OutDir,filesep,a(j).name])
        end
    end
    disp('done !')
end

nshock=estim_params_.nvx + estim_params_.nvn + estim_params_.ncx + estim_params_.ncn;
npar=estim_params_.np;
if ~isempty(options_.mode_file)
    load(options_.mode_file,'xparam1')
end
if options_.opt_gsa.ppost
    c=load([M_.dname filesep 'Output' filesep fname_,'_mean.mat'],'xparam1');
    xparam1_mean=c.xparam1;
    clear c
elseif ~isempty(options_.mode_file) && exist([M_.dname filesep 'Output' filesep fname_,'_mean.mat'],'file')==2
    c=load([M_.dname filesep 'Output' filesep fname_,'_mean.mat'],'xparam1');
    xparam1_mean=c.xparam1;
    clear c
end

if options_.opt_gsa.ppost
    fnamtmp=[fname_,'_post'];
    DirectoryName = CheckPath('metropolis',M_.dname);
else
    if options_.opt_gsa.pprior
        fnamtmp=[fname_,'_prior'];
        DirectoryName = CheckPath(['gsa' filesep 'prior'],M_.dname);
    else
        fnamtmp=[fname_,'_mc'];
        DirectoryName = CheckPath(['gsa' filesep 'mc'],M_.dname);
    end
end
if loadSA
    tmplist =load([OutDir,filesep,fnamtmp, '.mat'],'vvarvecm');
    if isempty(fieldnames(tmplist))
        disp('WARNING: cannot load results since the list of variables used is not present in the mat file')
        loadSA=0;
    elseif ~isequal(tmplist.vvarvecm,vvarvecm)
        disp('WARNING: cannot load results since the list of variables in the mat file differs from the one requested.')
        loadSA=0;
    end
end
if ~loadSA
    Y = transpose(dataset_.data);
    gend = dataset_.nobs;
    data_index = dataset_info.missing.aindex;
    missing_value = dataset_info.missing.state;
    filfilt = dir([DirectoryName filesep M_.fname '_filter_step_ahead*.mat']);
    filupdate = dir([DirectoryName filesep M_.fname '_update*.mat']);
    filparam = dir([DirectoryName filesep M_.fname '_param*.mat']);
    x=[];
    logpo2=[];
    sto_ys=[];
    for j=1:length(filparam)
        if isempty(strmatch([M_.fname '_param_irf'],filparam(j).name))
            temp=load([DirectoryName filesep filparam(j).name]);
            x=[x; temp.stock];
            logpo2=[logpo2; temp.stock_logpo];
            sto_ys=[sto_ys; temp.stock_ys];
            clear temp;
        end
    end
    nruns=size(x,1);
    nfilt=floor(pfilt*nruns);
    if options_.opt_gsa.ppost || (options_.opt_gsa.ppost==0 && options_.opt_gsa.lik_only==0)
        skipline()
        disp('Computing RMSE''s...')
        jxj=NaN(length(vvarvecm),1);
        js=NaN(length(vvarvecm),1);
        yss=NaN(length(vvarvecm),gend,size(sto_ys,1));
        for i = 1:length(vvarvecm)
            vj = vvarvecm{i};
            jxj(i) = strmatch(vj, M_.endo_names(oo_.dr.order_var), 'exact');
            js(i) = strmatch(vj, M_.endo_names, 'exact');
            yss(i,:,:)=repmat(sto_ys(:,js(i))',[gend,1]);
        end
        if exist('xparam1','var')
            [~,~,~,ahat,~,~,aK] = DsgeSmoother(xparam1,gend,Y,data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_);
            y0 = reshape( squeeze(aK(1,jxj,1:gend)),[gend length(jxj)]);
            yobs = transpose( ahat(jxj,:));
            rmse_mode = sqrt(mean((yobs(istart:end,:)-y0(istart:end,:)).^2));
            r2_mode = 1-sum((yobs(istart:end,:)-y0(istart:end,:)).^2)./sum(yobs(istart:end,:).^2);
        end
        
        y0=-yss;
        nbb=0;
        for j=1:length(filfilt)
            temp=load([DirectoryName filesep M_.fname '_filter_step_ahead',num2str(j),'.mat']);
            nb = size(temp.stock,4);
            y0(:,:,nbb+1:nbb+nb)=y0(:,:,nbb+1:nbb+nb)+reshape(temp.stock(1,js,1:gend,:),[length(js) gend nb]);
            nbb=nbb+nb;
            clear temp;
        end
        yobs=-yss;
        nbb=0;
        for j=1:length(filupdate)
            temp=load([DirectoryName filesep M_.fname '_update',num2str(j),'.mat']);
            nb = size(temp.stock,3);
            yobs(:,:,nbb+1:nbb+nb)=yobs(:,:,nbb+1:nbb+nb)+reshape(temp.stock(js,1:gend,:),[length(js) gend nb]);
            nbb=nbb+nb;
            clear temp;
        end
        rmse_MC=zeros(nruns,length(js));
        r2_MC=zeros(nruns,length(js));
        for j=1:nruns
            rmse_MC(j,:) = sqrt(mean((yobs(:,istart:end,j)'-y0(:,istart:end,j)').^2));
            r2_MC(j,:) = 1-mean((yobs(:,istart:end,j)'-y0(:,istart:end,j)').^2)./mean((yobs(:,istart:end,j)').^2);
        end
        if exist('xparam1_mean','var')
            [~,~,~,ahat,~,~,aK] = DsgeSmoother(xparam1_mean,gend,Y,data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_);
            y0 = reshape( squeeze(aK(1,jxj,1:gend)),[gend length(jxj)]);
            yobs = transpose( ahat(jxj,:));
            rmse_pmean = sqrt(mean((yobs(istart:end,:)-y0(istart:end,:)).^2));
            r2_pmean = 1-mean((yobs(istart:end,:)-y0(istart:end,:)).^2)./mean(yobs(istart:end,:).^2);
        end
        clear stock_filter;
    end
    lnprior=NaN(nruns,1);
    for j=1:nruns
        lnprior(j,1) = priordens(x(j,:)',bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
    end
    likelihood=logpo2(:)-lnprior(:);
    disp('... done!')

    if options_.opt_gsa.ppost
        save([OutDir,filesep,fnamtmp,'.mat'], 'x', 'logpo2', 'likelihood', 'rmse_MC', 'r2_MC', 'vvarvecm')
        if exist('xparam1_mean','var')
            save([OutDir,filesep,fnamtmp, '.mat'], 'rmse_pmean', 'r2_pmean','-append')
        end
        if exist('xparam1','var')
            save([OutDir,filesep,fnamtmp,'.mat'], 'rmse_mode', 'r2_mode','-append')
        end
    else
        if options_.opt_gsa.lik_only
            save([OutDir,filesep,fnamtmp, '.mat'], 'x', 'logpo2','likelihood', '-append')
        else
            save([OutDir,filesep,fnamtmp, '.mat'], 'x', 'logpo2','likelihood', 'rmse_MC', 'r2_MC', 'vvarvecm','-append')
            if exist('xparam1_mean','var')
                save([OutDir,filesep,fnamtmp, '.mat'], 'rmse_pmean', 'r2_pmean','-append')
            end
            if exist('xparam1','var')
                save([OutDir,filesep,fnamtmp,'.mat'], 'rmse_mode', 'r2_mode','-append')
            end
        end
    end
else % loadSA
    if options_.opt_gsa.lik_only && options_.opt_gsa.ppost==0
        load([OutDir,filesep,fnamtmp, '.mat'],'x','logpo2','likelihood');
    else
        load([OutDir,filesep,fnamtmp, '.mat'],'x','logpo2','likelihood','rmse_MC','rmse_mode','rmse_pmean', 'r2_MC', 'vvarvecm', 'r2_mode','r2_pmean');
    end
    lnprior=logpo2(:)-likelihood(:);
    nruns=size(x,1);
    nfilt=floor(pfilt*nruns);
end
% Smirnov tests
nfilt0 = nfilt*ones(length(vvarvecm), 1);
logpo2=logpo2(:);
if ~options_.opt_gsa.ppost
    [~, ipost]=sort(-logpo2);
    [~, ilik]=sort(-likelihood);
end

% visual scatter analysis!
if options_.opt_gsa.ppost
    tmp_title='R2 Posterior:';
    atitle='R2 Posterior:';
    asname='r2_post';
else
    if options_.opt_gsa.pprior
        tmp_title='R2 Prior:';
        atitle='R2 Prior:';
        asname='r2_prior';
    else
        tmp_title='R2 MC:';
        atitle='R2 MC:';
        asname='r2_mc';
    end
end
options_scatter.param_names = vvarvecm;
options_scatter.param_names_tex = vvarvecm_tex;
options_scatter.fname_ = fname_;
options_scatter.OutputDirectoryName = OutDir;
options_scatter.amcf_name = asname;
options_scatter.amcf_title = atitle;
options_scatter.title = tmp_title;
scatter_analysis(r2_MC, x,options_scatter, options_);
% end of visual scatter analysis

if ~options_.opt_gsa.ppost && options_.opt_gsa.lik_only
    if options_.opt_gsa.pprior
        anam='rmse_prior_post';
        atitle='RMSE prior: Log Posterior Kernel';
    else
        anam='rmse_mc_post';
        atitle='RMSE MC: Log Posterior Kernel';
    end
    options_mcf.pvalue_ks = alpha;
    options_mcf.pvalue_corr = pvalue;
    options_mcf.alpha2 = alpha2;
    if options_.TeX
        [pnames,pnames_tex]=get_LaTeX_parameter_names(M_,options_,estim_params_,bayestopt_);
        options_mcf.param_names = pnames;
        options_mcf.param_names_tex = pnames_tex;
    else
        [pnames]=get_LaTeX_parameter_names(M_,options_,estim_params_,bayestopt_);
        options_mcf.param_names = pnames;
        options_mcf.param_names_tex = {};
    end
    options_mcf.fname_ = fname_;
    options_mcf.OutputDirectoryName = OutDir;
    options_mcf.amcf_name = anam;
    options_mcf.amcf_title = atitle;
    options_mcf.title = atitle;
    options_mcf.beha_title = 'better posterior kernel';
    options_mcf.nobeha_title = 'worse posterior kernel';
    mcf_analysis(x, ipost(1:nfilt), ipost(nfilt+1:end), options_mcf, M_, options_, bayestopt_, estim_params_);
    if options_.opt_gsa.pprior
        anam = 'rmse_prior_lik';
        atitle = 'RMSE prior: Log Likelihood Kernel';
    else
        anam='rmse_mc_lik';
        atitle = 'RMSE MC: Log Likelihood Kernel';
    end
    options_mcf.amcf_name = anam;
    options_mcf.amcf_title = atitle;
    options_mcf.title = atitle;
    options_mcf.beha_title = 'better likelihood';
    options_mcf.nobeha_title = 'worse likelihood';
    mcf_analysis(x, ilik(1:nfilt), ilik(nfilt+1:end), options_mcf, M_, options_, bayestopt_, estim_params_);

else
    if options_.opt_gsa.ppost
        rmse_txt=rmse_pmean;
        r2_txt=r2_pmean;
    else
        if options_.opt_gsa.pprior || ~exist('rmse_pmean','var')
            if exist('rmse_mode','var')
                rmse_txt=rmse_mode;
                r2_txt=r2_mode;
            else
                rmse_txt=NaN(1,size(rmse_MC,2));
                r2_txt=NaN(1,size(r2_MC,2));
            end
        else
            rmse_txt=rmse_pmean;
            r2_txt=r2_pmean;
        end
    end
    ixx=NaN(size(rmse_MC,1),length(vvarvecm));
    for i = 1:length(vvarvecm)
        [~, ixx(:,i)] = sort(rmse_MC(:,i));
    end
    PP = ones(npar+nshock, length(vvarvecm));
    PPV = ones(length(vvarvecm), length(vvarvecm), npar+nshock);
    SS = zeros(npar+nshock, length(vvarvecm));
    for j = 1:npar+nshock
        for i = 1:length(vvarvecm)
            [~, P] = smirnov(x(ixx(nfilt0(i)+1:end,i),j),x(ixx(1:nfilt0(i),i),j), alpha);
            [H1] = smirnov(x(ixx(nfilt0(i)+1:end,i),j),x(ixx(1:nfilt0(i),i),j),alpha,1);
            [H2] = smirnov(x(ixx(nfilt0(i)+1:end,i),j),x(ixx(1:nfilt0(i),i),j),alpha,-1);
            if H1==0 && H2==0
                SS(j,i)=1;
            elseif H1==0
                SS(j,i)=-1;
            else
                SS(j,i)=0;
            end
            PP(j,i)=P;
        end
        for i = 1:length(vvarvecm)
            for l = 1:length(vvarvecm)
                if l~=i && PP(j,i)<alpha && PP(j,l)<alpha
                    [~,P] = smirnov(x(ixx(1:nfilt0(i),i),j),x(ixx(1:nfilt0(l),l),j), alpha);
                    PPV(i,l,j) = P;
                elseif l==i
                    PPV(i,l,j) = PP(j,i);
                end
            end
        end
    end
    if ~options_.nograph
        ifig=0;
        for i=1:length(vvarvecm)
            if options_.opt_gsa.ppost
                temp_name='RMSE Posterior: Log Prior';
            else
                if options_.opt_gsa.pprior
                    temp_name='RMSE Prior: Log Prior';
                else
                    temp_name='RMSE MC: Log Prior';
                end
            end
            if mod(i,9)==1
                ifig=ifig+1;
                hh_fig=dyn_figure(options_.nodisplay,'name',[temp_name,' ',int2str(ifig)]);
            end
            subplot(3,3,i-9*(ifig-1))
            h=cumplot(lnprior(ixx(1:nfilt0(i),i)));
            set(h,'color','blue','linewidth',2)
            hold on, h=cumplot(lnprior);
            set(h,'color','k','linewidth',1)
            h=cumplot(lnprior(ixx(nfilt0(i)+1:end,i)));
            set(h,'color','red','linewidth',2)
            title(vvarvecm{i},'interpreter','none')
            if mod(i,9)==0 || i==length(vvarvecm)
                if ~isoctave
                    annotation('textbox', [0.1,0,0.35,0.05],'String', 'Log-prior for BETTER R2','Color','Blue','horizontalalignment','center');
                    annotation('textbox', [0.55,0,0.35,0.05],'String', 'Log-prior for WORSE R2', 'Color','Red','horizontalalignment','center');
                end
                if options_.opt_gsa.ppost
                    dyn_saveas(hh_fig,[OutDir filesep fname_ '_rmse_post_lnprior',int2str(ifig)],options_.nodisplay,options_.graph_format);
                    if options_.TeX
                        create_TeX_loader(options_,[OutDir '/' fname_ '_rmse_post_lnprior',int2str(ifig)],ifig,[temp_name,' ',int2str(ifig)],'rmse_post_lnprior',options_.figures.textwidth*min((i-9*(ifig-1))/3,1))
                    end
                else
                    if options_.opt_gsa.pprior
                        dyn_saveas(hh_fig,[OutDir filesep fname_ '_rmse_prior_lnprior',int2str(ifig) ],options_.nodisplay,options_.graph_format);
                        if options_.TeX
                            create_TeX_loader(options_,[OutDir '/' fname_ '_rmse_prior_lnprior',int2str(ifig)],ifig,[temp_name,' ',int2str(ifig)],'rmse_prior_lnprior',options_.figures.textwidth*min((i-9*(ifig-1))/3,1))
                        end
                    else
                        dyn_saveas(hh_fig,[OutDir filesep fname_ '_rmse_mc_lnprior',int2str(ifig) ],options_.nodisplay,options_.graph_format);
                        if options_.TeX
                            create_TeX_loader(options_,[OutDir '/' fname_ '_rmse_mc_lnprior',int2str(ifig)],ifig,[temp_name,' ',int2str(ifig)],'rmse_mc_lnprior',options_.figures.textwidth*min((i-9*(ifig-1))/3,1))
                        end
                    end
                end
            end
        end
        ifig=0;
        for i=1:length(vvarvecm)
            if options_.opt_gsa.ppost
                temp_name='RMSE Posterior: Log Likelihood';
            else
                if options_.opt_gsa.pprior
                    temp_name='RMSE Prior: Log Likelihood';
                else
                    temp_name='RMSE MC: Log Likelihood';
                end
            end
            if mod(i,9)==1
                ifig=ifig+1;
                hh_fig = dyn_figure(options_.nodisplay,'Name',[temp_name,' ',int2str(ifig)]);
            end
            subplot(3,3,i-9*(ifig-1))
            h=cumplot(likelihood(ixx(1:nfilt0(i),i)));
            set(h,'color','blue','linewidth',2)
            hold on, h=cumplot(likelihood);
            set(h,'color','k','linewidth',1)
            h=cumplot(likelihood(ixx(nfilt0(i)+1:end,i)));
            set(h,'color','red','linewidth',2)
            title(vvarvecm{i},'interpreter','none')
            if options_.opt_gsa.ppost==0
                set(gca,'xlim',[min( likelihood(ixx(1:nfilt0(i),i)) ) max( likelihood(ixx(1:nfilt0(i),i)) )])
            end
            if mod(i,9)==0 || i==length(vvarvecm)
                if ~isoctave
                    annotation('textbox', [0.1,0,0.35,0.05],'String', 'Log-likelihood for BETTER R2','Color','Blue','horizontalalignment','center');
                    annotation('textbox', [0.55,0,0.35,0.05],'String', 'Log-likelihood for WORSE R2', 'Color','Red','horizontalalignment','center');
                end
                if options_.opt_gsa.ppost
                    dyn_saveas(hh_fig,[OutDir filesep fname_ '_rmse_post_lnlik',int2str(ifig) ],options_.nodisplay,options_.graph_format);
                    if options_.TeX
                        create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_post_lnprior',int2str(ifig)],ifig,[temp_name,' ',int2str(ifig)],'rmse_post_lnprior',options_.figures.textwidth*min((i-9*(ifig-1))/3,1));
                    end
                else
                    if options_.opt_gsa.pprior
                        dyn_saveas(hh_fig,[OutDir filesep fname_ '_rmse_prior_lnlik',int2str(ifig)],options_.nodisplay,options_.graph_format);
                        if options_.TeX
                            create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_prior_lnlik',int2str(ifig)],ifig,[temp_name,' ',int2str(ifig)],'rmse_prior_lnlik',options_.figures.textwidth*min((i-9*(ifig-1))/3,1));
                        end
                    else
                        dyn_saveas(hh_fig,[OutDir filesep fname_ '_rmse_mc_lnlik',int2str(ifig) ],options_.nodisplay,options_.graph_format);
                        if options_.TeX
                            create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_mc_lnlik',int2str(ifig) ],ifig,[temp_name,' ',int2str(ifig)],'rmse_mc_lnlik',options_.figures.textwidth*min((i-9*(ifig-1))/3,1));
                        end
                    end
                end
            end
        end
        ifig=0;
        for i=1:length(vvarvecm)
            if options_.opt_gsa.ppost
                temp_name='RMSE Posterior: Log Posterior';
            else
                if options_.opt_gsa.pprior
                    temp_name='RMSE Prior: Log Posterior';
                else
                    temp_name='RMSE MC: Log Posterior';
                end
            end
            if mod(i,9)==1
                ifig=ifig+1;
                hh_fig = dyn_figure(options_.nodisplay,'Name',[temp_name,' ',int2str(ifig)]);
            end
            subplot(3,3,i-9*(ifig-1))
            h=cumplot(logpo2(ixx(1:nfilt0(i),i)));
            set(h,'color','blue','linewidth',2)
            hold on, h=cumplot(logpo2);
            set(h,'color','k','linewidth',1)
            h=cumplot(logpo2(ixx(nfilt0(i)+1:end,i)));
            set(h,'color','red','linewidth',2)
            title(vvarvecm{i},'interpreter','none')
            if options_.opt_gsa.ppost==0
                set(gca,'xlim',[min( logpo2(ixx(1:nfilt0(i),i)) ) max( logpo2(ixx(1:nfilt0(i),i)) )])
            end
            if mod(i,9)==0 || i==length(vvarvecm)
                if ~isoctave
                    annotation('textbox', [0.1,0,0.35,0.05],'String', 'Log-posterior for BETTER R2','Color','Blue','horizontalalignment','center');
                    annotation('textbox', [0.55,0,0.35,0.05],'String', 'Log-posterior for WORSE R2', 'Color','Red','horizontalalignment','center');
                end
                if options_.opt_gsa.ppost
                    dyn_saveas(hh_fig,[OutDir filesep fname_ '_rmse_post_lnpost',int2str(ifig) ],options_.nodisplay,options_.graph_format);
                    if options_.TeX
                        create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_post_lnpost',int2str(ifig) ],ifig,[temp_name,' ',int2str(ifig)],'rmse_post_lnpost',options_.figures.textwidth*min((i-9*(ifig-1))/3,1));
                    end
                else
                    if options_.opt_gsa.pprior
                        dyn_saveas(hh_fig,[OutDir filesep fname_ '_rmse_prior_lnpost',int2str(ifig)],options_.nodisplay,options_.graph_format);
                        if options_.TeX
                            create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_prior_lnpost',int2str(ifig)],ifig,[temp_name,' ',int2str(ifig)],'rmse_prior_lnpost',options_.figures.textwidth*min((i-9*(ifig-1))/3,1));
                        end
                    else
                        dyn_saveas(hh_fig,[OutDir filesep fname_ '_rmse_mc_lnpost',int2str(ifig)],options_.nodisplay,options_.graph_format);
                        if options_.TeX
                            create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_mc_lnpost',int2str(ifig)],ifig,[temp_name,' ',int2str(ifig)],'rmse_mc_lnpost',options_.figures.textwidth*min((i-9*(ifig-1))/3,1));
                        end
                    end
                end
            end
        end
    end
    if options_.TeX
        [pnames,pnames_tex]=get_LaTeX_parameter_names(M_,options_,estim_params_,bayestopt_);
        param_names = pnames;
        param_names_tex = pnames_tex;
    else
        [pnames]=get_LaTeX_parameter_names(M_,options_,estim_params_,bayestopt_);
        param_names = pnames;
        param_names_tex = {};
    end
    skipline()
    title_string='RMSE over the MC sample:';
    data_mat=[min(rmse_MC)' max(rmse_MC)'];
    headers={'Variable'; 'min yr RMSE'; 'max yr RMSE'};
    dyntable(options_, title_string, headers, vvarvecm, data_mat, 0, 15, 5);
    if options_.TeX
        headers_tex = {'\text{Variable}'; '\text{min yr RMSE}'; '\text{max yr RMSE}'};
        dyn_latex_table(M_, options_, title_string, 'RMSE_MC', headers_tex, vvarvecm_tex, data_mat, 0, 15, 5);
    end
    invar = find( std(rmse_MC)./mean(rmse_MC)<=0.0001 );
    if ~isempty(invar)
        skipline(2)
        disp('RMSE is not varying significantly over the MC sample for the following variables:')
        disp(vvarvecm{invar})
        disp('These variables are excluded from SA')
        disp('[Unless you treat these series as exogenous, there is something wrong in your estimation !]')
    end
    vvarvecm0=vvarvecm;
    ivar = find( std(rmse_MC)./mean(rmse_MC)>0.0001 );
    vvarvecm = vvarvecm(ivar);
    rmse_MC = rmse_MC(:,ivar);
    skipline()
    disp(['Sample filtered the ',num2str(pfilt*100),'% best RMSE''s for each observed series ...' ])
    skipline(2)
    disp('RMSE ranges after filtering:')
    title_string='RMSE ranges after filtering:';
    if options_.opt_gsa.ppost==0 && options_.opt_gsa.pprior
        headers = {'Variable'; 'min'; 'max'; 'min'; 'max'; 'posterior mode'};
        headers_tex = {'\text{Variable}'; '\text{min}'; '\text{max}'; '\text{min}'; '\text{max}'; '\text{posterior mode}'};
    else
        headers = {'Variable'; 'min'; 'max'; 'min'; 'max'; 'posterior mean'};
        headers_tex = {'\text{Variable}'; '\text{min}'; '\text{max}'; '\text{min}'; '\text{max}'; '\text{posterior mean}'};
    end
    data_mat=NaN(length(vvarvecm),5);
    for j = 1:length(vvarvecm)
        data_mat(j,:)=[min(rmse_MC(ixx(1:nfilt0(j),j),j)) ...
                       max(rmse_MC(ixx(1:nfilt0(j),j),j))  ...
                       min(rmse_MC(ixx(nfilt0(j)+1:end,j),j)) ...
                       max(rmse_MC(ixx(nfilt0(j)+1:end,j),j)) ...
                       rmse_txt(j)];
    end
    %get formatting for additional header line
    val_width = 15;
    val_precis = 5;
    label_width = max(cellofchararraymaxlength(vertcat(headers{1}, vvarvecm))+2, 0);
    label_format_leftbound  = sprintf('%%-%ds', label_width);
    if all(~isfinite(data_mat))
        values_length = 4;
    else
        values_length = max(ceil(max(max(log10(abs(data_mat(isfinite(data_mat))))))),1)+val_precis+1;
    end
    if any(data_mat < 0) %add one character for minus sign
        values_length = values_length+1;
    end
    headers_length = cellofchararraymaxlength(headers(2:end));
    if ~isempty(val_width)
        val_width = max(max(headers_length,values_length)+2, val_width);
    else
        val_width = max(headers_length, values_length)+2;
    end
    header_string_format  = sprintf('%%%ds',val_width);
    if options_.opt_gsa.ppost==0 && options_.opt_gsa.pprior
        optional_header=sprintf([label_format_leftbound,header_string_format,header_string_format,header_string_format,header_string_format],'','',['best ',num2str(pfilt*100),'% filtered'],'','remaining 90%');
    else
        optional_header=sprintf([label_format_leftbound,header_string_format,header_string_format,header_string_format,header_string_format],'','','best  filtered','','remaining');
    end
    dyntable(options_, title_string, headers, vvarvecm, data_mat, 0, val_width, val_precis,optional_header);
    if options_.TeX
        if options_.opt_gsa.ppost==0 && options_.opt_gsa.pprior
            optional_header={[' & \multicolumn{2}{c}{best ',num2str(pfilt*100),' filtered} & \multicolumn{2}{c}{remaining 90\%}\\']};
        else
            optional_header={' & \multicolumn{2}{c}{best filtered} & \multicolumn{2}{c}{remaining}\\'};
        end
        dyn_latex_table(M_, options_, title_string, 'RMSE_ranges_after_filtering', headers_tex, vvarvecm_tex, data_mat, 0, val_width, val_precis, optional_header);
    end
    % R2 table
    vvarvecm=vvarvecm0;
    skipline()
    title_string='R2 over the MC sample:';
    data_mat=[min(r2_MC)' max(r2_MC)'];
    headers = {'Variable'; 'min yr R2'; 'max yr R2'};
    dyntable(options_, title_string, headers, vvarvecm, data_mat, 0, 15, 5);
    if options_.TeX
        headers_tex = {'\text{Variable}'; '\text{min yr R2}'; '\text{max yr R2}'};
        dyn_latex_table(M_, options_, title_string, 'R2_MC', headers_tex, vvarvecm_tex, data_mat, 0, 15, 5);
    end
    r2_MC=r2_MC(:,ivar);
    vvarvecm=vvarvecm(ivar);
    skipline()
    disp(['Sample filtered the ',num2str(pfilt*100),'% best R2''s for each observed series ...' ])
    skipline()
    disp('R2 ranges after filtering:')
    title_string='R2 ranges after filtering:';
    if options_.opt_gsa.ppost==0 && options_.opt_gsa.pprior
        headers = {'Variable'; 'min'; 'max'; 'min'; 'max'; 'posterior mode'};
        headers_tex = {'\text{Variable}'; '\text{min}'; '\text{max}'; '\text{min}'; '\text{max}'; '\text{posterior mode}'};
    else
        headers = {'Variable'; 'min'; 'max'; 'min'; 'max'; 'posterior mean'};
        headers_tex = {'\text{Variable}'; '\text{min}'; '\text{max}'; '\text{min}'; '\text{max}'; '\text{posterior mean}'};
    end
    data_mat=NaN(length(vvarvecm),5);
    for j = 1:length(vvarvecm)
        data_mat(j,:)=[min(r2_MC(ixx(1:nfilt0(j),j),j)) ...
                       max(r2_MC(ixx(1:nfilt0(j),j),j))  ...
                       min(r2_MC(ixx(nfilt0(j)+1:end,j),j)) ...
                       max(r2_MC(ixx(nfilt0(j)+1:end,j),j)) ...
                       r2_txt(j)];
    end
    %get formatting for additional header line
    val_width = 15;
    val_precis = 5;
    label_width = max(cellofchararraymaxlength(vertcat(headers{1}, vvarvecm))+2, 0);
    label_format_leftbound  = sprintf('%%-%ds', label_width);
    if all(~isfinite(data_mat))
        values_length = 4;
    else
        values_length = max(ceil(max(max(log10(abs(data_mat(isfinite(data_mat))))))),1)+val_precis+1;
    end
    if any(data_mat < 0) %add one character for minus sign
        values_length = values_length+1;
    end
    headers_length = cellofchararraymaxlength(headers(2:end));
    if ~isempty(val_width)
        val_width = max(max(headers_length, values_length)+2, val_width);
    else
        val_width = max(headers_length, values_length)+2;
    end
    header_string_format  = sprintf('%%%ds',val_width);

    if options_.opt_gsa.ppost==0 && options_.opt_gsa.pprior
        optional_header = sprintf([label_format_leftbound,header_string_format,header_string_format,header_string_format,header_string_format],'','',['best ',num2str(pfilt*100),'% filtered'],'','remaining 90%');
    else
        optional_header = sprintf([label_format_leftbound,header_string_format,header_string_format,header_string_format,header_string_format],'','','best  filtered','','remaining');
    end
    dyntable(options_, title_string, headers, vvarvecm, data_mat, 0, val_width, val_precis, optional_header);
    if options_.TeX
        if ~options_.opt_gsa.ppost && options_.opt_gsa.pprior
            optional_header = {[' & \multicolumn{2}{c}{best ',num2str(pfilt*100),' filtered} & \multicolumn{2}{c}{remaining 90\%}\\']};
        else
            optional_header = {' & \multicolumn{2}{c}{best filtered} & \multicolumn{2}{c}{remaining}\\'};
        end
        dyn_latex_table(M_, options_, title_string, 'R2_ranges_after_filtering', headers_tex, vvarvecm_tex, data_mat, 0, val_width, val_precis, optional_header);
    end
    %  R2 table
    SP = zeros(npar+nshock, length(vvarvecm));
    for j = 1:length(vvarvecm)
        ns=find(PP(:,j)<alpha);
        SP(ns,j)=ones(size(ns));
        SS(:,j)=SS(:,j).*SP(:,j);
    end
    nsp=NaN(npar+nshock,1);
    for j=1:npar+nshock
        nsp(j)=length(find(SP(j,:)));
    end
    snam0=param_names(nsp==0);
    snam1=param_names(nsp==1);
    snam2=param_names(nsp>1);
    nsnam=(find(nsp>1));
    skipline(2)
    disp('These parameters do not affect significantly the fit of ANY observed series:')
    disp(char(snam0))
    skipline()
    disp('These parameters affect ONE single observed series:')
    disp(char(snam1))
    skipline()
    disp('These parameters affect MORE THAN ONE observed series: trade off exists!')
    disp(char(snam2))
    pnam=bayestopt_.name;
    % plot trade-offs
    if ~options_.nograph
        a00=jet(length(vvarvecm));
        if options_.opt_gsa.ppost
            temp_name='RMSE Posterior Tradeoffs:';
            atitle='RMSE Posterior Map:';
            asname='rmse_post';
        else
            if options_.opt_gsa.pprior
                temp_name='RMSE Prior Tradeoffs:';
                atitle='RMSE Prior Map:';
                asname='rmse_prior';
            else
                temp_name='RMSE MC Tradeoffs:';
                atitle='RMSE MC Map:';
                asname='rmse_mc';
            end
        end
        % now I plot by observed variables
        options_mcf.pvalue_ks = alpha;
        options_mcf.pvalue_corr = pvalue;
        options_mcf.alpha2 = alpha2;
        options_mcf.param_names = param_names;
        options_mcf.param_names_tex = param_names_tex;
        options_mcf.fname_ = fname_;
        options_mcf.OutputDirectoryName = OutDir;
        for iy = 1:length(vvarvecm)
            options_mcf.amcf_name = [asname '_' vvarvecm{iy} '_map' ];
            options_mcf.amcf_title = [atitle ' ' vvarvecm{iy}];
            if options_.TeX
                options_mcf.beha_title = ['better fit of ' vvarvecm_tex{iy}];
                options_mcf.nobeha_title = ['worse fit of ' vvarvecm_tex{iy}];
            else
                options_mcf.beha_title = ['better fit of ' vvarvecm{iy}];
                options_mcf.nobeha_title = ['worse fit of ' vvarvecm{iy}];
            end
            options_mcf.title = ['the fit of ' vvarvecm{iy}];
            mcf_analysis(x, ixx(1:nfilt0(iy),iy), ixx(nfilt0(iy)+1:end,iy), options_mcf, M_, options_, bayestopt_, estim_params_);
        end
        for iy = 1:length(vvarvecm)
            ipar = find(any(squeeze(PPV(iy,:,:))<alpha));
            for ix=1:ceil(length(ipar)/5)
                hh_fig = dyn_figure(options_.nodisplay,'name',[temp_name,' observed variable ', vvarvecm{iy}]);
                for j=1+5*(ix-1):min(length(ipar),5*ix)
                    subplot(2,3,j-5*(ix-1))
                    h0=cumplot(x(:,ipar(j)));
                    set(h0,'color',[0 0 0])
                    hold on,
                    iobs=find(squeeze(PPV(iy,:,ipar(j)))<alpha);
                    for i = 1:length(vvarvecm)
                        if any(iobs==i) || i==iy
                            h0=cumplot(x(ixx(1:nfilt0(i),i),ipar(j)));
                            if ~isoctave
                                hcmenu = uicontextmenu;
                                uimenu(hcmenu,'Label',vvarvecm{i});
                                set(h0,'uicontextmenu',hcmenu)
                            end
                        else
                            h0=cumplot(x(ixx(1:nfilt0(i),i),ipar(j))*NaN);
                        end
                        set(h0,'color',a00(i,:),'linewidth',2)
                    end
                    ydum=get(gca,'ylim');
                    if exist('xparam1','var')
                        xdum=xparam1(ipar(j));
                        h1=plot([xdum xdum],ydum);
                        set(h1,'color',[0.85 0.85 0.85],'linewidth',2)
                    end
                    xlabel('')
                    title([pnam{ipar(j)}],'interpreter','none')
                end
                if isoctave
                    legend(vertcat('base',vvarvecm),'location','eastoutside');
                else
                    h0=legend(vertcat('base',vvarvecm));
                    set(h0,'fontsize',6,'position',[0.7 0.1 0.2 0.3],'interpreter','none');
                end
                if options_.opt_gsa.ppost
                    dyn_saveas(hh_fig,[ OutDir filesep fname_ '_rmse_post_' vvarvecm{iy} '_' int2str(ix)],options_.nodisplay,options_.graph_format);
                    if options_.TeX
                        create_TeX_loader(options_,[ OutDir filesep fname_ '_rmse_post_' vvarvecm{iy} '_' int2str(ix)],ix,[temp_name,' observed variable $',vvarvecm_tex{iy} '$'],['rmse_post_' vvarvecm{iy}],1)
                    end
                else
                    if options_.opt_gsa.pprior
                        dyn_saveas(hh_fig,[OutDir filesep fname_ '_rmse_prior_' vvarvecm{iy} '_' int2str(ix) ],options_.nodisplay,options_.graph_format);
                        if options_.TeX
                            create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_prior_' vvarvecm{iy} '_' int2str(ix) ],ix,[temp_name,' observed variable $',vvarvecm_tex{iy} '$'],['rmse_prior_' vvarvecm{iy}],1)
                        end
                    else
                        dyn_saveas(hh_fig,[OutDir filesep fname_ '_rmse_mc_' vvarvecm{iy} '_' int2str(ix)],options_.nodisplay,options_.graph_format);
                        if options_.TeX
                            create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_mc_' vvarvecm{iy} '_' int2str(ix)],ix,[temp_name,' observed variable $',vvarvecm_tex{iy} '$'],['rmse_mc_' vvarvecm{iy}],1)
                        end
                    end
                end
            end
        end
        % now I plot by individual parameters
        for ix=1:ceil(length(nsnam)/5)
            hh_fig = dyn_figure(options_.nodisplay,'name',[temp_name,' estimated params and shocks ',int2str(ix)]);
            for j=1+5*(ix-1):min(size(snam2,1),5*ix)
                subplot(2,3,j-5*(ix-1))
                h0=cumplot(x(:,nsnam(j)));
                set(h0,'color',[0 0 0])
                hold on,
                npx=find(SP(nsnam(j),:)==0);
                for i = 1:length(vvarvecm)
                    if any(npx==i)
                        h0=cumplot(x(ixx(1:nfilt0(i),i),nsnam(j))*NaN);
                    else
                        h0=cumplot(x(ixx(1:nfilt0(i),i),nsnam(j)));
                        if ~isoctave
                            hcmenu = uicontextmenu;
                            uimenu(hcmenu,'Label', vvarvecm{i});
                            set(h0,'uicontextmenu',hcmenu)
                        end
                    end
                    set(h0,'color',a00(i,:),'linewidth',2)
                end
                ydum=get(gca,'ylim');
                if exist('xparam1','var')
                    xdum=xparam1(nsnam(j));
                    h1=plot([xdum xdum],ydum);
                    set(h1,'color',[0.85 0.85 0.85],'linewidth',2)
                end
                xlabel('')
                title([pnam{nsnam(j)}],'interpreter','none')
            end
            %subplot(3,2,6)
            if isoctave
                legend(vertcat('base',vvarvecm),'location','eastoutside');
            else
                h0=legend(vertcat('base',vvarvecm));
                set(h0,'fontsize',6,'position',[0.7 0.1 0.2 0.3],'interpreter','none');
            end
            if options_.opt_gsa.ppost
                dyn_saveas(hh_fig,[ OutDir filesep fname_ '_rmse_post_params_' int2str(ix)],options_.nodisplay,options_.graph_format);
                if options_.TeX
                    create_TeX_loader(options_,[ OutDir filesep fname_ '_rmse_post_params_' int2str(ix)],ix,[temp_name,' estimated params and shocks ',int2str(ix)],'rmse_post_params',1)
                end
            else
                if options_.opt_gsa.pprior
                    dyn_saveas(hh_fig,[OutDir filesep fname_ '_rmse_prior_params_' int2str(ix) ],options_.nodisplay,options_.graph_format);
                    if options_.TeX
                        create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_prior_params_' int2str(ix) ],ix,[temp_name,' estimated params and shocks ',int2str(ix)],'rmse_prior_params',1)
                    end
                else
                    dyn_saveas(hh_fig,[OutDir filesep fname_ '_rmse_mc_params_' int2str(ix)],options_.nodisplay,options_.graph_format);
                    if options_.TeX
                        create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_mc_params_' int2str(ix)],ix,[temp_name,' estimated params and shocks ',int2str(ix)],'rmse_mc_params',1)
                    end
                end
            end
        end
    end
end

function []=create_TeX_loader(options_,figpath,label_number,caption,label_name,scale_factor)
if nargin<6
    scale_factor=1;
end
if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fidTeX = fopen([figpath '.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by filt_mc_.m (Dynare).\n');
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

function [pnames,pnames_tex]=get_LaTeX_parameter_names(M_,options_,estim_params_,bayestopt_)
np=size(bayestopt_.name,1);
pnames=cell(np,1);
pnames_tex=cell(np,1);
for ii=1:length(bayestopt_.name)
    if options_.TeX
        [param_name_temp, param_name_tex_temp]= get_the_name(ii,options_.TeX,M_,estim_params_,options_);
        pnames_tex{ii,1} = param_name_tex_temp;
        pnames{ii,1} = param_name_temp;
    else
        param_name_temp = get_the_name(ii,options_.TeX,M_,estim_params_,options_);
        pnames{ii,1} = param_name_temp;
    end
end
