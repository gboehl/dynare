function map_identification(OutputDirectoryName,opt_gsa,M_,oo_,options_,estim_params_,bayestopt_)
% map_identification(OutputDirectoryName,opt_gsa,M_,oo_,options_,estim_params_,bayestopt_)
% Inputs
%  - OutputDirectoryName [string]    name of the output directory
%  - opt_gsa             [structure]     GSA options structure
%  - M_                  [structure]     Matlab's structure describing the model
%  - oo_                 [structure]     Matlab's structure describing the results
%  - options_            [structure]     Matlab's structure describing the current options
%  - estim_params_       [structure]     characterizing parameters to be estimated
%  - bayestopt_          [structure]     describing the priors

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

fname_ = M_.fname;
dr=oo_.dr;
nliv   = opt_gsa.morris_nliv;
itrans = opt_gsa.trans_ident;

np = size(estim_params_.param_vals,1);

pnames = M_.param_names(estim_params_.param_vals(:,1));
if opt_gsa.pprior
    filetoload=[OutputDirectoryName '/' fname_ '_prior'];
else
    filetoload=[OutputDirectoryName '/' fname_ '_mc'];
end
load(filetoload,'lpmat','lpmat0','istable','T','yys')
if ~isempty(lpmat0)
    lpmatx=lpmat0(istable,:);
else
    lpmatx=[];
end
Nsam = size(lpmat,1);
nshock = size(lpmat0,2);
npT = np+nshock;

fname_ = M_.fname;

if opt_gsa.load_ident_files==0
    mss = yys(bayestopt_.mfys,:);
    mss = gsa.teff(mss(:,istable),Nsam,istable);
    yys = gsa.teff(yys(dr.order_var,istable),Nsam,istable);
    if exist('T','var')
        [vdec, cc, ac] = gsa.monte_carlo_moments(T, lpmatx, dr, M_, options_, estim_params_);
    else
        return
    end

    if opt_gsa.morris==2
        pdraws = identification.run(M_,oo_,options_,bayestopt_,estim_params_,options_.options_ident,[lpmatx lpmat(istable,:)]);
        if ~isempty(pdraws) && max(max(abs(pdraws-[lpmatx lpmat(istable,:)])))==0
            disp(['Sample check OK. Largest difference: ', num2str(max(max(abs(pdraws-[lpmatx lpmat(istable,:)]))))]),
            clear pdraws;
        end
        clear GAM gas
    end
    if opt_gsa.morris~=1 && M_.exo_nbr>1
        ifig=0;
        for j=1:M_.exo_nbr
            if mod(j,6)==1
                hh_fig=dyn_figure(options_.nodisplay,'name','Variance decomposition shocks');
                ifig=ifig+1;
                iplo=0;
            end
            iplo=iplo+1;
            subplot(2,3,iplo)
            gsa.boxplot(squeeze(vdec(:,j,:))',[],'.',[],10)
            set(gca,'xticklabel',' ','fontsize',10,'xtick',1:size(options_.varobs,1))
            set(gca,'xlim',[0.5 size(options_.varobs,1)+0.5])
            set(gca,'ylim',[-2 102])
            for ip=1:size(options_.varobs,1)
                if options_.TeX
                    text(ip,-4,deblank(opt_gsa.varobs_tex(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','latex')
                else
                    text(ip,-4,deblank(options_.varobs(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
                end
            end
            xlabel(' ')
            ylabel(' ')
            title(M_.exo_names{j},'interpreter','none')
            if mod(j,6)==0 || j==M_.exo_nbr
                dyn_saveas(hh_fig,[OutputDirectoryName,'/',fname_,'_vdec_exo_',int2str(ifig)],options_.nodisplay,options_.graph_format);
                create_TeX_loader(options_,[OutputDirectoryName,'/',fname_,'_vdec_exo_',int2str(ifig)],ifig,'Variance decomposition shocks','vdec_exo',options_.figures.textwidth*min(iplo/3,1))
            end
        end
    end
    for j=1:size(cc,1)
        cc(j,j,:)=gsa.standardize_columns(squeeze(log(cc(j,j,:))))./2;
    end
    [vdec, ~, ir_vdec, ic_vdec] = gsa.teff(vdec,Nsam,istable);
    [cc, ~, ir_cc, ic_cc] = gsa.teff(cc,Nsam,istable);
    [ac, ~, ir_ac, ic_ac] = gsa.teff(ac,Nsam,istable);

    nc1= size(T,2);
    endo_nbr = M_.endo_nbr;
    nstatic = M_.nstatic;
    nspred = M_.nspred;
    iv = (1:endo_nbr)';
    ic = [ nstatic+(1:nspred) endo_nbr+(1:size(dr.ghx,2)-nspred) ]';

    dr.ghx = T(:, 1:(nc1-M_.exo_nbr),1);
    dr.ghu = T(:, (nc1-M_.exo_nbr+1):end, 1);
    [Aa,Bb] = kalman_transition_matrix(dr,iv,ic);
    A = zeros(size(Aa,1),size(Aa,2)+size(Aa,1),length(istable));
    if ~isempty(lpmatx)
        M_=gsa.set_shocks_param(M_,estim_params_,lpmatx(1,:));
    end
    A(:,:,1)=[Aa, triu(Bb*M_.Sigma_e*Bb')];
    for j=2:length(istable)
        dr.ghx = T(:, 1:(nc1-M_.exo_nbr),j);
        dr.ghu = T(:, (nc1-M_.exo_nbr+1):end, j);
        [Aa,Bb] = kalman_transition_matrix(dr, iv, ic);
        if ~isempty(lpmatx)
            M_=gsa.set_shocks_param(M_,estim_params_,lpmatx(j,:));
        end
        A(:,:,j)=[Aa, triu(Bb*M_.Sigma_e*Bb')];
    end
    clear T
    clear lpmatx

    [yt, j0]=gsa.teff(A,Nsam,istable);
    yt = [yys yt];
    if opt_gsa.morris==2
        clear TAU A
    else
        clear A
    end
    save([OutputDirectoryName,'/',fname_,'_main_eff.mat'],'ac','cc','vdec','yt','mss')
else %load identification files
    load([OutputDirectoryName,'/',fname_,'_main_eff.mat'],'ac','cc','vdec','yt','mss')
end

if opt_gsa.morris==1
    if ~isempty(vdec)
        if opt_gsa.load_ident_files==0
            SAMorris=NaN(npT,3,size(vdec,2));
            for i=1:size(vdec,2)
                [~, SAMorris(:,:,i)] = gsa.Morris_Measure_Groups(npT, [lpmat0 lpmat], vdec(:,i),nliv);
            end
            SAvdec = squeeze(SAMorris(:,1,:))';
            save([OutputDirectoryName,'/',fname_,'_morris_IDE.mat'],'SAvdec','vdec','ir_vdec','ic_vdec')
        else
            load([OutputDirectoryName,'/',fname_,'_morris_IDE.mat'],'SAvdec')
        end

        hh_fig = dyn_figure(options_.nodisplay,'name','Screening identification: variance decomposition');
        gsa.boxplot(SAvdec,[],'.',[],10)
        set(gca,'xticklabel',' ','fontsize',10,'xtick',1:npT)
        set(gca,'xlim',[0.5 npT+0.5])
        ydum = get(gca,'ylim');
        set(gca,'ylim',[0 ydum(2)])
        set(gca,'position',[0.13 0.2 0.775 0.7])
        for ip=1:npT
            if options_.TeX
                [~, param_name_tex_temp]= get_the_name(ip,options_.TeX,M_,estim_params_,options_.varobs);
                text(ip,-2,param_name_tex_temp,'rotation',90,'HorizontalAlignment','right','interpreter','latex')
            else
                text(ip,-2,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
            end
        end
        xlabel(' ')
        title('Elementary effects variance decomposition')
        dyn_saveas(hh_fig,[OutputDirectoryName,'/',fname_,'_morris_vdec'],options_.nodisplay,options_.graph_format);
        create_TeX_loader(options_,[OutputDirectoryName,'/',fname_,'_morris_vdec'],1,'Screening identification: variance decomposition','morris_vdec',1)
    else
        save([OutputDirectoryName,'/',fname_,'_morris_IDE.mat'],'vdec')
    end

    if opt_gsa.load_ident_files==0
        ccac = [mss cc ac];
        SAMorris=NaN(npT,3,size(ccac,2));
        for i=1:size(ccac,2)
            [~, SAMorris(:,:,i)] = gsa.Morris_Measure_Groups(npT, [lpmat0 lpmat], [ccac(:,i)],nliv);
        end
        SAcc = squeeze(SAMorris(:,1,:))';
        SAcc = SAcc./(max(SAcc,[],2)*ones(1,npT));
        save([OutputDirectoryName,'/',fname_,'_morris_IDE.mat'],'SAcc','cc','ir_cc','ic_cc','-append')
        save([OutputDirectoryName,'/',fname_,'_morris_IDE.mat'],'ac','ir_ac','ic_ac','-append')
    else
        load([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAcc','cc','ir_cc','ic_cc')
        load([OutputDirectoryName,'/',fname_,'_morris_IDE'],'ac','ir_ac','ic_ac')
    end

    hh_fig=dyn_figure(options_.nodisplay,'name','Screening identification: theoretical moments');
    gsa.boxplot(SAcc,[],'.',[],10)
    set(gca,'xticklabel',' ','fontsize',10,'xtick',1:npT)
    set(gca,'xlim',[0.5 npT+0.5])
    set(gca,'ylim',[0 1])
    set(gca,'position',[0.13 0.2 0.775 0.7])
    for ip=1:npT
        if options_.TeX
            [~, param_name_tex_temp]= get_the_name(ip,options_.TeX,M_,estim_params_,options_.varobs);
            text(ip,-0.02,param_name_tex_temp,'rotation',90,'HorizontalAlignment','right','interpreter','latex')
        else
            text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
        end
    end
    xlabel(' ')
    title('Elementary effects in the moments')
    dyn_saveas(hh_fig,[OutputDirectoryName,'/',fname_,'_morris_moments'],options_.nodisplay,options_.graph_format);
    create_TeX_loader(options_,[OutputDirectoryName,'/',fname_,'_morris_moments'],1,'Screening identification: theoretical moments','morris_moments',1)

    if opt_gsa.load_ident_files==0
        SAMorris=NaN(npT,3,j0);
        for j=1:j0
            [~, SAMorris(:,:,j)] = gsa.Morris_Measure_Groups(npT, [lpmat0 lpmat], yt(:,j),nliv);
        end

        SAM = squeeze(SAMorris(1:end,1,:));
        SAnorm=NaN(npT,j0);
        irex=NaN(j0);
        for j=1:j0
            SAnorm(:,j)=SAM(:,j)./max(SAM(:,j));
            irex(j)=length(find(SAnorm(:,j)>0.01));
        end

        SAMmu = squeeze(SAMorris(1:end,2,:));
        SAmunorm=NaN(npT,j0);
        for j=1:j0
            SAmunorm(:,j)=SAMmu(:,j)./max(SAM(:,j));  % normalised w.r.t. mu*
        end
        SAMsig = squeeze(SAMorris(1:end,3,:));
        SAsignorm=NaN(npT,j0);
        for j=1:j0
            SAsignorm(:,j)=SAMsig(:,j)./max(SAMsig(:,j));
        end
        save([OutputDirectoryName,'/',fname_,'_morris_IDE.mat'],'SAnorm','SAmunorm','SAsignorm','-append')
    else
        load([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAnorm')
    end
    hh_fig=dyn_figure(options_.nodisplay,'name','Screening identification: model');
    gsa.boxplot(SAnorm',[],'.',[],10)
    set(gca,'xticklabel',' ','fontsize',10,'xtick',1:npT)
    set(gca,'xlim',[0.5 npT+0.5])
    set(gca,'ylim',[0 1])
    set(gca,'position',[0.13 0.2 0.775 0.7])
    xlabel(' ')
    for ip=1:npT
        if options_.TeX
            [~, param_name_tex_temp]= get_the_name(ip,options_.TeX,M_,estim_params_,options_.varobs);
            text(ip,-0.02,param_name_tex_temp,'rotation',90,'HorizontalAlignment','right','interpreter','latex')
        else
            text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
        end      
    end
    xlabel(' ')
    title('Elementary effects in the model')
    dyn_saveas(hh_fig,[OutputDirectoryName,'/',fname_,'_morris_par'],options_.nodisplay,options_.graph_format);
    create_TeX_loader(options_,[OutputDirectoryName,'/',fname_,'_morris_par'],1,'Screening identification: model','morris_par',1)

elseif opt_gsa.morris==3
    return
elseif opt_gsa.morris==2   % ISKREV stuff
    return
else  % main effects analysis
    if itrans==0
        fsuffix = '';
    elseif itrans==1
        fsuffix = '_log';
    else
        fsuffix = '_rank';
    end

    imap=1:npT;

    if isempty(lpmat0)
        x0=lpmat(istable,:);
    else
        x0=[lpmat0(istable,:), lpmat(istable,:)];
    end
    nrun=length(istable);
    nest=min(250,nrun);

    if opt_gsa.load_ident_files==0
        try
            EET=load([OutputDirectoryName,'/SCREEN/',fname_,'_morris_IDE'],'SAcc','ir_cc','ic_cc');
        catch
            EET=[];
        end
        ccac = gsa.standardize_columns([mss cc ac]);
        [pcc, dd] = eig(cov(ccac(istable,:)));
        [latent, isort] = sort(-diag(dd));
        latent = -latent;
        figure, bar(latent)
        title('Eigenvalues in PCA')
        pcc=pcc(:,isort);
        ccac = ccac*pcc;
        npca = max(find(cumsum(latent)./length(latent)<0.99))+1;
        siPCA = (EET.SAcc'*abs(pcc'))';
        siPCA = siPCA./(max(siPCA,[],2)*ones(1,npT));
        SAcc=zeros(size(ccac,2),npT);
        gsa_=NaN(npca);
        for j=1:npca %size(ccac,2),
            if itrans==0
                y0 = ccac(istable,j);
            elseif itrans==1
                y0 = gsa.log_transform(ccac(istable,j));
            else
                y0 = trank(ccac(istable,j));
            end
            if ~isempty(EET)
                imap=find(siPCA(j,:) >= (0.1.*max(siPCA(j,:))) );
            end
            gsa_(j) = gsa_sdp(y0(1:nest), x0(1:nest,imap), ...
                              2, [],[],[],0,[OutputDirectoryName,'/map_cc',fsuffix,int2str(j)], pnames);
            SAcc(j,imap)=gsa_(j).si;
            imap_cc{j}=imap;
        end
        save([OutputDirectoryName,'/map_cc',fsuffix,'.mat'],'gsa_')
        save([OutputDirectoryName,'/',fname_,'_main_eff.mat'],'imap_cc','SAcc','ccac','-append')
    else
        load([OutputDirectoryName,'/',fname_,'_main_eff'],'SAcc')
    end

    hh_fig=dyn_figure(options_.nodisplay,'Name',['Identifiability indices in the ',fsuffix,' moments.']);
    bar(sum(SAcc))
    set(gca,'xticklabel',' ','fontsize',10,'xtick',1:npT)
    set(gca,'xlim',[0.5 npT+0.5])
    ydum = get(gca,'ylim');
    set(gca,'ylim',[0 ydum(2)])
    set(gca,'position',[0.13 0.2 0.775 0.7])
    for ip=1:npT
        text(ip,-0.02*(ydum(2)),bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
    xlabel(' ')
    title(['Identifiability indices in the ',fsuffix,' moments.'],'interpreter','none')
    dyn_saveas(hh_fig,[OutputDirectoryName,'/',fname_,'_ident_ALL',fsuffix],options_.nodisplay,options_.graph_format);
    create_TeX_loader(options_,[OutputDirectoryName,'/',fname_,'_ident_ALL',fsuffix],1,['Identifiability indices in the ',fsuffix,' moments.'],['ident_ALL',fsuffix]',1)
end

function []=create_TeX_loader(options_,figpath,ifig_number,caption,label_name,scale_factor)
if nargin<6
    scale_factor=1;
end
if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fidTeX = fopen([figpath '.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by map_ident_.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s}\n',scale_factor,strrep(figpath,'\','/'));
    fprintf(fidTeX,'\\caption{%s.}',caption);
    fprintf(fidTeX,'\\label{Fig:%s:%u}\n',label_name,ifig_number);
    fprintf(fidTeX,'\\end{figure}\n\n');
    fprintf(fidTeX,'%% End Of TeX file. \n');
    fclose(fidTeX);
end


function yr = trank(y)
% yr is the rank transformation of y
yr=NaN(size(y));
[nr, nc] = size(y);
for j=1:nc
    [~, is]=sort(y(:,j));
    yr(is,j)=[1:nr]'./nr;
end
