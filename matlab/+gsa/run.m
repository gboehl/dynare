function x0=run(M_,oo_,options_,bayestopt_,estim_params_,options_gsa)
% x0=run(M_,oo_,options_,bayestopt_,estim_params_,options_gsa)
% Frontend to the Sensitivity Analysis Toolbox for DYNARE
% Inputs:
%  - M_                     [structure]     Matlab's structure describing the model
%  - oo_                    [structure]     Matlab's structure describing the results
%  - options_               [structure]     Matlab's structure describing the current options
%  - bayestopt_             [structure]     describing the priors
%  - estim_params_          [structure]     characterizing parameters to be estimated
%  - options_gsa            [structure]     Matlab's structure describing the GSA options
%
% Reference:
% M. Ratto (2008), Analysing DSGE Models with Global Sensitivity Analysis, 
% Computational Economics (2008), 31, pp. 115–139

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

if options_.dsge_var
    error('Identification does not support DSGE-VARs at the current stage')
end

fname_ = M_.fname;
x0=[];

% check user defined options
if isfield(options_gsa,'neighborhood_width') && options_gsa.neighborhood_width
    if isfield(options_gsa,'pprior') && options_gsa.pprior
        error('sensitivity:: neighborhood_width is incompatible with prior sampling')
    end
    if isfield(options_gsa,'ppost') && options_gsa.ppost
        error('sensitivity:: neighborhood_width is incompatible with posterior sampling')
    end
end

if isfield(options_gsa,'morris') && options_gsa.morris==1
    if isfield(options_gsa,'identification') && options_gsa.identification==0
    end
    if isfield(options_gsa,'ppost') && options_gsa.ppost
        error('sensitivity:: Morris is incompatible with posterior sampling')
    elseif isfield(options_gsa,'pprior') && options_gsa.pprior==0
        if ~(isfield(options_gsa,'neighborhood_width') && options_gsa.neighborhood_width)
            error('sensitivity:: Morris is incompatible with MC sampling with correlation matrix')
        end
    end
    if isfield(options_gsa,'rmse') && options_gsa.rmse
        error('sensitivity:: Morris is incompatible with rmse analysis')
    end
    if (isfield(options_gsa,'alpha2_stab') && options_gsa.alpha2_stab<1) || ...
            (isfield(options_gsa,'pvalue_ks') && options_gsa.pvalue_ks) || ...
            (isfield(options_gsa,'pvalue_corr') && options_gsa.pvalue_corr)
        error('sensitivity:: Morris is incompatible with Monte Carlo filtering')
    end
end

% end check user defined options
options_gsa = set_default_option(options_gsa,'datafile',[]);
options_gsa = set_default_option(options_gsa,'rmse',0);
options_gsa = set_default_option(options_gsa,'useautocorr',0);

options_gsa = set_default_option(options_gsa,'moment_calibration',options_.endogenous_prior_restrictions.moment);
options_gsa = set_default_option(options_gsa,'irf_calibration',options_.endogenous_prior_restrictions.irf);
if isfield(options_gsa,'nograph')
    options_.nograph=options_gsa.nograph;
end
if isfield(options_gsa,'nodisplay')
    options_.nodisplay=options_gsa.nodisplay;
end
if isfield(options_gsa,'graph_format')
    options_.graph_format=options_gsa.graph_format;
end
if isfield(options_gsa,'mode_file')
    options_.mode_file=options_gsa.mode_file;
elseif isfield(options_gsa,'neighborhood_width') && options_gsa.neighborhood_width>0
    options_.mode_file='';
end

if options_.order~=1
    warning('dynare_sensitivity: dynare_sensitivity does only support order=1, resetting to order=1.')
    options_.order = 1;
end

if ~isempty(options_gsa.datafile) || isempty(bayestopt_) || options_gsa.rmse
    if isempty(options_gsa.datafile) && options_gsa.rmse
        disp('The data file and all relevant estimation options ')
        disp('[first_obs, nobs, presample, prefilter, loglinear, lik_init, kalman_algo, ...]')
        disp('must be specified for RMSE analysis!');
        error('Sensitivity anaysis error!')
    end
    if isfield(options_gsa,'nobs')
        options_.nobs=options_gsa.nobs;
    end
    if ~isempty(options_.nobs) && length(options_.nobs)~=1
        error('dynare_sensitivity does not support recursive estimation. Please specify nobs as a scalar, not a vector.')
    end
    options_.datafile = options_gsa.datafile;
    if isfield(options_gsa,'first_obs')
        options_.first_obs=options_gsa.first_obs;
    end
    if isfield(options_gsa,'presample')
        options_.presample=options_gsa.presample;
    end
    if isfield(options_gsa,'prefilter')
        options_.prefilter=options_gsa.prefilter;
    end
    if isfield(options_gsa,'loglinear')
        options_.loglinear=options_gsa.loglinear;
    end
    if isfield(options_gsa,'lik_init')
        options_.lik_init=options_gsa.lik_init;
    end
    if isfield(options_gsa,'diffuse_filter')
        options_.diffuse_filter=options_gsa.diffuse_filter;
    end
    if isfield(options_gsa,'kalman_algo')
        options_.kalman_algo=options_gsa.kalman_algo;
    end
    options_.mode_compute = 0;
    options_.filtered_vars = 1;
    options_.plot_priors = 0;
    [dataset_,dataset_info,~,~, M_, options_, oo_, estim_params_, bayestopt_] = ...
        dynare_estimation_init(M_.endo_names, fname_, 1, M_, options_, oo_, estim_params_, bayestopt_);
    % computes a first linear solution to set up various variables
else
    if isempty(options_.qz_criterium)
        options_.qz_criterium = 1+1e-6;
    end
    if options_.prior_trunc==0
        options_.prior_trunc=1.e-10;
    end
end

if M_.exo_nbr==0
    error('dynare_sensitivity does not support having no varexo in the model. As a workaround you could define a dummy exogenous variable.')
end

[~,~,~,~,oo_.dr,M_.params] = dynare_resolve(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);

options_gsa = set_default_option(options_gsa,'identification',0);
if options_gsa.identification
    options_gsa.redform=0;
    options_gsa = set_default_option(options_gsa,'morris',1);
    options_gsa = set_default_option(options_gsa,'trans_ident',0);
    options_gsa = set_default_option(options_gsa,'load_ident_files',0);
    options_gsa = set_default_option(options_gsa,'ar',1);
    options_.ar = options_gsa.ar;
    if options_gsa.morris==0
        disp('The option morris = 0 is no longer supported (Type I errors)')
        disp('This option is reset at morris = 2 (local identification analysis).')
        options_gsa.morris=2;
    end
    if options_gsa.morris==2
        if isfield(options_,'options_ident')
            options_.options_ident.load_ident_files = options_gsa.load_ident_files;
            options_.options_ident.useautocorr = options_gsa.useautocorr;
            options_.options_ident.ar = options_gsa.ar;
        else
            options_ident=[];
            options_ident = set_default_option(options_ident,'load_ident_files',options_gsa.load_ident_files);
            options_ident = set_default_option(options_ident,'useautocorr',options_gsa.useautocorr);
            options_ident = set_default_option(options_ident,'ar',options_gsa.ar);
            options_.options_ident = options_ident;
        end
    end
end

% map stability
options_gsa = set_default_option(options_gsa,'stab',1);
options_gsa = set_default_option(options_gsa,'pprior',1);
options_gsa = set_default_option(options_gsa,'prior_range',1);
options_gsa = set_default_option(options_gsa,'ppost',0);
options_gsa = set_default_option(options_gsa,'neighborhood_width',0);
options_gsa = set_default_option(options_gsa,'ilptau',1);
options_gsa = set_default_option(options_gsa,'morris',0);
options_gsa = set_default_option(options_gsa,'glue',0);
options_gsa = set_default_option(options_gsa,'morris_nliv',6);
options_gsa = set_default_option(options_gsa,'morris_ntra',20);
options_gsa = set_default_option(options_gsa,'Nsam',2048);
options_gsa = set_default_option(options_gsa,'load_stab',0);
options_gsa = set_default_option(options_gsa,'alpha2_stab',0);
options_gsa = set_default_option(options_gsa,'pvalue_ks',0.001);
options_gsa = set_default_option(options_gsa,'pvalue_corr',1.e-5);
% REDFORM mapping
options_gsa = set_default_option(options_gsa,'redform',0);
options_gsa = set_default_option(options_gsa,'load_redform',0);
options_gsa = set_default_option(options_gsa,'logtrans_redform',0);
options_gsa = set_default_option(options_gsa,'threshold_redform',[]);
options_gsa = set_default_option(options_gsa,'ksstat_redform',0.001);
options_gsa = set_default_option(options_gsa,'alpha2_redform',1.e-5);
options_gsa = set_default_option(options_gsa,'namendo',{});
options_gsa = set_default_option(options_gsa,'namlagendo',{});
options_gsa = set_default_option(options_gsa,'namexo',{});
options_gsa = set_default_option(options_gsa,'namendo_tex',{});
options_gsa = set_default_option(options_gsa,'namlagendo_tex',{});
options_gsa = set_default_option(options_gsa,'namexo_tex',{});
if strmatch(':',options_gsa.namendo,'exact')
    options_gsa.namendo = M_.endo_names(1:M_.orig_endo_nbr);
end
if strmatch(':',options_gsa.namexo,'exact')
    options_gsa.namexo = M_.exo_names;
end
if strmatch(':',options_gsa.namlagendo,'exact')
    options_gsa.namlagendo = M_.endo_names(1:M_.orig_endo_nbr);
end

if options_.TeX
    [~,Locb]=ismember(options_gsa.namendo,M_.endo_names);
    options_gsa.namendo_tex=cellfun(@(x) horzcat('$', x, '$'), M_.endo_names_tex(Locb), 'UniformOutput', false);
    [~,Locb]=ismember(options_gsa.namlagendo,M_.endo_names);
    options_gsa.namlagendo_tex=cellfun(@(x) horzcat('$', x, '$'), M_.endo_names_tex(Locb), 'UniformOutput', false);
    [~,Locb]=ismember(options_gsa.namexo,M_.exo_names);
    options_gsa.namexo_tex=cellfun(@(x) horzcat('$', x, '$'), M_.exo_names_tex(Locb), 'UniformOutput', false);
end
% RMSE mapping
options_gsa = set_default_option(options_gsa,'load_rmse',0);
options_gsa = set_default_option(options_gsa,'lik_only',0);
options_gsa = set_default_option(options_gsa,'var_rmse', options_.varobs);
%get corresponding TeX-names;
options_gsa.var_rmse_tex={};
for ii=1:length(options_gsa.var_rmse)
    temp_name = M_.endo_names_tex{strmatch(options_gsa.var_rmse{ii}, M_.endo_names, 'exact')};
    options_gsa.var_rmse_tex = vertcat(options_gsa.var_rmse_tex, ['$' temp_name '$']);
end
options_gsa.varobs_tex = cellfun(@(x) horzcat('$', x, '$'), M_.endo_names_tex(options_.varobs_id), 'UniformOutput', false);
options_gsa = set_default_option(options_gsa,'pfilt_rmse', 0.1);
options_gsa = set_default_option(options_gsa,'istart_rmse', options_.presample+1);
options_gsa = set_default_option(options_gsa,'alpha_rmse', 0.001);
options_gsa = set_default_option(options_gsa,'alpha2_rmse', 1.e-5);

if options_gsa.neighborhood_width
    options_gsa.pprior=0;
    options_gsa.ppost=0;
end

if options_gsa.redform && options_gsa.neighborhood_width==0 && isempty(options_gsa.threshold_redform)
    options_gsa.pprior=1;
    options_gsa.ppost=0;
end

if options_gsa.morris>2
    disp('The option morris = 3 is no longer supported')
    disp('the option is reset at morris = 1 .')
    options_gsa.morris=1;
end

if options_gsa.morris==1
    if ~options_gsa.identification
        options_gsa.redform=1;
    end
    if options_gsa.neighborhood_width
        options_gsa.pprior=0;
    else
        options_gsa.pprior=1;
    end
    options_gsa.ppost=0;
    options_gsa.glue=0;
    options_gsa.rmse=0;
    options_gsa.load_rmse=0;
    options_gsa.alpha2_stab=1;
    options_gsa.pvalue_ks=0;
    options_gsa.pvalue_corr=0;
    OutputDirectoryName = CheckPath('gsa/screen',M_.dname);
else
    OutputDirectoryName = CheckPath('gsa',M_.dname);
end

if (options_gsa.load_stab || options_gsa.load_rmse || options_gsa.load_redform) && options_gsa.pprior
    filetoload=[OutputDirectoryName '/' fname_ '_prior.mat'];
    if ~exist(filetoload,'file')
        disp([filetoload,' not found!'])
        disp('You asked to load a non existent analysis')
        return
    else
        if isempty(strmatch('bkpprior',who('-file', filetoload),'exact'))
            disp('Warning! Missing prior info for saved sample') % trap for files previous
            disp('The saved files are generated with previous version of GSA package') % trap for files previous
        else
            load(filetoload,'bkpprior')
            if any(bayestopt_.pshape~=bkpprior.pshape) || ...
                    any(bayestopt_.p1~=bkpprior.p1) || ...
                    any(bayestopt_.p2~=bkpprior.p2) || ...
                    any(bayestopt_.p3(~isnan(bayestopt_.p3))~=bkpprior.p3(~isnan(bkpprior.p3))) || ...
                    any(bayestopt_.p4(~isnan(bayestopt_.p4))~=bkpprior.p4(~isnan(bkpprior.p4)))
                disp('WARNING!')
                disp('The saved sample has different priors w.r.t. to current ones!!')
                skipline()
                disp('Press ENTER to continue')
                pause
            end
        end
    end
end

if options_gsa.stab && ~options_gsa.ppost
    x0 = gsa.stability_mapping(OutputDirectoryName,options_gsa,M_,oo_,options_,bayestopt_,estim_params_);
    if isempty(x0)
        skipline()
        disp('Sensitivity computations stopped: no parameter set provided a unique solution')
        return
    end
end

options_.opt_gsa = options_gsa;
if ~isempty(options_gsa.moment_calibration) || ~isempty(options_gsa.irf_calibration)
    gsa.map_calibration(OutputDirectoryName, M_, options_, oo_, estim_params_,bayestopt_);
end

if options_gsa.identification
    gsa.map_identification(OutputDirectoryName,options_gsa,M_,oo_,options_,estim_params_,bayestopt_);
end

if options_gsa.redform && ~isempty(options_gsa.namendo)
    if options_gsa.ppost
        filnam = dir([M_.dname filesep 'metropolis' filesep '*param_irf*.mat']);
        lpmat=[];
        for j=1:length(filnam)
            load ([M_.dname filesep 'metropolis' filesep M_.fname '_param_irf' int2str(j) '.mat'],'stock')
            lpmat=[lpmat; stock];
        end
        clear stock
        nshock = estim_params_.nvx;
        nshock = nshock + estim_params_.nvn;
        nshock = nshock + estim_params_.ncx;
        nshock = nshock + estim_params_.ncn;

        lpmat0=lpmat(:,1:nshock);
        lpmat=lpmat(:,nshock+1:end);
        istable=(1:size(lpmat,1));
        iunstable=[];
        iwrong=[];
        iindeterm=[];
        save([OutputDirectoryName filesep M_.fname '_mc.mat'],'lpmat','lpmat0','istable','iunstable','iwrong','iindeterm')
        options_gsa.load_stab=1;

        x0 = gsa.stability_mapping(OutputDirectoryName,options_gsa,M_,oo_,options_,bayestopt_,estim_params_);
    end
    if options_gsa.morris==1
        gsa.reduced_form_screening(OutputDirectoryName,options_gsa, estim_params_, M_, oo_.dr, options_, bayestopt_);
    else
        % check existence of the SS_ANOVA toolbox
        if isempty(options_gsa.threshold_redform) && ~(exist('gsa_sdp','file')==6 || exist('gsa_sdp','file')==2)
            fprintf('\nThe "SS-ANOVA-R: MATLAB Toolbox for the estimation of Smoothing Spline ANOVA models with Recursive algorithms" is missing.\n')
            fprintf('To obtain it, go to:\n\n')
            fprintf('https://ec.europa.eu/jrc/en/macro-econometric-statistical-software/ss-anova-r-downloads \n\n')
            fprintf('and follow the instructions there.\n')
            fprintf('After obtaining the files, you need to unpack them and set a Matlab Path to those files.\n')
            error('SS-ANOVA-R Toolbox missing!')
        end
        gsa.reduced_form_mapping(OutputDirectoryName,options_gsa,M_,estim_params_,options_,bayestopt_,oo_);
    end
end
% RMSE mapping
options_.opt_gsa = options_gsa;
if options_gsa.rmse
    if ~options_gsa.ppost
        if options_gsa.pprior
            a=whos('-file',[OutputDirectoryName,'/',fname_,'_prior'],'logpo2');
        else
            a=whos('-file',[OutputDirectoryName,'/',fname_,'_mc'],'logpo2');
        end
        if isoctave()
            aflag=0;
            for ja=1:length(a)
                aflag=aflag+strcmp('logpo2',a(ja).name);
            end
            if aflag==0
                a=[];
            else
                a=1;
            end
        end
        if isempty(a)
            if options_gsa.lik_only
                options_.smoother=0;
                options_.filter_step_ahead=[];
                options_.forecast=0;
                options_.filtered_vars=0;
            end
            if options_gsa.pprior
                TmpDirectoryName = ([M_.dname filesep 'gsa' filesep 'prior']);
            else
                TmpDirectoryName = ([M_.dname filesep 'gsa' filesep 'mc']);
            end
            if exist(TmpDirectoryName,'dir')
                mydelete([M_.fname '_filter_step_ahead*.mat'],[TmpDirectoryName filesep]);
                mydelete([M_.fname '_inno*.mat'],[TmpDirectoryName filesep]);
                mydelete([M_.fname '_smooth*.mat'],[TmpDirectoryName filesep]);
                mydelete([M_.fname '_update*.mat'],[TmpDirectoryName filesep]);
                filparam = dir([TmpDirectoryName filesep M_.fname '_param*.mat']);
                for j=1:length(filparam)
                    if isempty(strmatch([M_.fname '_param_irf'],filparam(j).name))
                        delete([TmpDirectoryName filesep filparam(j).name]);
                    end
                end
            end
            oo_=prior_posterior_statistics('gsa',dataset_, dataset_info,M_,oo_,options_,estim_params_,bayestopt_,'gsa::mcmc');
            if options_.bayesian_irf
                oo_=PosteriorIRF('gsa',options_,estim_params_,oo_,M_,bayestopt_,dataset_,dataset_info,'gsa::mcmc');
            end
            options_gsa.load_rmse=0;
        end
    end
    clear a;
    gsa.monte_carlo_filtering(OutputDirectoryName,options_gsa,dataset_,dataset_info,M_,oo_,options_,bayestopt_,estim_params_);
end
options_.opt_gsa = options_gsa;

if options_gsa.glue
    dr_ = oo_.dr;
    if options_gsa.ppost
        load([OutputDirectoryName,'/',fname_,'_post']);
        DirectoryName = CheckPath('metropolis',M_.dname);
    else
        if options_gsa.pprior
            load([OutputDirectoryName,'/',fname_,'_prior']);
        else
            load([OutputDirectoryName,'/',fname_,'_mc']);
        end
    end
    if ~exist('x','var')
        disp('No RMSE analysis is available for current options')
        disp('No GLUE file prepared')
        return,
    end
    gend = options_.nobs;
    rawdata = read_variables(options_.datafile,options_.varobs,[],options_.xls_sheet,options_.xls_range);
    rawdata = rawdata(options_.first_obs:options_.first_obs+gend-1,:);
    if options_.loglinear
        rawdata = log(rawdata);
    end
    if options_.prefilter == 1
        data = transpose(rawdata-ones(gend,1)*mean(rawdata,1));
    else
        data = transpose(rawdata);
    end

    Obs.data = data;
    Obs.time = 1:gend;
    Obs.num  = gend;
    for j=1:length(options_.varobs)
        Obs.name{j} = options_.varobs{j};
        vj = options_.varobs{j};

        jxj = strmatch(vj,M_.endo_names(dr_.order_var),'exact');
        if ~options_gsa.ppost
            Out(j).data = jxj;
            Out(j).time = [pwd,'/',OutputDirectoryName];
        else
            Out(j).data = jxj;
            Out(j).time = [pwd,'/',DirectoryName];
        end
        Out(j).name = vj;
        Out(j).ini  = 'yes';
        Lik(j).name = ['rmse_',vj];
        Lik(j).ini  = 'yes';
        Lik(j).isam = 1;
        Lik(j).data = rmse_MC(:,j)';

        Out1=Out;
        ismoo(j)=jxj;

    end
    jsmoo = length(options_.varobs);
    for j=1:M_.endo_nbr
        if ~ismember(j,ismoo)
            jsmoo=jsmoo+1;
            vj = M_.endo_names{dr_.order_var(j)};
            if ~options_gsa.ppost
                Out1(jsmoo).data = j;
                Out1(jsmoo).time = [pwd,'/',OutputDirectoryName];
            else
                Out1(jsmoo).data = j;
                Out1(jsmoo).time = [pwd,'/',DirectoryName];
            end
            Out1(jsmoo).name = vj;
            Out1(jsmoo).ini  = 'yes';
        end
    end
    tit = M_.exo_names;
    for j=1:M_.exo_nbr
        Exo(j).name = tit{j};
    end
    if ~options_gsa.ppost
        Lik(length(options_.varobs)+1).name = 'logpo';
        Lik(length(options_.varobs)+1).ini  = 'yes';
        Lik(length(options_.varobs)+1).isam = 1;
        Lik(length(options_.varobs)+1).data = -logpo2;
    end
    Sam.name = bayestopt_.name;
    Sam.dim  = [size(x) 0];
    Sam.data = x;

    Rem.id = 'Original';
    Rem.ind= 1:size(x,1);

    Info.dynare=M_.fname;
    Info.order_var=dr_.order_var;
    Out=Out1;
    if options_gsa.ppost
        Info.TypeofSample='post';
        save([OutputDirectoryName,'/',fname_,'_glue_post.mat'], 'Out', 'Sam', 'Lik', 'Obs', 'Rem','Info', 'Exo')
    else
        if options_gsa.pprior
            Info.TypeofSample='prior';
            save([OutputDirectoryName,'/',fname_,'_glue_prior.mat'], 'Out', 'Sam', 'Lik', 'Obs', 'Rem','Info', 'Exo')
        else
            Info.TypeofSample='mc';
            save([OutputDirectoryName,'/',fname_,'_glue_mc.mat'], 'Out', 'Sam', 'Lik', 'Obs', 'Rem','Info', 'Exo')
        end
    end
end