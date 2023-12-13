function oo_=prior_posterior_statistics(type,dataset_,dataset_info,M_,oo_,options_,estim_params_,bayestopt_,dispString)
% oo_=prior_posterior_statistics(type,dataset_,dataset_info,M_,oo_,options_,estim_params_,bayestopt_,dispString))
% Computes Monte Carlo filter smoother and forecasts
%
% INPUTS
%    type           [string]        posterior, prior, or gsa
%   o dataset_      [structure]     storing the dataset
%   o dataset_info  [structure]     Various information about the dataset
%   o M_            [structure]     storing the model information
%   o oo_           [structure]     storing the results
%   o options_      [structure]     storing the options
%   o estim_params_ [structure]     storing information about estimated parameters
%   o bayestopt_    [structure]     storing information about priors
%    dispString:   	[string] 		display info in the command window
% OUTPUTS
%    oo_:           [structure]     storing the results
%
% SPECIAL REQUIREMENTS
%    none
%
% PARALLEL CONTEXT
% See the comments in the posterior_sampler.m funtion.


% Copyright © 2005-2023 Dynare Team
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

if nargin < 9
    dispString = 'prior_posterior_statistics';
end
localVars=[];

Y = transpose(dataset_.data);
gend = dataset_.nobs;

nvn  = estim_params_.nvn;
npar = estim_params_.nvx+nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.np;
naK = length(options_.filter_step_ahead);

MaxNumberOfBytes=options_.MaxNumberOfBytes;
endo_nbr=M_.endo_nbr;
exo_nbr=M_.exo_nbr;
meas_err_nbr=length(M_.Correlation_matrix_ME);
horizon = options_.forecast;
if horizon
    i_last_obs = gend+(1-M_.maximum_endo_lag:0);
end
maxlag = M_.maximum_endo_lag;

if strcmpi(type,'posterior')
    DirectoryName = CheckPath('metropolis',M_.dname);
    B = options_.sub_draws;
elseif strcmpi(type,'gsa')
    RootDirectoryName = CheckPath('gsa',M_.dname);
    if options_.opt_gsa.pprior
        DirectoryName = CheckPath(['gsa',filesep,'prior'],M_.dname);
        load([ RootDirectoryName filesep  M_.fname '_prior.mat'],'lpmat0','lpmat','istable')
    else
        DirectoryName = CheckPath(['gsa',filesep,'mc'],M_.dname);
        load([ RootDirectoryName filesep  M_.fname '_mc.mat'],'lpmat0','lpmat','istable')
    end
    if ~isempty(lpmat0)
        x=[lpmat0(istable,:) lpmat(istable,:)];
    else
        x=lpmat(istable,:);
    end
    clear lpmat lpmat0 istable
    B = size(x,1);
elseif strcmpi(type,'prior')
    DirectoryName = CheckPath('prior',M_.dname);
    B = options_.prior_draws;
end

MAX_nruns = min(B,ceil(MaxNumberOfBytes/(npar+2)/8));
MAX_nsmoo = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*gend)/8));
MAX_n_smoothed_constant = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*gend)/8));
MAX_n_smoothed_trend = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*gend)/8));
MAX_n_trend_coeff = min(B,ceil(MaxNumberOfBytes/endo_nbr/8));
MAX_ninno = min(B,ceil(MaxNumberOfBytes/(exo_nbr*gend)/8));
MAX_nerro = min(B,ceil(MaxNumberOfBytes/(size(options_.varobs,1)*gend)/8));

if naK
    MAX_naK   = min(B,ceil(MaxNumberOfBytes/(endo_nbr* ...
                                             length(options_.filter_step_ahead)*(gend+max(options_.filter_step_ahead)))/8));
end

if horizon
    MAX_nforc1 = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*(horizon+maxlag))/8));
    MAX_nforc2 = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*(horizon+maxlag))/8));
    if ~isequal(M_.H,0)
        MAX_nforc_ME = min(B,ceil(MaxNumberOfBytes/((size(options_.varobs,1))*(horizon+maxlag))/8));
    end
end
MAX_momentsno = min(B,ceil(MaxNumberOfBytes/(get_moments_size(options_)*8)));

if options_.filter_covariance
    MAX_filter_covariance = min(B,ceil(MaxNumberOfBytes/(endo_nbr^2*(gend+1))/8));
end

if options_.smoothed_state_uncertainty
    MAX_n_smoothed_state_uncertainty = min(B,ceil(MaxNumberOfBytes/((endo_nbr*endo_nbr)*gend)/8));
end

varlist = options_.varlist;
if isempty(varlist)
    varlist = sort(M_.endo_names(1:M_.orig_endo_nbr));
end

n_variables_to_fill=13;

irun = ones(n_variables_to_fill,1);
ifil = zeros(n_variables_to_fill,1);

run_smoother = 0;
if options_.smoother || options_.forecast || ~isempty(options_.filter_step_ahead) || options_.smoothed_state_uncertainty
    run_smoother = 1;
    if options_.loglinear
        oo_.Smoother.loglinear = true;
    else
        oo_.Smoother.loglinear = false;
    end
end

filter_covariance=0;
if options_.filter_covariance
    filter_covariance=1;
end

smoothed_state_uncertainty=0;
if options_.smoothed_state_uncertainty
    smoothed_state_uncertainty=1;
end

% Store the variable mandatory for local/remote parallel computing.

localVars.type=type;
localVars.run_smoother=run_smoother;
localVars.filter_covariance=filter_covariance;
localVars.smoothed_state_uncertainty=smoothed_state_uncertainty;
localVars.gend=gend;
localVars.Y=Y;
localVars.data_index=dataset_info.missing.aindex;
localVars.missing_value=dataset_info.missing.state;
localVars.varobs=options_.varobs;
localVars.mean_varobs=dataset_info.descriptive.mean;
localVars.irun=irun;
localVars.endo_nbr=endo_nbr;
localVars.nvn=nvn;
localVars.naK=naK;
localVars.horizon=horizon;
localVars.iendo=1:endo_nbr;
localVars.IdObs=bayestopt_.mfys;
if horizon
    localVars.i_last_obs=i_last_obs;
    localVars.MAX_nforc1=MAX_nforc1;
    localVars.MAX_nforc2=MAX_nforc2;
    if ~isequal(M_.H,0)
        localVars.MAX_nforc_ME = MAX_nforc_ME;
    end
end
localVars.exo_nbr=exo_nbr;
localVars.maxlag=maxlag;
localVars.MAX_nsmoo=MAX_nsmoo;
localVars.MAX_ninno=MAX_ninno;
localVars.MAX_nerro = MAX_nerro;
if naK
    localVars.MAX_naK=MAX_naK;
end
if options_.filter_covariance
    localVars.MAX_filter_covariance = MAX_filter_covariance;
end
if options_.smoothed_state_uncertainty
    localVars.MAX_n_smoothed_state_uncertainty = MAX_n_smoothed_state_uncertainty ;
end
localVars.MAX_n_smoothed_constant=MAX_n_smoothed_constant;
localVars.MAX_n_smoothed_trend=MAX_n_smoothed_trend;
localVars.MAX_n_trend_coeff=MAX_n_trend_coeff;
localVars.MAX_nruns=MAX_nruns;
localVars.MAX_momentsno = MAX_momentsno;
localVars.ifil=ifil;
localVars.DirectoryName = DirectoryName;

localVars.M_=M_;
localVars.oo_=oo_;
localVars.options_=options_;
localVars.estim_params_=estim_params_;
localVars.bayestopt_=bayestopt_;


if strcmpi(type,'posterior')
    record=load_last_mh_history_file(DirectoryName, M_.fname);
    [nblck, npar] = size(record.LastParameters);
    FirstMhFile = record.KeepedDraws.FirstMhFile;
    FirstLine = record.KeepedDraws.FirstLine;
    TotalNumberOfMhFiles = sum(record.MhDraws(:,2));
    TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
    NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
    mh_nblck = options_.mh_nblck;
    if B==NumberOfDraws*mh_nblck
        % we load all retained MH runs !
        logpost=GetAllPosteriorDraws(options_, M_.dname, M_.fname, 0, FirstMhFile, FirstLine, TotalNumberOfMhFiles, NumberOfDraws, nblck);
        for column=1:npar
            x(:,column) = GetAllPosteriorDraws(options_, M_.dname, M_.fname, column, FirstMhFile, FirstLine, TotalNumberOfMhFiles, NumberOfDraws, nblck);
        end
    else
        logpost=NaN(B,1);
        for b=1:B
            [x(b,:), logpost(b)] = GetOneDraw(type,M_,estim_params_,oo_,options_,bayestopt_);
        end
    end
    localVars.logpost=logpost;
end

if ~strcmpi(type,'prior')
    localVars.x=x;
end

% Like sequential execution!
if isnumeric(options_.parallel)
    [fout] = prior_posterior_statistics_core(localVars,1,B,0);
    % Parallel execution!
else
    [~, totCPU, nBlockPerCPU] = distributeJobs(options_.parallel, 1, B);
    ifil=zeros(n_variables_to_fill,totCPU);
    for j=1:totCPU-1
        if run_smoother
            nfiles = ceil(nBlockPerCPU(j)/MAX_nsmoo);
            ifil(1,j+1) =ifil(1,j)+nfiles;
            nfiles = ceil(nBlockPerCPU(j)/MAX_ninno);
            ifil(2,j+1) =ifil(2,j)+nfiles;
            nfiles = ceil(nBlockPerCPU(j)/MAX_nerro);
            ifil(3,j+1) =ifil(3,j)+nfiles;
        end
        if naK
            nfiles = ceil(nBlockPerCPU(j)/MAX_naK);
            ifil(4,j+1) =ifil(4,j)+nfiles;
        end
        nfiles = ceil(nBlockPerCPU(j)/MAX_nruns);
        ifil(5,j+1) =ifil(5,j)+nfiles;
        if horizon
            nfiles = ceil(nBlockPerCPU(j)/MAX_nforc1);
            ifil(6,j+1) =ifil(6,j)+nfiles;
            nfiles = ceil(nBlockPerCPU(j)/MAX_nforc2);
            ifil(7,j+1) =ifil(7,j)+nfiles;
            if ~isequal(M_.H,0)
                nfiles = ceil(nBlockPerCPU(j)/MAX_nforc_ME);
                ifil(12,j+1) =ifil(12,j)+nfiles;
            end
        end
        if options_.filter_covariance
            nfiles = ceil(nBlockPerCPU(j)/MAX_filter_covariance);
            ifil(8,j+1) =ifil(8,j)+nfiles;
        end
        if run_smoother
            nfiles = ceil(nBlockPerCPU(j)/MAX_n_trend_coeff);
            ifil(9,j+1) =ifil(9,j)+nfiles;
            nfiles = ceil(nBlockPerCPU(j)/MAX_n_smoothed_constant);
            ifil(10,j+1) =ifil(10,j)+nfiles;
            nfiles = ceil(nBlockPerCPU(j)/MAX_n_smoothed_trend);
            ifil(11,j+1) =ifil(11,j)+nfiles;
            if smoothed_state_uncertainty
                nfiles = ceil(nBlockPerCPU(j)/MAX_n_smoothed_state_uncertainty);
                ifil(13,j+1) =ifil(13,j)+nfiles;
            end
        end
    end
    localVars.ifil = ifil;
    globalVars = [];
    % which files have to be copied to run remotely
    NamFileInput(1,:) = {'',[M_.fname '.static.m']};
    NamFileInput(2,:) = {'',[M_.fname '.dynamic.m']};
    if M_.set_auxiliary_variables
        NamFileInput(3,:) = {'',[M_.fname '.set_auxiliary_variables.m']};
    end
    if options_.steadystate_flag
        if options_.steadystate_flag == 1
            NamFileInput(length(NamFileInput)+1,:)={'',[M_.fname '_steadystate.m']};
        else
            NamFileInput(length(NamFileInput)+1,:)={'',[M_.fname '.steadystate.m']};
        end
    end
    [fout] = masterParallel(options_.parallel, 1, B,NamFileInput,'prior_posterior_statistics_core', localVars,globalVars, options_.parallel_info);

end
ifil = fout(end).ifil;

stock_gend=gend;
stock_data=Y;
save([DirectoryName '/' M_.fname '_data.mat'],'stock_gend','stock_data');

if strcmpi(type,'gsa')
    return
end

if ~isnumeric(options_.parallel)
    leaveSlaveOpen = options_.parallel_info.leaveSlaveOpen;
    if options_.parallel_info.leaveSlaveOpen == 0
        % Commenting for testing!!!
        options_.parallel_info.leaveSlaveOpen = 1; % Force locally to leave open remote matlab sessions (repeated pm3 calls)
    end
end

if options_.smoother
    oo_=pm3(M_,options_,oo_,endo_nbr,gend,ifil(1),B,'Smoothed variables',...
        varlist, M_.endo_names_tex,M_.endo_names,...
        varlist,'SmoothedVariables',DirectoryName,'_smooth',dispString);
    oo_=pm3(M_,options_,oo_,exo_nbr,gend,ifil(2),B,'Smoothed shocks',...
        M_.exo_names,M_.exo_names_tex,M_.exo_names,...
        M_.exo_names,'SmoothedShocks',DirectoryName,'_inno',dispString);
    oo_=pm3(M_,options_,oo_,endo_nbr,1,ifil(9),B,'Trend_coefficients',...
        varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'TrendCoeff',DirectoryName,'_trend_coeff',dispString);
    oo_=pm3(M_,options_,oo_,endo_nbr,gend,ifil(10),B,'Smoothed constant',...
        varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'Constant',DirectoryName,'_smoothed_constant',dispString);
    oo_=pm3(M_,options_,oo_,endo_nbr,gend,ifil(11),B,'Smoothed trend',...
        varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'Trend',DirectoryName,'_smoothed_trend',dispString);
    oo_=pm3(M_,options_,oo_,endo_nbr,gend,ifil(1),B,'Updated Variables',...
        varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'UpdatedVariables',DirectoryName, ...
        '_update',dispString);
    if smoothed_state_uncertainty
        oo_=pm3(M_,options_,oo_,endo_nbr,endo_nbr,ifil(13),B,'State Uncertainty',...
            varlist,M_.endo_names_tex,M_.endo_names,...
            varlist,'StateUncertainty',DirectoryName,'_state_uncert',dispString);
    end

    if nvn
        for obs_iter=1:length(options_.varobs)
            meas_error_names{obs_iter,1}=['SE_EOBS_' M_.endo_names{strmatch(options_.varobs{obs_iter},M_.endo_names,'exact')}];
            texnames{obs_iter,1}=['\sigma^{ME}_' M_.endo_names_tex{strmatch(options_.varobs{obs_iter},M_.endo_names,'exact')}];
        end
        oo_=pm3(M_,options_,oo_,meas_err_nbr,gend,ifil(3),B,'Smoothed measurement errors',...
            meas_error_names,texnames,meas_error_names,...
            meas_error_names,'SmoothedMeasurementErrors',DirectoryName,'_error',dispString);
    end
end

if options_.filtered_vars
    oo_=pm3(M_,options_,oo_,endo_nbr,gend,ifil(4),B,'One step ahead forecast (filtered variables)',...
        varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'FilteredVariables',DirectoryName,'_filter_step_ahead',dispString);
end

if options_.forecast
    oo_=pm3(M_,options_,oo_,endo_nbr,horizon,ifil(6),B,'Forecasted variables (mean)',...
        varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'MeanForecast',DirectoryName,'_forc_mean',dispString);
    oo_=pm3(M_,options_,oo_,endo_nbr,horizon,ifil(7),B,'Forecasted variables (point)',...
        varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'PointForecast',DirectoryName,'_forc_point',dispString);
    if ~isequal(M_.H,0) && ~isempty(intersect(options_.varobs,varlist))
        texnames = cell(length(options_.varobs), 1);
        obs_names = cell(length(options_.varobs), 1);
        for obs_iter=1:length(options_.varobs)
            obs_names{obs_iter}=M_.endo_names{strmatch(options_.varobs{obs_iter},M_.endo_names,'exact')};
            texnames{obs_iter}=M_.endo_names_tex{strmatch(options_.varobs{obs_iter},M_.endo_names,'exact')};
        end
        varlist_forecast_ME=intersect(options_.varobs,varlist);
        oo_=pm3(M_,options_,oo_,meas_err_nbr,horizon,ifil(12),B,'Forecasted variables (point) with ME',...
            varlist_forecast_ME,texnames,obs_names,...
            varlist_forecast_ME,'PointForecastME',DirectoryName,'_forc_point_ME',dispString);
    end
end

if options_.filter_covariance
    oo_=pm3(M_,options_,oo_,endo_nbr,endo_nbr,ifil(8),B,'Filtered covariances',...
        varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'FilterCovariance',DirectoryName,'_filter_covar',dispString);
end


if ~isnumeric(options_.parallel)
    options_.parallel_info.leaveSlaveOpen = leaveSlaveOpen;
    if leaveSlaveOpen == 0
        closeSlave(options_.parallel,options_.parallel_info.RemoteTmpFolder),
    end
end
