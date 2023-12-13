function myoutput=PosteriorIRF_core1(myinputs,fpar,B,whoiam, ThisMatlab)
%myoutput=PosteriorIRF_core1(myinputs,fpar,B,whoiam, ThisMatlab)
%   Generates and stores Posterior IRFs
%   PARALLEL CONTEXT
%   This function perfoms in parallel execution a portion of the PosteriorIRF.m code.
%   This is a special kind of parallel function. Unlike of other parallel functions,
%   that running in parallel a 'for' cycle, this function run in parallel a
%   'while' loop! The parallelization of 'while' loop (when possible) is a more
%   sophisticated procedure.
%
%   See also the comment in posterior_sampler_core.m funtion.
%
% INPUTS
%   See the comment in posterior_sampler_core.m funtion.
%
% OUTPUTS
% o myoutput  [struc]
%  Contained:
%  OutputFileName_dsge, OutputFileName_param and OutputFileName_bvardsge.
%
% ALGORITHM
%   Portion of PosteriorIRF.m function. Specifically the 'while' cycle.
%
% SPECIAL REQUIREMENTS.
%   None.
%
% Copyright © 2006-2023 Dynare Team
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

if nargin<4
    whoiam=0;
end

% Reshape 'myinputs' for local computation.
% In order to avoid confusion in the name space, the instruction struct2local(myinputs) is replaced by:

options_= myinputs.options_;
estim_params_= myinputs.estim_params_;
M_= myinputs.M_;
oo_= myinputs.oo_;
bayestopt_= myinputs.bayestopt_;
dataset_= myinputs.dataset_;
dataset_info= myinputs.dataset_info;

IRUN = myinputs.IRUN;
irun =myinputs.irun;
irun2=myinputs.irun2;
npar=myinputs.npar;
type=myinputs.type;
if ~strcmpi(type,'prior')
    x=myinputs.x;
end

nvar=myinputs.nvar;
IndxVariables=myinputs.IndxVariables;
MAX_nirfs_dsgevar=myinputs.MAX_nirfs_dsgevar;
MAX_nirfs_dsge=myinputs.MAX_nirfs_dsge;
MAX_nruns=myinputs.MAX_nruns;

NumberOfIRFfiles_dsge=myinputs.NumberOfIRFfiles_dsge;
NumberOfIRFfiles_dsgevar=myinputs.NumberOfIRFfiles_dsgevar;
ifil2=myinputs.ifil2;

if options_.dsge_var
    nvobs=myinputs.nvobs;
    NumberOfParametersPerEquation = myinputs.NumberOfParametersPerEquation;
    NumberOfLagsTimesNvobs = myinputs.NumberOfLagsTimesNvobs;
    Companion_matrix = myinputs.Companion_matrix;
    stock_irf_bvardsge = zeros(options_.irf,nvobs,M_.exo_nbr,MAX_nirfs_dsgevar);
    bounds = prior_bounds(bayestopt_,options_.prior_trunc);
end

if whoiam
    Parallel=myinputs.Parallel;
end

if strcmpi(type,'posterior')
    MhDirectoryName = CheckPath('metropolis',M_.dname);
elseif strcmpi(type,'gsa')
    if options_.opt_gsa.pprior
        MhDirectoryName = CheckPath(['gsa' filesep 'prior'],M_.dname);
    else
        MhDirectoryName = CheckPath(['gsa' filesep 'mc'],M_.dname);
    end
else
    MhDirectoryName = CheckPath('prior',M_.dname);
end

RemoteFlag = 0;

if whoiam
    if Parallel(ThisMatlab).Local==0
        RemoteFlag =1;
    end
    prct0={0,whoiam,Parallel(ThisMatlab)};
else
    prct0=0;
end
if strcmpi(type,'posterior')
    h = dyn_waitbar(prct0,'Bayesian (posterior) IRFs...');
elseif strcmpi(type,'gsa')
    h = dyn_waitbar(prct0,'GSA (prior) IRFs...');
else
    h = dyn_waitbar(prct0,'Bayesian (prior) IRFs...');
end

OutputFileName_bvardsge = {};
OutputFileName_dsge = {};
OutputFileName_param = {};

fpar = fpar-1;
fpar0=fpar;
nosaddle=0;

if whoiam
    ifil2=ifil2(whoiam);
    NumberOfIRFfiles_dsge=NumberOfIRFfiles_dsge(whoiam);
    NumberOfIRFfiles_dsgevar=NumberOfIRFfiles_dsgevar(whoiam);
end

% Parallel 'while' very good!!!
stock_param=zeros(MAX_nruns,npar);
stock_irf_dsge=zeros(options_.irf,nvar,M_.exo_nbr,MAX_nirfs_dsge);

if strcmp(type, 'prior')
    Prior = dprior(bayestopt_, options_.prior_trunc);
end

while fpar<B
    fpar = fpar + 1;
    irun = irun+1;
    irun2 = irun2+1;
    if strcmpi(type,'prior')
        deep = Prior.draw();
    else
        deep = x(fpar,:);
    end
    stock_param(irun2,:) = deep;
    M_ = set_parameters_locally(M_, deep);
    [oo_.dr,info,M_.params] =compute_decision_rules(M_,options_,oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    if info(1)
        nosaddle = nosaddle + 1;
        fpar = fpar - 1;
        irun = irun-1;
        irun2 = irun2-1;
        if info(1) == 1
            errordef = 'Static variables are not uniquely defined';
        elseif info(1) == 2
            errordef = 'Dll problem';
        elseif info(1) == 3
            errordef = 'No stable trajectory';
        elseif info(1) == 4
            errordef = 'Indeterminacy';
        elseif info(1) == 5
            errordef = 'Rank condition  is not satisfied';
        else
            errordef = get_error_message(info, options_);
        end
        if strcmpi(type,'prior')
            disp(['PosteriorIRF :: Dynare is unable to solve the model (' errordef ')'])
            continue
        else
            error(['PosteriorIRF :: Dynare is unable to solve the model (' errordef ') with sample ' type])
        end
    end
    SS = get_lower_cholesky_covariance(M_.Sigma_e,options_.add_tiny_number_to_cholesky);
    irf_shocks_indx = getIrfShocksIndx(M_, options_);
    for i=irf_shocks_indx
        if SS(i,i) > 5e-7
            if options_.order>1 && options_.relative_irf % normalize shock to 0.01 before IRF generation for GIRFs; multiply with 100 later
                y=irf(M_,options_,oo_.dr,SS(:,i)./SS(i,i)/100, options_.irf, options_.drop,options_.replic,options_.order);
            else
                y=irf(M_,options_,oo_.dr,SS(:,i), options_.irf, options_.drop,options_.replic,options_.order);
            end
            if options_.relative_irf && options_.order==1 %multiply with 100 for backward compatibility
                y = 100*y/SS(i,i);
            end
            for j = 1:nvar
                if max(y(IndxVariables(j),:)) - min(y(IndxVariables(j),:)) > 1e-12
                    stock_irf_dsge(:,j,i,irun) = transpose(y(IndxVariables(j),:));
                end
            end
        end
    end
    if MAX_nirfs_dsgevar
        IRUN = IRUN+1;
        [~,~,~,~,~,~,~,PHI,SIGMAu,iXX] =  dsge_var_likelihood(deep',dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
        dsge_prior_weight = M_.params(strmatch('dsge_prior_weight', M_.param_names));
        DSGE_PRIOR_WEIGHT = floor(dataset_.nobs*(1+dsge_prior_weight));
        SIGMA_inv_upper_chol = chol(inv(SIGMAu*dataset_.nobs*(dsge_prior_weight+1)));
        explosive_var  = 1;
        while explosive_var
            % draw from the marginal posterior of SIGMA
            SIGMAu_draw = rand_inverse_wishart(dataset_.vobs, DSGE_PRIOR_WEIGHT-NumberOfParametersPerEquation, ...
                                               SIGMA_inv_upper_chol);
            % draw from the conditional posterior of PHI
            PHI_draw = rand_matrix_normal(NumberOfParametersPerEquation,dataset_.vobs, PHI, ...
                                          chol(SIGMAu_draw)', chol(iXX)');
            Companion_matrix(1:dataset_.vobs,:) = transpose(PHI_draw(1:NumberOfLagsTimesNvobs,:));
            % Check for stationarity
            explosive_var = any(abs(eig(Companion_matrix))>1.000000001);
        end
        % Get the mean
        mu = zeros(1,dataset_.vobs);
        % Get rotation
        if dsge_prior_weight > 0
            Atheta(oo_.dr.order_var,:) = oo_.dr.ghu*sqrt(M_.Sigma_e);
            A0 = Atheta(bayestopt_.mfys,:);
            OMEGAstar = qr2(A0');
        end
        SIGMAu_chol = chol(SIGMAu_draw)';
        SIGMAtrOMEGA = SIGMAu_chol*OMEGAstar';
        PHIpower = eye(NumberOfLagsTimesNvobs);
        irfs = zeros (options_.irf,dataset_.vobs*M_.exo_nbr);
        tmp3 = PHIpower(1:dataset_.vobs,1:dataset_.vobs)*SIGMAtrOMEGA;
        irfs(1,:) = tmp3(:)';
        for t = 2:options_.irf
            PHIpower = Companion_matrix*PHIpower;
            tmp3 = PHIpower(1:dataset_.vobs,1:dataset_.vobs)*SIGMAtrOMEGA;
            irfs(t,:)  = tmp3(:)'+kron(ones(1,M_.exo_nbr),mu);
        end
        tmp_dsgevar = kron(ones(options_.irf,1),mu);
        for j = 1:(dataset_.vobs*M_.exo_nbr)
            if max(irfs(:,j)) - min(irfs(:,j)) > 1e-10
                tmp_dsgevar(:,j) = (irfs(:,j));
            end
        end
        if IRUN < MAX_nirfs_dsgevar
            stock_irf_bvardsge(:,:,:,IRUN) = reshape(tmp_dsgevar,options_.irf,dataset_.vobs,M_.exo_nbr);
        else
            stock_irf_bvardsge(:,:,:,IRUN) = reshape(tmp_dsgevar,options_.irf,dataset_.vobs,M_.exo_nbr);
            save([MhDirectoryName '/' M_.fname '_irf_bvardsge' int2str(NumberOfIRFfiles_dsgevar) '.mat'], 'stock_irf_bvardsge');
            if RemoteFlag==1
                OutputFileName_bvardsge = [OutputFileName_bvardsge; {[MhDirectoryName filesep], [M_.fname '_irf_bvardsge' int2str(NumberOfIRFfiles_dsgevar) '.mat']}];
            end
            NumberOfIRFfiles_dsgevar = NumberOfIRFfiles_dsgevar+1;
            IRUN =0;
        end
    end
    if irun == MAX_nirfs_dsge || irun == B || fpar == B
        if fpar == B
            stock_irf_dsge = stock_irf_dsge(:,:,:,1:irun);
            if MAX_nirfs_dsgevar && (fpar == B || IRUN == B)
                stock_irf_bvardsge = stock_irf_bvardsge(:,:,:,1:IRUN);
                save([MhDirectoryName '/' M_.fname '_irf_bvardsge' int2str(NumberOfIRFfiles_dsgevar) '.mat'], 'stock_irf_bvardsge');
                NumberOfIRFfiles_dsgevar = NumberOfIRFfiles_dsgevar+1;
                if RemoteFlag==1
                    OutputFileName_bvardsge = [OutputFileName_bvardsge; {[MhDirectoryName filesep], [M_.fname '_irf_bvardsge' int2str(NumberOfIRFfiles_dsgevar) '.mat']}];
                end
                irun = 0;
            end
        end
        save([MhDirectoryName '/' M_.fname '_irf_dsge' int2str(NumberOfIRFfiles_dsge) '.mat'],'stock_irf_dsge');
        if RemoteFlag==1
            OutputFileName_dsge = [OutputFileName_dsge; {[MhDirectoryName filesep], [M_.fname '_irf_dsge' int2str(NumberOfIRFfiles_dsge) '.mat']}];
        end
        NumberOfIRFfiles_dsge = NumberOfIRFfiles_dsge+1;
        irun = 0;
    end
    if irun2 == MAX_nruns || fpar == B
        if fpar == B
            stock_param = stock_param(1:irun2,:);
        end
        stock = stock_param;
        save([MhDirectoryName '/' M_.fname '_param_irf' int2str(ifil2) '.mat'],'stock');
        if RemoteFlag==1
            OutputFileName_param = [OutputFileName_param; {[MhDirectoryName filesep], [M_.fname '_param_irf' int2str(ifil2) '.mat']}];
        end
        ifil2 = ifil2 + 1;
        irun2 = 0;
    end
    dyn_waitbar((fpar-fpar0)/(B-fpar0),h);
end

dyn_waitbar_close(h);

if whoiam==0
    if nosaddle
        disp(['PosteriorIRF :: Percentage of discarded posterior draws = ' num2str(nosaddle/(B+nosaddle))])
    end
end

% Copy the rusults of computation on the call machine (specifically in the
% directory on call machine that contain the model).

myoutput.OutputFileName = [OutputFileName_dsge;
                           OutputFileName_param;
                           OutputFileName_bvardsge];

myoutput.nosaddle = nosaddle;
