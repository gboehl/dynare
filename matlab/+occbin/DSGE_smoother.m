function [alphahat,etahat,epsilonhat,ahat0,SteadyState,trend_coeff,aKK,T0,R0,P,PKK,decomp,Trend,state_uncertainty,oo_,bayestopt_,alphahat0,state_uncertainty0] = DSGE_smoother(xparam1,gend,Y,data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_,dataset_, dataset_info)
% [alphahat,etahat,epsilonhat,ahat0,SteadyState,trend_coeff,aKK,T0,R0,P,PKK,decomp,Trend,state_uncertainty,oo_,bayestopt_,alphahat0,state_uncertainty0] = DSGE_smoother(xparam1,gend,Y,data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_,dataset_, dataset_info)
% Runs a DSGE smoother with occasionally binding constraints
%
% INPUTS
% - xparam1       [double]        (p*1) vector of (estimated) parameters.
% - gend          [integer]       scalar specifying the number of observations
% - Y             [double]        (n*T) matrix of data.
% - data_index    [cell]          1*smpl cell of column vectors of indices.
% - missing_value [boolean]       1 if missing values, 0 otherwise
% - M_            [structure]     Matlab's structure describing the model (M_).
% - oo_           [structure]     Matlab's structure containing the results (oo_).
% - options_      [structure]     Matlab's structure describing the current options (options_).
% - bayestopt_    [structure]     describing the priors
% - estim_params_ [structure]     characterizing parameters to be estimated
% - dataset_      [structure]     the dataset after required transformation
% - dataset_info  [structure]     Various informations about the dataset (descriptive statistics and missing observations)
%
% OUTPUTS
% - alphahat      [double]  (m*T) matrix, smoothed endogenous variables (a_{t|T})  (decision-rule order)
% - etahat        [double]  (r*T) matrix, smoothed structural shocks (r>=n is the number of shocks).
% - epsilonhat    [double]  (n*T) matrix, smoothed measurement errors.
% - ahat0         [double]  (m*T) matrix, updated (endogenous) variables (a_{t|t}) (decision-rule order)
% - SteadyState   [double]  (m*1) vector specifying the steady state level of each endogenous variable (declaration order)
% - trend_coeff   [double]  (n*1) vector, parameters specifying the slope of the trend associated to each observed variable.
% - aKK           [double]  (K,n,T+K) array, k (k=1,...,K) steps ahead
%                                   filtered (endogenous) variables  (decision-rule order)
% - T0 and R0     [double]  Matrices defining the state equation (T is the (m*m) transition matrix).
% - P:            (m*m*(T+1)) 3D array of one-step ahead forecast error variance
%                       matrices (decision-rule order)
% - PKK           (K*m*m*(T+K)) 4D array of k-step ahead forecast error variance
%                       matrices (meaningless for periods 1:d) (decision-rule order)
% - decomp        (K*m*r*(T+K)) 4D array of shock decomposition of k-step ahead
%                       filtered variables (decision-rule order)
% - Trend         [double] (n*T) pure trend component; stored in options_.varobs order
% - state_uncertainty [double] (K,K,T) array, storing the uncertainty
%                                   about the smoothed state (decision-rule order)
% - M_            [structure] decribing the model
% - oo_           [structure] storing the results
% - options_      [structure] describing the options
% - bayestopt_    [structure] describing the priors
% - alphahat0     [double]  (m*1) array, smoothed endogenous variables in period 0 (a_{0|T})  (decision-rule order)
% - state_uncertainty0 [double] (K,K,1) array, storing the uncertainty in period 0

% Copyright © 2021-2023 Dynare Team
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

smoother_field_list = {'SmoothedVariables', 'UpdatedVariables', 'SmoothedShocks'};

regime_history=[];
if  options_.occbin.smoother.linear_smoother && nargin==12
    %% linear smoother
    options_.occbin.smoother.status=false;
    [alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK,T0,R0,P,PK,decomp,Trend,state_uncertainty,oo_,bayestopt_] = DsgeSmoother(xparam1,gend,Y,data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_);
    tmp_smoother=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK,P,PK,decomp,Trend,state_uncertainty);
    for jf=1:length(smoother_field_list)
        oo_.occbin.linear_smoother.(smoother_field_list{jf}) = tmp_smoother.(smoother_field_list{jf});
    end
    oo_.occbin.linear_smoother.alphahat=alphahat;
    oo_.occbin.linear_smoother.etahat=etahat;
    oo_.occbin.linear_smoother.epsilonhat=epsilonhat;
    oo_.occbin.linear_smoother.ahat=ahat;
    oo_.occbin.linear_smoother.SteadyState=SteadyState;
    oo_.occbin.linear_smoother.trend_coeff=trend_coeff;
    oo_.occbin.linear_smoother.aK=aK;
    oo_.occbin.linear_smoother.T0=T0;
    oo_.occbin.linear_smoother.R0=R0;
    oo_.occbin.linear_smoother.decomp=decomp;
    
    fprintf('\nOccbin: linear smoother done.\n')
    options_.occbin.smoother.status=true;
end
% if init_mode
%% keep these commented lines for the moment, should some nasty problem in the future call back for this option
%     shocks1     = etahat(:,2:end)';
%     clear regime_history ;
%         
%     for ii = 1:gend-1
% 
%         fprintf('Period number %d out of %d \n',ii,gend-1)
%         opts_simul.SHOCKS = shocks1(ii,:);
%         opts_simul.endo_init = alphahat(oo_.dr.inv_order_var,ii);
%         out =  occbin.runsim_fn(opts_simul, M_, oo_, options_); 
%         regime_history(ii) = out.regime_history;
%     end    
%  
%     [TT,RR] = dynare_resolve(M_,options_,oo_);
%     CC = zeros([size(RR,1),gend ]);
% else

opts_simul = options_.occbin.simul;
opts_simul.curb_retrench = options_.occbin.smoother.curb_retrench;
opts_simul.maxit = options_.occbin.smoother.maxit;
opts_simul.waitbar = options_.occbin.smoother.waitbar;
opts_simul.periods = options_.occbin.smoother.periods;
opts_simul.check_ahead_periods = options_.occbin.smoother.check_ahead_periods;
opts_simul.max_check_ahead_periods = options_.occbin.smoother.max_check_ahead_periods;
opts_simul.periodic_solution = options_.occbin.smoother.periodic_solution;
opts_simul.full_output = options_.occbin.smoother.full_output;
opts_simul.piecewise_only = options_.occbin.smoother.piecewise_only;
occbin_options = struct();

occbin_options.first_period_occbin_update = options_.occbin.smoother.first_period_occbin_update;
occbin_options.opts_regime = opts_simul; % this builds the opts_simul options field needed by occbin.solver
occbin_options.opts_regime.binding_indicator = options_.occbin.likelihood.init_binding_indicator;
occbin_options.opts_regime.regime_history=options_.occbin.likelihood.init_regime_history;
[alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK,T0,R0,P,PK,decomp,Trend,state_uncertainty,oo_,bayestopt_,alphahat0,state_uncertainty0] = DsgeSmoother(xparam1,gend,Y,data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_,occbin_options);%     T1=TT;

oo_.occbin.smoother.realtime_regime_history = oo_.occbin.smoother.regime_history;
regime_history = oo_.occbin.smoother.regime_history;
opts_regime.regime_history = oo_.occbin.smoother.regime_history;

ahat0 = ahat;
aKK=aK;
PKK=PK;
clear aK PK;

occbin_options.first_period_occbin_update = inf;

opts_regime.binding_indicator=[];
regime_history0 = regime_history;

fprintf('Occbin smoother iteration 1.\n')
opts_simul.SHOCKS = [etahat(:,1:end)'; zeros(1,M_.exo_nbr)];
opts_simul.exo_pos = 1:M_.exo_nbr;
opts_simul.endo_init = alphahat0(oo_.dr.inv_order_var,1);
opts_simul.init_regime=regime_history; % use realtime regime for guess, to avoid multiple solution issues!
opts_simul.periods = size(opts_simul.SHOCKS,1);
options_.occbin.simul=opts_simul;
options_.noprint = true;
[~, out, ss] = occbin.solver(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
if out.error_flag
    fprintf('Occbin smoother:: simulation within smoother did not converge.\n')
    print_info(out.error_flag, options_.noprint, options_)
    oo_.occbin.smoother.error_flag=321;
    return;
end
regime_history = out.regime_history;
if options_.smoother_redux
    occbin_options.opts_simul.restrict_state_space =1;  
    [T0,R0] = dynare_resolve(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
    oo_.occbin.linear_smoother.T0=T0;
    oo_.occbin.linear_smoother.R0=R0;
end
TT = ss.T(oo_.dr.order_var,oo_.dr.order_var,:);
RR = ss.R(oo_.dr.order_var,:,:);
CC = ss.C(oo_.dr.order_var,:);

opts_regime.regime_history = regime_history;
opts_regime.binding_indicator = [];
[TT, RR, CC, regime_history] = occbin.check_regimes(TT, RR, CC, opts_regime, M_, options_ , oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
is_changed = ~isequal(regime_history0,regime_history);
if isempty(regime_history0)
    regime_history0 = regime_history;
end    
iter=1;

is_periodic = 0;
is_changed_start = 0;
maxiter = options_.occbin.smoother.max_number_of_iterations;
occbin_smoother_fast = options_.occbin.smoother.fast;
occbin_smoother_debug=options_.occbin.smoother.debug;

sto_alphahat=alphahat;
sto_etahat={etahat};
sto_CC = CC;
sto_RR = RR;
sto_TT = TT;
sto_eee=NaN(size(TT,1),size(TT,3));
for k=1:size(TT,3)
    sto_eee(:,k) = eig(TT(:,:,k));
end


while is_changed && maxiter>iter && ~is_periodic
    iter=iter+1;
    fprintf('Occbin smoother iteration %u.\n', iter)
    occbin_options.opts_regime.regime_history=regime_history;
    [alphahat,etahat,epsilonhat,~,SteadyState,trend_coeff,~,T0,R0,P,~,decomp,Trend,state_uncertainty,oo_,bayestopt_,alphahat0,state_uncertainty0]...
        = DsgeSmoother(xparam1,gend,Y,data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_,occbin_options,TT,RR,CC);
    sto_etahat(iter)={etahat};
    regime_history0(iter,:) = regime_history;
    if occbin_smoother_debug
        save('Occbin_smoother_debug_regime_history','regime_history0');
    end
    
    sto_CC = CC;
    sto_RR = RR;
    sto_TT = TT;
    
    opts_simul.SHOCKS = [etahat(:,1:end)'; zeros(1,M_.exo_nbr)];
    opts_simul.endo_init = alphahat0(oo_.dr.inv_order_var,1);
    options_.occbin.simul=opts_simul;
    [~, out, ss] = occbin.solver(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
    if out.error_flag
        fprintf('Occbin smoother:: simulation within smoother did not converge.\n')
        print_info(out.error_flag, false, options_)
        oo_.occbin.smoother.error_flag=321;
        return;
    end
    regime_history = out.regime_history;
    TT = ss.T(oo_.dr.order_var,oo_.dr.order_var,:);
    RR = ss.R(oo_.dr.order_var,:,:);
    CC = ss.C(oo_.dr.order_var,:);

    opts_regime.regime_history = regime_history;
    [TT, RR, CC, regime_history] = occbin.check_regimes(TT, RR, CC, opts_regime, M_, options_ , oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    is_changed = ~isequal(regime_history0(iter,:),regime_history);
    isdiff_regime=NaN(size(regime_history0,2),M_.occbin.constraint_nbr);
    isdiff_start=NaN(size(isdiff_regime));
    isdiff_=NaN(size(isdiff_regime));
    if M_.occbin.constraint_nbr==2
        for k=1:size(regime_history0,2)
            isdiff_regime(k,1) = ~isequal(regime_history0(end,k).regime1,regime_history(k).regime1);
            isdiff_start(k,1) = ~isequal(regime_history0(end,k).regimestart1,regime_history(k).regimestart1);
            isdiff_(k,1) =  isdiff_regime(k,1) ||  isdiff_start(k,1);
            isdiff_regime(k,2) = ~isequal(regime_history0(end,k).regime2,regime_history(k).regime2);
            isdiff_start(k,2) = ~isequal(regime_history0(end,k).regimestart2,regime_history(k).regimestart2);
            isdiff_(k,2) =  isdiff_regime(k,2) ||  isdiff_start(k,2);
        end
        is_changed_regime = any(isdiff_regime(:,1)) || any(isdiff_regime(:,2));
        is_changed_start = any(isdiff_start(:,1)) || any(isdiff_start(:,2));
    else
        for k=1:size(regime_history0,2)
            isdiff_regime(k,1) = ~isequal(regime_history0(end,k).regime,regime_history(k).regime);
            isdiff_start(k,1) = ~isequal(regime_history0(end,k).regimestart,regime_history(k).regimestart);
            isdiff_(k,1) =  isdiff_regime(k,1) ||  isdiff_start(k,1);
        end
        is_changed_regime = any(isdiff_regime(:,1));
        is_changed_start = any(isdiff_start(:,1));
    end
    if occbin_smoother_fast
        is_changed = is_changed_regime;
    end
    if iter>1
        for kiter=1:iter-1
            is_tmp = isequal(regime_history0(kiter,:),regime_history);
            if is_tmp
                break
            end
        end
        is_periodic = (is_changed && is_tmp);
    end
    
    if is_changed
        eee=NaN(size(TT,1),size(TT,3));
        for k=1:size(TT,3)
            eee(:,k) = eig(TT(:,:,k));
        end
        if options_.debug
            err_eig(iter-1) = max(max(abs(sort(eee)-sort(sto_eee))));
            err_alphahat(iter-1) = max(max(max(abs(alphahat-sto_alphahat))));
            err_etahat(iter-1) = max(max(max(abs(etahat-sto_etahat{iter-1}))));
            err_CC(iter-1) = max(max(max(abs(CC-sto_CC))));
            err_RR(iter-1) = max(max(max(abs(RR-sto_RR))));
            err_TT(iter-1) = max(max(max(abs(TT-sto_TT))));
        end
    end

    if occbin_smoother_debug || is_periodic
        regime_ = cell(0);
        regime_new = regime_;
        start_ = regime_;
        start_new = regime_;
        if M_.occbin.constraint_nbr==2
            indx_init_1 = find(isdiff_(:,1));
            if ~isempty(indx_init_1)
                qq={regime_history0(end,indx_init_1).regime1}';
                for j=1:length(qq), regime_(j,1) = {int2str(qq{j})}; end
                qq={regime_history0(end,indx_init_1).regimestart1}';
                for j=1:length(qq), start_(j,1) = {int2str(qq{j})}; end
                qq={regime_history(indx_init_1).regime1}';
                for j=1:length(qq), regime_new(j,1) = {int2str(qq{j})}; end
                qq={regime_history(indx_init_1).regimestart1}';
                for j=1:length(qq), start_new(j,1) = {int2str(qq{j})}; end
                disp('Time points where regime 1 differs')
                if ~isoctave
                    disp(table(indx_init_1, regime_, start_, regime_new, start_new))
                else % The table() function is not implemented in Octave, print something more or less equivalent (though much less readable)
                    disp(vertcat({'indx_init_1', 'regime_', 'start_', 'regime_new', 'start_new'}, ...
                                 horzcat(num2cell(indx_init_1), regime_, start_, regime_new, start_new)))
                end
            end
            
            indx_init_2 = find(isdiff_(:,2));
            if ~isempty(indx_init_2)
                regime_ = cell(0);
                regime_new = regime_;
                start_ = regime_;
                start_new = regime_;
                qq={regime_history0(end,indx_init_2).regime2}';
                for j=1:length(qq), regime_(j,1) = {int2str(qq{j})}; end
                qq={regime_history0(end,indx_init_2).regimestart2}';
                for j=1:length(qq), start_(j,1) = {int2str(qq{j})}; end
                qq={regime_history(indx_init_2).regime2}';
                for j=1:length(qq), regime_new(j,1) = {int2str(qq{j})}; end
                qq={regime_history(indx_init_2).regimestart2}';
                for j=1:length(qq), start_new(j,1) = {int2str(qq{j})}; end
                disp('Time points where regime 2 differs ')
                if ~isoctave
                    disp(table(indx_init_2, regime_, start_, regime_new, start_new))
                else % The table() function is not implemented in Octave, print something more or less equivalent (though much less readable)
                    disp(vertcat({'indx_init_2', 'regime_', 'start_', 'regime_new', 'start_new'}, ...
                                 horzcat(num2cell(indx_init_2), regime_, start_, regime_new, start_new)))
                end
            end
        else
            indx_init_1 = find(isdiff_(:,1));
            if ~isempty(indx_init_1)
                qq={regime_history0(end,indx_init_1).regime}';
                for j=1:length(qq), regime_(j,1) = {int2str(qq{j})}; end
                qq={regime_history0(end,indx_init_1).regimestart}';
                for j=1:length(qq), start_(j,1) = {int2str(qq{j})}; end
                qq={regime_history(indx_init_1).regime}';
                for j=1:length(qq), regime_new(j,1) = {int2str(qq{j})}; end
                qq={regime_history(indx_init_1).regimestart}';
                for j=1:length(qq), start_new(j,1) = {int2str(qq{j})}; end
                disp('Time points where regime differs')
                if ~isoctave
                    disp(table(indx_init_1, regime_, start_, regime_new, start_new))
                else % The table() function is not implemented in Octave, print something more or less equivalent (though much less readable)
                    disp(vertcat({'indx_init_1', 'regime_', 'start_', 'regime_new', 'start_new'}, ...
                                 horzcat(num2cell(indx_init_1), regime_, start_, regime_new, start_new)))
                end
            end

        end
    end
    if is_changed
        sto_alphahat=alphahat;
        sto_eee = eee;
    end
end
regime_history0(max(iter+1,1),:) = regime_history;
oo_.occbin.smoother.regime_history=regime_history0(end,:);
oo_.occbin.smoother.regime_history_iter=regime_history0;
if occbin_smoother_debug
    save('Occbin_smoother_debug_regime_history','regime_history0')
end

if (maxiter==iter && is_changed) || is_periodic
    disp('occbin.DSGE_smoother: smoother did not converge.')
    fprintf('occbin.DSGE_smoother: The algorithm did not reach a fixed point for the smoothed regimes.\n')
    if is_periodic
        oo_.occbin.smoother.error_flag=0;
        fprintf('occbin.DSGE_smoother: For the periods indicated above, regimes loops between the "regime_" and the "regime_new_" pattern displayed above.\n')
        fprintf('occbin.DSGE_smoother: We provide smoothed shocks consistent with "regime_" in oo_.\n')
    else
        fprintf('occbin.DSGE_smoother: The respective fields in oo_ will be left empty.\n')
        oo_.occbin.smoother=[];
        oo_.occbin.smoother.error_flag=322;
    end
else
    disp('occbin.DSGE_smoother: smoother converged.')
    oo_.occbin.smoother.error_flag=0;
    if occbin_smoother_fast && is_changed_start
        disp('occbin.DSGE_smoother: WARNING: fast algo is used, regime duration was not forced to converge')
    end
end
if (~is_changed || occbin_smoother_debug) && nargin==12
    if is_changed
        CC = sto_CC;
        RR = sto_RR;
        TT = sto_TT;
        oo_.occbin.smoother.regime_history=regime_history0(end-1,:);
    end
    tmp_smoother=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,alphahat,etahat,epsilonhat,ahat0,SteadyState,trend_coeff,aKK,P,PKK,decomp,Trend,state_uncertainty);
    for jf=1:length(smoother_field_list)
        oo_.occbin.smoother.(smoother_field_list{jf}) = tmp_smoother.(smoother_field_list{jf});
    end
    oo_.occbin.smoother.alphahat=alphahat;
    oo_.occbin.smoother.etahat=etahat;
    oo_.occbin.smoother.epsilonhat=epsilonhat;
    oo_.occbin.smoother.ahat=ahat0;
    oo_.occbin.smoother.SteadyState=SteadyState;
    oo_.occbin.smoother.trend_coeff=trend_coeff;
    oo_.occbin.smoother.aK=aKK;
    oo_.occbin.smoother.T0=TT;
    oo_.occbin.smoother.R0=RR;
    oo_.occbin.smoother.C0=CC;
    oo_.occbin.smoother.simul.piecewise = out.piecewise(1:end-1,:);
    if ~options_.occbin.simul.piecewise_only
        oo_.occbin.smoother.simul.linear = out.linear(1:end-1,:);
    end        
    if options_.occbin.smoother.plot
        GraphDirectoryName = CheckPath('graphs',M_.fname);
        latexFolder = CheckPath('latex',M_.dname);

        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fidTeX = fopen([latexFolder filesep M_.fname '_OccBin_smoother_plots.tex'],'w');
            fprintf(fidTeX,'%% TeX eps-loader file generated by occbin.DSGE_smoother.m (Dynare).\n');
            fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
            fprintf(fidTeX,' \n');
        end
        j1=0;
        ifig=0;
        for j=1:M_.exo_nbr
            if max(abs(oo_.occbin.smoother.etahat(j,:)))>1.e-8
                j1=j1+1;
                if mod(j1,9)==1
                    hh_fig = dyn_figure(options_.nodisplay,'name','Occbin smoothed shocks');
                    ifig=ifig+1;
                    isub=0;
                end
                isub=isub+1;
                subplot(3,3,isub)
                if  options_.occbin.smoother.linear_smoother
                    plot(oo_.occbin.linear_smoother.etahat(j,:)','linewidth',2)
                    hold on,
                end
                plot(oo_.occbin.smoother.etahat(j,:)','r--','linewidth',2)
                hold on, plot([0 options_.nobs],[0 0],'k--')
                set(gca,'xlim',[0 options_.nobs])
                if options_.TeX
                    title(['$' M_.exo_names_tex{j,:} '$'],'interpreter','latex')
                else
                    title(M_.exo_names{j,:},'interpreter','none')
                end

                if mod(j1,9)==0
                    if  options_.occbin.smoother.linear_smoother
                        annotation('textbox', [0.1,0,0.35,0.05],'String', 'Linear','Color','Blue','horizontalalignment','center','interpreter','none');
                    end
                    annotation('textbox', [0.55,0,0.35,0.05],'String', 'Piecewise','Color','Red','horizontalalignment','center','interpreter','none');
                    dyn_saveas(gcf,[GraphDirectoryName filesep M_.fname,'_smoothedshocks_occbin',int2str(ifig)],options_.nodisplay,options_.graph_format);
                    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                        % TeX eps loader file
                        fprintf(fidTeX,'\\begin{figure}[H]\n');
                        fprintf(fidTeX,'\\centering \n');
                        fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_smoothedshocks_occbin%s}\n',options_.figures.textwidth*min(j1/3,1),[GraphDirectoryName '/' M_.fname],int2str(ifig)); % don't use filesep as it will create issues with LaTeX on Windows
                        fprintf(fidTeX,'\\caption{Check plots.}');
                        fprintf(fidTeX,'\\label{Fig:smoothedshocks_occbin:%s}\n',int2str(ifig));
                        fprintf(fidTeX,'\\end{figure}\n');
                        fprintf(fidTeX,' \n');
                    end
                end
            end
        end

        if mod(j1,9)~=0 && j==M_.exo_nbr
            annotation('textbox', [0.1,0,0.35,0.05],'String', 'Linear','Color','Blue','horizontalalignment','center','interpreter','none');
            annotation('textbox', [0.55,0,0.35,0.05],'String', 'Piecewise','Color','Red','horizontalalignment','center','interpreter','none');
            dyn_saveas(hh_fig,[GraphDirectoryName filesep M_.fname,'_smoothedshocks_occbin',int2str(ifig)],options_.nodisplay,options_.graph_format);
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                % TeX eps loader file
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_smoothedshocks_occbin%s}\n',options_.figures.textwidth*min(j1/3,1),[GraphDirectoryName '/' M_.fname],int2str(ifig)); % don't use filesep as it will create issues with LaTeX on Windows
                fprintf(fidTeX,'\\caption{Check plots.}');
                fprintf(fidTeX,'\\label{Fig:smoothedshocks_occbin:%s}\n',int2str(ifig));
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,' \n');
            end
        end
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fclose(fidTeX);
        end
    end
end
