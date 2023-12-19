function [oo_, error_flag] = forecast(options_,M_,oo_,forecast) %,hist_period)
%function oo_ = forecast(options_,M_,oo_,forecast)
% forecast

opts = options_.occbin.forecast;

options_.occbin.simul.maxit = opts.maxit;
options_.occbin.simul.check_ahead_periods = opts.check_ahead_periods;
options_.occbin.simul.periods = forecast;
SHOCKS0 = opts.SHOCKS0;
if ~isempty(SHOCKS0)
    if iscell(SHOCKS0)
        for j=1:length(SHOCKS0)
            sname = SHOCKS0{j}{1};
            inds(j)=strmatch(sname,M_.exo_names);
            SHOCKS1(j,:)=SHOCKS0{j}{2};
        end
    elseif isreal(SHOCKS0)
        SHOCKS1=SHOCKS0;
        inds = 1:M_.exo_nbr;
    end
end


if opts.replic
    h = dyn_waitbar(0,'Please wait occbin forecast replic ...');
    ishock = find(sqrt(diag((M_.Sigma_e))));
    effective_exo_nbr=  length(ishock);
    effective_exo_names = M_.exo_names(ishock);
    effective_Sigma_e = M_.Sigma_e(ishock,ishock);
    [U,S] = svd(effective_Sigma_e);
    if opts.qmc
        opts.replic =2^(round(log2(opts.replic+1)))-1;
        SHOCKS_ant =   qmc_sequence(forecast*effective_exo_nbr, int64(1), 1, opts.replic)';
    else
        SHOCKS_ant = randn(forecast*effective_exo_nbr,opts.replic)';
    end
    zlin0=zeros(forecast,M_.endo_nbr,opts.replic);
    zpiece0=zeros(forecast,M_.endo_nbr,opts.replic);
    for iter=1:opts.replic

        if ~isempty(SHOCKS0)
            for j=1:length(SHOCKS0)
                SHOCKS(:,inds(j))=SHOCKS1(j,:);
            end
        end

        error_flagx=1;
        % while error_flagx,
        % SHOCKS=transpose(sqrt(diag(diag(effective_Sigma_e)))*(reshape(SHOCKS_ant(iter,:),forecast,effective_exo_nbr))');
        SHOCKS=transpose(U*sqrt(S)*(reshape(SHOCKS_ant(iter,:),forecast,effective_exo_nbr))');
        %             SHOCKS=transpose(U*sqrt(S)*randn(forecast,M_.exo_nbr)');   %realized shocks
        options_.occbin.simul.endo_init = M_.endo_histval(:,1)-oo_.steady_state;
        options_.occbin.simul.init_regime = opts.frcst_regimes;
        options_.occbin.simul.init_binding_indicator = [];
        options_.occbin.simul.exo_pos=ishock;
        options_.occbin.simul.SHOCKS = SHOCKS;
        options_.occbin.simul.waitbar=0;
        [~, out] = occbin.solver(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
        zlin0(:,:,iter)=out.linear;
        zpiece0(:,:,iter)=out.piecewise;
        ys=out.ys;
        frcst_regime_history(iter,:)=out.regime_history;
        error_flag(iter)=out.error_flag;
        error_flagx = error_flag(iter);
        % end
        simul_SHOCKS(:,:,iter) = SHOCKS;

        if error_flag(iter) && debug_flag
            %             display('no solution')

%             keyboard;
            save no_solution SHOCKS zlin0 zpiece0 iter frcst_regime_history
        end
        dyn_waitbar(iter/opts.replic,h,['OccBin MC forecast replic ',int2str(iter),'/',int2str(opts.replic)])
    end
    dyn_waitbar_close(h);
    save temp zlin0 zpiece0 simul_SHOCKS error_flag
    inx=find(error_flag==0);
    zlin0=zlin0(:,:,inx);
    zpiece0=zpiece0(:,:,inx);
    zlin   = mean(zlin0,3);
    zpiece = mean(zpiece0,3);
    zpiecemin = min(zpiece0,[],3);
    zpiecemax = max(zpiece0,[],3);
    zlinmin = min(zlin0,[],3);
    zlinmax = max(zlin0,[],3);

    for i=1:M_.endo_nbr
        for j=1:forecast
            [post_mean(j,1), post_median(j,1), post_var(j,1), hpd_interval(j,:), post_deciles(j,:)] = posterior_moments(squeeze(zlin0(j,i,:)),options_.forecasts.conf_sig);
        end
        oo_.occbin.linear_forecast.Mean.(M_.endo_names{i})=post_mean;
        oo_.occbin.linear_forecast.Median.(M_.endo_names{i})=post_median;
        oo_.occbin.linear_forecast.Var.(M_.endo_names{i})=post_var;
        oo_.occbin.linear_forecast.HPDinf.(M_.endo_names{i})=hpd_interval(:,1);
        oo_.occbin.linear_forecast.HPDsup.(M_.endo_names{i})=hpd_interval(:,2);
        oo_.occbin.linear_forecast.Deciles.(M_.endo_names{i})=post_deciles;
        oo_.occbin.linear_forecast.Min.(M_.endo_names{i})=zlinmin(:,i);
        oo_.occbin.linear_forecast.Max.(M_.endo_names{i})=zlinmax(:,i);
        for j=1:forecast
            [post_mean(j,1), post_median(j,1), post_var(j,1), hpd_interval(j,:), post_deciles(j,:)] = posterior_moments(squeeze(zpiece0(j,i,:)),options_.forecasts.conf_sig);
        end
        oo_.occbin.forecast.Mean.(M_.endo_names{i})=post_mean;
        oo_.occbin.forecast.Median.(M_.endo_names{i})=post_median;
        oo_.occbin.forecast.Var.(M_.endo_names{i})=post_var;
        oo_.occbin.forecast.HPDinf.(M_.endo_names{i})=hpd_interval(:,1);
        oo_.occbin.forecast.HPDsup.(M_.endo_names{i})=hpd_interval(:,2);
        oo_.occbin.forecast.Deciles.(M_.endo_names{i})=post_deciles;
        oo_.occbin.forecast.Min.(M_.endo_names{i})=zpiecemin(:,i);
        oo_.occbin.forecast.Max.(M_.endo_names{i})=zpiecemax(:,i);
        %   eval([M_.endo_names{i},'_ss=zdatass(i);']);
    end

else
    SHOCKS = zeros(forecast,M_.exo_nbr);
    if ~isempty(SHOCKS0)
        for j=1:length(SHOCKS0)
            SHOCKS(:,inds(j))=SHOCKS1(j,:);
        end
    end
    options_.occbin.simul.endo_init = M_.endo_histval(:,1)-oo_.steady_state;
    options_.occbin.simul.init_regime = opts.frcst_regimes;
    options_.occbin.simul.init_violvecbool = [];
    options_.occbin.simul.irfshock = M_.exo_names;
    options_.occbin.simul.SHOCKS = SHOCKS;
    [~, out] = occbin.solver(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);

    zlin=out.linear;
    zpiece=out.piecewise;
    frcst_regime_history=out.regime_history;
    error_flag=out.error_flag;
    for i=1:M_.endo_nbr
        oo_.occbin.linear_forecast.Mean.(M_.endo_names{i})= zlin(:,i);
        oo_.occbin.forecast.Mean.(M_.endo_names{i})= zpiece(:,i);
        oo_.occbin.forecast.HPDinf.(M_.endo_names{i})= nan;
        oo_.occbin.forecast.HPDsup.(M_.endo_names{i})= nan;
    end

end

oo_.occbin.forecast.regimes=frcst_regime_history;

