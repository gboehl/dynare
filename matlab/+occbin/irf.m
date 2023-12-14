function [oo_] = irf(M_,oo_,options_)

shocknames = options_.occbin.irf.exo_names;
shocksigns = options_.occbin.irf.shocksigns;
shocksize = options_.occbin.irf.shocksize;
t0 = options_.occbin.irf.t0;
options_.occbin.simul.init_regime = options_.occbin.irf.init_regime;
options_.occbin.simul.check_ahead_periods = options_.occbin.irf.check_ahead_periods;
options_.occbin.simul.maxit = options_.occbin.irf.maxit;
options_.occbin.simul.periods = options_.irf;

% Run inital conditions + other shocks
if t0 == 0
    shocks0= zeros(options_.occbin.simul.periods+1,M_.exo_nbr);
    options_.occbin.simul.endo_init = [];
else
    % girf conditional to smoothed states in t0 and shocks in t0+1
    shocks0= [oo_.occbin.smoother.etahat(:,t0+1)'; zeros(options_.occbin.simul.periods,M_.exo_nbr)];
    options_.occbin.simul.SHOCKS=shocks0;
    options_.occbin.simul.endo_init = oo_.occbin.smoother.alphahat(oo_.dr.inv_order_var,t0);
end
options_.occbin.simul.SHOCKS=shocks0;
[~, out0] = occbin.solver(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
zlin0 = out0.linear;
zpiece0 = out0.piecewise;

% Select shocks of interest
jexo_all =zeros(size(shocknames,1),1);

for i=1:length(shocknames)
    jexo_all(i) = strmatch(shocknames{i},M_.exo_names,'exact');
end

oo_.occbin.linear_irfs = struct();
oo_.occbin.irfs = struct();
% Set shock size
if isempty(shocksize)
    %         if isfield(oo_.posterior_mode.shocks_std,M_.exo_names{jexo})
    shocksize = sqrt(diag(M_.Sigma_e(jexo_all,jexo_all))); %oo_.posterior_mode.shocks_std.(M_.exo_names{jexo});
    if any(shocksize < 1.e-9)
        shocksize(shocksize < 1.e-9) = 0.01;
    end
end
if numel(shocksize)==1
    shocksize=repmat(shocksize,[length(shocknames),1]);
end

% Run IRFs
for counter = 1:length(jexo_all)
  
    jexo = jexo_all(counter);
    
    if ~options_.noprint
        fprintf('Producing GIRFs for shock %s. Simulation %d out of %d. \n',M_.exo_names{jexo},counter,size(jexo_all,1));
    end 
    
    if ismember('pos',shocksigns)
        % (+) shock
        shocks1=shocks0;
        shocks1(1,jexo)=shocks0(1,jexo)+shocksize(counter);
        if t0 == 0
            options_.occbin.simul.SHOCKS=shocks1;
            options_.occbin.simul.endo_init = [];
            [~, out_pos] = occbin.solver(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
        else
            options_.occbin.simul.SHOCKS=shocks1;
            options_.occbin.simul.endo_init = oo_.occbin.smoother.alphahat(oo_.dr.inv_order_var,t0);
            [~, out_pos] = occbin.solver(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
        end
        if out_pos.error_flag
            warning('Occbin error.')
            return
        end
        zlin_pos = out_pos.linear;
        zpiece_pos = out_pos.piecewise;
        % Substract inital conditions + other shocks
        zlin_pos_diff       = zlin_pos-zlin0;
        zpiece_pos_diff      = zpiece_pos-zpiece0;
    end
    
    if ismember('neg',shocksigns)
        % (-) shock
        shocks_1=shocks0;
        shocks_1(1,jexo)=shocks0(1,jexo)-shocksize(counter);
        if t0 == 0
            options_.occbin.simul.SHOCKS=shocks_1;
            options_.occbin.simul.endo_init = [];
            [~, out_neg] = occbin.solver(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
        else
            options_.occbin.simul.SHOCKS=shocks_1;
            options_.occbin.simul.endo_init = oo_.occbin.smoother.alphahat(oo_.dr.inv_order_var,t0);
            [~, out_neg] = occbin.solver(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
        end
        if out_neg.error_flag
            warning('Occbin error.')
            return
        end
        zlin_neg    = out_neg.linear;
        zpiece_neg  = out_neg.piecewise;
        zlin_neg_diff    = zlin_neg-zlin0;
        zpiece_neg_diff  = zpiece_neg-zpiece0;
    end

    % Save
    if ~isfield(oo_.occbin,'linear_irfs')
        oo_.occbin.linear_irfs = struct();
    end
    if ~isfield(oo_.occbin,'irfs')
        oo_.occbin.irfs = struct();
    end

    for jendo=1:M_.endo_nbr
%         oo_.occbin.irfs.([M_.endo_names{jendo} '_' M_.exo_names{jexo} '1'])   = zpiece_pos (:,jendo);
%         oo_.occbin.irfs.([M_.endo_names{jendo} '_' M_.exo_names{jexo} '_1'])  = zpiece_neg (:,jendo);
%         oo_.occbin.linear_irfs.([M_.endo_names{jendo} '_' M_.exo_names{jexo} '1'])   = zlin_pos (:,jendo);
%         oo_.occbin.linear_irfs.([M_.endo_names{jendo} '_' M_.exo_names{jexo} '_1'])  = zlin_neg(:,jendo);
   
        if ismember('pos',shocksigns)
        oo_.occbin.irfs.([M_.endo_names{jendo} '_' M_.exo_names{jexo} '_pos'])   = zpiece_pos_diff (:,jendo);
        oo_.occbin.linear_irfs.([M_.endo_names{jendo} '_' M_.exo_names{jexo} '_pos'])   = zlin_pos_diff (:,jendo);
        end
        
        if ismember('neg',shocksigns)
            oo_.occbin.irfs.([M_.endo_names{jendo} '_' M_.exo_names{jexo} '_neg'])  = zpiece_neg_diff (:,jendo);
            oo_.occbin.linear_irfs.([M_.endo_names{jendo} '_' M_.exo_names{jexo} '_neg'])  = zlin_neg_diff (:,jendo);
        end
        
% %         
%         oo_.occbin.irfs0.([M_.endo_names{jendo} '_' M_.exo_names{jexo} '1'])   = zpiece0(:,jendo);
%         oo_.occbin.linear_irfs0.([M_.endo_names{jendo} '_' M_.exo_names{jexo} '1'])  = zlin0(:,jendo);
%         oo_.occbin.irfs0.([M_.endo_names{jendo} '_' M_.exo_names{jexo} '_1'])   = zpiece0(:,jendo);
%         oo_.occbin.linear_irfs0.([M_.endo_names{jendo} '_' M_.exo_names{jexo} '_1'])  = zlin0(:,jendo);

    end
end



end

