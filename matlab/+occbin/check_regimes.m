function [TT, RR, CC, regime_history] = check_regimes(TT, RR, CC, opts_regime, M_, options_, dr ,steady_state, exo_steady_state, exo_det_steady_state)
%function [TT, RR, CC, regime_history] = check_regimes(TT, RR, CC, opts_regime, M_, options_ dr ,steady_state, exo_steady_state, exo_det_steady_state)
%
% INPUTS
% - TT            [N by N]          transition matrix of state space
% - RR            [N by N_exo]      shock impact matrix of state space
% - CC            [N by 1]          constant of state space
% - opts_regime_  [structure]       structure describing the regime
% - M_            [structure]       Matlab's structure describing the model
% - options_      [structure]       Matlab's structure describing the current options
% - dr                  [structure]     Reduced form model.
% - endo_steady_state   [vector]        steady state value for endogenous variables
% - exo_steady_state    [vector]        steady state value for exogenous variables
% - exo_det_steady_state    [vector]    steady state value for exogenous deterministic variables                                    
%
% OUTPUTS
% - TT              [N by N]            transition matrix of state space for each period
% - RR              [N by N_exo by T]   shock impact matrix of state space for each period
% - CC              [N by N_exo by T]   constant of state space for each period
% - regime_history  [structure]         contains the regime history


% Copyright © 2021 Dynare Team
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

binding_indicator = opts_regime.binding_indicator;
regime_history = opts_regime.regime_history;
gend = size(binding_indicator,1);
if gend ==0
    gend = length(regime_history)+1;
end

if M_.occbin.constraint_nbr==1
    base_regime.regime = 0;
    base_regime.regimestart = 1;
else
    base_regime.regime1 = 0;
    base_regime.regimestart1 = 1;
    base_regime.regime2 = 0;
    base_regime.regimestart2 = 1;
end

if isempty(regime_history)
    % set default unconstrained regimes
    regime_history = base_regime;
    regime_history(2:gend-1) = base_regime;
end
init_=0;
if ndim(TT)==2
    CC = zeros([size(RR,1),gend ]);
    TT = repmat(TT, [1 1 gend]);
    RR = repmat(RR, [1 1 gend]);
    init_=1;
end
opts_simul = options_.occbin.simul;
opts_simul.full_output=0;
opts_simul.piecewise_only=1;

for tp=1:gend-1
    change_it = 0;
    if ~isempty(binding_indicator)
        if M_.occbin.constraint_nbr==1
            if regime_history(tp).regime(1) > binding_indicator(tp+1,1)
                regime_history(tp).regime = [0 regime_history(tp).regime];
                regime_history(tp).regimestart = [1 regime_history(tp).regimestart+1];
                change_it = 1;
            end
            if regime_history(tp).regime(1) < binding_indicator(tp+1,1)
                regime_history(tp).regime = [1 regime_history(tp).regime];
                regime_history(tp).regimestart = [1 regime_history(tp).regimestart+1];
                change_it = 1;
            end
            
        else
            if regime_history(tp).regime1(1) > binding_indicator(tp+1,1)
                regime_history(tp).regime1 = [0 regime_history(tp).regime1];
                regime_history(tp).regimestart1 = [1 regime_history(tp).regimestart1+1];
                change_it = 1;
            end
            if regime_history(tp).regime1(1) < binding_indicator(tp+1,1)
                regime_history(tp).regime1 = [1 regime_history(tp).regime1];
                regime_history(tp).regimestart1 = [1 regime_history(tp).regimestart1+1];
                change_it = 1;
            end
            if regime_history(tp).regime2(1) > binding_indicator(tp+1,2)
                regime_history(tp).regime2 = [0 regime_history(tp).regime2];
                regime_history(tp).regimestart2 = [1 regime_history(tp).regimestart2+1];
                change_it = 1;
            end
            if regime_history(tp).regime2(1) < binding_indicator(tp+1,2)
                regime_history(tp).regime2 = [1 regime_history(tp).regime2];
                regime_history(tp).regimestart2 = [1 regime_history(tp).regimestart2+1];
                change_it = 1;
            end
        end
    end
    if change_it || init_
        check_it = 0;
        jt=1;
        while jt<=tp-1 && check_it==0
            check_it = isequal(regime_history(tp),regime_history(jt));
            jt = jt+1;
        end
        if ~check_it && isequal(regime_history(tp),base_regime) && init_
            check_it = true;
            jt=1;
        end
        
        if check_it %|| tp==1 % the matrices for regime have been already computed
            TT(:,:,tp+1) = TT(:,:,jt);
            RR(:,:,tp+1) = RR(:,:,jt);
            CC(:,tp+1) = CC(:,jt);
        else
            if M_.occbin.constraint_nbr==1
                nperi = max(regime_history(tp).regimestart);
            else
                nperi = max([regime_history(tp).regimestart1 regime_history(tp).regimestart2]);
            end
            
            opts_simul.endo_init_=zeros(size(RR,1),1);
            opts_simul.init_binding_indicator=[];
            opts_simul.init_regime=regime_history(tp);
            opts_simul.SHOCKS=zeros(1,size(RR,2));
            opts_simul.maxit=1;
            opts_simul.periods=nperi;
            options_.occbin.simul=opts_simul;
            [~, ~, ss] = occbin.solver(M_,options_,dr,steady_state,exo_steady_state,exo_det_steady_state);
            TT(:,:,tp+1) = ss.T(dr.order_var,dr.order_var);
            RR(:,:,tp+1) = ss.R(dr.order_var,:);
            CC(:,tp+1) = ss.C(dr.order_var);
        end
        
    end
end
