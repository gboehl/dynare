function [ax, a1x, Px, P1x, vx, Tx, Rx, Cx, regx, info, M_, likx, etahat, alphahat, V, Fix, Kix, TTx,RRx,CCx] = ...
                    kalman_update_engine(a0,a1,P0,P1,t,data_index,Z,vv,Y,H,Qt,T0,R0,TT,RR,CC,regimes_,base_regime,d_index,M_,...
                    dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options, Fi,Ki,kalman_tol,nk)
% [ax, a1x, Px, P1x, vx, Tx, Rx, Cx, regx, info, M_, likx, etahat, alphahat, V, Fix, Kix, TTx,RRx,CCx] = kalman_update_engine(
%                                       a0,a1,P0,P1,t,data_index,Z,vv,Y,H,Qt,T0,R0,TT,RR,CC,regimes_,base_regime,d_index,M_,
%                                       dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options, Fi,Ki,kalman_tol,nk)

% Copyright © 2023 Dynare Team
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

is_multivariate = true;
Fix=[];
Kix=[];
TTx=[];
RRx=[];
CCx=[];
V=[];
if nargin>26
    is_multivariate = false;
end

use_relaxation = false;

if is_multivariate
    [ax, a1x, Px, P1x, vx, Tx, Rx, Cx, regx, info, M_, likx, etahat, alphahat, V] = occbin.kalman_update_algo_1(a0,a1,P0,P1,data_index,Z,vv,Y,H,Qt,T0,R0,TT,RR,CC,struct(),M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
else
    [ax, a1x, Px, P1x, vx, Fix, Kix, Tx, Rx, Cx, regx, info, M_, likx, alphahat, etahat,TTx,RRx,CCx] = occbin.kalman_update_algo_3(a0,a1,P0,P1,data_index,Z,vv,Fi,Ki,Y,H,Qt,T0,R0,TT,RR,CC,struct(),M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options,kalman_tol,nk);
end
info0=info;
if info
    if ~isequal(regimes_(1:2),[base_regime base_regime])
        if is_multivariate
            [ax, a1x, Px, P1x, vx, Tx, Rx, Cx, regx, info, M_, likx, etahat, alphahat, V] = occbin.kalman_update_algo_1(a0,a1,P0,P1,data_index,Z,vv,Y,H,Qt,T0,R0,TT,RR,CC,regimes_(1:2),M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
        else
            [ax, a1x, Px, P1x, vx, Fix, Kix, Tx, Rx, Cx, regx, info, M_, likx, alphahat, etahat,TTx,RRx,CCx] = occbin.kalman_update_algo_3(a0,a1,P0,P1,data_index,Z,vv,Fi,Ki,Y,H,Qt,T0,R0,TT,RR,CC,regimes_(1:2),M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options,kalman_tol,nk);
        end
    end
    info1=info;
else
    if ~isequal(regimes_(1:2),[base_regime base_regime])
        if is_multivariate
            [ax1, a1x1, Px1, P1x1, vx1, Tx1, Rx1, Cx1, regx1, info1, M_1, likx1, etahat1, alphahat1, V1] = occbin.kalman_update_algo_1(a0,a1,P0,P1,data_index,Z,vv,Y,H,Qt,T0,R0,TT,RR,CC,regimes_(1:2),M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
        else
            [ax1, a1x1, Px1, P1x1, vx1, Fix1, Kix1, Tx1, Rx1, Cx1, regx1, info1, M_1, likx1, alphahat1, etahat1,TTx1,RRx1,CCx1] = occbin.kalman_update_algo_3(a0,a1,P0,P1,data_index,Z,vv,Fi,Ki,Y,H,Qt,T0,R0,TT,RR,CC,regimes_(1:2),M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options,kalman_tol,nk);
        end
        if info1==0 && likx1<likx
            ax=ax1;
            a1x=a1x1;
            Px=Px1;
            P1x=P1x1;
            vx=vx1;
            Tx=Tx1;
            Rx=Rx1;
            Cx=Cx1;
            regx=regx1;
            info=info1;
            M_= M_1;
            likx=likx1;
            etahat=etahat1;
            alphahat=alphahat1;
            if is_multivariate
            V=V1;
            else
                Fix = Fix1;
                Kix = Kix1;
                TTx = TTx1;
                RRx = RRx1;
                CCx = CCx1;
            end
        end
    else
        if t>options_.occbin.likelihood.number_of_initial_periods_with_extra_regime_guess
            info1=0;
        else
        % may help in first 2 periods to try some other guess regime, due to
        % larger state uncertainty
            info1=1;
            options_.occbin.likelihood.brute_force_regime_guess   = true;
            options_.occbin.likelihood.loss_function_regime_guess = true;
        end

    end
end

diffstart=0;
if info==0
if M_.occbin.constraint_nbr==1
    oldstart = regimes_(1).regimestart(end);
    newstart = regx(1).regimestart(end);
    diffstart = newstart-oldstart;
else
    newstart1 = regx(1).regimestart1(end);
    newstart2 = regx(1).regimestart2(end);
    oldstart1 = regimes_(1).regimestart1(end);
    oldstart2 = regimes_(1).regimestart2(end);
    diffstart = max(newstart1-oldstart1,newstart2-oldstart2);
end
end

if options_.occbin.filter.use_relaxation && diffstart>2
    if info0==0
        % make sure we match criteria to enter further solution attempts
        info1=1;
    end
    options_.occbin.likelihood.brute_force_regime_guess   = true;
    options_.occbin.likelihood.loss_function_regime_guess = true;
    use_relaxation = true;
end

if options_.occbin.likelihood.brute_force_regime_guess && (info0 || info1) %|| (info==0 &&  ~isequal(regx(1),base_regime))

    guess_regime = [base_regime base_regime];
    options_.occbin.filter.guess_regime = true;

    use_index = 0;
    if M_.occbin.constraint_nbr==1
        for k=1:5
            guess_regime(1).regimestart=[1 5 5+4*k];
            guess_regime(1).regime=[0 1 0];
            if is_multivariate
                [ax2{1}, a1x2{1}, Px2{1}, P1x2{1}, vx2{1}, Tx2{1}, Rx2{1}, Cx2{1}, regx2{1}, info2, M_2{1}, likx2{1}, etahat2{1}, alphahat2{1}, V2{1}] = occbin.kalman_update_algo_1(a0,a1,P0,P1,data_index,Z,vv,Y,H,Qt,T0,R0,TT,RR,CC,guess_regime,M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
            else
                [ax2{1}, a1x2{1}, Px2{1}, P1x2{1}, vx2{1}, Fix2{1}, Kix2{1}, Tx2{1}, Rx2{1}, Cx2{1}, regx2{1}, info2, M_2{1}, likx2{1}, alphahat2{1}, etahat2{1},TTx2{1},RRx2{1},CCx2{1}] = occbin.kalman_update_algo_3(a0,a1,P0,P1,data_index,Z,vv,Fi,Ki,Y,H,Qt,T0,R0,TT,RR,CC,guess_regime,M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options,kalman_tol,nk);
            end
            if info2==0
                use_index= 1;
                if not(info==0 && isequal(regx2{1},regx)) && not(use_relaxation && likx2{1}>=likx)
                    % found a solution, different from previous or
                    % use_relaxation and likelihood is better
                    break
                end
            end

            guess_regime(1).regimestart=[1 1+4*k];
            guess_regime(1).regime=[1 0];
            if is_multivariate
                [ax2{2}, a1x2{2}, Px2{2}, P1x2{2}, vx2{2}, Tx2{2}, Rx2{2}, Cx2{2}, regx2{2}, info2, M_2{2}, likx2{2}, etahat2{2}, alphahat2{2}, V2{2}] = occbin.kalman_update_algo_1(a0,a1,P0,P1,data_index,Z,vv,Y,H,Qt,T0,R0,TT,RR,CC,guess_regime,M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
            else
                [ax2{2}, a1x2{2}, Px2{2}, P1x2{2}, vx2{2}, Fix2{2}, Kix2{2}, Tx2{2}, Rx2{2}, Cx2{2}, regx2{2}, info2, M_2{2}, likx2{2}, alphahat2{2}, etahat2{2},TTx2{2},RRx2{2},CCx2{2}] = occbin.kalman_update_algo_3(a0,a1,P0,P1,data_index,Z,vv,Fi,Ki,Y,H,Qt,T0,R0,TT,RR,CC,guess_regime,M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options,kalman_tol,nk);
            end
            if info2==0
                use_index = 2;
                % if use_relaxation and we are here, previous guess did not
                % improve solution, so we test for this one
            end

            if use_index
                % in case the second guess does not find a solution!
                info2=0;
                % a solution was found
                break
            end
        end
    end

    if M_.occbin.constraint_nbr==2
        for jk=0:1 % loop over other regime duration. this loop is shorter for parsimony. one may add an option ...
            for k=1:5 % loop over current regime duration
                gindex = 0;
                for jr=1:2 % loop over current regime 1 or 2
                    if jr==1
                        regstart1 = 'regimestart1';
                        reg1 = 'regime1';
                        regstart2 = 'regimestart2';
                        reg2 = 'regime2';
                    else
                        regstart1 = 'regimestart2';
                        reg1 = 'regime2';
                        regstart2 = 'regimestart1';
                        reg2 = 'regime1';
                    end
                    for kk=1:2 % loop over current regime binding in expectation vs binding in current period
                        if kk==1
                            guess_regime(1).(regstart1)=[1 5 5+4*k];
                            guess_regime(1).(reg1)=[0 1 0];
                        else
                            guess_regime(1).(regstart1)=[1 1+4*k];
                            guess_regime(1).(reg1)=[1 0];
                        end
                        for kj=1:1+1*(jk>0)
                            % loop over other regime slack or binding in current period or binding in
                            % expectation
                            if jk==0
                                % other regime is slack
                                guess_regime(1).(regstart2) = 1;
                                guess_regime(1).(reg2) = 0;
                            else % jk>0
                                if kj==1
                                    % other regime binding in current period
                                    guess_regime(1).(regstart2)=[1 1+4*jk];
                                    guess_regime(1).(reg2) = [1 0];
                                else
                                    % other regime binding in expectation
                                    guess_regime(1).(regstart2)=[1 5 5+4*jk];
                                    guess_regime(1).(reg2) = [0 1 0];
                                end
                            end
                            gindex = gindex+1;
                            if is_multivariate
                                [ax2{gindex}, a1x2{gindex}, Px2{gindex}, P1x2{gindex}, vx2{gindex}, Tx2{gindex}, Rx2{gindex}, Cx2{gindex}, regx2{gindex}, info2, M_2{gindex}, likx2{gindex}, etahat2{gindex}, alphahat2{gindex}, V2{gindex}] = occbin.kalman_update_algo_1(a0,a1,P0,P1,data_index,Z,vv,Y,H,Qt,T0,R0,TT,RR,CC,guess_regime,M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
                            else
                                [ax2{gindex}, a1x2{gindex}, Px2{gindex}, P1x2{gindex}, vx2{gindex}, Fix2{gindex}, Kix2{gindex}, Tx2{gindex}, Rx2{gindex}, Cx2{gindex}, regx2{gindex}, info2, M_2{gindex}, likx2{gindex}, alphahat2{gindex}, etahat2{gindex},TTx2{gindex},RRx2{gindex},CCx2{gindex}] = occbin.kalman_update_algo_3(a0,a1,P0,P1,data_index,Z,vv,Fi,Ki,Y,H,Qt,T0,R0,TT,RR,CC,guess_regime,M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options,kalman_tol,nk);
                            end
                            if info2==0
                                use_index= gindex;
                                if not(info==0 && isequal(regx2{gindex},regx)) && not(use_relaxation && likx2{gindex}>=likx)
                                    % found a solution, different from previous one
                                    % use_relaxation and likelihood improves
                                    break
                                end
                            end
                        end % loop over other regime slack, binding in expectation or binding in current period

                        if info2==0
                            if not(info==0 && isequal(regx2{gindex},regx)) && not(use_relaxation && likx2{gindex}>=likx)
                                % found a solution, different from previous one
                                % use_relaxation and likelihood improves
                                break
                            end
                        end

                    end % loop over current regime binding in expectation vs binding in current period

                    if info2==0
                        if not(info==0 && isequal(regx2{gindex},regx)) && not(use_relaxation && likx2{gindex}>=likx)
                            % found a solution, different from previous one
                            % use_relaxation and likelihood improves
                            break
                        end
                    end

                end % loop over current regime 1 or 2

                if use_index
                    info2=0;
                    break
                end
            end % loop over current regime duration

            if use_index
                break
            end
        end % loop over other regime duration
    end % 2 constraints


    if info2==0
        % so that we DO NOT enter IVF step
        info0=0;
        info1=0;
    end
    if info2==0 && likx2{use_index}<likx
        ax=ax2{use_index};
        a1x=a1x2{use_index};
        Px=Px2{use_index};
        P1x=P1x2{use_index};
        vx=vx2{use_index};
        Tx=Tx2{use_index};
        Rx=Rx2{use_index};
        Cx=Cx2{use_index};
        regx=regx2{use_index};
        info=info2;
        M_= M_2{use_index};
        likx=likx2{use_index};
        etahat=etahat2{use_index};
        alphahat=alphahat2{use_index};
        if is_multivariate
            V=V2{use_index};
        else
            Fix = Fix2{use_index};
            Kix = Kix2{use_index};
            TTx = TTx2{use_index};
            RRx = RRx2{use_index};
            CCx = CCx2{use_index};
        end
    end
    options_.occbin.filter.guess_regime = false;
end
if options_.occbin.likelihood.loss_function_regime_guess && (info0 || info1) %|| (info==0 &&  ~isequal(regx(1),base_regime))
    [~, out] = occbin.findmin(d_index, a0, P1, Qt, Y, Z, occbin_options.opts_simul,M_, dr,endo_steady_state,exo_steady_state,exo_det_steady_state, options_);
    if out.error_flag==0
        options_.occbin.filter.guess_regime = true;
        guess_regime=out.regime_history;
        guess_regime = [guess_regime base_regime];
        if is_multivariate
            [ax2, a1x2, Px2, P1x2, vx2, Tx2, Rx2, Cx2, regx2, info2, M_2, likx2, etahat2, alphahat2, V2] = occbin.kalman_update_algo_1(a0,a1,P0,P1,data_index,Z,vv,Y,H,Qt,T0,R0,TT,RR,CC,guess_regime,M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
        else
            [ax2, a1x2, Px2, P1x2, vx2, Fix2, Kix2, Tx2, Rx2, Cx2, regx2, info2, M_2, likx2, alphahat2, etahat2,TTx2,RRx2,CCx2] = occbin.kalman_update_algo_3(a0,a1,P0,P1,data_index,Z,vv,Fi,Ki,Y,H,Qt,T0,R0,TT,RR,CC,guess_regime,M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options,kalman_tol,nk);
        end
        options_.occbin.filter.guess_regime = false;
        if info2==0 && likx2<likx
            ax=ax2;
            a1x=a1x2;
            Px=Px2;
            P1x=P1x2;
            vx=vx2;
            Tx=Tx2;
            Rx=Rx2;
            Cx=Cx2;
            regx=regx2;
            info=info2;
            likx=likx2;
            M_= M_2;
            etahat=etahat2;
            alphahat=alphahat2;
            if is_multivariate
                V=V2;
            else
                Fix = Fix2;
                Kix = Kix2;
                TTx = TTx2;
                RRx = RRx2;
                CCx = CCx2;
            end
        end
    end
end
end
