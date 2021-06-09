% DSGE model based on replication files of
% Andreasen, Fernandez-Villaverde, Rubio-Ramirez (2018), The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications, Review of Economic Studies, 85, p. 1-49
% Original code by Martin M. Andreasen, Jan 2016
% Adapted for Dynare by Willi Mutschler (@wmutschl, willi@mutschler.eu), Jan 2021
% =========================================================================
% Copyright (C) 2021 Dynare Team
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
% =========================================================================

%--------------------------------------------------------------------------
% Variable declaration
%--------------------------------------------------------------------------
var
ln_k
ln_s
ln_a
ln_g
ln_d

ln_c
ln_r
ln_pai
ln_h
ln_q
ln_evf
ln_iv
ln_x2
ln_la
ln_goy
ln_Esdf

xhr20
xhr40
Exhr

@#for i in 1:40
ln_p@{i}
@#endfor

Obs_Gr_C
Obs_Gr_I
Obs_Infl
Obs_r1
Obs_r40
Obs_xhr40
Obs_GoY
Obs_hours
;

predetermined_variables ln_k ln_s;

varobs Obs_Gr_C Obs_Gr_I Obs_Infl Obs_r1 Obs_r40 Obs_xhr40 Obs_GoY Obs_hours;

%--------------------------------------------------------------------------
% Exogenous shocks
%--------------------------------------------------------------------------
varexo 
eps_a
eps_d
eps_g
;

%--------------------------------------------------------------------------
% Parameter declaration
%--------------------------------------------------------------------------
parameters 
BETTA
B
INHABIT
H
PHI1
PHI2
RRA
PHI4
KAPAone
KAPAtwo
DELTA
THETA
ETA
ALFA
CHI
RHOR
BETTAPAI
BETTAY
MYYPS
MYZ
RHOA
%STDA
RHOG
%STDG
RHOD
%STDD
CONSxhr40
BETTAxhr
BETTAxhr40
CONSxhr20
PAI
GAMA
GoY

%auxiliary
PHIzero
AA
PHI3
negVf
;


%--------------------------------------------------------------------------
% Model equations
%--------------------------------------------------------------------------
% Based on DSGE_model_NegVf_yieldCurve.m and DSGE_model_PosVf_yieldCurve.m
% Note that we include an auxiliary parameter negVf to distinguish whether
% the steady state value function is positive (negVf=0) or negative (negVf=1).
% This parameter is endogenously determined in the steady_state_model block.

model;
%--------------------------------------------------------------------------
% Auxiliary expressions
%--------------------------------------------------------------------------
% do exp transform such that variables are logged variables
@#for var in [ "k", "s", "c", "r", "a", "g", "d", "pai", "h", "q", "evf", "iv", "x2", "la", "goy", "Esdf" ]
#@{var}_ba1  = exp(ln_@{var}(-1));
#@{var}_cu   = exp(ln_@{var});
#@{var}_cup  = exp(ln_@{var}(+1));
@#endfor
@#for i in 1:40
#p@{i}_cu   = exp(ln_p@{i});
#p@{i}_cup  = exp(ln_p@{i}(+1));
@#endfor
% these variables are not transformed
#xhr20_cu  = xhr20;
#xhr20_cup = xhr20(+1);
#xhr40_cu  = xhr40;
#xhr40_cup = xhr40(+1);
#Exhr_cu   = Exhr;
#Exhr_cup  = Exhr(+1);

% auxiliary steady state variables
#K  = exp(steady_state(ln_k));
#IV = exp(steady_state(ln_iv));
#C  = exp(steady_state(ln_c));
#Y  = (C + IV)/(1-GoY);
#R  = exp(steady_state(ln_r));
#G  = Y-C-IV;

#removeMeanXhr = 1;

% The atemporal relations if possible
% No stochastic trend in investment specific shocks
#myyps_cu  = MYYPS;
#myyps_cup = MYYPS;

% No stochastic trend in non-stationary technology shocks
#myz_cu    = MYZ;
#myz_cup   = MYZ;

% Defining myzstar
#MYZSTAR    = MYYPS^(THETA/(1-THETA))*MYZ;
#myzstar_cu = myyps_cu ^(THETA/(1-THETA))*myz_cu;
#myzstar_cup= myyps_cup^(THETA/(1-THETA))*myz_cup;

% The expression for the value function (only valid for deterministic trends!)
% Note that we make use of auxiliary parameter negVf to switch signs
#mvf_cup   = -negVf*(d_cup/(1-PHI2)*((c_cup-B*c_cu*MYZSTAR^-1)^(1-PHI2)-1) + d_cup*PHIzero/(1-PHI1)*(1-h_cup)^(1-PHI1) - negVf* BETTA*MYZSTAR^((1-PHI4)*(1-PHI2))*AA*evf_cup^(1/(1-PHI3)));

% The growth rate in lambda
#myla_cup  = (la_cup/la_cu)*(AA*evf_cu^(1/(1-PHI3))/mvf_cup)^PHI3*myzstar_cup^(-PHI2*(1-PHI4)-PHI4);

% The relation between the optimal price for the firms and the pris and inflation 
%ptil_cu   = ((1-ALFA*(pai_ba1^CHI/pai_cu )^(1-ETA))/(1-ALFA))^(1/(1-ETA));
%ptil_cup  = ((1-ALFA*(pai_cu ^CHI/pai_cup)^(1-ETA))/(1-ALFA))^(1/(1-ETA));
#ptil_cu   = ((1-ALFA*(1/pai_cu )^(1-ETA))/(1-ALFA))^(1/(1-ETA));
#ptil_cup  = ((1-ALFA*(1/pai_cup)^(1-ETA))/(1-ALFA))^(1/(1-ETA));
        
% From the households' FOC for labor
#w_cu      = d_cu*PHIzero*(1-h_cu )^(-PHI1)/la_cu;
#w_cup     = d_cu*PHIzero*(1-h_cup)^(-PHI1)/la_cup;
% Shouldn't w_cup include d_cup? Let's stick to the original (wrong) code in the replication files as results don't change dramatically... [@wmutschl]

% The firms' FOC for labor 
#mc_cu     = w_cu /((1-THETA)*a_cu *myyps_cu ^(-THETA/(1-THETA))*myz_cu ^-THETA *k_cu ^THETA*h_cu ^(-THETA));
#mc_cup    = w_cup/((1-THETA)*a_cup*myyps_cup^(-THETA/(1-THETA))*myz_cup^-THETA *k_cup^THETA*h_cup^(-THETA));

% The firms' FOC for capital 
#rk_cu     = mc_cu *THETA* a_cu *myyps_cu *myz_cu ^(1-THETA)*k_cu ^(THETA-1)*h_cu ^(1-THETA);
#rk_cup    = mc_cup*THETA* a_cup*myyps_cup*myz_cup^(1-THETA)*k_cup^(THETA-1)*h_cup^(1-THETA);

% The income identity
#y_cu = c_cu + iv_cu + g_cu;

%--------------------------------------------------------------------------
% Actual model equations
%--------------------------------------------------------------------------

[name='Expected value of the value function']
0 = -evf_cu + (mvf_cup/AA)^(1-PHI3);

[name='Households FOC for capital']
0 = -q_cu+BETTA*myla_cup/myyps_cup*(rk_cup+q_cup*(1-DELTA) -q_cup*KAPAtwo/2*(iv_cup/k_cup*myyps_cup*myzstar_cup - IV/K*MYYPS*MYZSTAR)^2 +q_cup*KAPAtwo*(iv_cup/k_cup*myyps_cup*myzstar_cup - IV/K*MYYPS*MYZSTAR)*iv_cup/k_cup*myyps_cup*myzstar_cup);

[name='Households FOC for investments']
0 = -1+q_cu*(1-KAPAone/2*(iv_cu/IV-1)^2-iv_cu/IV*KAPAone*(iv_cu/IV-1)-KAPAtwo*(iv_cu/k_cu*myyps_cu*myzstar_cu - IV/K*MYYPS*MYZSTAR));

[name='Euler equation for consumption']
0 = -1+BETTA*r_cu*exp(CONSxhr40*xhr40_cu + CONSxhr20*xhr20_cu)*myla_cup/pai_cup;

[name='Households FOC for consumption']
0 = -la_cu  + d_cu*(c_cu -B*c_ba1*myzstar_cu^-1)^(-PHI2) -INHABIT*B*BETTA*d_cup*(AA*evf_cu^(1/(1-PHI3))/mvf_cup)^PHI3*(c_cup -B*c_cu*myzstar_cup^-1)^(-PHI2)*myzstar_cup^(-PHI2*(1-PHI4)-PHI4);

[name='Nonlinear pricing, relation for x1 = (ETA-1)/ETA*x2']
0= -(ETA-1)/ETA*x2_cu+y_cu*mc_cu*ptil_cu^(-ETA-1) +ALFA*BETTA*myla_cup*(ptil_cu/ptil_cup)^(-ETA-1)*(1/pai_cup)^(-ETA)*(ETA-1)/ETA*x2_cup*myzstar_cup;
  
[name='Nonlinear pricing, relation for x2']
0=-x2_cu+y_cu*ptil_cu^-ETA +ALFA*BETTA*myla_cup*(ptil_cu/ptil_cup)^(-ETA)*(1/pai_cup)^(1-ETA)*x2_cup*myzstar_cup;

[name='Nonlinear pricing, relation for s']
0= -s_cup+(1-ALFA)*ptil_cu^(-ETA)+ALFA*(pai_cu/1)^ETA*s_cu;
  
[name='Interest rate rule']
0 = -log(r_cu/R)+RHOR*log(r_ba1/R)+(1-RHOR)*(BETTAPAI*log(pai_cu/PAI)+BETTAY*log(y_cu/Y) + BETTAxhr*(BETTAxhr40*xhr40_cu - removeMeanXhr*Exhr_cu));

[name='Production function']
0 = -y_cu*s_cup + a_cu *(k_cu *myyps_cu ^(-1/(1-THETA))*myz_cu ^-1)^THETA*h_cu ^(1-THETA);

[name='Relation for physical capital stock']
0= -k_cup + (1-DELTA)*k_cu*(myyps_cu*myzstar_cu)^-1 + iv_cu - iv_cu*KAPAone/2*(iv_cu/IV-1)^2 - k_cu*(myyps_cu*myzstar_cu)^-1*KAPAtwo/2*(iv_cu/k_cu*myyps_cu*myzstar_cu - IV/K*MYYPS*MYZSTAR)^2;

[name='Goverment spending over output']
0=-goy_cu + g_cu/y_cu;

[name='The yield curve: p1']
0= -p1_cu + 1/r_cu;

@#for i in 2:40
[name='The yield curve: p@{i}']
0= -p@{i}_cu + BETTA*myla_cup/pai_cup*p@{i-1}_cup;
@#endfor

[name='Stochastic discount factor']
0= -Esdf_cu+ BETTA*myla_cup/pai_cup;

[name='Expected 5 year excess holding period return']
0= -xhr20_cu+ log(p19_cup) - log(p20_cu) - log(r_cu);

[name='Expected 10 year excess holding period return']
0= -xhr40_cu+ log(p39_cup) - log(p40_cu) - log(r_cu);

[name='Mean of expected excess holding period return in Taylor rule']
0= -Exhr_cu + (1-GAMA)*(BETTAxhr40*xhr40_cu) + GAMA*Exhr_cup;

[name='Exogenous process for productivity']
0 = -log(a_cu)+RHOA*log(a_ba1) + eps_a;

[name='Exogenous process for government spending']
0 = -log(g_cu/G)+RHOG*log(g_ba1/G) + eps_g;

[name='Exogenous process for discount factor shifter']
0 = -log(d_cu)+RHOD*log(d_ba1) + eps_d;

[name='Observable annualized consumption growth']
Obs_Gr_C  = 4*( ln_c -ln_c(-1) + log(MYZSTAR));

[name='Observable annualized investment growth']
Obs_Gr_I  = 4*( ln_iv - ln_iv(-1) + log(MYZSTAR)+log(MYYPS));

[name='Observable annualized inflation']
Obs_Infl  = 4*( ln_pai);

[name='Observable annualized one-quarter nominal yield']
Obs_r1    = 4*( ln_r);

[name='Observable annualized 10-year nominal yield']
Obs_r40   = 4*( -1/40*ln_p40);

[name='Observable annualized 10-year ex post excess holding period return']
Obs_xhr40 = 4*( ln_p39 - ln_p40(-1) - ln_r(-1) );

[name='Observable annualized log ratio of government spending to GDP']
Obs_GoY   = 4*( 1/4*ln_goy);

[name='Observable annualized log of hours']
Obs_hours = 4*( 1/100*ln_h);
end;


%--------------------------------------------------------------------------
% Steady State Computations
%--------------------------------------------------------------------------
% Based on DSGE_model_yieldCurve_ss.m, getPHI3.m, ObjectGMM.m
% Note that we include an auxiliary parameter negVf to distinguish whether
% the steady state value function is positive (negVf=0) or negative (negVf=1).
% This parameter is endogenously determined in the steady_state_model block.


steady_state_model;

% The growth rate in the firms' fixed costs
MYZSTARBAR       = MYYPS^(THETA/(1-THETA))*MYZ;

% The growth rate for lampda
MYLABAR          = MYZSTARBAR^(-PHI2*(1-PHI4)-PHI4);

% The relative optimal price for firms
PTILBAR 		  = ((1-ALFA*PAI^((CHI-1)*(1-ETA)))/(1-ALFA))^(1/(1-ETA)); 

% The state variable s for distortions between output and produktion
SBAR             = ((1-ALFA)*PTILBAR^(-ETA))/(1-ALFA*PAI^((1-CHI)*ETA));

% The 1-period interest rate
RBAR 			  = PAI/(BETTA*MYLABAR);

% The market price of capital
QBAR             = 1;

% The real price of renting capital
RKBAR            = QBAR*(MYYPS/(BETTA*MYLABAR)-(1-DELTA));

% The marginal costs in the firms
MCBAR            = (1-ALFA*BETTA*MYLABAR*PAI^((1-CHI)*ETA)*MYZSTARBAR)*(ETA-1)/ETA*PTILBAR/(1-ALFA*BETTA*MYLABAR*PAI^((CHI-1)*(1-ETA))*MYZSTARBAR);
            
% The capital stock 
KBAR             = H*(RKBAR/(MCBAR*THETA*MYYPS*MYZ^(1-THETA)))^(1/(THETA-1));

% The wage level 
WBAR             = MCBAR*(1-THETA)*MYYPS^(-THETA/(1-THETA))*MYZ^-THETA*(KBAR/H)^THETA;

% The level of investment
IVBAR            = KBAR - (1-DELTA)*KBAR*MYYPS^(-1/(1-THETA))*MYZ^-1; 

% The consumption level
CBAR              = ((1-GoY)*(KBAR*MYYPS^(-1/(1-THETA))*MYZ^-1)^THETA*H^(1-THETA))/SBAR-IVBAR;

% The output level 
YBAR             = (CBAR + IVBAR)/(1-GoY); 

% The value of lambda
LABAR            = (CBAR-B*CBAR*MYZSTARBAR^-1)^-PHI2 - INHABIT*B*BETTA*(CBAR-B*CBAR*MYZSTARBAR^-1)^-PHI2*MYZSTARBAR^(-PHI2*(1-PHI4)-PHI4);

% The value of PHIzero
PHIzero       = LABAR*WBAR*(1-H)^PHI1;

% The level of the value function 
VFBAR            = 1/(1-BETTA*MYZSTARBAR^((1-PHI4)*(1-PHI2)))*(1/(1-PHI2)*((CBAR-B*CBAR*MYZSTARBAR^-1)^(1-PHI2)-1)+PHIzero/(1-PHI1)*(1-H)^(1-PHI1));
UBAR             = 1/(1-PHI2)*((CBAR-B*CBAR*MYZSTARBAR^-1)^(1-PHI2)-1)+PHIzero/(1-PHI1)*(1-H)^(1-PHI1);
[AA, EVFBAR, PHI3, negVf, info]= AFVRR_steady_helper(VFBAR,RBAR,IVBAR,CBAR,KBAR,LABAR,QBAR,YBAR,  BETTA,B,PAI,H,PHIzero,PHI1,PHI2,THETA,MYYPS,MYZ,INHABIT,RRA,CONSxhr40);
% The value of X2
X2BAR            = YBAR*PTILBAR^(-ETA)/(1-BETTA*ALFA*MYLABAR*PAI^((CHI-1)*(1-ETA))*MYZSTARBAR);

% Government spending
GBAR             = GoY*YBAR;
%**************************************************************************        

% map into model variables
ln_k = log(KBAR);
ln_s = log(SBAR);
ln_c_ba1 = log(CBAR);
ln_r_ba1 = log(RBAR);
ln_a = log(1);
ln_g = log(GBAR);
ln_d = log(1);

ln_c = log(CBAR);
ln_r = log(RBAR);
ln_pai = log(PAI);
ln_h = log(H);
ln_q = log(QBAR);
ln_evf = log(EVFBAR);
ln_iv = log(IVBAR);
ln_x2 = log(X2BAR);
ln_la = log(LABAR);
ln_goy = log(GoY);
ln_Esdf = log(1/RBAR);
xhr20 = 0;
xhr40 = 0;
Exhr = 0;
% The yield curve
ln_p1 = log((1/RBAR)^1);
ln_p2 = log((1/RBAR)^2);
ln_p3 = log((1/RBAR)^3);
ln_p4 = log((1/RBAR)^4);
ln_p5 = log((1/RBAR)^5);
ln_p6 = log((1/RBAR)^6);
ln_p7 = log((1/RBAR)^7);
ln_p8 = log((1/RBAR)^8);
ln_p9 = log((1/RBAR)^9);
ln_p10 = log((1/RBAR)^10);
ln_p11 = log((1/RBAR)^11);
ln_p12 = log((1/RBAR)^12);
ln_p13 = log((1/RBAR)^13);
ln_p14 = log((1/RBAR)^14);
ln_p15 = log((1/RBAR)^15);
ln_p16 = log((1/RBAR)^16);
ln_p17 = log((1/RBAR)^17);
ln_p18 = log((1/RBAR)^18);
ln_p19 = log((1/RBAR)^19);
ln_p20 = log((1/RBAR)^20);
ln_p21 = log((1/RBAR)^21);
ln_p22 = log((1/RBAR)^22);
ln_p23 = log((1/RBAR)^23);
ln_p24 = log((1/RBAR)^24);
ln_p25 = log((1/RBAR)^25);
ln_p26 = log((1/RBAR)^26);
ln_p27 = log((1/RBAR)^27);
ln_p28 = log((1/RBAR)^28);
ln_p29 = log((1/RBAR)^29);
ln_p30 = log((1/RBAR)^30);
ln_p31 = log((1/RBAR)^31);
ln_p32 = log((1/RBAR)^32);
ln_p33 = log((1/RBAR)^33);
ln_p34 = log((1/RBAR)^34);
ln_p35 = log((1/RBAR)^35);
ln_p36 = log((1/RBAR)^36);
ln_p37 = log((1/RBAR)^37);
ln_p38 = log((1/RBAR)^38);
ln_p39 = log((1/RBAR)^39);
ln_p40 = log((1/RBAR)^40);

Obs_Gr_C  = 4*( log(MYZSTARBAR) );
Obs_Gr_I  = 4*( log(MYZSTARBAR)+log(MYYPS) );
Obs_Infl  = 4*( ln_pai );
Obs_r1    = 4*( ln_r );
Obs_r40   = 4*( -1/40*ln_p40 );
Obs_xhr40 = 4*( xhr40 );
Obs_GoY   = 4*( 1/4*ln_goy );
Obs_hours = 4*( 1/100*ln_h );
end;

%--------------------------------------------------------------------------
% Declare moments to use in estimation
%--------------------------------------------------------------------------
% These are the moments used in the paper; corresponds to momentSet=2 in the replication files

matched_moments;
%mean
Obs_Gr_C;
Obs_Gr_I;
Obs_Infl;
Obs_r1;
Obs_r40;
Obs_xhr40;
Obs_GoY;
Obs_hours;

% all variances
Obs_Gr_C*Obs_Gr_C;
Obs_Gr_I*Obs_Gr_I;
Obs_Infl*Obs_Infl;
Obs_r1*Obs_r1;
Obs_r40*Obs_r40;
Obs_xhr40*Obs_xhr40;
Obs_GoY*Obs_GoY;
Obs_hours*Obs_hours;

% covariance excluding GoY and hours
Obs_Gr_C*Obs_Gr_I;
Obs_Gr_C*Obs_Infl;
Obs_Gr_C*Obs_r1;
Obs_Gr_C*Obs_r40;
Obs_Gr_C*Obs_xhr40;
%Obs_Gr_C*Obs_GoY;
%Obs_Gr_C*Obs_hours;

Obs_Gr_I*Obs_Infl;
Obs_Gr_I*Obs_r1;
Obs_Gr_I*Obs_r40;
Obs_Gr_I*Obs_xhr40;
%Obs_Gr_I*Obs_GoY;
%Obs_Gr_I*Obs_hours;

Obs_Infl*Obs_r1;
Obs_Infl*Obs_r40;
Obs_Infl*Obs_xhr40;
%Obs_Infl*Obs_GoY;
%Obs_Infl*Obs_hours;

Obs_r1*Obs_r40;
Obs_r1*Obs_xhr40;
%Obs_r1*Obs_GoY;
%Obs_r1*Obs_hours;

Obs_r40*Obs_xhr40;
%Obs_r40*Obs_GoY;
%Obs_r40*Obs_hours;

%Obs_xhr40*Obs_GoY;
%Obs_xhr40*Obs_hours;

%Obs_GoY*Obs_hours;

%first autocovariance
Obs_Gr_C*Obs_Gr_C(-1);
Obs_Gr_I*Obs_Gr_I(-1);
Obs_Infl*Obs_Infl(-1);
Obs_r1*Obs_r1(-1);
Obs_r40*Obs_r40(-1);
Obs_xhr40*Obs_xhr40(-1);
Obs_GoY*Obs_GoY(-1);
Obs_hours*Obs_hours(-1);
end;

%--------------------------------------------------------------------------
% Create Data
%--------------------------------------------------------------------------
@#ifdef CreateData
verbatim;
% From 1961Q3 to 2007Q4
DataUS = xlsread('Data_PruningPaper_v5.xlsx','Data_used','E3:M188');
%                                         ANNUALIZED (except for hours and GoY)
%           1   2      3     4     5        6      7    8     9        
% Lables: Date Gr_C   Gr_I  GoY   hours   Infl_C   r1   r40  xhr40   
%label_data = {'Gr_C   ', 'Gr_I   ','Infl  ', 'r1    ', 'r40   ', 'xhr40  ','GoY    ', 'hours  '};
%DataUS     = [DataUS(:,2:3) DataUS(:,6:8)  DataUS(:,9) log(DataUS(:,4)) 4*log(DataUS(:,5))/100];
Obs_Gr_C  = DataUS(:,2);
Obs_Gr_I  = DataUS(:,3);
Obs_Infl  = DataUS(:,6);
Obs_r1    = DataUS(:,7);
Obs_r40   = DataUS(:,8);
Obs_xhr40 = DataUS(:,9);
Obs_GoY   = log(DataUS(:,4));
Obs_hours = 4*log(DataUS(:,5))/100;

save('AFVRR_data.mat','Obs_Gr_C','Obs_Gr_I','Obs_Infl','Obs_r1','Obs_r40','Obs_xhr40','Obs_GoY','Obs_hours');
pause(1);
end;
@#endif