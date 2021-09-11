/*
 * This file implements the baseline New Keynesian model of Jordi Gal� (2015): Monetary Policy, Inflation,
 * and the Business Cycle, Princeton University Press, Second Edition, Chapter 3
 *
 * Note that this mod-file implements the non-linear first order conditions and that the IRFs show the log-deviations
 * from steady state.
 *
 * THIS MOD-FILE REQUIRES DYNARE 4.5 OR HIGHER
 *
 * Notes:
 *  - in the LOM for the discount rate shock z the shock enters with a minus sign in this mod-file to generate the 
 *      IRF to a -0.5% shock
 *
 * This implementation was written by Johannes Pfeifer. In case you spot mistakes,
 * email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model. 
 */

/*
 * Copyright (C) 2016-20 Johannes Pfeifer
 * Copyright (C) 2020 Dynare Team
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * It is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For a copy of the GNU General Public License,
 * see <https://www.gnu.org/licenses/>.
 */

var C               ${C}$           (long_name='Consumption')
    W_real          ${\frac{W}{P}}$ (long_name='Real Wage')
    Pi              ${\Pi}$         (long_name='inflation')
    A               ${A}$           (long_name='AR(1) technology process')
    N               ${N}$           (long_name='Hours worked')
    R               ${R^n}$         (long_name='Nominal Interest Rate') 
    realinterest    ${R^{r}}$       (long_name='Real Interest Rate')
    Y               ${Y}$           (long_name='Output') 
    Q               ${Q}$           (long_name='Bond price')
    Z               ${Z}$           (long_name='AR(1) preference shock process')
    S               ${S}$           (long_name='Price dispersion')
    Pi_star         ${\Pi^*}$       (long_name='Optimal reset price')
    x_aux_1         ${x_1}$         (long_name='aux. var. 1 recursive price setting')
    x_aux_2         ${x_2}$         (long_name='aux. var. 2 recursive price setting')
    MC              ${mc}$          (long_name='real marginal costs')
    M_real          ${M/P}$         (long_name='real money stock')
    i_ann           ${i^{ann}}$     (long_name='annualized nominal interest rate')
    pi_ann          ${\pi^{ann}}$   (long_name='annualized inflation rate')
    r_real_ann      ${r^{r,ann}}$   (long_name='annualized real interest rate')
    P               ${P}$           (long_name='price level')
    log_m_nominal   ${log(M)}$      (long_name='log nominal money stock')
    log_y           ${log(Y)}$      (long_name='log output')
    log_W_real      ${log(W/P)}$    (long_name='log real wage')
    log_N           ${log(N)}$      (long_name='log hours')
    log_P           ${log(P)}$      (long_name='log price level')
    log_A           ${log(A)}$      (long_name='log technology level')
    log_Z           ${log(Z)}$      (long_name='log preference shock')
    y_hat
    pi
    ;

varexo eps_a        ${\varepsilon_a}$   (long_name='technology shock')
       eps_z        ${\varepsilon_z}$   (long_name='preference shock')
       ;   

parameters alppha       ${\alpha}$      (long_name='capital share')
    betta               ${\beta}$       (long_name='discount factor')
    rho_a               ${\rho_a}$      (long_name='autocorrelation technology shock')
    rho_z               ${\rho_{z}}$    (long_name='autocorrelation monetary demand shock')
    siggma              ${\sigma}$      (long_name='inverse EIS')
    varphi              ${\varphi}$     (long_name='inverse Frisch elasticity')
    phi_pi              ${\phi_{\pi}}$  (long_name='inflation feedback Taylor Rule')
    phi_y               ${\phi_{y}}$    (long_name='output feedback Taylor Rule')
    eta                 ${\eta}$        (long_name='semi-elasticity of money demand')
    epsilon             ${\epsilon}$    (long_name='demand elasticity')
    theta               ${\theta}$      (long_name='Calvo parameter')
    tau                 ${\tau}$      (long_name='labor subsidy')
    ;
    
%----------------------------------------------------------------
% Parametrization, p. 67  and p. 113-115
%----------------------------------------------------------------
siggma = 1;
varphi=5;
phi_pi = 1.5;
phi_y  = 0.125;
theta=3/4;
rho_z  = 0.5;
rho_a  = 0.9;
betta  = 0.99;
eta  =3.77; %footnote 11, p. 115
alppha=1/4;
epsilon=9;
tau=0; //1/epsilon;

%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------

model;
    [name='FOC Wages, eq. (2)']
    W_real=C^siggma*N^varphi;
    [name='Euler equation eq. (3)']
    Q=betta*(C(+1)/C)^(-siggma)*(Z(+1)/Z)/Pi(+1);
    [name='Definition nominal interest rate), p. 22 top']
    R=1/Q;
    [name='Aggregate output, above eq. (14)']
    Y=A*(N/S)^(1-alppha);
    [name='Definition Real interest rate']
    R=realinterest*Pi(+1);
%     @#if money_growth_rule==0
%         [name='Monetary Policy Rule, p. 26 bottom/eq. (22)']
%         R=1/betta*Pi^phi_pi*(Y/steady_state(Y))^phi_y;
%     @#endif
    [name='Market Clearing, eq. (15)']
    C=Y;
    [name='Technology Shock, eq. (6)']
    log(A)=rho_a*log(A(-1))+eps_a;
    [name='Preference Shock, p.54']
    log(Z)=rho_z*log(Z(-1))-eps_z;
    [name='Definition marginal cost']
    MC=W_real/((1-alppha)*Y/N*S);
    [name='LOM prices, eq. (7)']
    1=theta*Pi^(epsilon-1)+(1-theta)*(Pi_star)^(1-epsilon);
    [name='LOM price dispersion']
    S=(1-theta)*Pi_star^(-epsilon/(1-alppha))+theta*Pi^(epsilon/(1-alppha))*S(-1);
    [name='FOC price setting']
    Pi_star^(1+epsilon*(alppha/(1-alppha)))=x_aux_1/x_aux_2*(1-tau)*epsilon/(epsilon-1);
    [name='Auxiliary price setting recursion 1']
    x_aux_1=Z*C^(-siggma)*Y*MC+betta*theta*Pi(+1)^(epsilon+alppha*epsilon/(1-alppha))*x_aux_1(+1);
    [name='Auxiliary price setting recursion 2']
    x_aux_2=Z*C^(-siggma)*Y+betta*theta*Pi(+1)^(epsilon-1)*x_aux_2(+1);
    [name='Definition log output']
    log_y = log(Y);
    [name='Definition log real wage']
    log_W_real=log(W_real);
    [name='Definition log hours']
    log_N=log(N);
    [name='Annualized inflation']
    pi_ann=4*log(Pi);
    [name='Annualized nominal interest rate']
    i_ann=4*log(R);
    [name='Annualized real interest rate']
    r_real_ann=4*log(realinterest);
    [name='Real money demand, eq. (4)']
    M_real=Y/R^eta;
    [name='definition nominal money stock']
    log_m_nominal=log(M_real*P);
    [name='Definition price level']
    Pi=P/P(-1);
    [name='Definition log price level']
    log_P=log(P);
    [name='Definition log TFP']
    log_A=log(A);
    [name='Definition log preference']
    log_Z=log(Z);
    [mcp='a']
    y_hat=log(Y)-STEADY_STATE(log(Y));
    
    pi=log(Pi)-STEADY_STATE(log(Pi));
end;

%----------------------------------------------------------------
%  Steady state values
%---------------------------------------------------------------

steady_state_model;
A=1;
Z=1;
S=1;
Pi_star=1;
P=1;
MC=(epsilon-1)/epsilon/(1-tau);
R=1/betta;
Pi=1;
Q=1/R;
realinterest=R;
N=((1-alppha)*MC)^(1/((1-siggma)*alppha+varphi+siggma));
C=A*N^(1-alppha);
W_real=C^siggma*N^varphi;
Y=C;
money_growth=0;
money_growth_ann=0;
nu=0;
x_aux_1=C^(-siggma)*Y*MC/(1-betta*theta*Pi^(epsilon/(1-alppha)));
x_aux_2=C^(-siggma)*Y/(1-betta*theta*Pi^(epsilon-1));
log_y = log(Y);
log_W_real=log(W_real);
log_N=log(N);
pi_ann=4*log(Pi);
i_ann=4*log(R);
r_real_ann=4*log(realinterest);
M_real=Y/R^eta;
log_m_nominal=log(M_real*P);
log_P=log(P);
log_A=0;
log_Z=0;
end;

%----------------------------------------------------------------
%  define shock variances
%---------------------------------------------------------------

shocks;
var eps_a  = 0.5^2; //unit shock to preferences 
end;

% steady;
% check;
% stoch_simul;
planner_objective 0.5*((siggma+(varphi+alppha)/(1-alppha))*y_hat^2+epsilon/0.0215*pi^2)/100;
discretionary_policy(order=1,instruments=(R),irf=20,planner_discount=betta, periods=0) y_hat pi_ann log_y log_N log_W_real log_P;

temp=load(['Gali_2015_chapter_3' filesep 'Output' filesep 'Gali_2015_chapter_3_results.mat']);
if abs(oo_.planner_objective_value.unconditional-temp.oo_.planner_objective_value.unconditional)>1e-6 ... 
    || abs(oo_.planner_objective_value.conditional.zero_initial_multiplier-temp.oo_.planner_objective_value.conditional.zero_initial_multiplier)>1e-6 ... 
    ||  abs(oo_.planner_objective_value.conditional.steady_initial_multiplier-temp.oo_.planner_objective_value.conditional.steady_initial_multiplier)>1e-6 ... 
   warning('Planner objective does not match linear model')
end
if max(max(abs([temp.oo_.irfs.y_eps_a; temp.oo_.irfs.w_real_eps_a; temp.oo_.irfs.n_eps_a; temp.oo_.irfs.pi_ann_eps_a]-...
        [oo_.irfs.log_y_eps_a; oo_.irfs.log_W_real_eps_a; oo_.irfs.log_N_eps_a; oo_.irfs.pi_ann_eps_a])))>1e-6
    error('Policy is different')
end