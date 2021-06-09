% This mod file compares the functionality of Dynare's pruned_state_space.m with the
% external Dynare pruning toolbox of Andreasen, Fernández-Villaverde and Rubio-Ramírez (2018):
% "The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications",
% Review of Economic Studies, Volume 85, Issue 1, Pages 1–49.
% The model under study is taken from An and Schorfheide (2007): "Bayesian Analysis of DSGE Models",
% Econometric Reviews, Volume 26, Issue 2-4, Pages 113-172.
% Note that we use version 2 of the toolbox, i.e. the one which is called
% "Third-order GMM estimate package for DSGE models (version 2)" and can be
% downloaded from https://sites.google.com/site/mandreasendk/home-1
%
% Created by @wmutschl (Willi Mutschler, willi@mutschler.eu)
%
% =========================================================================
% Copyright (C) 2020 Dynare Team
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

% set this to 1 if you want to recompute using the Andreasen et al toolbox
% otherwise the results are loaded from Andreasen_et_al_2018_Dynare44Pruning_v2.mat
@#define Andreasen_et_al_toolbox = 0

var YGR INFL INT 
    c p R g y z; %if ordering of var is changed comparison code below needs to be adapted
varexo e_r e_g e_z;
parameters tau nu kap cyst psi1 psi2 rhor rhog rhoz rrst pist gamst;

tau   = 2;
nu    = 0.1;
kap   = 0.33;
cyst  = 0.85;
psi1  = 1.5;
psi2  = 0.125;
rhor  = 0.75;
rhog  = 0.95;
rhoz  = 0.9;
rrst  = 1;
pist  = 3.2;
gamst = 0.55;
sig_r = .2;
sig_g = .6;
sig_z = .3;

model;
#pist2 = exp(pist/400);
#rrst2 = exp(rrst/400);
#bet   = 1/rrst2;
#phi   = tau*(1-nu)/nu/kap/pist2^2;
#gst   = 1/cyst;
#cst   = (1-nu)^(1/tau);
#yst   = cst*gst;
#dy    = y-y(-1);
1 = exp(-tau*c(+1)+tau*c+R-z(+1)-p(+1));
(1-nu)/nu/phi/(pist2^2)*(exp(tau*c)-1) = (exp(p)-1)*((1-1/2/nu)*exp(p)+1/2/nu) - bet*(exp(p(+1))-1)*exp(-tau*c(+1)+tau*c+y(+1)-y+p(+1));
exp(c-y) = exp(-g) - phi*pist2^2*gst/2*(exp(p)-1)^2;
R = rhor*R(-1) + (1-rhor)*psi1*p + (1-rhor)*psi2*(y-g) + e_r;
g = rhog*g(-1) + e_g;
z = rhoz*z(-1) + e_z;
YGR = gamst+100*(dy+z);
INFL = pist+400*p;
INT = pist+rrst+4*gamst+400*R;
end;

shocks;
var e_r = sig_r^2;
var e_g = sig_g^2;
var e_z = sig_z^2;
end;

steady_state_model;
y    = 0;
R    = 0;
g    = 0;
z    = 0;
c    = 0;
p    = 0;
YGR  = gamst;
INFL = pist;
INT  = pist + rrst + 4*gamst;
end;

steady; check; model_diagnostics;

@#for orderApp in [1, 2, 3]
    stoch_simul(order=@{orderApp},pruning,irf=0,periods=0);
    pruned_state_space.order_@{orderApp} = pruned_state_space_system(M_, options_, oo_.dr, [], options_.ar, 1, 0);                                                    
    @#if Andreasen_et_al_toolbox
        addpath('Dynare44Pruning_v2/simAndMoments3order'); %provide path to toolbox
        optPruning.orderApp   = @{orderApp};
        outAndreasenetal.order_@{orderApp} = RunDynarePruning(optPruning,oo_,M_,[oo_.dr.ghx oo_.dr.ghu]);
        rmpath('Dynare44Pruning_v2/simAndMoments3order');
        close all;
    @#endif
@#endfor

@#if Andreasen_et_al_toolbox
    delete Andreasen_et_al_2018_Dynare44Pruning_v2.mat;
    pause(3);
    save('Andreasen_et_al_2018_Dynare44Pruning_v2.mat', 'outAndreasenetal')
    pause(3);
@#endif

load('Andreasen_et_al_2018_Dynare44Pruning_v2.mat')

% Make comparisons only at orders 1 and 2
for iorder = 1:3
    fprintf('ORDER %d:\n',iorder);
    pruned       = pruned_state_space.(sprintf('order_%d',iorder));
    outAndreasen = outAndreasenetal.(sprintf('order_%d',iorder));
    %make sure variable ordering is correct
    if ~isequal(M_.endo_names,[outAndreasen.label_y; outAndreasen.label_v(1:M_.nspred)])
        error('variable ordering is not the same, change declaration order');
    end
    norm_E_yx   = norm(pruned.E_y(oo_.dr.inv_order_var)    - [outAndreasen.Mean_y; outAndreasen.Mean_v(1:M_.nspred)] , Inf);
    fprintf('max(sum(abs(E[y;x]''))): %d\n',norm_E_yx);
    norm_Var_y = norm(pruned.Var_y(oo_.dr.inv_order_var(1:(M_.endo_nbr-M_.nspred)),oo_.dr.inv_order_var(1:(M_.endo_nbr-M_.nspred))) - outAndreasen.Var_y , Inf);
    fprintf('max(sum(abs(Var[y]''))):: %d\n',norm_Var_y);
    norm_Var_x = norm(pruned.Var_y(M_.nstatic+(1:M_.nspred),M_.nstatic+(1:M_.nspred)) - outAndreasen.Var_v(1:M_.nspred,1:M_.nspred) , Inf);
    fprintf('max(sum(abs(Var[x]''))): %d\n',norm_Var_x);
    norm_Corr_yi1 = norm(pruned.Corr_yi(oo_.dr.inv_order_var(1:(M_.endo_nbr-M_.nspred)),oo_.dr.inv_order_var(1:(M_.endo_nbr-M_.nspred)),1) - outAndreasen.Corr_y(:,:,1) , Inf);
    fprintf('max(sum(abs(Corr[y,y(-1)]''))): %d\n',norm_Corr_yi1);
    norm_Corr_yi2 = norm(pruned.Corr_yi(oo_.dr.inv_order_var(1:(M_.endo_nbr-M_.nspred)),oo_.dr.inv_order_var(1:(M_.endo_nbr-M_.nspred)),2) - outAndreasen.Corr_y(:,:,2) , Inf);
    fprintf('max(sum(abs(Corr[y,y(-2)]''))): %d\n',norm_Corr_yi2);
    norm_Corr_xi1 = norm(pruned.Corr_yi(M_.nstatic+(1:M_.nspred),M_.nstatic+(1:M_.nspred),1) - outAndreasen.Corr_v(1:M_.nspred,1:M_.nspred,1) , Inf);
    fprintf('max(sum(abs(Corr[x,x(-1)]''))): %d\n',norm_Corr_xi1);
    norm_Corr_xi2 = norm(pruned.Corr_yi(M_.nstatic+(1:M_.nspred),M_.nstatic+(1:M_.nspred),2) - outAndreasen.Corr_v(1:M_.nspred,1:M_.nspred,2) , Inf);
    fprintf('max(sum(abs(Corr[x,x(-2)]''))): %d\n',norm_Corr_xi2);
    
    if iorder < 3 && any([norm_E_yx norm_Var_y norm_Var_x norm_Corr_yi1 norm_Corr_yi2 norm_Corr_xi1 norm_Corr_xi2] > 1e-5)
        error('Something wrong with pruned_state_space.m compared to Andreasen et al 2018 Toolbox v2 at order %d.',iorder);
    end
end
% skipline();
% fprintf('Note that at third order, there is an error in the computation of Var_z in Andreasen et al (2018)''s toolbox, @wmutschl is in contact to clarify this.\n');
% fprintf('EXAMPLE:\n')
% fprintf('  Consider Var[kron(kron(xf,xf),xf)] = E[kron(kron(kron(kron(kron(xf,xf),xf),xf),xf),xf)] - E[kron(kron(xf,xf),xf)]*E[kron(kron(xf,xf),xf)].''\n');
% fprintf('  Now note that xf=hx*xf(-1)+hu*u is Gaussian, that is E[kron(kron(xf,xf),xf)]=0, and Var[kron(kron(xf,xf),xf)] are the sixth-order product moments\n');
% fprintf('  which can be computed using the prodmom.m function by providing E[xf*xf''] as covariance matrix.\n');
% fprintf('  In order to replicate this you have to change UnconditionalMoments_3rd_Lyap.m to also output Var_z.\n')
% 
% dynare_nx    = M_.nspred;
% dynare_E_xf2 = pruned_state_space.order_3.Var_z(1:dynare_nx,1:dynare_nx);
% dynare_E_xf6 = pruned_state_space.order_3.Var_z((end-dynare_nx^3+1):end,(end-dynare_nx^3+1):end);
% dynare_E_xf6 = dynare_E_xf6(:);
%         
% Andreasen_nx    = M_.nspred+M_.exo_nbr;
% Andreasen_E_xf2 = outAndreasenetal.order_3.Var_z(1:Andreasen_nx,1:Andreasen_nx);
% Andreasen_E_xf6 = outAndreasenetal.order_3.Var_z((end-Andreasen_nx^3+1):end,(end-Andreasen_nx^3+1):end);
% Andreasen_E_xf6 = Andreasen_E_xf6(:);
% 
% fprintf('Second-order product moments of xf and u are the same:\n')
% norm_E_xf2 = norm(dynare_E_xf2-Andreasen_E_xf2(1:M_.nspred,1:M_.nspred),Inf)
% norm_E_uu  = norm(M_.Sigma_e-Andreasen_E_xf2(M_.nspred+(1:M_.exo_nbr),M_.nspred+(1:M_.exo_nbr)),Inf)
% 
% % Compute unique sixth-order product moments of xf, i.e. unique(E[kron(kron(kron(kron(kron(xf,xf),xf),xf),xf),xf)],'stable')
% dynare_nx6     = dynare_nx*(dynare_nx+1)/2*(dynare_nx+2)/3*(dynare_nx+3)/4*(dynare_nx+4)/5*(dynare_nx+5)/6;
% dynare_Q6Px    = Q6_plication(dynare_nx);
% dynare_COMBOS6 = flipud(allVL1(dynare_nx, 6)); %all possible (unique) combinations of powers that sum up to six
% dynare_true_E_xf6 = zeros(dynare_nx6,1); %only unique entries
% for j6 = 1:size(dynare_COMBOS6,1)
%     dynare_true_E_xf6(j6) = prodmom(dynare_E_xf2, 1:dynare_nx, dynare_COMBOS6(j6,:));
% end
% dynare_true_E_xf6 = dynare_Q6Px*dynare_true_E_xf6; %add duplicate entries
% norm_dynare_E_xf6 = norm(dynare_true_E_xf6 - dynare_E_xf6, Inf);
% 
% Andreasen_nx6     = Andreasen_nx*(Andreasen_nx+1)/2*(Andreasen_nx+2)/3*(Andreasen_nx+3)/4*(Andreasen_nx+4)/5*(Andreasen_nx+5)/6;
% Andreasen_Q6Px    = Q6_plication(Andreasen_nx);
% Andreasen_COMBOS6 = flipud(allVL1(Andreasen_nx, 6)); %all possible (unique) combinations of powers that sum up to six
% Andreasen_true_E_xf6   = zeros(Andreasen_nx6,1); %only unique entries
% for j6 = 1:size(Andreasen_COMBOS6,1)
%     Andreasen_true_E_xf6(j6) = prodmom(Andreasen_E_xf2, 1:Andreasen_nx, Andreasen_COMBOS6(j6,:));
% end
% Andreasen_true_E_xf6 = Andreasen_Q6Px*Andreasen_true_E_xf6; %add duplicate entries
% norm_Andreasen_E_xf6 = norm(Andreasen_true_E_xf6 - Andreasen_E_xf6, Inf);
% 
% fprintf('Sixth-order product moments of xf and u are not the same!\n');
% fprintf('    Dynare maximum absolute deviations of sixth-order product moments of xf: %d\n',norm_dynare_E_xf6)
% fprintf('    Andreasen et al maximum absolute deviations of sixth-order product moments of xf: %d\n',norm_Andreasen_E_xf6)
% skipline();
% fprintf('Note that the standard deviations are set quite high to make the numerical differences more apparent.\n');
