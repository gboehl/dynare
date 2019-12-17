% =========================================================================
% Copyright (C) 2019 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
% =========================================================================


/*
 Check the policy functions obtained by perturbation at a high approximation
 order, using the Burnside (1998, JEDC) model (for which the analytical form of
 the policy function is known).

 As shown by Burnside, the policy function for yₜ is:

  yₜ =  βⁱ exp[aᵢ+bᵢ(xₜ−xₛₛ)]


 where:
                      θ²   ⎛     2ρ             1−ρ²ⁱ⎞
 — aᵢ = iθxₛₛ + σ² ─────── ⎢i − ────(1−ρⁱ) + ρ² ─────⎥
                   2(1−ρ)² ⎝    1−ρ             1−ρ² ⎠

         θρ
 — bᵢ = ───(1−ρⁱ)
        1−ρ

 — xₛₛ is the steady state of x
 — σ is the standard deviation of e.

 With some algebra, it can be shown that the derivative of yₜ at the deterministic
 steady state is equal to:

     ∂ᵐ⁺ⁿ⁺²ᵖ yₜ       ∞              (2p)!
  ──────────────── =  ∑  βⁱ bᵢᵐ⁺ⁿ ρᵐ ───── cᵢᵖ exp(iθxₛₛ)
  ∂ᵐxₜ₋₁ ∂ⁿeₜ ∂²ᵖs   ⁱ⁼¹               p!

 where:
 — s is the stochastic scale factor

           θ²   ⎛     2ρ             1−ρ²ⁱ⎞
 — cᵢ = ─────── ⎢i − ────(1−ρⁱ) + ρ² ─────⎥
        2(1−ρ)² ⎝    1−ρ             1−ρ² ⎠

 Note that derivatives with respect to an odd order for s (i.e. ∂²ᵖ⁺¹s) are always
 equal to zero.

 The policy function as returned in the oo_.dr.g_* matrices has the following properties:
 — its elements are pre-multiplied by the Taylor coefficients;
 — derivatives w.r.t. the stochastic scale factor have already been summed up;
 — symmetric elements are folded (and they are not pre-multiplied by the number of repetitions).

 As a consequence, the element gₘₙ corresponding to the m-th derivative w.r.t.
 to xₜ₋₁ and the n-th derivative w.r.t. to eₜ is given by:

          1                  ∞               cᵢᵖ
  gₘₙ = ──────      ∑        ∑  βⁱ bᵢᵐ⁺ⁿ ρᵐ ──── exp(iθxₛₛ)
        (m+n)!  0≤2p≤k-m-n  ⁱ⁼¹              p!

 where k is the order of approximation.

 */
 
@#define ORDER = 3

var y x;
varobs y;
varexo e;

parameters beta theta rho xbar;
xbar = 0.0179;
rho = -0.139;
theta = -1.5;
theta = -10;
beta = 0.95;

model;
y = beta*exp(theta*x(+1))*(1+y(+1));
x = (1-rho)*xbar + rho*x(-1)+e;
end;

shocks;
var e; stderr 0.0348;
end;

steady_state_model;
x = xbar;
y = beta*exp(theta*xbar)/(1-beta*exp(theta*xbar));
end;

estimated_params;
stderr e, normal_pdf, 0.0348,0.01;
beta, normal_pdf, 0.95, 0.01;
theta, normal_pdf, -10, 0.01;
rho, normal_pdf, -0.139, 0.01;
xbar, normal_pdf, 0.0179, 0.01;
end;

steady;check;model_diagnostics;
stoch_simul(order=@{ORDER},k_order_solver,irf=0,drop=0,periods=0,nograph);
identification(order=@{ORDER},nograph,no_identification_strength);

%make sure everything is computed at prior mean
xparam_prior = set_prior(estim_params_,M_,options_);
M_ = set_all_parameters(xparam_prior,estim_params_,M_);
[~,info,M_,options_,oo_] = resol(0,M_, options_, oo_);

indpmodel = estim_params_.param_vals(:,1);
indpstderr = estim_params_.var_exo(:,1);
indpcorr = estim_params_.corrx(:,1:2);
totparam_nbr = length(indpmodel) + length(indpstderr) + size(indpcorr,1);

%% Verify that the policy function coefficients are correct
i = 1:800;
SE_e=sqrt(M_.Sigma_e);
aux1 = rho*(1-rho.^i)/(1-rho);
aux2 = (1-rho.^(2*i))/(1-rho^2);
aux3 = 1/((1-rho)^2);
aux4 = (i-2*aux1+rho^2*aux2);
aux5=aux3*aux4;
b = theta*aux1;
c = 1/2*theta^2*SE_e^2*aux5;

%derivatives wrt to rho only
daux1_drho = zeros(1,length(i));
daux2_drho = zeros(1,length(i));
daux3_drho = 2/((1-rho)^3);
daux4_drho = zeros(1,length(i));
daux5_drho = zeros(1,length(i));
for ii = 1:length(i)
    if ii == 1
        daux1_drho(ii) = 1;
        daux2_drho(ii) = 0;
    else
        daux1_drho(ii) = rho/(rho^2 - 2*rho + 1) - 1/(rho - 1) + ((ii+1)*rho^ii)/(rho - 1) - rho^(ii+1)/(rho^2 - 2*rho + 1);
        daux2_drho(ii) = (2*rho)/(rho^4 - 2*rho^2 + 1) + (2*ii*rho^(2*ii-1))/(rho^2 - 1) - (2*rho^(2*ii+1))/(rho^4 - 2*rho^2 + 1);
    end
    daux4_drho(ii) = -2*daux1_drho(ii) + 2*rho*aux2(ii) + rho^2*daux2_drho(ii);
    daux5_drho(ii) = daux3_drho*aux4(ii) + aux3*daux4_drho(ii);
end
%derivatives of b and c wrt to all parameters
db = zeros(size(b,1),size(b,2),M_.exo_nbr+M_.param_nbr);
db(:,:,3) = aux1;%wrt theta
db(:,:,4) = theta*daux1_drho;%wrt rho
dc = zeros(size(c,1),size(c,2),M_.exo_nbr+M_.param_nbr);
dc(:,:,1) = theta^2*SE_e*aux3*aux4;%wrt SE_e
dc(:,:,3) = theta*SE_e^2*aux3*aux4;%wrt theta
dc(:,:,4) = 1/2*theta^2*SE_e^2*daux5_drho; %wrt rho

d2flag=0;
g_0 = 1/2*oo_.dr.ghs2; if ~isequal(g_0,oo_.dr.g_0); error('something wrong'); end
g_1 = [oo_.dr.ghx oo_.dr.ghu] +3/6*[oo_.dr.ghxss oo_.dr.ghuss]; if ~isequal(g_1,oo_.dr.g_1); error('something wrong'); end
g_2 = 1/2*[oo_.dr.ghxx oo_.dr.ghxu oo_.dr.ghuu]; if ~isequal(g_2,oo_.dr.g_2); error('something wrong'); end
g_3 = 1/6*[oo_.dr.ghxxx oo_.dr.ghxxu oo_.dr.ghxuu oo_.dr.ghuuu]; if ~isequal(g_3,oo_.dr.g_3); error('something wrong'); end

tols = [1e-4 1e-4 1e-12 1e-12];
KRONFLAGS = [-1 -2 0 1];
for k = 1:length(KRONFLAGS)
    fprintf('KRONFLAG=%d\n',KRONFLAGS(k));
    options_.analytic_derivation_mode = KRONFLAGS(k);
    DERIVS = get_perturbation_params_derivs(M_, options_, estim_params_, oo_, indpmodel, indpstderr, indpcorr, d2flag);
    oo_.dr.dg_0 = permute(1/2*DERIVS.dghs2,[1 3 2]);
    oo_.dr.dg_1 = cat(2,DERIVS.dghx,DERIVS.dghu) + 3/6*cat(2,DERIVS.dghxss,DERIVS.dghuss);
    oo_.dr.dg_2 = 1/2*cat(2,DERIVS.dghxx,DERIVS.dghxu,DERIVS.dghuu);
    oo_.dr.dg_3 = 1/6*[DERIVS.dghxxx DERIVS.dghxxu DERIVS.dghxuu DERIVS.dghuuu];
    
    for ord = 0:@{ORDER}
        g = oo_.dr.(['g_' num2str(ord)])(2,:); % Retrieve computed policy function for variable y
        dg = oo_.dr.(['dg_' num2str(ord)])(2,:,:);    
        for m = 0:ord % m is the derivation order with respect to x(-1)
            v = 0;
            dv = zeros(1,M_.exo_nbr + M_.param_nbr);
            for p = 0:floor((@{ORDER}-ord)/2) % 2p is the derivation order with respect to s
                if ord+2*p > 0 % Skip the deterministic steady state constant
                    v = v + sum(beta.^i.*exp(theta*xbar*i).*b.^ord.*rho^m.*c.^p)/factorial(ord)/factorial(p);
                    %derivatives
                    dv(:,1) = dv(:,1) + sum( beta.^i.*exp(theta*xbar*i).*ord.*b.^(ord-1).*db(:,:,1).*rho^m.*c.^p...
                                            +beta.^i.*exp(theta*xbar*i).*b.^ord.*rho^m.*p.*c.^(p-1).*dc(:,:,1)...
                                           )/factorial(ord)/factorial(p);%wrt SE_E
                    dv(:,2) = dv(:,2) + sum( i.*beta.^(i-1).*exp(theta*xbar*i).*b.^ord.*rho^m.*c.^p...
                                            +beta.^i.*exp(theta*xbar*i).*ord.*b.^(ord-1).*db(:,:,2).*rho^m.*c.^p...
                                            +beta.^i.*exp(theta*xbar*i).*b.^ord.*rho^m.*p.*c.^(p-1).*dc(:,:,2)...
                                           )/factorial(ord)/factorial(p);%wrt beta
                    dv(:,3) = dv(:,3) + sum( beta.^i.*exp(theta*xbar*i).*xbar.*i.*b.^ord.*rho^m.*c.^p...
                                            +beta.^i.*exp(theta*xbar*i).*ord.*b.^(ord-1).*db(:,:,3).*rho^m.*c.^p...
                                            +beta.^i.*exp(theta*xbar*i).*b.^ord.*rho^m.*p.*c.^(p-1).*dc(:,:,3)...
                                           )/factorial(ord)/factorial(p);%wrt theta
                    dv(:,4) = dv(:,4) + sum( beta.^i.*exp(theta*xbar*i).*b.^ord.*m.*rho^(m-1).*c.^p...
                                            +beta.^i.*exp(theta*xbar*i).*ord.*b.^(ord-1).*db(:,:,4).*rho^m.*c.^p...
                                            +beta.^i.*exp(theta*xbar*i).*b.^ord.*rho^m.*p.*c.^(p-1).*dc(:,:,4)...
                                           )/factorial(ord)/factorial(p);%wrt rho
                    dv(:,5) = dv(:,5) + sum( beta.^i.*exp(theta*xbar*i).*theta.*i.*b.^ord.*rho^m.*c.^p...
                                            +beta.^i.*exp(theta*xbar*i).*ord.*b.^(ord-1).*db(:,:,5).*rho^m.*c.^p...
                                            +beta.^i.*exp(theta*xbar*i).*b.^ord.*rho^m.*p.*c.^(p-1).*dc(:,:,5)...
                                           )/factorial(ord)/factorial(p);%wrt xbar
                end                
            end
            if abs(v-g(ord+1-m)) > 1e-14
                error(['Error in matrix oo_.dr.g_' num2str(ord)])
            end
            chk_dg = squeeze(dg(:,ord+1-m,:))';
            if isempty(indpstderr)
                chk_dv = dv(:,M_.exo_nbr+indpmodel);
            elseif isempty(indpmodel)
                chk_dv = dv(:,1:M_.exo_nbr);
            else
                chk_dv = dv;
            end
            fprintf('Max absolute deviation for dg_%d(2,%d,:): %e\n',ord,ord+1-m,norm( chk_dv - chk_dg, Inf));
            if norm( chk_dv - chk_dg, Inf) > tols(k)
                error(['Error in matrix dg_' num2str(ord)])
                chk_dv
                chk_dg
            end
        end
    end
    fprintf('\n');
end