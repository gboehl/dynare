var labobs robs pinfobs dy dc dinve dw eta_w_ma eta_p_ma zcapf rkf kf pkf cf invef yf labf wf rrf mc zcap rk k_s q c i y l pinf w r eps_a eps_b eps_g eps_i eps_r eps_p eps_w kpf k;

varexo eta_a eta_b eta_g  eta_i eta_r eta_p eta_w;

parameters curv_w rho_ga curv_p l_bar pi_bar beta_const Mu_w Mu_p alpha psi phi delta sigma_c lambda phi_p iota_w xi_w iota_p xi_p sigma_l phi_w r_pi r_dy r_y rho
	rho_a rho_b rho_g rho_i rho_r rho_p rho_w gamma_bar G;

delta = 0.025;
phi_w = 1.5;
G = 0.18;
curv_p = 10;
curv_w = 10;
phi = 6.3325;
sigma_c = 1.2312;
lambda = 0.7205;
xi_w = 0.7937;
sigma_l = 2.8401;
xi_p = 0.7813;
iota_w = 0.4425;
iota_p = 0.3291;
psi = 0.2648;
phi_p = 1.4672;
r_pi = 1.7985;
rho = 0.8258;
r_y = 0.0893;
r_dy = 0.2239;
pi_bar = 0.7;
beta_const = 0.7420;
l_bar = 1.2918;
gamma_bar = 0.3982;
alpha = 0.24;
rho_a = .9676;
rho_b = .2703;
rho_g = .9930;
rho_i = .5724;
rho_r = .3;
rho_p = .8692;
rho_w = .9546;
Mu_p = .7652;
Mu_w = .8936;
rho_ga = 0.05;

model;
	# PI_star = 1 + pi_bar/100;
	# gamma = 1 + gamma_bar/100 ;
	# beta = 1/(1 + beta_const/100);
	# beta_bar = beta*gamma^(-sigma_c);
	# Rk = (beta^(-1)) * (gamma^sigma_c) - (1-delta);
	# W = (alpha^alpha*(1-alpha)^(1-alpha)/(phi_p*Rk^alpha))^(1/(1-alpha));
	# I_K_bar = (1-(1-delta)/gamma);
	# I_K = (1-(1-delta)/gamma)*gamma;
	# L_K = ((1-alpha)/alpha)*(Rk/W);
	# K_Y = phi_p*(L_K)^(alpha-1);
	# I_Y = I_K * K_Y;
	# C_Y = 1 - G - I_K*K_Y;
	# Z_Y = Rk*K_Y;
	# WL_C = (1/phi_w)*(1-alpha)/alpha*Rk*K_Y/C_Y;
	# r_bar=((PI_star/(beta*gamma^(-sigma_c)))-1)*100;
	eps_a = alpha * rkf + (1-alpha)*wf;
	zcapf = (1/(psi/(1-psi))) * rkf;
	rkf = (wf)+labf-kf;
	kf = kpf(-1) + zcapf;
	invef = (1/(1+beta_bar*gamma))*(invef(-1) + beta_bar*gamma*invef(1)+(1/(gamma^2*phi))*pkf) + eps_i;
	pkf = ((1-delta)/(Rk+(1-delta))) * pkf(+1) + (Rk/(Rk+(1-delta))) * rkf(+1) + (-rrf) + (1/((1-lambda/gamma)/(sigma_c*(1+lambda/gamma)))) * eps_b  ;
	cf = (lambda/gamma)/(1+lambda/gamma)*cf(-1) + (1/(1+lambda/gamma))*cf(+1) +((sigma_c-1)*WL_C/(sigma_c*(1+lambda/gamma)))*(labf-labf(+1)) - (1-lambda/gamma)/(sigma_c*(1+lambda/gamma))*(rrf+0*eps_b) + eps_b ;
	yf = C_Y*cf+I_Y*invef+eps_g  +  Z_Y*zcapf;
	yf = phi_p*( alpha*kf+(1-alpha)*labf +eps_a );
	wf = sigma_l*labf 	+(1/(1-lambda/gamma))*cf - (lambda/gamma)/(1-lambda/gamma)*cf(-1) ;
	kpf =  (1-I_K_bar)*kpf(-1)+(I_K_bar)*invef + (I_K_bar)*(gamma^2*phi)*eps_i ;
	mc = alpha*rk + (1-alpha)*w - eps_a;
	zcap =  ((1 - psi)/psi) * rk;
	rk =  w + l - k_s;
	k_s =  k(-1) + zcap;
	i = (1/(1 + beta_bar*gamma)) * (i(-1) + (beta_bar * gamma) * i(1) + (1/(gamma^2*phi)) * q) + eps_i;
	q = ((1-delta)/(Rk+(1-delta)))*q(1) + (Rk/(Rk+(1-delta))) * rk(1) - r + pinf(+1) + (1/((1-lambda/gamma)/(sigma_c*(1+lambda/gamma)))) * eps_b ;
	c = (lambda/gamma)/(1+lambda/gamma) * c(-1) + (1/(1+lambda/gamma)) * c(+1) + ((sigma_c-1)*WL_C/(sigma_c*(1+lambda/gamma))) * (l - l(+1)) - (1-lambda/gamma)/(sigma_c*(1+lambda/gamma)) * (r - pinf(+1)) + eps_b;
	y = C_Y * c + I_Y * i + eps_g + Z_Y * zcap;
	y = phi_p * (alpha * k_s + (1-alpha) * l + eps_a);
	pinf = (1/(1+beta_bar*gamma*iota_p)) * (beta_bar*gamma*pinf(+1) + iota_p * pinf(-1) + ((1-xi_p)*(1-beta_bar*gamma*xi_p)/xi_p)/((phi_p-1)*curv_p+1) * mc) + eps_p ;
	w =  (1/(1+beta_bar*gamma))*w(-1)
	   +(beta_bar*gamma/(1+beta_bar*gamma))*w(1)
	   +(iota_w/(1+beta_bar*gamma))*pinf(-1)
	   -(1+beta_bar*gamma*iota_w)/(1+beta_bar*gamma)*pinf
	   +(beta_bar*gamma)/(1+beta_bar*gamma)*pinf(1)
	   +(1-xi_w)*(1-beta_bar*gamma*xi_w)/((1+beta_bar*gamma)*xi_w)*(1/((phi_w-1)*curv_w+1))*
	   (sigma_l*l + (1/(1-lambda/gamma))*c - ((lambda/gamma)/(1-lambda/gamma))*c(-1) -w)
	   + 1*eps_w ;
	r =  r_pi * (1-rho) * pinf + r_y * (1-rho) * (y-yf) + r_dy * ( y - yf - (y(-1) - yf(-1))) + rho * r(-1) + eps_r;
	eps_a = rho_a * eps_a(-1) + eta_a;
	eps_b = rho_b * eps_b(-1) + eta_b;
	eps_g = rho_g * eps_g(-1) + eta_g + rho_ga * eta_a;
	eps_i = rho_i * eps_i(-1) + eta_i;
	eps_r = rho_r * eps_r(-1) + eta_r;
	eps_p = rho_p * eps_p(-1) + eta_p_ma - Mu_p * eta_p_ma(-1);
	eta_p_ma = eta_p;
	eps_w = rho_w * eps_w(-1) + eta_w_ma - Mu_w * eta_w_ma(-1);
	eta_w_ma = eta_w;
	k = (1-I_K_bar) * k(-1) + I_K_bar * i + I_K_bar*gamma^2*phi*eps_i;
	dy = y - y(-1) + gamma_bar;
	dc = c - c(-1) + gamma_bar;
	dinve = i - i(-1) + gamma_bar;
	dw = w - w(-1) + gamma_bar;
	pinfobs = pinf + pi_bar;
	robs = r + r_bar;
	labobs = l + l_bar;
end;

shocks;
    var eta_a;
    periods 100;
    values -.5;
end;

steady;

check;

simul(periods=300);
endo_simul_0 = oo_.endo_simul;
if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

simul(linear_approximation,periods=300);
endo_simul_1 = oo_.endo_simul;
if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

if max(abs(endo_simul_0(:)-endo_simul_1(:)))>.01*options_.dynatol.f
    error('Something is wrong!')
end
