var pi_h y_hat i r_nat a y_star nu s y y_nat z z_star z_e dY_o dY_star_o dP_o I_o p_h e i_star de_o;

varexo eps_a eps_nu eps_star eps_z eps_z_star eps_z_e;

parameters alpha varphi beta theta gamma sigma eta rho_nu rho_a rho_star phi_pi phi_y sig_a sig_nu sig_star sig_z sig_z_e gY gY_star sig_z_star gP gi ge;

gY = 0.3783;
gY_star = 0.3783;
gP = 0.3845;
gi = 2.5197;
ge = 0;
beta = 0.99;
alpha = 0.53;
sigma = 0.7;
gamma = 1.6;
eta = 1.5;
varphi = 5;
theta = 3/4;
rho_nu = 0.5;
rho_a = 0.9;
rho_star = 0.5;
phi_pi = 1.5;
phi_y = 0.5/4;
sig_a = 0.0405;
sig_nu = 0.0031;
sig_star = 0.1;

sig_z = 0.0109;
sig_z_star = 0.0109;
sig_z_e = 0.0109;


model(linear);
    #omega = sigma*gamma + (1-alpha)*(sigma*eta-1);
    #lambda = (1-beta*theta)*(1-theta)/theta;
    #sigma_alpha = sigma/(1+alpha*(omega-1));
    #Theta = omega-1;
    #Gamma_a = (1+varphi)/(sigma_alpha+varphi);
    #Gamma_star = -alpha*Theta*sigma_alpha/(sigma_alpha+varphi);
    #kappa_alpha = lambda*(sigma_alpha*varphi);


    pi_h = beta*pi_h(+1) + kappa_alpha*y_hat; //NKPC
    y_hat = y_hat(+1) - 1/sigma_alpha*(i-pi_h(+1)-r_nat);   //DIS
    r_nat = -sigma_alpha*Gamma_a*(1-rho_a)*a+alpha*Theta*sigma_alpha*varphi/(sigma_alpha+varphi)*(y_star(+1)-y_star);   //Natural interest rate
    a = rho_a*a(-1)+eps_a;  //Technology
    i = phi_pi*pi_h + phi_y*y_hat + nu;  //Taylor rule
    nu = rho_nu*nu(-1) + eps_nu;    //Monetary shock
    s = sigma_alpha*(y-y_star);     //Equation 29
    y_hat = y - y_nat;  //Output gap
    y_nat = Gamma_a*a + Gamma_star*y_star;  //Natural output
    y_star = rho_star*y_star(-1)+eps_star;  //Foreign GDP growth
    e = s + p_h; //PPP
    //e-e(-1) = s - s(-1) + pi_h; //PPP Modified
    pi_h = p_h - 0.9999999*p_h(-1); //Definition of inflation. If coef = 1, problem of singularity and estimation doesn't work
    i = i_star + e(+1) - e; //UIP

    //Observables
    dY_o = gY + 100*(y - y(-1) + z) ;
    z = eps_z;
    dY_star_o = gY_star + 100*(y_star - y_star(-1) + z_star) ;
    z_star = eps_z_star;
    dP_o = gP + 100*pi_h ;
    I_o  = gi + 400*i ;
    de_o = ge + 100*(e - e(-1) + z_e);
    z_e = eps_z_e;
    
end;


shocks;
    var eps_a=sig_a^2;
    var eps_nu=sig_nu^2;
    var eps_star=sig_star^2;
    var eps_z = sig_z^2;
    var eps_z_star = sig_z_star^2;
    var eps_z_e = sig_z_e^2;
end;

steady;

varobs dY_o dY_star_o dP_o I_o de_o;

estimated_params;

    gY,  normal_pdf, 0.3783, 0.10 ;
    gY_star,  normal_pdf, 0.3783, 0.10 ;
    gP, normal_pdf, 0.3845, 0.10;
    gi,  normal_pdf, 2.5197, 0.50;
    ge, normal_pdf, 0, 0.1;
    sigma, normal_pdf, 0.7, 0.1 ;
    gamma, normal_pdf, 1.6, 0.3 ;
    eta, normal_pdf, 1.5, 0.3;
    varphi, normal_pdf, 2.05, 0.6;
    theta, beta_pdf, 3/4, 0.05;
    rho_nu, beta_pdf, 0.5, 0.1;
    rho_a, beta_pdf, 0.7, 0.1;
    rho_star, beta_pdf, 0.5, 0.1;
    phi_pi, normal_pdf, 1.5, 0.1;
    phi_y, normal_pdf, 0.125, 0.05;
    stderr eps_a, uniform_pdf, , , 0, 0.1;
    stderr eps_z, uniform_pdf, , , 0, 0.2;
    stderr eps_nu, uniform_pdf, , , 0, 0.1;
    stderr eps_star, uniform_pdf, , , 0, 0.1;
    stderr eps_z_star, uniform_pdf, , , 0, 0.1;
    stderr eps_z_e, uniform_pdf, , , 0, 0.1;

end;

options_.diffuse_filter=1;
stoch_simul(nograph);
calib_smoother(datafile=SOE_data_file,smoothed_state_uncertainty) dY_o dY_star_o dP_o I_o de_o;