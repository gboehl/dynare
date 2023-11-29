@#include "ls2003_model.inc"

estimated_params;
psi1 , gamma_pdf,1.5,0.5;
psi2 , gamma_pdf,0.25,0.125;
psi3 , gamma_pdf,0.25,0.125;
rho_R ,beta_pdf,0.5,0.2;
alpha ,beta_pdf,0.3,0.1;
rr ,gamma_pdf,2.5,1;
k , gamma_pdf,0.5,0.25;
tau ,gamma_pdf,0.5,0.2;
rho_q ,beta_pdf,0.4,0.2;
rho_A ,beta_pdf,0.5,0.2;
rho_ys ,beta_pdf,0.8,0.1;
rho_pies,beta_pdf,0.7,0.15;
stderr e_R,inv_gamma_pdf,1.2533,0.6551;
stderr e_q,inv_gamma_pdf,2.5066,1.3103;
stderr e_A,inv_gamma_pdf,1.2533,0.6551;
stderr e_ys,inv_gamma_pdf,1.2533,0.6551;
stderr e_pies,inv_gamma_pdf,1.88,0.9827;
end;

  
disp('CREATE SCREENING SAMPLE, CHECK FOR STABILITY AND PERFORM SENSITIVITY ANALYSIS');
% disp('PRESS ENTER TO CONTUNUE');
% pause;

dynare_sensitivity(morris=1, morris_nliv=6, morris_ntra=20, redform=1, 
  namendo=(:), namexo=(:), namlagendo=(:));

stoch_simul(order=1,irf=40);
