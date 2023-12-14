var c, y_ro, psi, s, tt, pi_h ${\pi_h}$, mc_ro, i_ro, pi_ro ${\pi_\rho}$, pi_eu, pi_f, y_eu, mc_eu, i_eu, g, a, cp, bnr, gs, as ;
varexo eps_a, eps_g, eps_cp, eps_s, eps_m $\varepsilon_m$, eps_as, eps_gs, eps_ms;
parameters alfa, niu, delt, teta_h, teta_f, bet, fi, sigm, h, ro, psi_pi $\psi_\pi$, tetas, rhos, psi_pis, psi_y, psi_deltay, rho_g, rho_a, rho_cp, rho_s, rho_gs, rho_as ; 

alfa = 0.425 ;
bet = 0.99 ; 

niu = 0.6 ; 
delt = 0.3 ;
teta_h = 0.7 ;
teta_f = 0.7 ;
fi = 2 ; 
sigm = 0.5 ;
h = 0.8 ;
ro = 0.5 ;
psi_pi = 1.5 ;
tetas = 0.7 ;
rhos = 0.5 ;
psi_pis = 1.5 ;
psi_y = 0.25 ;
psi_deltay = 0.25 ;
rho_g = 0.5 ;
rho_a = 0.5 ;
rho_cp = 0.5 ;
rho_s = 0.5 ;
rho_gs = 0.5 ;
rho_as = 0.5 ;

model (linear) ; 

// (1) mk clearing condition
(1-alfa) * c = y_ro - alfa * niu* (2-alfa)* tt - alfa * niu * psi - alfa * y_eu ;
// (2) terms of trade
tt - tt(-1) = pi_f - pi_h ;
// (3) terms of trade, real exchange rate, law of one price gap
s = psi + (1 - alfa) * tt;
// (4) domestic producers price setting
pi_h  - delt * pi_h(-1) = (1 - teta_h) * (1- bet * teta_h) / teta_h * mc_ro + bet * ( pi_h(+1)  - delt * pi_h) ;
// (5) domestic real mg costs
mc_ro = fi * y_ro - (1 + fi) * a + alfa * tt + sigm / (1-h) * (c - h* c(-1)) ;
// (6) importers price setting 
pi_f  - delt * pi_f(-1) = (1 - teta_f) * (1 - bet * teta_f) / teta_f * psi + bet * (pi_f (+1)  - delt * pi_f) + cp;
// (7) complete mk assumption
c - h * c(-1) = y_eu - h * y_eu(-1) + (1 - h) / sigm * (psi + (1 - alfa) * tt) + g ;
// (8) domestic inflation
pi_ro = pi_h + alfa * (tt - tt(-1)) ;
// (9) UIP
i_ro - pi_ro(+1) - i_eu + pi_eu(+1) = s(+1) - s + bnr ;
// (10) domestic monetary policy rule
i_ro = ro * i_ro(-1) + (1 - ro) * (psi_pi * pi_ro + psi_y * y_ro + psi_deltay * (y_ro - y_ro(-1)) ) + eps_m ;
// (11) foreign euler equation
y_eu - h * y_eu(-1) = y_eu(+1) - h * y_eu - (1 - h) / sigm * (i_eu - pi_eu(+1)) + gs - gs(+1) ;
// (12) foreign producers price setting
pi_eu = (1 - tetas) * (1 - bet * tetas) / tetas * mc_eu + bet * pi_eu(+1) ;  
// (13) foreign real mg costs
mc_eu = fi * y_eu - (1 + fi) * as + alfa * tt + sigm / (1-h) * (y_eu - h * y_eu(-1)) ;
// (14) foreign monetary policy rule
i_eu = rhos * i_eu(-1) + (1 - rhos) * (psi_pis * pi_eu + psi_y * y_eu) + eps_ms ;

// exogenous processes
// (15) domestic shock in preferences (demand sh)
g = rho_g * g(-1) + eps_g ;
// (16) domestic technological shock (supply sh)
a = rho_a * a(-1) + eps_a ;
// (17) exchange rate shock 
bnr = rho_s * bnr(-1) + eps_s ;
// (18) foreign shock in preferences
gs = rho_gs * gs(-1) + eps_gs ;
// (19) foreign technological shock
as = rho_as * as(-1) + eps_as ;
// (20) cost-push shock
cp = rho_cp * cp(-1) + eps_cp ;

end ;

shocks;
var eps_g; stderr 1 ;
var eps_a; stderr 1 ;
var eps_cp; stderr 1 ;
var eps_s; stderr 1 ;
var eps_m ; stderr 1 ;
var eps_gs; stderr 1 ;
var eps_as; stderr 1 ;
var eps_ms ; stderr 1 ; 
end ; 

/* steady ; 
*/

check ;

estimated_params ;
niu, 0.6, 0.01, 0.999 ; 
delt, 0.3, 0.01, 0.999 ;
teta_h, 0.7, 0.01, 0.999 ;
teta_f, 0.7, 0.01, 0.999 ;
fi, 2, 0.01, 5 ; 
sigm, 0.5, 0.01, 0.999 ;
h, 0.8, 0.01, 0.999 ;
ro, 0.5, 0.01, 0.999 ;
psi_pi, 1.4, 0.01, 8 ;
tetas, 0.6, 0.01, 0.999 ;
rhos, 0.5, 0.01, 0.999 ;
psi_pis, 1.4, 0.01, 8 ;
psi_y, 0.25, 0.01, 0.999 ;
psi_deltay, 0.25, 0.01, 0.999 ;
rho_g, 0.3, 0.01, 0.999 ;
rho_a, 0.3, 0.01, 0.999 ;
rho_cp, 0.3, 0.01, 0.999 ;
rho_s, 0.3, 0.01, 0.999 ;
rho_gs, 0.3, 0.01, 0.999 ;
rho_as, 0.3, 0.01, 0.999 ;

stderr eps_g, 1, 0.01, 10 ;
stderr eps_a, 1, 0.01, 10 ;
stderr eps_cp, 1, 0.01, 10 ;
stderr eps_s, 1, 0.01, 10 ;
stderr eps_gs, 1, 0.01, 10 ;
stderr eps_as, 1, 0.01, 10 ;
% corr eps_g, eps_a, 0, -1, 1;
end;

varobs y_ro, pi_ro, i_ro, s, y_eu, pi_eu, i_eu, tt ;
options_.TeX    =1;
sensitivity (identification=1, nsam = 500, lik_only = 1, morris=2) ;

stoch_simul(order=2,irf=20)  y_ro pi_ro i_ro s ;
//order 1 - impulse response functions are simply the algebraic forward iteration of the model's policy 
//order 2 - impulse response functions will be the result of actual Monte Carlo simulations of future shocks  
//specify sucient periods in IRFs to see the graphs return to steady state 
 
check ;

/* shock_decomposition y_ro ;
*/
