@#define countries = 1:2
@#define assets = 1:2
var
@#for c in countries
  yk_@{c} yl_@{c} c_@{c}
@#endfor
a ze_1 r_1
;

varexo
@#for c in countries
  ek_@{c} el_@{c}
@#endfor
;

parameters rho eta omega sig_k sig_l ykbar ylbar;

rho = 1;
psi = 0.7;
eta = 0.9;
g = 0.25;
ykbar = 1;
ylbar = 1;
cbar = ykbar+ylbar;
nu = 0.5;
omega = 0.75;
sig_k = 0.02;
sig_l = 0.01;
sig_m = 0;
betabar = omega*cbar^(-eta);


model;
@#for c in countries
  yk_@{c} = log(ykbar) + sig_k*ek_@{c};
  yl_@{c} = log(ylbar) + sig_l*el_@{c};
  exp(c_@{c})^(eta-rho) = omega*exp(c_@{c}(+1))^(-rho)*r_1(+1);
@#endfor

a = a(-1)*r_1 + exp(yk_1) + exp(yl_1) - exp(c_1);
-a = -a(-1)*r_1 + exp(yk_2) + exp(yl_2) - exp(c_2);
//  r_1+ze_1  = ze_1+exp(yk_1-ze_1(-1));
  r_1  = exp(yk_1-ze_1(-1))+ze_1-ze_1;


end;

shocks;
@#for c in countries
  var ek_@{c} = 1;
  var el_@{c} = 1;
@#endfor
end;

initval;
@#for c in countries
  yk_@{c} = log(ykbar);
  yl_@{c} = log(ylbar);
  c_@{c} = log(ykbar+ylbar);
@#endfor
  r_1 = 1/betabar;
  ze_1 = log(betabar)+yk_1;
end;

resid;
steady;
model_diagnostics;
check;
stoch_simul(order=2,irf=0);

o1 = load(['ds1' filesep 'Output' filesep 'ds1_results'],'oo_');
oo1 = o1.oo_.dr;
oo2 = oo_.dr;
if any(abs(oo1.ghxx-oo2.ghxx) > 1e-14);error('ds1 with missing variable doesn''t reproduce ds2');end;
if any(abs(oo1.ghuu-oo2.ghuu) > 1e-14);error('ds1 with missing variable doesn''t reproduce ds2');end;
if any(abs(oo1.ghxu-oo2.ghxu) > 1e-14);error('ds1 with missing variable doesn''t reproduce ds2');end;
if any(abs(oo1.ghs2-oo2.ghs2) > 1e-14);error('ds1 with missing variable doesn''t reproduce ds2');end;
