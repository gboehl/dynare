var c q m b u R uaux ca;

varexo eeps nnu;

parameters bbetaL ies crra gbar bbar eeta Rbar Rpers;

bbetaL = 0.98;
ies = 1.5;
crra = 12.0;
gbar = .01;
bbar = -.04;
eeta = 2.0;
Rbar = 1.03;
Rpers = .5;

model;

    q = m(+1)*R(+1);
    uaux = u(+1)^(1-crra);
    m = bbetaL*exp(-gbar/ies)*uaux(-1)^((ies-1)/(ies*(1-crra))-1)*u^(1/ies-crra)*c^(-1/ies)/c(-1)^(-1/ies);
    b = bbar*q^eeta*exp(nnu);
    c + q*b*exp(gbar) = 1 + R*b(-1);
    u = ( (1-bbetaL)*c^((ies-1)/ies) + bbetaL*exp(gbar*(ies-1)/ies)*uaux^((ies-1)/(ies*(1-crra))))^(ies/(ies-1));
    R = (1-Rpers)*Rbar + Rpers*R(-1) + eeps;
    ca = b*exp(gbar)-b(-1);
  
end;

steady_state_model;
R = Rbar;
m = bbetaL*exp(-gbar/ies);
q = bbetaL*exp(-gbar/ies)*Rbar;
b = bbar*(bbetaL*exp(-gbar/ies)*Rbar)^eeta;
c = 1 + Rbar*bbar*(bbetaL*exp(-gbar/ies)*Rbar)^eeta - exp(gbar)*bbar*(bbetaL*exp(-gbar/ies)*Rbar)^eeta*bbetaL*exp(-gbar/ies)*Rbar;
u = ((1-bbetaL)/(1-bbetaL*exp(gbar*(ies-1)/ies )) )^(ies/(ies-1))*(1 + Rbar*bbar*(bbetaL*exp(-gbar/ies)*Rbar)^eeta - exp(gbar)*bbar*(bbetaL*exp(-gbar/ies)*Rbar)^eeta*bbetaL*exp(-gbar/ies)*Rbar);
uaux = u^(1-crra);
ca = (exp(gbar)-1)*bbar*(bbetaL*exp(-gbar/ies)*Rbar)^eeta;
end;

resid;
steady(solve_algo=3);

estimated_params;
ies, normal_pdf, 1.5, .2;
end;

options_.plot_priors=0;
