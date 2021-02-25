var U k z c W;

varexo e;

parameters beta gamma alpha delta rho s;

beta = 0.987;
gamma = 1;
delta = 0.012;
alpha = 0.4;
rho = 0.95;
s = 0.007;

model;
c^(-gamma)=beta*c(+1)^(-gamma)*(alpha*exp(z(+1))*k^(alpha-1)+1-delta);
W = U + beta*W(+1);
k=exp(z)*k(-1)^(alpha)-c+(1-delta)*k(-1);
z=rho*z(-1)+s*e;
U=ln(c);
end;

steady_state_model;
k = ((1/beta-(1-delta))/alpha)^(1/(alpha-1));
c = k^alpha-delta*k;
z = 0;
U = ln(c);
W = U/(1-beta);
end;
