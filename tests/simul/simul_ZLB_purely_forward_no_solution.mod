var c, h, pi, w, R, r_e, y, gdp, gdp_hat, k, u, g, c_hat, w_hat, y_hat, h_hat;
varexo d, z, eta; 

parameters alpha, beta, sigma, gamma, theta, ni, tau_w, phi_p, phi_y;

beta  = 0.997;
sigma = 2;
gamma = 25;
theta = 7.67;
tau_w = 0.2;
ni    = 0.28;
phi_p = 1.5;
phi_y = 0.125;
alpha = 0.064;

model;

// log-deviation of _ from its steady state value
gdp_hat=log(gdp)-log(steady_state(gdp));
c_hat=log(c)-log(steady_state(c));
w_hat=log(w)-log(steady_state(w));
y_hat=log(y)-log(steady_state(y));
h_hat=log(h)-log(steady_state(h));

// real interest rate
r_e=1/(beta*d(+1))-1;

//FOC labor
c^sigma*h^ni=steady_state(w)*(1-tau_w);

//Euler equation 1
1=beta*d(+1)*(1+R)/(1+pi(+1))*(c/c(+1))^sigma;

//Euler equation 2
0=1/(1-alpha)*(steady_state(w)/z)*h^alpha-1-gamma/theta*pi*(1+pi)+beta*d(+1)*(c/c(+1))^sigma * y(+1)/y*gamma/theta*pi(+1)*(1+pi(+1));

// Taylor rule with ZLB
R=max(0,r_e+phi_p*pi+phi_y*gdp_hat);

//output
y=z*h^(1-alpha);

//aggregate resource constraint
c=(1-k-eta)*y;

// resource cost of price adjustment
k=(gamma/2)*(pi^2);

//government purchases
g=eta*y;

// GDP
gdp=(1-k)*y;

//utility
u=(c^(1-sigma))/(1-sigma)-(h^(1+ni))/(1+ni);
end;

initval;
z=1;
d=1;
pi=0;
k=(gamma/2)*(pi^2);
r_e=1/(beta*d)-1;
eta=0.2;
h=1;
y=z*h;
g=eta*y;
c=(1-k-eta)*y;
w=z;
gdp=(1-k)*y;
R=r_e;
end;

steady;
check;

shocks;
//5% preference shock
var d;
periods 1:10;
values 1.05;

//technology shock
var z;
periods 1:10;
values 1.05;
end;

simul(periods=40,maxit=1000);

rplot gdp_hat;
rplot R;

if oo_.deterministic_simulation.status==1
    error('This model has no solution');
end