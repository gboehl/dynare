// --+ options: nostrict +--
/*
** Same as ex1.mod without the first equation -> h is an exogenous variable, psi is not a model object (because it
** is not used in the remaining equations), and theta is an exogenous equation (see the last equation).
*/

model;
//c*theta|p*h|e^(1+psi|p)=(1-alpha)*y;
k|e = beta|p*(((exp(b)*c)/(exp(b(+1))*c(+1)))
    *(exp(b(+1))*alpha|p*y(+1)+(1-delta)*k));
y|e = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c|e)+(1-delta|p)*k(-1);
a|e = rho|p*a(-1)+tau*b(-1) + e|x;
b|e = tau|p*a(-1)+rho*b(-1) + u|x + .0*theta; // The last term is here just for testing purpose.
end;

alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
// If the following line is uncommented, this will trigger an error from the preprocessor
// because it is not allowed to give a value to an exogenous variable.
//theta = 2.95;
phi   = 0.1;

if ~isequal(length(intersect(M_.endo_names, {'c'; 'y'; 'k'; 'b'; 'a'})), 5)
   error('Endogenous variables are wrong.')
end

if ~isequal(length(intersect(M_.param_names, {'alpha'; 'beta'; 'delta'; 'rho'; 'tau'})), 5)
   error('Parameters are wrong.')
end

if ~isequal(length(intersect(M_.exo_names, {'e'; 'u'; 'h'; 'theta'})), 4)
   error('Exogenous variables are wrong.')
end