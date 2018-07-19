// --+ options: nostrict +--
/*
** In the following example, the types (exogenous variable, endogenous variable, parameter) of the objects appearing in the model are declared
** on the fly in the equations. With the nostrict option, an object is by default an exogenous variable (if the type of an object is not declared
** in one of the equations, it will be interpreted as an exogenous variable). Without the nostrict option, Dynare will raise an error if the model
** contains untyped objects.
**
** An object followed by |e is an endogenous variable,
**                       |x is an exogenous variable,
**                       |p is a parameter.
**
** Example. If the first equation (consumption/leisure arbitrage) is removed from the following model block, then h (hours)
** will be interpreted as an exogenous variable.
*/

model;
c*theta|p*h|e^(1+psi|p)=(1-alpha)*y;
k|e = beta|p*(((exp(b)*c)/(exp(b(+1))*c(+1)))
    *(exp(b(+1))*alpha|p*y(+1)+(1-delta)*k));
y|e = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c|e)+(1-delta|p)*k(-1);
a|e = rho|p*a(-1)+tau*b(-1) + e|x;
b|e = tau|p*a(-1)+rho*b(-1) + u|x;
end;

alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;
phi   = 0.1;

if ~isequal(length(intersect(M_.endo_names, {'c'; 'h'; 'y'; 'k'; 'b'; 'a'})), 6)
   error('Endogenous variables are wrong.')
end

if ~isequal(length(intersect(M_.param_names, {'theta'; 'psi'; 'alpha'; 'beta'; 'delta'; 'rho'; 'tau'})), 7)
   error('Parameters are wrong.')
end

if ~isequal(length(intersect(M_.exo_names, {'e'; 'u'})), 2)
   error('Exogenous variables are wrong.')
end