// Declaration of some endogenous variables.
var c;

/*
** REMARKS
**
** Some declarations for the endogenous variables are missing. They are declared below in
** the model block.
**
*/

// Declaration of some of the parameters.
parameters beta, rho, theta, psi, tau;

rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
psi   = 0;
theta = 2.95;
phi   = 0.1;

/*
** REMARKS
**
** Same remark, some parameters are defined in the below in the model block (alpha and delta). These
** variables must be calibrated after the model block (ie after the calibration).
**
** We do not declare the exogenous variables in the preamble. The innovations are decalred in the model
** block below.
**
*/

model;
c*theta*h|e^(1+psi)=(1-alpha)*y;
k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
    *(exp(b(+1))*alpha|p*y(+1)+(1-delta)*k));
y|e = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k|e = exp(b)*(y-c)+(1-delta|p)*k(-1);
a|e = rho*a(-1) + e|x;
b|e = rho*b(-1) + u|x;
end;

/*
** REMARKS
**
** On the fly declaration works as follows:
**
**  - A symbol followed by |e implicitely declares an endogenous variable.
**  - A symbol followed by |x implicitely declares an exogenous variable.
**  - A symbol followed by |p implicitely declares a parameter.
**
** A parameter declared on the fly has to be calibrated after its declaration, ie after the model
** block. If the user tries to give a value to a parameter before its declaration, then Dynare
** interpret the calibration as a matlab statement, because Dynare doesn't know that the symbol is
** a parameter.
**
** By default an undeclared symbol is interpreted as an exogenous variable. Consequently, if the user
** removes an equation where an endogenous variable is declared, the status of the variable changes.
** For instance, if one comments the last equation, the law of motion of b, symbol b becomes an
** exogenous variable (whose value is zero by default). A more interesting use case, is to remove
** the first equation, the arbitrage condition between leisure and consumption. Labour supply is
** then supposed to be exogenous.
*/

alpha = 0.36;
delta = 0.025;

initval;
y = 1.08068253095672;
c = 0.80359242014163;
h = 0.29175631001732;
k = 11.08360443260358;
a = 0;
b = 0;
e = 0;
u = 0;
end;

shocks;
  var e; stderr 0.009;
  var u; stderr 0.009;
  var e, u = phi*0.009*0.009;
end;

stoch_simul(order=1);
