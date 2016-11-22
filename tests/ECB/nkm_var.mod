var i, pi, y;
varexo ee, nu;
parameters bet, gam, sig, phi, psi;


bet = 0.99;
gam = 0.1275;
sig = 1;
phi = 1.5;
psi = 0.125;

var_model(model_name=my_var_est, order=1) pi y;

model;

% Taylor Rule
  i = phi*pi + psi*y + ee;

% Phillips Curve
  pi = bet*pi(1) + gam*y;

% IS
  y = var_expectation(y, model_name=my_var_est) + y(-1) - 1/sig*(i - var_expectation(pi, model_name=my_var_est)) + nu;

end;

shocks;
var ee; stderr 0.0205;
var nu; stderr 0.0305;
end;

stoch_simul(order=1, periods=200);

