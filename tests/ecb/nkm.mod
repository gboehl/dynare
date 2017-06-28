var i, pi, y;
varexo ee, nu;
parameters bet, bet1, gam, sig, phi, psi;


bet = 0.99;
bet1 = 0.99;
gam = 0.1275;
sig = 1;
phi = 1.5;
psi = 0.125;


model;

% Taylor Rule
  i = phi*pi + psi*y + ee;

% Phillips Curve
  pi = bet*pi(1) + gam*y;

% IS
  y = y(1) + y(-1) - 1/sig*(i - pi(1)) + nu;

end;

shocks;
var ee; stderr 0.0205;
var nu; stderr 0.898;
end;

stoch_simul(order=1, periods=200);

