/*
 * This file is a modified version of 'fs2000.mod'.
 *
 * The difference is that, here, the equations are written in non-stationary form, and we test if 
 * we are able to identify the trends.
 *
 */

/*
 * Copyright (C) 2019 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

var gM M;
var gA A k c y;
var P;             // follows M(-1)/A
var W l d;         // follows M(-1)
var R n;

varexo e_a e_m;

parameters alp bet gam mst rho psi del;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

model;
A = gA*A(-1);
M = gM*M(-1);
gA = exp(gam+e_a);
log(gM) = (1-rho)*log(mst) + rho*log(gM(-1))+e_m;
c+k = k(-1)^alp*(A*n)^(1-alp)+(1-del)*k(-1);
P*c = M;
P/(c(+1)*P(+1))=bet*P(+1)*(alp*k^(alp-1)*(A(+1)*n(+1))^(1-alp)+(1-del))/(c(+2)*P(+2));
(psi/(1-psi))*(c*P/(1-n))=W;
R = P*(1-alp)*k(-1)^alp*A^(1-alp)*n^(-alp)/W;
W = l/n;
M-M(-1)+d = l;
1/(c*P)=bet*R/(c(+1)*P(+1));
y = k(-1)^alp*(A*n)^(1-alp);
end;

verbatim;

    bgp.write(M_);
    options = optimoptions('fsolve','Display','iter','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',1000000,'MaxIterations',100000,'SpecifyObjectiveGradient',true,'FunctionTolerance',1e-6,'StepTolerance',1e-6);
    y = 1+(rand(M_.endo_nbr,1)-.5)*.5;
    g = ones(M_.endo_nbr,1);% 1+(rand(M_.endo_nbr,1)-.5)*.1;
    [y, fval, exitflag] = fsolve(@fs2000.bgpfun, [y;g], options);
    y(1:M_.orig_endo_nbr)
    y(M_.endo_nbr+(1:M_.orig_endo_nbr))

end;