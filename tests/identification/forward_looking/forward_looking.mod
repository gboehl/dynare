% Forward-looking example model from Koop, Pesaran, Smith (2013, JBES)
% created by Willi Mutschler (@wmutschl, willi@mutschler.eu)
% =========================================================================
% Copyright (C) 2020 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
% =========================================================================
var r x p;
varexo e_M e_D e_S;
varobs r x p;

parameters PSI TAU BETA KAPPA;

PSI=1.1;
TAU=2;
BETA=0.9;
KAPPA=0.6;

model;
r = PSI*p + e_M;
x = x(+1) - 1/TAU*(r-p(+1)) + e_D;
p = BETA*p(+1) + KAPPA*x + e_S;
end;

shocks;
var e_M = 1;
var e_D = 1;
var e_S = 1;
end;

steady;
check;

estimated_params;
PSI,   1.1;
TAU,   2;
BETA,  0.9;
KAPPA, 0.6;
end;

identification; %this triggers sylvester3a with empty ghx

% as a side note, we have the true solution:
% [r;x;p] = TRUE_SOLUTION*[e_M;e_D;e_S] (ghx is empty)
A = [1 0 -PSI; 1/TAU 1 0; 0 -KAPPA 1];
TRUE_SOLUTION1 = inv(A);
TRUE_SOLUTION2 = 1/(KAPPA*PSI/TAU +1)*[1         KAPPA*PSI PSI;
                                      -1/TAU     1         -PSI/TAU;
                                      -KAPPA/TAU KAPPA     1];
% note that BETA drops out from the solution

if max(max(abs(TRUE_SOLUTION1 - oo_.dr.ghu))) > 1e-15
    error('Something wrong with perturbation');
end
if max(max(abs(TRUE_SOLUTION2 - oo_.dr.ghu))) > 1e-15
    error('Something wrong with perturbation');
end