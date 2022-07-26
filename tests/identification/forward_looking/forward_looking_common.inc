% Forward-looking example model from Koop, Pesaran, Smith (2013, JBES)
% created by Willi Mutschler (@wmutschl, willi@mutschler.eu)
% =========================================================================
% Copyright © 2020-2022 Dynare Team
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
var
  r  ${r}$    (long_name='interest rate')
  x  ${x}$    (long_name='output-gap')
  p  ${\pi}$  (long_name='inflation')
;
varexo
  e_M  ${\varepsilon^M}$  (long_name='monetary policy shock')
  e_D  ${\varepsilon^D}$  (long_name='demand shock')
  e_S  ${\varepsilon^S}$  (long_name='supply shock')
;
parameters
  PSI    ${\psi}$    (long_name='inflation elasticity Taylor rule')
  TAU    ${\tau}$    (long_name='intertemporal elasticity of subsitution')
  BETA   ${\beta}$   (long_name='discount factor')
  KAPPA  ${\kappa}$  (long_name='slope Phillips curve')
;

PSI = 1.1;
TAU = 2;
BETA = 0.9;
KAPPA = 0.6;

model;
[name='Taylor rule']
r = PSI*p + e_M;
[name='New Keynesian IS curve']
x = x(+1) - 1/TAU*(r-p(+1)) + e_D;
[name='New Keynesian Phillips curve']
p = BETA*p(+1) + KAPPA*x + e_S;
end;

shocks;
  var e_M = 1;
  var e_D = 1;
  var e_S = 1;
end;
