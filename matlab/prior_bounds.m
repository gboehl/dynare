function bounds = prior_bounds(bayestopt_, priortrunc)

% computes bounds for prior density.
%
% INPUTS
% - bayestopt   [struct]  characterizing priors (shape, mean, p1..p4)
% - priortrunc  [double]  scalar, probability mass in the tails to be removed
%
% OUTPUTS
% - bounds     [struct]  prior bounds (lb, lower bounds, and ub, upper bounds, fields are n×1 vectors)

% Copyright © 2003-2023 Dynare Team
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

if nargin<2, priortrunc = 0.0; end

assert(priortrunc>=0 && priortrunc<=1, 'Second input argument must be non negative and not larger than one.')

pshape = bayestopt_.pshape;
p3 = bayestopt_.p3;
p4 = bayestopt_.p4;
p6 = bayestopt_.p6;
p7 = bayestopt_.p7;

bounds.lb = zeros(size(p6));
bounds.ub = zeros(size(p6));

for i=1:length(p6)
    switch pshape(i)
      case 1
        if priortrunc==0
            bounds.lb(i) = p3(i);
            bounds.ub(i) = p4(i);
        else
            bounds.lb(i) = betainv(priortrunc, p6(i), p7(i))*(p4(i)-p3(i))+p3(i);
            bounds.ub(i) = betainv(1.0-priortrunc, p6(i), p7(i))*(p4(i)-p3(i))+p3(i);
        end
      case 2
        if priortrunc==0
            bounds.lb(i) = p3(i);
            bounds.ub(i) = Inf;
        else
            bounds.lb(i) = gaminv(priortrunc, p6(i), p7(i))+p3(i);
            bounds.ub(i) = gaminv(1.0-priortrunc, p6(i), p7(i))+p3(i);
        end
      case 3
        if priortrunc == 0
            bounds.lb(i) = max(-Inf, p3(i));
            bounds.ub(i) = min(Inf, p4(i));
        else
            bounds.lb(i) = max(norminv(priortrunc, p6(i), p7(i)), p3(i));
            bounds.ub(i) = min(norminv(1-priortrunc, p6(i), p7(i)), p4(i));
        end
      case 4
        if priortrunc==0
            bounds.lb(i) = p3(i);
            bounds.ub(i) = Inf;
        else
            bounds.lb(i) = 1.0/sqrt(gaminv(1.0-priortrunc, p7(i)/2.0, 2.0/p6(i)))+p3(i);
            bounds.ub(i) = 1.0/sqrt(gaminv(priortrunc, p7(i)/2.0, 2.0/p6(i)))+p3(i);
        end
      case 5
        if priortrunc == 0
            bounds.lb(i) = p6(i);
            bounds.ub(i) = p7(i);
        else
            bounds.lb(i) = p6(i)+(p7(i)-p6(i))*priortrunc;
            bounds.ub(i) = p7(i)-(p7(i)-p6(i))*priortrunc;
        end
      case 6
        if priortrunc == 0
            bounds.lb(i) = p3(i);
            bounds.ub(i) = Inf;
        else
            bounds.lb(i) = 1.0/gaminv(1.0-priortrunc, p7(i)/2.0, 2.0/p6(i))+p3(i);
            bounds.ub(i) = 1.0/gaminv(priortrunc, p7(i)/2.0, 2.0/p6(i))+ p3(i);
        end
      case 8
        if priortrunc == 0
            bounds.lb(i) = p3(i);
            bounds.ub(i) = Inf;
        else
            bounds.lb(i) = p3(i)+wblinv(priortrunc, p6(i), p7(i));
            bounds.ub(i) = p3(i)+wblinv(1.0-priortrunc, p6(i), p7(i));
        end
      otherwise
        error('prior_bounds: unknown distribution shape (index %d, type %d)', i, pshape(i));
    end
end
