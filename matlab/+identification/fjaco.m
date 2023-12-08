function fjac = fjaco(f,x,varargin)
% fjac = fjaco(f,x,varargin)
% FJACO Computes two-sided finite difference Jacobian
% USAGE
%   fjac = fjaco(f,x,P1,P2,...)
% INPUTS
%   f         : name of function of form fval = f(x)
%   x         : evaluation point
%   P1,P2,... : additional arguments for f (optional)
% OUTPUT
%   fjac      : finite difference Jacobian
%
% Copyright © 2010-2017,2019-2023 Dynare Team
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

ff=feval(f,x,varargin{:});

tol = eps.^(1/3); %some default value
if strcmp(func2str(f),'identification.get_perturbation_params_derivs_numerical_objective') || strcmp(func2str(f),'identification.numerical_objective')
    tol= varargin{4}.dynatol.x;
end
h = tol.*max(abs(x),1);
xh1=x+h; xh0=x-h;
h=xh1-xh0;
fjac = NaN(length(ff),length(x));
for j=1:length(x)
    xx = x;
    xx(j) = xh1(j); f1=feval(f,xx,varargin{:});
    if isempty(f1) && (strcmp(func2str(f),'identification.get_perturbation_params_derivs_numerical_objective') || strcmp(func2str(f),'identification.numerical_objective') )
        [~,info]=feval(f,xx,varargin{:});
        disp_info_error_identification_perturbation(info,j);
    end
    xx(j) = xh0(j); f0=feval(f,xx,varargin{:});
    if isempty(f0) && (strcmp(func2str(f),'identification.get_perturbation_params_derivs_numerical_objective') || strcmp(func2str(f),'identification.numerical_objective') )
        [~,info]=feval(f,xx,varargin{:});
        disp_info_error_identification_perturbation(info,j)
    end
    fjac(:,j) = (f1-f0)/h(j);
end

feval(f,x,varargin{:});

%Auxiliary functions
function disp_info_error_identification_perturbation(info,j)
    % there are errors in the solution algorithm
    probl_par = get_the_name(j,varargin{4}.TeX,varargin{3},varargin{2},varargin{4}.varobs);
    skipline()
    message = get_error_message(info,varargin{4});
    fprintf('Parameter error in numerical two-sided difference method:\n')
    fprintf('Cannot solve the model for %s (info = %d, %s)\n', probl_par, info(1), message);
    fprintf('Possible solutions:\n')
    fprintf('  -- check your mod file, calibration and steady state computations carefully\n');
    fprintf('  -- use analytic derivatives, i.e. set analytic_derivation_mode=0\n');
    fprintf('  -- use an estimated_params block without %s or change its value\n', probl_par);
    fprintf('  -- change numerical tolerance level in fjaco.m (you can tune ''options_.dynatol.x'' or change fjaco.m function directly)\n');
    error('fjaco.m: numerical two-sided difference method yields errors in solution algorithm');
end

end %main function end