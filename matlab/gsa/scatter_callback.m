function  scatter_callback(K, type)
%
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu
%

% Copyright © 2017 European Commission
% Copyright © 2018 Dynare Team
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

global oo_ M_ options_ bayestopt_ estim_params_

x=get(gcf,'userdata');
r2=x{1};
x=x{2};

xparam1=x(K,:)';

switch type
  case 'save'
    save(['my_params_' int2str(K)],'xparam1')
  case 'eval'
    disp('Evaluating smoother ...')
    [oo_, M_]=evaluate_smoother(xparam1, M_.endo_names, M_, oo_, options_, bayestopt_, estim_params_);
end
