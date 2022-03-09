function [tlogpostkern,loglik] = tempered_likelihood(TargetFun,xparam1,lambda,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_)
% function [tlogpostkern,loglik] = tempered_likelihood(TargetFun,xparam1,lambda,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_)

% Copyright © 2022 Dynare Team
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

    logpostkern = -feval(TargetFun,xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_);
    logprior = priordens(xparam1,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
    loglik = logpostkern-logprior ;
    tlogpostkern = lambda*loglik + logprior;
