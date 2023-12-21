function bayestopt_ = transform_prior_to_laplace_prior(bayestopt_)
% bayestopt_ = transform_prior_to_laplace_prior(bayestopt_)
% -------------------------------------------------------------------------
% Transforms the prior specification to a Laplace type of approximation:
% only the prior mean and standard deviation are relevant, all other shape
% information, except for the parameter bounds, is ignored.
% -------------------------------------------------------------------------
% INPUTS
% bayestopt_    [structure] prior information
% -------------------------------------------------------------------------
% OUTPUT
% bayestopt_    [structure] Laplace prior information
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run
% -------------------------------------------------------------------------

% Copyright © 2023 Dynare Team
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

if any(setdiff([0;bayestopt_.pshape],[0,3]))
    fprintf('\nNon-normal priors specified. Penalized estimation uses a Laplace type of approximation:');
    fprintf('\nOnly the prior mean and standard deviation are relevant, all other shape information, except for the parameter bounds, is ignored.\n\n');
    non_normal_priors = (bayestopt_.pshape~=3);
    bayestopt_.pshape(non_normal_priors) = 3;
    bayestopt_.p3(non_normal_priors) = -Inf*ones(sum(non_normal_priors),1);
    bayestopt_.p4(non_normal_priors) = Inf*ones(sum(non_normal_priors),1);
    bayestopt_.p6(non_normal_priors) = bayestopt_.p1(non_normal_priors);
    bayestopt_.p7(non_normal_priors) = bayestopt_.p2(non_normal_priors);
    bayestopt_.p5(non_normal_priors) = bayestopt_.p1(non_normal_priors);
end
if any(isinf(bayestopt_.p2)) % find infinite variance priors
    inf_var_pars = bayestopt_.name(isinf(bayestopt_.p2));
    disp_string = [inf_var_pars{1,:}];
    for ii = 2:size(inf_var_pars,1)
        disp_string = [disp_string,', ',inf_var_pars{ii,:}];
    end
    fprintf('The parameter(s) %s have infinite prior variance. This implies a flat prior.\n',disp_string);
    fprintf('Dynare disables the matrix singularity warning\n');
    if isoctave
        warning('off','Octave:singular-matrix');
    else
        warning('off','MATLAB:singularMatrix');
    end
end