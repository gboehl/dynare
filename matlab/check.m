function [eigenvalues_,result,info] = check(M_, options_, oo_)
%[eigenvalues_,result,info] = check(M_, options_, oo_)
% Checks determinacy conditions by computing the generalized eigenvalues.
%
% INPUTS
% - M_            [structure]     Matlab's structure describing the model
% - options_      [structure]     Matlab's structure describing the current options
% - oo_           [structure]     Matlab's structure containing the results
%
% OUTPUTS
% - eigenvalues_  [double]        vector, eigenvalues.
% - result        [integer]       scalar, equal to 1 if Blanchard and Kahn conditions are satisfied, zero otherwise.
% - info          [integer]       scalar or vector, error code as returned by resol routine.

% Copyright © 2001-2023 Dynare Team
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


if ~options_.initval_file && M_.exo_nbr > 1
    oo_.exo_simul = ones(M_.maximum_lead+M_.maximum_lag+1,1)*oo_.exo_steady_state';
end

options_.order = 1;

if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end

oo_.dr=set_state_space(oo_.dr,M_);

[dr,info] = resol(1,M_,options_,oo_.dr ,oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);

if info(1) ~= 0 && info(1) ~= 3 && info(1) ~= 4
    print_info(info, 0, options_);
end

eigenvalues_ = dr.eigval;
[m_lambda,i]=sort(abs(eigenvalues_));

result = 0;
if (M_.nsfwrd == dr.edim) && (dr.full_rank)
    result = 1;
end

if ~options_.noprint
    skipline()
    disp('EIGENVALUES:')
    fprintf('%16s %16s %16s\n','Modulus','Real','Imaginary');
    z=[m_lambda real(eigenvalues_(i)) imag(eigenvalues_(i))]';
    fprintf('%16.4g %16.4g %16.4g\n',z);
    fprintf('\nThere are %d eigenvalue(s) larger than 1 in modulus ', dr.edim);
    fprintf('for %d forward-looking variable(s)', M_.nsfwrd);
    skipline()
    if result
        disp('The rank condition is verified.')
    else
        disp('The rank condition ISN''T verified!')
    end
    skipline()
end
