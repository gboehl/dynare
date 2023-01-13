function [eigenvalues_,result,info] = check(M, options, oo)

% Checks determinacy conditions by computing the generalized eigenvalues.
%
% INPUTS
% - M             [structure]     Matlab's structure describing the model (M_).
% - options       [structure]     Matlab's structure describing the current options (options_).
% - oo            [structure]     Matlab's structure containing the results (oo_).
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


if ~options.initval_file && M.exo_nbr > 1
    oo.exo_simul = ones(M.maximum_lead+M.maximum_lag+1,1)*oo.exo_steady_state';
end

options.order = 1;

if isempty(options.qz_criterium)
    options.qz_criterium = 1+1e-6;
end

oo.dr=set_state_space(oo.dr,M,options);

[dr,info,M,~] = resol(1,M,options,oo);

if info(1) ~= 0 && info(1) ~= 3 && info(1) ~= 4
    print_info(info, 0, options);
end

eigenvalues_ = dr.eigval;
[m_lambda,i]=sort(abs(eigenvalues_));

result = 0;
if (M.nsfwrd == dr.edim) && (dr.full_rank)
    result = 1;
end

if ~options.noprint
    skipline()
    disp('EIGENVALUES:')
    disp(sprintf('%16s %16s %16s\n','Modulus','Real','Imaginary'))
    z=[m_lambda real(eigenvalues_(i)) imag(eigenvalues_(i))]';
    disp(sprintf('%16.4g %16.4g %16.4g\n',z))
    disp(sprintf('\nThere are %d eigenvalue(s) larger than 1 in modulus ', dr.edim));
    disp(sprintf('for %d forward-looking variable(s)', M.nsfwrd));
    skipline()
    if result
        disp('The rank condition is verified.')
    else
        disp('The rank condition ISN''T verified!')
    end
    skipline()
end
