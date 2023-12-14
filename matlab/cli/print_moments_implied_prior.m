function print_moments_implied_prior(M_, mm, vm, mv, vv)
%function print_moments_implied_prior(M_, mm, vm, mv, vv)
% This routine prints in the command window some descriptive statistics
% about the endogenous variables implied prior moments.
% Inputs:
%   - M_            [structure]             Dynare's model structure
%   - mm            [endo_nbr*1]            mean first moments of the endogenous
%                                           variables
%   - vm            [endo_nbr*1]            variance of the first moments of the
%                                           endogenous variables
%   - mv            [endo_nbr*endo_nbr]     mean first moments of the endogenous
%                                           variables
%   - vv            [endo_nbr]              variance of the first moments of the
%                                           endogenous variables


% Copyright © 2016-2023 Dynare Team
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

% First order moments.

disp('(Implied) Prior moments of the endogenous variables'' expectation')
disp(printline(64, '-'))

T1 = 'VARIABLE ';
T2 = sprintf('Prior mean \t Prior st. dev.');

for i=1:M_.orig_endo_nbr
    Name = M_.endo_names{i};
    T1 = strvcat(T1, Name);
    str = sprintf(' %6.4f \t %6.4f', mm(i), sqrt(vm(i)));
    T2 = strvcat(T2, str);
end

T0 = repmat('  ', M_.orig_endo_nbr+1, 1);

TT = [T1, T0, T2];
l0 = printline(size(TT, 2)+1, '-');
TT = strvcat(l0, TT(1,:), l0, TT(2:end,:), l0);

skipline(2)
disp(TT)
skipline(2)

disp('(Implied) Prior moments of the endogenous variables'' (co)variance')
disp(printline(61, '-'))

T1a = 'VARIABLE-1';
T1b = 'VARIABLE-2';
T2a = 'Prior mean';
T2b = 'Prior st.dev.';

for i=1:M_.orig_endo_nbr
    for j=i:M_.orig_endo_nbr
        Name1 = M_.endo_names{i};
        Name2 = M_.endo_names{j};
        T1a = strvcat(T1a, Name1);
        T1b = strvcat(T1b, Name2);
        sta = sprintf('%12.8f', mv(i,j));
        stb = sprintf('%12.8f', vv(i,j));
        T2a = strvcat(T2a, sta);
        T2b = strvcat(T2b, stb);
    end
end

T0 = repmat('  ', M_.orig_endo_nbr*(M_.orig_endo_nbr+1)/2+1, 1);

TT = [T1a, T0, T1b, T0, T2a, T0, T2b];
l0 = printline(size(TT, 2)+1, '-');
TT = strvcat(l0, TT(1,:), l0, TT(2:end,:), l0);

skipline(2)
disp(TT)
skipline(2)