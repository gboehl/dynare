function errorcode = check(eqname, errorflag)

% Checks that error correction term is well defined in PAC equation.
%
% INPUTS
% - eqname       [string]    Name of the pac equation.
% OUTPUTS
% - errorcode    [integer]   Error code, positive if there is an issue with the PAC equation
%                            0  ->  No error,
%                            1  ->  eqname has to PAC expectation term,
%                            2  ->  LHS is not a variable in difference,
%                            3  ->  Possible calibration issue on the error correction term (should be positive),
%                            4  ->  Error correction term is missing.

% Copyright (C) 2018 Dynare Team
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

global M_

% Set default for second input argument.
if nargin<2
    errorflag = true;
end

% Set default for the returned argument.
errorcode = 0;

% Get the original equation to be estimated
[LHS, RHS] = get_lhs_and_rhs(eqname, M_, true);

% Check that the equation has a PAC expectation term.
if ~contains(RHS, 'pac_expectation', 'IgnoreCase', true)
    errorcode = 1;
    if errorflag
        skipline()
        error('This is not a PAC equation.')
    end
end

% Get the name of the PAC model.
pattern = '(\(model_name\s*=\s*)(?<name>\w+)\)';
pacmodl = regexp(RHS, pattern, 'names');
pacmodl = pacmodl.name;

% Get the index of the lhs variable from M_.pac.
lhsid = M_.pac.(pacmodl).lhs_var;

% Get info about lhs variable.
auxinfo = M_.aux_vars(get_aux_variable_id(lhsid));

% Check that the LHS variable is a variable in difference.
if ~isequal(auxinfo.type, 8)
    errorcode = 2;
    if errorflag
        skipline()
        error('LHS variable in %s has to be a variable in difference.', eqname)
    end
end

% Check that the level of the lhs variable appear on the RHS.
if ~ismember(auxinfo.orig_index, M_.pac.(pacmodl).ec.vars) || M_.params(M_.pac.(pacmodl).ec.params)<=0
    if errorflag
        msg = sprintf('\nPAC equation %s must have an error correction term of the form\n\n', eqname);
        msg = sprintf('%s  ... +  a*(x(-1)-y(-1)) + ... \n\n', msg);
        msg = sprintf('%swhere a (%s) is a positive parameter, x (%s) is a the trend/target,\n', ...
                      msg, M_.param_names{M_.pac.(pacmodl).ec.params}, M_.endo_names{M_.pac.(pacmodl).ec.vars(1)});
        msg = sprintf('%sand y (%s) is the level of the LHS variable.\n\n', ...
                      msg, M_.endo_names{auxinfo.orig_index});
        if M_.params(M_.pac.(pacmodl).ec.params)<=0
            msg = sprintf('%sPlease change the calibration.\n', msg);
        else
            msg = sprintf('%sPlease change the PAC equation.\n', msg);
        end
        skipline()
        error(msg)
    else
        if M_.params(M_.pac.(pacmodl).ec.params)<=0
            errorcode = 3;
        else
            errorcode = 4;
        end
    end
end