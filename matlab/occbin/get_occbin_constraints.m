function [i_base,i_alt,c_base,c_alt] = get_occbin_constraints(M,steady_state,ramsey_policy)

% Copyright (C) 2015 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

nrow = 1;
if ramsey_policy
    if isfield(M,'ramsey_model_constraints')
        rc = M.ramsey_model_constraints;
        for i = 1:length(rc)
            switch rc{i}{2}
              case {'>','>='}
                ivar(nrow) = rc{i}{1};
                ieq(nrow) = rc{i}{1};
                lb(nrow) = eval(rc{i}{3});
              case {'<','<='}
                ivar(nrow) = rc{i}{1};
                ieq(nrow) = rc{i}{1};
                ub(nrow) = eval(rc{i}{3});
              otherwise
                error('Wrong operator in get_complementarity_conditions')
            end
            nrow = nrow + 1;
        end
    end
end

i_base = {};
i_alt = {};
etags = M.equations_tags;
m = 1;
base = true;
for i=1:size(etags,1)
    [iv,boundary,operator] = parse_constraint(etags{i,3},M.endo_names,M.params,M.param_names);
    if strcmp(etags{i,2},'OCCBIN')
        if base
            i_alt{m} = 1:M.eq_nbr;
            i_alt{m}(etags{i,1}) = [];
            c_base{m,1} = etags{i,1};
            c_base{m,2} = iv;
            c_base{m,3} = boundary - steady_state(iv);
            c_base{m,4} = operator;
            base = false;
        else
            i_base{m} = 1:M.eq_nbr;
            i_base{m}(etags{i,1}) = [];
            c_alt{m,1} = etags{i,1};
            c_alt{m,2} = iv;
            c_alt{m,3} = boundary - steady_state(iv);
            c_alt{m,4} = operator;
            base = true;
            m = m + 1;
        end
    end
end
if ~base
    error('OCCBIN: constraints must come by pair')
end

function [iv,boundary,operator] = parse_constraint(str,endo_names,params,param_names)
delim = {'<=','>=','<','>'};
[c,operator] = strsplit(str,delim);
operator = operator{1};
iv = strmatch(strtrim(c{1}),endo_names);
% try for a number
boundary = str2num(strtrim(c{2}));
% if not a number try for a parameter name
if isempty(boundary)
    k = strmatch(strtrim(c{2}),param_names);
    if isempty(k)
        error(['OCCBIN: illegal constraint ' str]);
    end
    boundary = params(k);
end
