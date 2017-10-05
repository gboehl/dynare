function [ivar,ieq,lb,ub] = get_occbin_complementarity_conditions(M,ramsey_policy)

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

etags = M.equations_tags;
for i=1:size(etags,1)
    if strcmp(etags{i,2},'mcp')
        str = etags{i,3};
        kop = strfind(etags{i,3},'<');
        if ~isempty(kop)
                k = find(strcmp(strtrim(str(1:kop-1)),cellstr(M.endo_names)));
                if isempty(k)
                    error(sprintf(['Complementarity condition %s: variable %s is ' ...
                                   'not recognized',etags{i,3},b{1}]))
                end
                ivar(nrow) = k;
                ieq(nrow) = etags{i,1};
                ub(nrow) = eval(str(kop+1:end));
        else
            kop = strfind(etags{i,3},'>');
            if ~isempty(kop)
                k = find(strcmp(strtrim(str(1:kop-1)),cellstr(M.endo_names)));
                if isempty(k)
                    error(sprintf(['Complementarity condition %s: variable %s is ' ...
                                   'not recognized',etags{i},b{1}]))
                end
                ivar(nrow) = k;
                ieq(nrow) = etags{i,1};
                lb(k) = eval(str(kop+1:end));
            else
                error(sprintf(['Complementarity condition %s can''t be ' ...
                               'parsed'],etags{i,3}))
            end
        end
        nrow = nrow + 1;
    end
end

