function M_ = get_info(M_)
%function M_ = get_info(M_)
% Parses constraint to clean spaces and provide information on endogenous variables 
% involved in constraints
%
% INPUTS
% - M_         [struct]     Definition of the model.
%
% OUTPUTS
% - M_         [struct]     Definition of the model.

% Copyright (C) 2021 Dynare Team
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

delimiters = char(',',';','(',')','+','-','^','*','/','>','<','=');
iwish=0;
indx_=[];
wish_list='';
M_.occbin.constraint_nbr=length(M_.occbin.constraint);
for k=1:M_.occbin.constraint_nbr
    %get rid of redundant spaces;
    M_.occbin.constraint(k).bind=strrep(M_.occbin.constraint(k).bind,' ','');
    M_.occbin.constraint(k).relax=strrep(M_.occbin.constraint(k).relax,' ','');
    bind_clean = occbin.tokenize(M_.occbin.constraint(k).bind,delimiters);
    relax_clean = occbin.tokenize(M_.occbin.constraint(k).relax,delimiters);
    for j=1:M_.endo_nbr
        tmp=ismember(M_.endo_names{j},bind_clean);
        tmp1=ismember(M_.endo_names{j},relax_clean);
        if tmp || tmp1
            iwish = iwish+1;
            indx_(iwish)=j;
            if iwish>1
                wish_list = char(wish_list,M_.endo_names{j});
            else
                wish_list = M_.endo_names{j};
            end
        end
    end
end

M_.occbin.wish_list.endo = wish_list;
M_.occbin.wish_list.iendo = indx_;
end