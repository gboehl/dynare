function [M_, oo_, options_] = initialize(M_,oo_,options_)
% function [M_, oo_, options_] = initialize(M_,oo_,options_)
% Initializes run of Occbin: sets default options, parses constraints, creates
% eval_difference-file
%
% INPUT: 
% - M_                  [structure]     Matlab's structure describing the model
% - oo_                 [structure]     Matlab's structure containing the results
% - options_            [structure]     Matlab's structure containing the options
%
% OUTPUT: 
% - M_                  [structure]     Matlab's structure describing the model
% - oo_                 [structure]     Matlab's structure containing the results
% - options_            [structure]     Matlab's structure containing the options

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

options_.occbin = struct();
options_.occbin = occbin.set_default_options(options_.occbin, M_);
   
oo_.dr=set_state_space(oo_.dr,M_,options_);

M_ = occbin.get_info(M_);

if M_.occbin.constraint_nbr>2
    error('Occbin: Only up to two constraints are supported')
end
for k=1:M_.occbin.constraint_nbr
    M_.occbin.constraint(k).bind_difference = occbin.process_constraint(M_.occbin.constraint(k).bind,'_difference',M_.endo_names,false,M_.param_names);
    M_.occbin.constraint(k).bind_error = occbin.process_error_constraint(M_.occbin.constraint(k).bind_difference);
    if isempty(M_.occbin.constraint(k).relax)
        M_.occbin.constraint(k).relax_difference = ['~binding.constraint_',int2str(k)];
    else
        M_.occbin.constraint(k).relax_difference = occbin.process_constraint(M_.occbin.constraint(k).relax,'_difference',M_.endo_names,false,M_.param_names);
    end
    M_.occbin.constraint(k).relax_error = occbin.process_error_constraint(M_.occbin.constraint(k).relax_difference);
    M_.occbin.constraint(k).pswitch_index = strmatch(M_.occbin.constraint(k).pswitch,M_.param_names);
end
nwishes_ = size(M_.occbin.wish_list.endo,1);
fid = fopen(['+' M_.fname '/eval_difference.m'],'w+');
fprintf(fid,'function [binding, relax, err]=eval_difference(zdatalinear_,M_,ys);\n');
for i_indx_=1:nwishes_
    fprintf(fid,'%s_difference=zdatalinear_(:,%i);\r\n',deblank(M_.occbin.wish_list.endo(i_indx_,:)),M_.occbin.wish_list.iendo(i_indx_));
end
for k=1:M_.occbin.constraint_nbr
    fprintf(fid,'binding.constraint_%i = %s;\r\n',k,M_.occbin.constraint(k).bind_difference);
    fprintf(fid,'relax.constraint_%i = %s;\r\n',k,M_.occbin.constraint(k).relax_difference);
%         fprintf(fid,'try\r\n');            

    fprintf(fid,'    err.binding_constraint_%i = abs(%s);\r\n',k,M_.occbin.constraint(k).bind_error);
    fprintf(fid,'    err.relax_constraint_%i = abs(%s);\r\n',k,M_.occbin.constraint(k).relax_error);

%         fprintf(fid,'catch\r\n');            
%         
%         fprintf(fid,'    err.newviolvecbool_%i = nan(size(err.newviolvecbool_%i));\r\n',k,k);
%         fprintf(fid,'    err.relaxconstraint_%i = nan(size(err.relaxconstraint_%i));\r\n',k,k);

%         fprintf(fid,'end\r\n');            
end
fclose(fid);
