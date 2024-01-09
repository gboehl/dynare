function model_info(options_model_info_)
%function model_info(options_model_info_)

% Copyright © 2008-2022 Dynare Team
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

global M_;

dynamic_ = isfield(options_model_info_, 'block_dynamic') && options_model_info_.block_dynamic;
static_ = isfield(options_model_info_, 'block_static') && options_model_info_.block_static;
incidence = isfield(options_model_info_, 'incidence') && options_model_info_.incidence;

if static_
    temp_string=sprintf('\nInformation about %s (static model)\n',M_.fname);
    fprintf(temp_string);
    block_structure_str = 'block_structure_stat';
    if ~isfield(M_,'block_structure_stat')
        fprintf('\nmodel_info: block information not present; skipping display.\n')
        return;
    else
        nb_leadlag = 1;
    end
else
    temp_string=sprintf('\nInformation about %s (dynamic model)\n',M_.fname);
    fprintf(temp_string);
    block_structure_str = 'block_structure';
    if dynamic_ && isfield(M_,'block_structure')
        fprintf('\nmodel_info: block information not present; skipping display.\n')        
        return;
    elseif dynamic_
        nb_leadlag = length([M_.(block_structure_str).incidence.lead_lag]);
    end
end


if dynamic_ || static_ || incidence %block information requested
    if static_
        block_structure = M_.block_structure_stat;
    else
        block_structure = M_.block_structure;
    end
    fprintf([char(ones(1,length(temp_string))*'='),'\n']);
    nb_blocks=length(block_structure.block);
    fprintf('The model has %d equations and is decomposed into %d blocks as follows:\n',M_.endo_nbr,nb_blocks);
    fprintf('================================================================================================================================\n');
    fprintf('| %10s | %10s | %30s | %31s | %31s |\n','Block no','Size','Block Type','Equation','Dependent variable');
    fprintf('|============|============|================================|=================================|=================================|\n');
    for i=1:nb_blocks
        size_block=length(block_structure.block(i).equation);
        if(i>1)
            fprintf('|------------|------------|--------------------------------|---------------------------------|---------------------------------|\n');
        end
        for j=1:size_block
            if(j==1)
                fprintf('| %10d | %10d | %30s | %-6d %24s | %-6d %24s |\n',i,size_block,Sym_type(block_structure.block(i).Simulation_Type),block_structure.block(i).equation(j),get_equation_name_by_number(block_structure.block(i).equation(j), M_),block_structure.block(i).variable(j),M_.endo_names{block_structure.block(i).variable(j)});
            else
                fprintf('| %10s | %10s | %30s | %-6d %24s | %-6d %24s |\n','','','',block_structure.block(i).equation(j),get_equation_name_by_number(block_structure.block(i).equation(j), M_),block_structure.block(i).variable(j),M_.endo_names{block_structure.block(i).variable(j)});
            end
        end
    end
    fprintf('================================================================================================================================\n');
    fprintf('\n');
    if static_
        fprintf('%-30s %s','The variable','is used contemporaneously in the following equations:');
        if(size(block_structure.incidence.sparse_IM,1)>0)
            IM=sortrows(block_structure.incidence.sparse_IM,2);
        else
            IM=[];
        end
        size_IM=size(IM,1);
        last=99999999;
        for i=1:size_IM
            if(last~=IM(i,2))
                fprintf('\n%-30s',M_.endo_names{IM(i,2)});
            end
            fprintf(' %5d',IM(i,1));
            last=IM(i,2);
        end
        fprintf('\n\n');
    else %dynamic model
        for k=1:M_.maximum_endo_lag+M_.maximum_endo_lead+1
            if(k==M_.maximum_endo_lag+1)
                fprintf('%-30s %s','The variable','is used in the following equations contemporaneously');
            elseif(k<M_.maximum_endo_lag+1)
                fprintf('%-30s %s %d','The variable','is used in the following equations with lag ',M_.maximum_endo_lag+1-k);
            else
                fprintf('%-30s %s %d','The variable','is used in equations with lead ',k-(M_.maximum_endo_lag+1));
            end
            if(size(block_structure.incidence(k).sparse_IM,1)>0)
                IM=sortrows(block_structure.incidence(k).sparse_IM,2);
            else
                IM=[];
            end
            size_IM=size(IM,1);
            last=99999999;
            for i=1:size_IM
                if(last~=IM(i,2))
                    fprintf('\n%-30s',M_.endo_names{IM(i,2)});
                end
                fprintf(' %5d',IM(i,1));
                last=IM(i,2);
            end
            fprintf('\n\n');
        end
    end
    if incidence

        %printing the gross incidence matrix
        IM_star = char([kron(ones(M_.endo_nbr, M_.endo_nbr-1), double(blanks(3))) double(blanks(M_.endo_nbr)')]);
        for i = 1:nb_leadlag
            n = size(block_structure.incidence(i).sparse_IM,1);
            for j = 1:n
                if ismember(block_structure.incidence(i).sparse_IM(j,2), M_.state_var)
                    IM_star(block_structure.incidence(i).sparse_IM(j,1), 3 * (block_structure.incidence(i).sparse_IM(j,2) - 1) + 1) = 'X';
                else
                    IM_star(block_structure.incidence(i).sparse_IM(j,1), 3 * (block_structure.incidence(i).sparse_IM(j,2) - 1) + 1) = '1';
                end
            end
        end
        seq = 1: M_.endo_nbr;
        blank = [ blanks(cellofchararraymaxlength(M_.endo_names)); blanks(cellofchararraymaxlength(M_.endo_names))];
        for i = 1:M_.endo_nbr
            if i == 1
                var_names = char(blank, M_.endo_names{i});
            else
                var_names = char(var_names, blank, M_.endo_names{i});
            end
        end
        topp = [char(kron(double(blanks(ceil(log10(M_.endo_nbr)))),ones(cellofchararraymaxlength(M_.endo_names),1))) var_names' ];
        bott = [int2str(seq') blanks(M_.endo_nbr)' blanks(M_.endo_nbr)' IM_star];
        fprintf('\nGross incidence matrix\n');
        fprintf('=======================\n');
        disp(topp);
        skipline;
        disp(bott);
        
        %printing the reordered incidence matrix
        IM_star_reordered = char([kron(ones(M_.endo_nbr, M_.endo_nbr-1), double(blanks(3))) double(blanks(M_.endo_nbr)')]);
        eq(block_structure.equation_reordered) = seq;
        va(block_structure.variable_reordered) = seq;
        barre_blank = [ barre(cellofchararraymaxlength(M_.endo_names)); blanks(cellofchararraymaxlength(M_.endo_names))];
        cur_block = 1;
        for i = 1:M_.endo_nbr
            past_block = cur_block;
            while ismember(block_structure.variable_reordered(i), block_structure.block(cur_block).variable) == 0
                cur_block = cur_block + 1;
            end
            if i == 1
                var_names = char(blank, M_.endo_names{block_structure.variable_reordered(i)});
            else
                if past_block ~= cur_block
                    var_names = char(var_names, barre_blank, M_.endo_names{block_structure.variable_reordered(i)});
                else
                    var_names = char(var_names, blank, M_.endo_names{block_structure.variable_reordered(i)});
                end
            end
        end
        topp = [char(kron(double(blanks(ceil(log10(M_.endo_nbr)))),ones(cellofchararraymaxlength(M_.endo_names),1))) var_names' ];
        n_state_var = length(M_.state_var);
        IM_state_var = zeros(n_state_var, n_state_var);
        inv_variable_reordered(block_structure.variable_reordered) = 1:M_.endo_nbr;
        state_equation = block_structure.equation_reordered(inv_variable_reordered(M_.state_var));
        for i = 1:nb_leadlag
            n = size(block_structure.incidence(i).sparse_IM,1);
            for j = 1:n
                [tf, loc] = ismember(block_structure.incidence(i).sparse_IM(j,2), M_.state_var);
                if tf
                    IM_star_reordered(eq(block_structure.incidence(i).sparse_IM(j,1)), 3 * (va(block_structure.incidence(i).sparse_IM(j,2)) - 1) + 1) = 'X';
                    [tfi, loci] = ismember(block_structure.incidence(i).sparse_IM(j,1), state_equation);
                    if tfi
                        IM_state_var(loci, loc) = 1;
                    end
                else
                    IM_star_reordered(eq(block_structure.incidence(i).sparse_IM(j,1)), 3 * (va(block_structure.incidence(i).sparse_IM(j,2)) - 1) + 1) = '1';
                end
            end
        end
        fprintf('\n1: non-null element, X: non-null element related to a state variable\n');
        
        cur_block = 1;
        i_last = 0;
        block = {};
        for i = 1:n_state_var
            past_block = cur_block;
            while ismember(M_.state_var(i), block_structure.block(cur_block).variable) == 0
                cur_block = cur_block + 1;
            end
            if (past_block ~= cur_block) || (past_block == cur_block && i == n_state_var)
                block(past_block).IM_state_var(1:(i - 1 - i_last), 1:i - 1) = IM_state_var(i_last+1:i - 1, 1:i - 1);
                i_last = i - 1;
            end
        end
        cur_block = 1;
        for i = 1:M_.endo_nbr
            past_block = cur_block;
            while ismember(block_structure.variable_reordered(i), block_structure.block(cur_block).variable) == 0
                cur_block = cur_block + 1;
            end
            if past_block ~= cur_block
                for j = 1:i-1
                    IM_star_reordered(j, 3 * (i - 1) - 1) = '|';
                end
            end
        end
        
        bott = [int2str(block_structure.equation_reordered') blanks(M_.endo_nbr)' blanks(M_.endo_nbr)' IM_star_reordered];
        fprintf('\nReordered incidence matrix\n');
        fprintf('==========================\n');
        disp(topp);
        skipline;
        disp(bott);
        fprintf('\n1: non-null element, X: non-null element related to a state variable\n');
    end
else %non-block information
    % print states
    if M_.maximum_endo_lag~=0
        lag_index=find(M_.lead_lag_incidence(1,:));
        fprintf('\nThe following variables appear with a lag and are therefore states:\n')
        for var_iter=1:length(lag_index)
            print_line(M_.endo_names,lag_index(var_iter),-1,M_)
        end
    else
        lag_index=[];
    end
    
    %print forward-looking variables
    if M_.maximum_endo_lead~=0
        lead_index = find(M_.lead_lag_incidence(M_.maximum_lag+2,:));
        fprintf('\nThe following variables appear with a lead and are therefore forward-looking variables:\n')
        for var_iter=1:length(lead_index)
            print_line(M_.endo_names,lead_index(var_iter),1,M_)
        end
    else
        lead_index=[];
    end
    
    %print purely static ones
    static_index = setdiff(1:M_.endo_nbr,union(lag_index,lead_index));
    if ~isempty(static_index)
        fprintf('\nThe following variables do not appear with a lead or lag and are therefore purely static variables:\n')
        for var_iter=1:length(static_index)
            print_line(M_.endo_names,static_index(var_iter),0,M_)
        end
    end
    skipline;
end

end

function print_line(names,var_index,lead_lag,M_)

    if var_index<=M_.orig_endo_nbr
        if iscell(names)
            if lead_lag==0
                fprintf('%s\n',names{var_index})
            else
                fprintf('%s(%d)\n',names{var_index},lead_lag)
            end
        else
            if lead_lag==0
                fprintf('%s\n',names(var_index))
            else
                fprintf('%s(%d)\n',names(var_index),lead_lag)
            end
        end
    else
        aux_index=find([M_.aux_vars(:).endo_index]==var_index);
        aux_type=M_.aux_vars(aux_index).type;
        if lead_lag==0
            str = subst_auxvar(var_index, [], M_);
        else
            str = subst_auxvar(var_index, lead_lag, M_);
        end
        aux_orig_expression=M_.aux_vars(aux_index).orig_expr;
        if isempty(aux_orig_expression)
            fprintf('%s\n',str);
        else
            fprintf('%s (original expression %s) \n',str,aux_orig_expression);
        end
    end

end

function ret=Sym_type(type)
    UNKNOWN=0;
    EVALUATE_FORWARD=1;
    EVALUATE_BACKWARD=2;
    SOLVE_FORWARD_SIMPLE=3;
    SOLVE_BACKWARD_SIMPLE=4;
    SOLVE_TWO_BOUNDARIES_SIMPLE=5;
    SOLVE_FORWARD_COMPLETE=6;
    SOLVE_BACKWARD_COMPLETE=7;
    SOLVE_TWO_BOUNDARIES_COMPLETE=8;
    EVALUATE_FORWARD_R=9;
    EVALUATE_BACKWARD_R=10;
    switch (type)
        case (UNKNOWN)
            ret='UNKNOWN                     ';
        case {EVALUATE_FORWARD,EVALUATE_FORWARD_R}
            ret='EVALUATE FORWARD            ';
        case {EVALUATE_BACKWARD,EVALUATE_BACKWARD_R}
            ret='EVALUATE BACKWARD            ';
        case SOLVE_FORWARD_SIMPLE
            ret='SOLVE FORWARD SIMPLE        ';
        case SOLVE_BACKWARD_SIMPLE
            ret='SOLVE BACKWARD SIMPLE        ';
        case SOLVE_TWO_BOUNDARIES_SIMPLE
            ret='SOLVE TWO BOUNDARIES SIMPLE  ';
        case SOLVE_FORWARD_COMPLETE
            ret='SOLVE FORWARD COMPLETE      ';
        case SOLVE_BACKWARD_COMPLETE
            ret='SOLVE BACKWARD COMPLETE      ';
        case SOLVE_TWO_BOUNDARIES_COMPLETE
            ret='SOLVE TWO BOUNDARIES COMPLETE';
    end
end

function ret = barre(n)
    s = [];
    for i=1:n
        s = [s '|'];
    end
    ret = s;
end
