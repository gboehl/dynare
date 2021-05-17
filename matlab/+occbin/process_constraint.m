function constraint_parsed = process_constraint(constraint,suffix,endo_names,invert_switch,param_names)
%function constraint_parsed = process_constraint(constraint,suffix,endo_names,invert_switch,param_names)
% Processes constraints for use in Occbin, appending endogenous variables
% with suffix and replacing parameters by their value in M_.params
%
% INPUTS
% - constraint      [char]     constraint to be parsed
% - suffix          [char]     suffix to be appended
% - endo_names      [cell]     names of endogenous variables
% - invert_switch   [bool]     if true, invert direction of the inequality constraint
% - param_names     [cell]     names of parameters
%
% OUTPUTS
% - constraint_parsed   [char]     parsed constraint 

% Original authors: Luca Guerrieri and Matteo Iacoviello 
% Original file downloaded from:
% https://www.matteoiacoviello.com/research_files/occbin_20140630.zip
% Adapted for Dynare by Dynare Team.
%
% This code is in the public domain and may be used freely.
% However the authors would appreciate acknowledgement of the source by
% citation of any of the following papers:
%
% Luca Guerrieri and Matteo Iacoviello (2015): "OccBin: A toolkit for solving
% dynamic models with occasionally binding constraints easily"
% Journal of Monetary Economics 70, 22-38

% create a list of delimiters that can separate parameters and endogenoous
% variables in the string that expresses the constraint
delimiters = char(',',';','(',')','+','-','^','*','/','>','<','=');

% split the string that holds the constraint into tokens
tokens = occbin.tokenize(constraint,delimiters);

ntokens = length(tokens);

endo_ss=strcat(endo_names,repmat('_ss',length(endo_names),1));

% search for tokens that match the list of endogenous variables
for i=1:ntokens
    valid_token=0;
    if ~isempty(find(strcmp(tokens(i),delimiters)))
        % if the invert_switch is true
        % reverse the direction of the inequality
        if invert_switch
            if  strcmp(tokens(i),cellstr('>'))
                tokens(i) = cellstr('<');
            elseif strcmp(tokens(i),cellstr('<'))
                tokens(i) = cellstr('>');
            end
        end
        valid_token=1;
        continue;
    end
        
    if ~isempty(find(strcmp(tokens(i),endo_names)))
        % when there is a match with an endogenous variable append the suffix
        tokens(i) = cellstr([char(tokens(i)),suffix]);
        valid_token=1;
        continue;
    end
    par_index=find(strcmp(tokens(i),param_names));
    if ~isempty(par_index)
        tokens(i) = {['M_.params(',num2str(par_index),')']};
        valid_token=1;
        continue;
    end
    ss_index=find(strcmp(tokens(i),endo_ss));
    if ~isempty(ss_index)
        tokens(i) = {['ys(',num2str(ss_index),')']};
        valid_token=1;
        continue;
    end
    if isnan(str2double(tokens(i)))
        error('Occbin: Constraint %s contains the uninterpretable token %s', constraint, tokens{i});
    end
end

% reassemble the tokens to create a string that expresses the constraint
constraint_parsed = regexprep(strjoin(tokens),' ','');