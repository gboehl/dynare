function tokens = tokenize(source,delimiter)
%function tokens = tokenize(source,delimiter)
% Breaks down strings into its components (tokens)
% INPUTS
% - source          [string]        string to be broken into tokens
% - delimiter       [char array]    single character delimiters
%
% OUTPUTS:
% -tokens           [cell]          array containing the tokens
 
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

posdelims = [];

% assumes that delimiter cannot be in the first position or the last position
ndelimiters = size(delimiter,1);
for i=1:ndelimiters
    newpositions = strfind(source,delimiter(i,:));
    if ~isempty(newpositions)
        posdelims =[posdelims, newpositions];
    end
end

% reorder posdelims in ascending order
posdelims = sort(posdelims);

if isempty(posdelims)
    tokens = cellstr(source);
else
    ndelims = length(posdelims);
    % build positions for substrings
    delims = zeros(ndelims+1,2);
    for i=1:ndelims+1
        if i==1
            if posdelims(1) == 1
                tokens = cellstr(source(1));
            else
                delims(i,:) = [1,posdelims(i)-1];
                tokens = cellstr(source([delims(i,1):delims(i,2)]));
                tokens = [tokens, source(posdelims(i))];
            end
        elseif  i==ndelims+1
            if (posdelims(i-1) < length(source))
                delims(i,:) = [posdelims(i-1)+1,length(source)];
                tokens = [tokens, cellstr(source([delims(i,1):delims(i,2)]))];
            end
        else
            if posdelims(i)>posdelims(i-1)+1
                delims(i,:) = [posdelims(i-1)+1,posdelims(i)-1];
                tokens = [tokens, cellstr(source([delims(i,1):delims(i,2)]))];
            end
            tokens = [tokens, source(posdelims(i))];
        end
    end    
end