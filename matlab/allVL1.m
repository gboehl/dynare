function v = allVL1(n, L1, L1ops, MaxNbSol)
% All integer permutations with sum criteria
%
% function v=allVL1(n, L1); OR
% v=allVL1(n, L1, L1opt);
% v=allVL1(n, L1, L1opt, MaxNbSol);
% 
% INPUT
%    n: length of the vector
%    L1: target L1 norm
%    L1ops: optional string ('==' or '<=' or '<')
%           default value is '=='
%    MaxNbSol: integer, returns at most MaxNbSol permutations.
%    When MaxNbSol is NaN, allVL1 returns the total number of all possible
%    permutations, which is useful to check the feasibility before getting
%    the permutations.
% OUTPUT:
%    v: (m x n) array such as: sum(v,2) == L1,
%       (or <= or < depending on L1ops)                            
%       all elements of v is naturel numbers {0,1,...}
%       v contains all (=m) possible combinations
%       v is sorted by sum (L1 norm), then by dictionnary sorting criteria
%    class(v) is same as class(L1) 
% Algorithm:
%    Recursive
% Remark:
%    allVL1(n,L1-n)+1 for natural numbers defined as {1,2,...}
% Example:
%    This function can be used to generate all orders of all
%    multivariable polynomials of degree p in R^n:
%         Order = allVL1(n, p)
% Author: Bruno Luong
% Original, 30/nov/2007
% Version 1.1, 30/apr/2008: Add H1 line as suggested by John D'Errico
%         1.2, 17/may/2009: Possibility to get the number of permutations
%                           alone (set fourth parameter MaxNbSol to NaN)
%         1.3, 16/Sep/2009: Correct bug for number of solution
%         1.4, 18/Dec/2010: + non-recursive engine

% Retrieved from https://www.mathworks.com/matlabcentral/fileexchange/17818-all-permutations-of-integers-with-sum-criteria

% Copyright (c) 2009 Bruno Luong
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

global MaxCounter;

if nargin<3 || isempty(L1ops)
    L1ops = '==';
end

n = floor(n); % make sure n is integer

if n<1
    v = [];
    return
end

if nargin<4  || isempty(MaxNbSol)
    MaxCounter = Inf;
else
    MaxCounter = MaxNbSol;
end
Counter(0);

switch L1ops
    case {'==' '='},
        if isnan(MaxCounter)
            % return the number of solutions
            v = nchoosek(n+L1-1,L1); % nchoosek(n+L1-1,n-1)
        else
            v = allVL1eq(n, L1);
        end
    case '<=', % call allVL1eq for various sum targets
        if isnan(MaxCounter)
            % return the number of solutions
            %v = nchoosek(n+L1,L1)*factorial(n-L1); BUG <- 16/Sep/2009: 
            v = 0;
            for j=0:L1
                v = v + nchoosek(n+j-1,j);
            end
            % See pascal's 11th identity, the sum doesn't seem to
            % simplify to a fix formula
        else
            v = cell2mat(arrayfun(@(j) allVL1eq(n, j), (0:L1)', ...
                         'UniformOutput', false));
        end
    case '<',
        v = allVL1(n, L1-1, '<=', MaxCounter);
    otherwise
        error('allVL1: unknown L1ops')
end

end % allVL1

%%
function v = allVL1eq(n, L1)

global MaxCounter;

n = feval(class(L1),n);
s = n+L1;
sd = double(n)+double(L1);
notoverflowed = double(s)==sd;
if isinf(MaxCounter) && notoverflowed
    v = allVL1nonrecurs(n, L1);
else
    v = allVL1recurs(n, L1);
end

end % allVL1eq

%% Recursive engine
function v = allVL1recurs(n, L1, head)
% function v=allVL1eq(n, L1);
% INPUT
%    n: length of the vector
%    L1: desired L1 norm
%    head: optional parameter to by concatenate in the first column
%          of the output
% OUTPUT:
%    if head is not defined
%      v: (m x n) array such as sum(v,2)==L1
%         all elements of v is naturel numbers {0,1,...}
%         v contains all (=m) possible combinations
%         v is (dictionnary) sorted
% Algorithm:
%    Recursive

global MaxCounter;

if n==1
    if Counter < MaxCounter
        v = L1;
    else
        v = zeros(0,1,class(L1));
    end
else % recursive call
    v = cell2mat(arrayfun(@(j) allVL1recurs(n-1, L1-j, j), (0:L1)', ...
                 'UniformOutput', false));
end

if nargin>=3 % add a head column
    v = [head+zeros(size(v,1),1,class(head)) v];
end

end % allVL1recurs

%%
function res=Counter(newval)
    persistent counter;
    if nargin>=1
        counter = newval;
        res = counter;
    else
        res = counter;
        counter = counter+1;
    end
end % Counter

%% Non-recursive engine
function v = allVL1nonrecurs(n, L1)
% function v=allVL1eq(n, L1);
% INPUT
%    n: length of the vector
%    L1: desired L1 norm
% OUTPUT:
%    if head is not defined
%      v: (m x n) array such as sum(v,2)==L1
%         all elements of v is naturel numbers {0,1,...}
%         v contains all (=m) possible combinations
%         v is (dictionnary) sorted
% Algorithm:
%    NonRecursive

% Chose (n-1) the splitting points of the array [0:(n+L1)]
s = nchoosek(1:n+L1-1,n-1);
m = size(s,1);

s1 = zeros(m,1,class(L1));
s2 = (n+L1)+s1;

v = diff([s1 s s2],1,2); % m x n
v = v-1;

end % allVL1nonrecurs
