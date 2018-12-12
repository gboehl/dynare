function rnd = gamrnd(a, b, method)

% This function produces independent random variates from the Gamma distribution.
%
%  INPUTS
%    a       [double]    n*1 vector of positive parameters.
%    b       [double]    n*1 vector of positive parameters.
%    method  [string]    'BawensLubranoRichard' or anything else (see below).
%
%  OUTPUT
%    rnd     [double]    n*1 vector of independent variates from the gamma(a,b) distribution.
%                        rnd(i) is gamma distributed with mean a(i)b(i) and variance a(i)b(i)^2.
%
%  ALGORITHMS
%    Described in Bauwens, Lubrano and Richard (1999, page 316) and Devroye (1986, chapter 9).

% Copyright (C) 2006-2017 Dynare Team
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

if (nargin < 2)
    error('gamrnd:: Two input arguments are needed!');
end

if nargin==2
    method= 'BauwensLubranoRichard';
    if any(a<1)
        method = 'Devroye';
        Devroye.small = 'Best'; % 'Weibull' , 'Johnk' , 'Berman' , 'GS' , 'Best'
                                % REMARK: The first algorithm (Weibull) is producing too much extreme values.
    end
    if ~strcmpi(method,'BauwensLubranoRichard')
        Devroye.big = 'Best'; % 'Cheng' , 'Best'
                              % REMARK 1: The first algorithm (Cheng) is still producing obviously wrong simulations.
                              % REMARK 2: The second algorithm seems slightly slower than the algorithm advocated by Bauwens,
                              %           Lubrano and Richard, but the comparison depends on the value of a (this should be
                              %           investigated further).
    end
else
    error('gamrnd:: Selection of method not yet implemented')
end

[ma,na] = size(a);
[mb,nb] = size(b);

if ma~=mb || na~=nb
    error('gamrnd:: Input arguments must have the same size!');
end

if na~=1
    error('gamrnd:: Input arguments must be column vectors');
end

if (any(a<0)) || (any(b<0)) || (any(a==Inf)) || (any(b==Inf))
    error('gamrnd:: Input arguments must be finite and positive!');
end

[~,integer_idx,double_idx] = isint(a);

number_of_integer_a = length(integer_idx);
number_of_double_a = length(double_idx);

rnd = NaN(ma,1);

if number_of_integer_a
    small_idx = find(a(integer_idx)<30);
    big_idx = find(a(integer_idx)>=30);
    number_of_small_a = length(small_idx);
    number_of_big_a = length(big_idx);
    if number_of_small_a
        % Exact sampling.
        for i=1:number_of_small_a
            rnd(integer_idx(small_idx(i))) = sum(exprnd(ones(a(integer_idx(small_idx(i))),1)))*b(integer_idx(small_idx(i)));
        end
    end
    if number_of_big_a
        % Gaussian approximation.
        rnd(integer_idx(big_idx)) = sqrt(a(integer_idx(big_idx))).* b(integer_idx(big_idx)) .* randn(number_of_big_a, 1) + a(integer_idx(big_idx)) .* b(integer_idx(big_idx));
    end
end


if number_of_double_a
    if strcmpi(method,'BauwensLubranoRichard')
        % Algorithm given in Bauwens, Lubrano & Richard (1999) page 316.
        rnd(double_idx) = gamrnd.knuth(a(double_idx),b(double_idx));
    else% Algorithm given in  Devroye (1986, chapter 9)
        small_idx = find(a(double_idx)<1);
        big_idx = find(a(double_idx)>1);
        number_of_small_a = length(small_idx);
        number_of_big_a = length(big_idx);
        if number_of_small_a
            if strcmpi(Devroye.small,'Weibull')
                % Algorithm given in Devroye (1986, page 415) [Rejection from the Weibull density]
                rnd(double_idx(small_idx)) = gamrnd.weibull_rejection(a(double_idx(small_idx)),b(double_idx(small_idx)));
            elseif strcmpi(Devroye.small,'Johnk')
                % Algorithm given in Devroye (1986, page 418) [Johnk's gamma generator]
                rnd(double_idx(small_idx)) = gamrnd.johnk(a(double_idx(small_idx)),b(double_idx(small_idx)));
            elseif strcmpi(Devroye.small,'Berman')
                % Algorithm given in Devroye (1986, page 418) [Berman's gamma generator]
                rnd(double_idx(small_idx)) = gamrnd.berman(a(double_idx(small_idx)),b(double_idx(small_idx)));
            elseif strcmpi(Devroye.small,'GS')
                % Algorithm given in Devroye (1986, page 425) [Ahrens and Dieter, 1974]
                rnd(double_idx(small_idx)) = gamrnd.ahrens_dieter(a(double_idx(small_idx)),b(double_idx(small_idx)));
            elseif strcmpi(Devroye.small,'Best')
                % Algorithm given in Devroye (1986, page 426) [Best, 1983]
                rnd(double_idx(small_idx)) = gamrnd.best_1983(a(double_idx(small_idx)),b(double_idx(small_idx)));
            end
        end
        if number_of_big_a
            if strcmpi(Devroye.big,'Cheng')
                % Algorithm given in Devroye (1986, page 413) [Cheng's rejection algorithm GB]
                rnd(double_idx(big_idx)) = gamrnd.cheng(a(double_idx(big_idx)),b(double_idx(big_idx)));
            elseif strcmpi(Devroye.big,'Best')
                % Algorithm given in Devroye (1986, page 410) [Best's rejection algorithm XG]
                rnd(double_idx(big_idx)) = gamrnd.best_1978(a(double_idx(big_idx)),b(double_idx(big_idx)));
            end
        end
    end
end