function indx = smc_resampling(weights,noise,number)
% function indx = smc_resampling(weights,noise,number)

% Copyright © 2022 Dynare Team
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

    indx = zeros(number,1);
    cumweights = cumsum(weights);
    randvec = (transpose(1:number)-1+noise(:))/number;
    j = 1;
    for i=1:number
        while (randvec(i)>cumweights(j))
            j = j+1;
        end
        indx(i) = j;
    end
