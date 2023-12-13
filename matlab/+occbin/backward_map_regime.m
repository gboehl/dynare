function [binding_indicator, A, regime_string] = backward_map_regime(regime, regime_start)
% [binding_indicator, A, regime_string] = backward_map_regime(regime, regime_start)
% Map regime information into regime indicator
%
% Inputs:
% - regime              [integer]   [1 by n_transitions] vector of regime number indices
% - regime_start        [integer]   [1 by n_transitions] vectors with period numbers in which regime starts
%
% Outputs:
% - binding_indicator   [integer]   [nperiods by 1] vector of regime indices
% - A                   [bin]       binary representation of binding indicator
% - error_flag          [boolean]   1 if regime never leaves 1 or is still there at the end of nperiods
%                                   0 otherwise

% Copyright © 2023 Dynare Team
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

if nargin ==1
% polymorphism
    if ~isstruct(regime)
        disp('error::backward_map_regime')
        disp('input arguments may be 1 structure with regime info')
        disp('or two arrays: regime and regimestart ')
        error('wrong input')
    end

    fnam = fieldnames(regime);
    if length(fnam) == 2
        [binding_indicator, A, regime_string] = occbin.backward_map_regime(regime.regime, regime.regimestart);
    else
        for k=1:2
            nperiods(k) = regime.(['regimestart' int2str(k)])(end);
            number_of_binary_tokens(k) = ceil((nperiods(k)-1)/50);
        end
        binding_indicator = false(max(nperiods),2);
        A = int64(zeros(max(number_of_binary_tokens),2));

        for k=1:2
            [binding_indicator(1:nperiods(k),k), A(1:number_of_binary_tokens(k),k), tmp{k}] = ...
                occbin.backward_map_regime(regime.(['regime' int2str(k)]), regime.(['regimestart' int2str(k)]));
        end
        regime_string = char(tmp{1},tmp{2});
    end

    return
else
    if isstruct(regime)
        disp('error::backward_map_regime')
        disp('input arguments may be ONE structure with regime info')
        disp('or TWO arrays: regime and regimestart ')
        error('wrong input')
    end
end

regime_string = char(mat2str(double(regime)),mat2str(regime_start));

nperiods_0 = regime_start(end);
number_of_binary_tokens = max(1,ceil((nperiods_0-1)/50));
A = int64(zeros(number_of_binary_tokens,1));
binding_indicator = false(nperiods_0,1);
if length(regime)>1
    for ir=1:length(regime)-1
        binding_indicator(regime_start(ir):regime_start(ir+1)-1,1) = regime(ir);
        for k=regime_start(ir):regime_start(ir+1)-1
            this_token = ceil(k/50);
            A(this_token) = int64(bitset(A(this_token),k-50*(this_token-1),regime(ir)));
        end
    end
end

binding_indicator = logical(binding_indicator);
% to convert regime in a readable string array
% a = dec2bin(A);
% bindicator = [a(end:-1:1) '0'];