function options_=set_dynare_seed_local_options(options_,a,b)
% options_=set_dynare_seed_local_options(options_,a,b)
% Set seeds depending on Matlab (octave) version
% Inputs:
%   o options_              options structure
%   o a                     first input argument, 
%                           for single argument input, either
%                               [number]       seed
%                               [string]        'default'
%                                               'reset'
%       	                                    'clock'
%                           for two inputs: 
%                               [string]        indicating feasible RNG algorithm (default: 'mt19937ar')
%   o b                     second input argument
%                               [number]        seed for algorithm
%                                               specified with first input

% Copyright © 2010-2023 Dynare Team
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

if nargin<2
    error('set_dynare_seed:: I need at least two input arguments!')
end

matlab_random_streams = ~(isoctave || options_.parallel_info.isHybridMatlabOctave);

if matlab_random_streams% Use new matlab interface.
    if nargin==2
        if ischar(a) && strcmpi(a,'default')
            options_.DynareRandomStreams.algo = 'mt19937ar';
            options_.DynareRandomStreams.seed = 0;
            s = RandStream(options_.DynareRandomStreams.algo,'Seed',options_.DynareRandomStreams.seed);
            reset(RandStream.setGlobalStream(s));
            return
        end
        if ischar(a) && strcmpi(a,'reset')
            s = RandStream(options_.DynareRandomStreams.algo,'Seed',options_.DynareRandomStreams.seed);
            reset(RandStream.setGlobalStream(s));
            return
        end
        if ~ischar(a) || (ischar(a) && strcmpi(a, 'clock'))
            options_.DynareRandomStreams.algo = 'mt19937ar';
            if ischar(a)
                options_.DynareRandomStreams.seed = rem(floor(now*24*60*60), 2^32);
            else
                options_.DynareRandomStreams.seed = a;
            end
            s = RandStream(options_.DynareRandomStreams.algo,'Seed',options_.DynareRandomStreams.seed);
            reset(RandStream.setGlobalStream(s));
            return
        end
        error('set_dynare_seed:: something is wrong in the calling sequence!')
    elseif nargin==3
        if ~ischar(a) || ~( strcmpi(a,'mcg16807') || ...
                            strcmpi(a,'mlfg6331_64') || ...
                            strcmpi(a,'mrg32k3a') || ...
                            strcmpi(a,'mt19937ar') || ...
                            strcmpi(a,'shr3cong') || ...
                            strcmpi(a,'swb2712') )
            disp('set_dynare_seed:: First argument must be string designing the uniform random number algorithm!')
            RandStream.list
            skipline()
            disp('set_dynare_seed:: Change the first input accordingly...')
            skipline()
            error(' ')
        end
        if ~isint(b)
            error('set_dynare_seed:: The second input argument must be an integer!')
        end
        options_.DynareRandomStreams.algo = a;
        options_.DynareRandomStreams.seed = b;
        s = RandStream(options_.DynareRandomStreams.algo,'Seed',options_.DynareRandomStreams.seed);
        reset(RandStream.setGlobalStream(s));
    end
else% Use old matlab interface.
    if nargin==2
        if ischar(a) && strcmpi(a,'default')
            if isoctave
                options_.DynareRandomStreams.algo = 'state';
            else
                options_.DynareRandomStreams.algo = 'twister';
            end
            options_.DynareRandomStreams.seed = 0;
            rand(options_.DynareRandomStreams.algo,options_.DynareRandomStreams.seed);
            randn('state',options_.DynareRandomStreams.seed);
            return
        end
        if ischar(a) && strcmpi(a,'reset')
            rand(options_.DynareRandomStreams.algo,options_.DynareRandomStreams.seed);
            randn('state',options_.DynareRandomStreams.seed);
            return
        end
        if (~ischar(a) && isint(a)) || (ischar(a) && strcmpi(a,'clock'))
            if ischar(a)
                options_.DynareRandomStreams.seed = floor(now*24*60*60);
            else
                options_.DynareRandomStreams.seed = a;
            end
            rand(options_.DynareRandomStreams.algo,options_.DynareRandomStreams.seed);
            randn('state',options_.DynareRandomStreams.seed);
            return
        end
        error('set_dynare_seed:: Something is wrong in the calling sequence!')
    else
        error('set_dynare_seed:: Cannot use more than one input argument with your version of Matlab/Octave!')
    end
end
