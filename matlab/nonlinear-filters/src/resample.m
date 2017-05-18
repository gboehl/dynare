function resampled_output = resample(particles,weights,ParticleOptions)
% Resamples particles.
% if particles = 0, returns the resampling index (except for smooth resampling)
% Otherwise, returns the resampled particles set.

%@info:
%! @deftypefn {Function File} {@var{indx} =} resample (@var{weights}, @var{method})
%! @anchor{particle/resample}
%! @sp 1
%! Resamples particles.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item weights
%! n*1 vector of doubles, particles' weights.
%! @item method
%! string equal to 'residual' or 'traditional'.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item indx
%! n*1 vector of intergers, indices.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{particle/sequantial_importance_particle_filter}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{residual_resampling}, @ref{traditional_resampling}
%! @sp 2
%! @end deftypefn
%@eod:

% Copyright (C) 2011-2014 Dynare Team
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

defaultmethod = 1; % For residual based method set this variable equal to 0.

if defaultmethod
    if ParticleOptions.resampling.method.kitagawa
        resampled_output = traditional_resampling(particles,weights,rand);
    elseif ParticleOptions.resampling.method.stratified
        resampled_output = traditional_resampling(particles,weights,rand(size(weights)));
    elseif ParticleOptions.resampling.method.smooth
        if particles==0
            error('Particle = 0 is incompatible with this resampling method!')
        end
        resampled_output = multivariate_smooth_resampling(particles,weights);
    else
        error('Unknown sampling method!')
    end
else
    if ParticleOptions.resampling.method.kitagawa
        resampled_output = residual_resampling(particles,weights,rand);
    elseif ParticleOptions.resampling.method.stratified
        resampled_output = residual_resampling(particles,weights,rand(size(weights)));
    else
        error('Unknown sampling method!')
    end
end
