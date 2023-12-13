function [info,nG,nKAP,nD,nGAM,nPHI] = cet_steady_helper(justCheck,Kw,J,V,U,  dolagZBAR,mu,thetaG,thetaKAP,thetaD,thetaGAM,thetaPHI)
% Based on cet_steadystate.m of the replication codes for
% Christiano, Eichenbaum, Trabandt (2016, Econometrica) - Unemployment and the Business Cycle;
% slightly modified such that most code is now in a steady_state_model block
% and only the two if statements are computed in this helper function.
% =========================================================================
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
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
% =========================================================================

info = 0;
if justCheck
    if (Kw<0)
        disp('could not compute steadystate: Kw<0');
    elseif (J<0)
        disp('could not compute steadystate: J<0');
    elseif (V-U<0)
        disp('could not compute steadystate: V-U<0')
        info=1;
    end
    return
else
    if dolagZBAR==1
        nG=mu^(-1/thetaG);
        nKAP=mu^(-1/thetaKAP);
        nD=mu^(-1/thetaD);
        nGAM=mu^(-1/thetaGAM);
        nPHI=mu^(-1/thetaPHI);
    else
        nG=(1/mu)^((1-thetaG)/thetaG);
        nKAP=(1/mu)^((1-thetaKAP)/thetaKAP);
        nD=(1/mu)^((1-thetaD)/thetaD);
        nGAM=(1/mu)^((1-thetaGAM)/thetaGAM);
        nPHI=(1/mu)^((1-thetaPHI)/thetaPHI);
    end
end