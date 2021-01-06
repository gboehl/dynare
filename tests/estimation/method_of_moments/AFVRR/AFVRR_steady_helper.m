% DSGE model based on replication files of
% Andreasen, Fernandez-Villaverde, Rubio-Ramirez (2018), The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications, Review of Economic Studies, 85, p. 1-49
% Adapted for Dynare by Willi Mutschler (@wmutschl, willi@mutschler.eu), Jan 2021
% =========================================================================
% Copyright (C) 2021 Dynare Team
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
% =========================================================================

% This is a helper function to compute steady state values and endogenous parameters
% Based on DSGE_model_yieldCurve_ss.m, getPHI3.m, ObjectGMM.m

function [AA, EVFBAR, PHI3, negVf, info]= AFVRR_steady_helper(VFBAR,RBAR,IVBAR,CBAR,KBAR,LABAR,QBAR,YBAR,  BETTA,B,PAI,H,PHIzero,PHI1,PHI2,THETA,MYYPS,MYZ,INHABIT,RRA,CONSxhr40)         
% We get nice values of EVF by setting AA app. equal to VF. 
% The value of the expected value function raised to the power 1-PHI3
% Also we check bounds on other variables
% % Adding PHI3 to params. Note that PHI3 only affects the value function in
% % steady state, hence the value we assign to PHI3 is irrelevant
%     PHI3 = -100;

info=0;
AA = NaN;
EVFBAR = NaN;
PHI3 = NaN;
negVf = NaN;

MYZSTAR = MYYPS^(THETA/(1-THETA))*MYZ;
% The wage level
WBAR      = PHIzero*(1-H)^(-PHI1)/LABAR;
RRAc   = RRA;
if INHABIT == 1
    PHI3 = (RRAc - PHI2/((1-B*MYZSTAR^-1)/(1-BETTA*B)+PHI2/PHI1*WBAR*(1-H)/CBAR))/((1-PHI2)/((1-B*MYZSTAR^-1)/(1-BETTA*B)-(CBAR-B*CBAR*MYZSTAR^-1)^PHI2/((1-BETTA*B)*CBAR)+WBAR*(1-H)/CBAR*(1-PHI2)/(1-PHI1)));
else    
    PHI3 = (RRAc - PHI2/(1-B*MYZSTAR^-1+PHI2/PHI1*WBAR*(1-H)/CBAR))/((1-PHI2)/(1-B*MYZSTAR^-1-(CBAR-B*CBAR*MYZSTAR^-1)^PHI2/((1-BETTA*B)*CBAR)+WBAR*(1-H)/CBAR*(1-PHI2)/(1-PHI1)));
end
if abs(PHI3) > 30000
    disp('abs of PHI3 exceeds 30000')
    info=1;
    return
end

if CONSxhr40 > 1
   info=1;
   return
end


if VFBAR < 0   
    AA        = -VFBAR; 
    EVFBAR    = (-VFBAR/AA)^(1-PHI3);
    negVf      = 1;    
else
    AA        = VFBAR;
    EVFBAR    = (VFBAR/AA)^(1-PHI3);
    negVf      = -1;
    disp('Positive Value Function');
end


if RBAR < 1 || IVBAR < 0 || CBAR < 0 || KBAR < 0 || PAI < 1 || H < 0 || H > 1 || QBAR < 0 || YBAR < 0
   info = 1;
end

end



