function [nCPU]= GiveCPUnumber (ComputerInformations, Environment)
% PARALLEL CONTEXT
% In a parallel context this function return the CPUs or cores numer avaiable
% on the computer used for run parallel code.
%
% INPUTS
% an array contained several fields that describe the hardaware
% software enviroments of a generic computer.
%
% OUTPUTS
% The CPUs or Cores numbers of computer.
%
% SPECIAL REQUIREMENTS
% none

% Copyright © 2010-2017 Dynare Team
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


nCPU='';

if nargin < 2
    % Determine a specific operating system or software version when necessary
    % for different command (sintax, name, ...).
    Environment=~ispc;
end

switch Environment
  case 0          %WINDOWS OPERATING SYSTEM

    OffSet=27;

    SringPosition=strfind(ComputerInformations, 'Processors:');
    nCPU=ComputerInformations(SringPosition+OffSet);

    % We check if there are Processors/Cores more than 9.


    t0=ComputerInformations(SringPosition+OffSet+1);
    t1=str2num(t0);
    t1=isempty(t1);

    % if t1 is 0 the machine have more than 9 CPU.

    if t1==0
        nCPU=strcat(nCPU,t0);
    end

    nCPU=str2num(nCPU);

    return
  case 1            %LIKE UNIX OPERATING SYSTEM

    nCPU=str2num(ComputerInformations);

  case 2            %MAC-OS OPERATING SYSTEM

    nCPU=str2num(ComputerInformations);

end
