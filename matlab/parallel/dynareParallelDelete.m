function dynareParallelDelete(fname,pname,Parallel)
% PARALLEL CONTEXT
% In a parallel context, this is a specialized version of delete() function.
%
% INPUTS
%  o fname      []   ...
%  o pname      []   ...
%  o Parallel   []   ...
%
%  OUTPUTS
%  None
%
%
% Copyright © 2009-2020 Dynare Team
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

if nargin ~= 3
    disp('dynareParallelDelete(fname,pname,Parallel)')
    return
end

if ~isempty(pname)
    pname=[pname,filesep];
end

for indPC=1:length(Parallel)
    if ~ispc || strcmpi('unix',Parallel(indPC).OperatingSystem)
        if ~isempty(Parallel(indPC).Port)
            ssh_token = ['-p ',Parallel(indPC).Port ' '];
        else
            ssh_token = ' ';
        end
        username = Parallel(indPC).UserName;
        if ~isempty(username)
            username = [username '@'];
        end
        directory = Parallel(indPC).RemoteDirectory;
        if ~isempty(directory)
            directory = [directory '/'];
        end
        system(['ssh ',ssh_token,username,Parallel(indPC).ComputerName,' ''/bin/bash --norc -c "rm -f ',directory,pname,fname,'"''']);
    else
        fname_temp=['\\',Parallel(indPC).ComputerName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteDirectory,'\',pname,fname];
        if exist(fname_temp,'file')
            delete(fname_temp);
        end
    end

end
