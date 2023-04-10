function draws = GetAllPosteriorDraws(options_, dname, fname, column, FirstMhFile, FirstLine, TotalNumberOfMhFile, NumberOfDraws, nblcks, blck)

% Gets all posterior draws.
%
% INPUTS
% - options_               [struct]   Dynare's options.
% - dname                  [char]     name of directory with results.
% - fname                  [char]     name of mod file.
% - column                 [integer]  scalar, parameter index.
% - FirstMhFile            [integer]  scalar, first MH file.
% - FirstLine              [integer]  scalar, first line in first MH file.
% - TotalNumberOfMhFile    [integer]  scalar, total number of MH file.
% - NumberOfDraws          [integer]  scalar, number of posterior draws.
% - nblcks                 [integer]  scalar, total number of blocks.
% - blck:                  [integer]  scalar, desired block to read.
%
% OUTPUTS
% - draws:                 [double]   NumberOfDraws×1 vector, draws from posterior distribution.
%
% REMARKS
% Only the first and third input arguments are required for SMC samplers.

% Copyright © 2005-2023 Dynare Team
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

if ishssmc(options_)
    % Load draws from the posterior distribution
    pfiles = dir(sprintf('%s/hssmc/particles-*.mat', dname));
    posterior = load(sprintf('%s/hssmc/particles-%u-%u.mat', dname, length(pfiles), length(pfiles)));
    draws = transpose(posterior.particles(column,:));
else
    iline = FirstLine;
    linee = 1;
    DirectoryName = CheckPath('metropolis',dname);
    if nblcks>1 && nargin<10
        draws = zeros(NumberOfDraws*nblcks,1);
        iline0=iline;
        if column>0
            for blck = 1:nblcks
                iline=iline0;
                for file = FirstMhFile:TotalNumberOfMhFile
                    load([DirectoryName '/'  fname '_mh' int2str(file) '_blck' int2str(blck)],'x2')
                    NumberOfLines = size(x2(iline:end,:),1);
                    draws(linee:linee+NumberOfLines-1) = x2(iline:end,column);
                    linee = linee+NumberOfLines;
                    iline = 1;
                end
            end
        else
            for blck = 1:nblcks
                iline=iline0;
                for file = FirstMhFile:TotalNumberOfMhFile
                    load([DirectoryName '/'  fname '_mh' int2str(file) '_blck' int2str(blck)],'logpo2')
                    NumberOfLines = size(logpo2(iline:end),1);
                    draws(linee:linee+NumberOfLines-1) = logpo2(iline:end);
                    linee = linee+NumberOfLines;
                    iline = 1;
                end
            end
        end
    else
        if nblcks==1
            blck=1;
        end
        if column>0
            for file = FirstMhFile:TotalNumberOfMhFile
                load([DirectoryName '/'  fname '_mh' int2str(file) '_blck' int2str(blck)],'x2')
                NumberOfLines = size(x2(iline:end,:),1);
                draws(linee:linee+NumberOfLines-1) = x2(iline:end,column);
                linee = linee+NumberOfLines;
                iline = 1;
            end
        else
            for file = FirstMhFile:TotalNumberOfMhFile
                load([DirectoryName '/'  fname '_mh' int2str(file) '_blck' int2str(blck)],'logpo2')
                NumberOfLines = size(logpo2(iline:end,:),1);
                draws(linee:linee+NumberOfLines-1) = logpo2(iline:end);
                linee = linee+NumberOfLines;
                iline = 1;
            end
        end
    end
end
