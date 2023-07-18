function Draws = GetAllPosteriorDraws(column, FirstMhFile, FirstLine, TotalNumberOfMhFile, NumberOfDraws, nblcks, blck)
% function Draws = GetAllPosteriorDraws(column, FirstMhFile, FirstLine, TotalNumberOfMhFile, NumberOfDraws, nblcks, blck)
% Gets all posterior draws
%
% INPUTS
%    column:               column of desired parameter in draw matrix
%    FirstMhFile:          first mh file
%    FirstLine:            first line
%    TotalNumberOfMhFile:  total number of mh file
%    NumberOfDraws:        number of draws
%    nblcks:               total number of blocks 
%    blck:                 desired block to read
%
% OUTPUTS
%    Draws:                draws from posterior distribution
%
% SPECIAL REQUIREMENTS
%    none

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

global M_ 

iline = FirstLine;
linee = 1;
DirectoryName = CheckPath('metropolis',M_.dname);

if nblcks>1 && nargin<7
    Draws = zeros(NumberOfDraws*nblcks,1);
    iline0=iline;
    if column>0
        for blck = 1:nblcks
            iline=iline0;
            for file = FirstMhFile:TotalNumberOfMhFile
                load([DirectoryName '/'  M_.fname '_mh' int2str(file) '_blck' int2str(blck)],'x2')
                NumberOfLines = size(x2(iline:end,:),1);
                Draws(linee:linee+NumberOfLines-1) = x2(iline:end,column);
                linee = linee+NumberOfLines;
                iline = 1;
            end
        end
    else
        for blck = 1:nblcks
            iline=iline0;
            for file = FirstMhFile:TotalNumberOfMhFile
                load([DirectoryName '/'  M_.fname '_mh' int2str(file) '_blck' int2str(blck)],'logpo2')
                NumberOfLines = size(logpo2(iline:end),1);
                Draws(linee:linee+NumberOfLines-1) = logpo2(iline:end);
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
            load([DirectoryName '/'  M_.fname '_mh' int2str(file) '_blck' int2str(blck)],'x2')
            NumberOfLines = size(x2(iline:end,:),1);
            Draws(linee:linee+NumberOfLines-1) = x2(iline:end,column);
            linee = linee+NumberOfLines;
            iline = 1;
        end
    else
        for file = FirstMhFile:TotalNumberOfMhFile
            load([DirectoryName '/'  M_.fname '_mh' int2str(file) '_blck' int2str(blck)],'logpo2')
            NumberOfLines = size(logpo2(iline:end,:),1);
            Draws(linee:linee+NumberOfLines-1) = logpo2(iline:end);
            linee = linee+NumberOfLines;
            iline = 1;
        end
    end
end