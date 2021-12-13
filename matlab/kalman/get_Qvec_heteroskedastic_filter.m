function Qvec=get_Qvec_heteroskedastic_filter(Q,smpl,Model)
% function Qvec=get_Qvec_heteroskedastic_filter(Q,smpl,Model)
% 
% INPUTS
%   Q:      baseline non-heteroskadastic covariance matrix of shocks
%   smpl:   scalar storing end of sample
%   Model:  structure storing the model information
% Outputs:
%   Qvec:   [n_exo by n_exo by smpl] array of covariance matrices

% Copyright (C) 2020-21 Dynare Team
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

isqdiag = all(all(abs(Q-diag(diag(Q)))<1e-14)); % ie, the covariance matrix is diagonal...
Qvec=repmat(Q,[1 1 smpl+1]);
for k=1:smpl
    inx = ~isnan(Model.heteroskedastic_shocks.Qvalue(:,k));
    if any(inx)
        if isqdiag
            Qvec(inx,inx,k)=diag(Model.heteroskedastic_shocks.Qvalue(inx,k));
        else
            inx = find(inx);
            for s=1:length(inx)
                if Q(inx(s),inx(s))>1.e-14
                    tmpscale = sqrt(Model.heteroskedastic_shocks.Qvalue(inx(s),k)./Q(inx(s),inx(s)));
                    Qvec(inx(s),:,k) = Qvec(inx(s),:,k).*tmpscale;
                    Qvec(:,inx(s),k) = Qvec(:,inx(s),k).*tmpscale;
                else
                    Qvec(inx(s),inx(s),k)=Model.heteroskedastic_shocks.Qvalue(inx(s),k);
                end
            end
        end
    end
    inx = ~isnan(Model.heteroskedastic_shocks.Qscale(:,k));
    if any(inx)
        if isqdiag
            Qvec(inx,inx,k)=Qvec(inx,inx,k).*diag(Model.heteroskedastic_shocks.Qscale(inx,k));
        else
            inx = find(inx);
            for s=1:length(inx)
                tmpscale = sqrt(Model.heteroskedastic_shocks.Qscale(inx(s),k));
                Qvec(inx(s),:,k) = Qvec(inx(s),:,k).*tmpscale;
                Qvec(:,inx(s),k) = Qvec(:,inx(s),k).*tmpscale;
            end
        end
    end
end