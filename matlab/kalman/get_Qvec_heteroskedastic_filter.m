function Qvec=get_Qvec_heteroskedastic_filter(Q,smpl,Model)
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