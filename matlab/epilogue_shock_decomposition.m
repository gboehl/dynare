function [zout, zss] = epilogue_shock_decomposition(zout ,M_, oo_)

ivar = varlist_indices(M_.epilogue_var_list_,M_.endo_names);
y = nan(length(M_.epilogue_names),size(zout,2),size(zout,3));
z=zout(ivar,:,:);
if isfield(oo_.dr,'ys')
    zss = oo_.dr.ys;
else
    zss = oo_.steady_state;
end

ztmp = dseries(zss(ivar),'',M_.epilogue_var_list_);
fname = M_.fname;
h_dynamic = str2func([fname '.epilogue_dynamic']);
h_static = str2func([fname '.epilogue_static']);
yss = extract(h_static(M_.params, ztmp),M_.epilogue_names{:});
yss = yss.data;
yss1 = repmat(yss(:),[1,1,size(y,3)]);
for k=1:size(z,2)
    ztmp = dseries(squeeze(z(:,k,:))+zss(ivar),'',M_.epilogue_var_list_);
    tmp = extract(h_dynamic(M_.params, ztmp),M_.epilogue_names{:});
    y(:,k,:) = tmp.data';
    y(:,k,:) = y(:,k,:)-yss1;
end

nterms = size(z,2);
for k=1:size(y,1)
    ytmp  = squeeze(y(k,:,:));
    yres = ytmp(end,:) - sum(ytmp(1:end-1,:));
    if ~isoctave && matlab_ver_less_than('9.1') % Automatic broadcasting was introduced in MATLAB R2016b
        w = bsxfun(@rdivide, abs(ytmp(1:end-1,:)), sum(abs(ytmp(1:end-1,:))))
    else
        w = abs(ytmp(1:end-1,:))./sum(abs(ytmp(1:end-1,:)));
    end
    %     ytmp(1:end-1,:) = ytmp(1:end-1,:) + repmat(yres,[nterms-1 1])/(nterms-1);
    ytmp(1:end-1,:) = ytmp(1:end-1,:) + repmat(yres,[nterms-1 1]).*w;
    y(k,:,:) = ytmp;
end

zout = cat(1,zout,y);
zss = [zss; yss(:)];
