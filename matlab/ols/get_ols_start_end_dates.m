function [fp, lp] = get_ols_start_end_dates(Y, lhssub, X, jsonmodel)
% function [fp, lp] = get_ols_start_end_dates(Y, lhssub, X, jsonmodel)

fp = max(Y.firstobservedperiod, X.firstobservedperiod);
lp = min(Y.lastobservedperiod, X.lastobservedperiod);
if ~isempty(lhssub)
    fp = max(fp, lhssub.firstobservedperiod);
    lp = min(lp, lhssub.lastobservedperiod);
end
if isfield(jsonmodel, 'tags') ...
        && isfield(jsonmodel.tags, 'sample') ...
        && ~isempty(jsonmodel.tags.sample)
    colon_idx = strfind(jsonmodel.tags.sample, ':');
    fsd = dates(jsonmodel.tags.sample(1:colon_idx-1));
    lsd = dates(jsonmodel.tags.sample(colon_idx+1:end));
    if fp > fsd
        warning(['The sample over which you want to estimate contains NaNs. '...
            'Adjusting estimation range to begin on: ' fp.char])
    else
        fp = fsd;
    end
    if lp < lsd
        warning(['The sample over which you want to estimate contains NaNs. '...
            'Adjusting estimation range to end on: ' lp.char])
    else
        lp = lsd;
    end
end

end
