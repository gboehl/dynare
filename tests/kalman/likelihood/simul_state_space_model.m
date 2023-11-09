function observed_data = simul_state_space_model(T,R,Q,mf,nobs,H)
    pp = length(mf);
    mm = length(T);
    rr = length(Q);

    upper_cholesky_Q = chol(Q);
    if nargin>5 && ~isempty(H)
        upper_cholesky_H = chol(H);
    end

    state_data = zeros(mm,1);

    if (nargin==5 || isempty(H))
        for t = 1:nobs
            state_data = T*state_data + R* upper_cholesky_Q * randn(rr,1);
            observed_data(:,t) = state_data(mf);
        end
    elseif (nargin==6 && ~isempty(H))
        for t = 1:nobs
            state_data = T*state_data + R* upper_cholesky_Q * randn(rr,1);
            observed_data(:,t) = state_data(mf);
        end
        observed_data = observed_data + upper_cholesky_H * randn(pp,nobs);
    else
        error('simul_state_space_model:: I don''t understand what you want!!!')
    end
