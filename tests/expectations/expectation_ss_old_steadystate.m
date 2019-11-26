function [ys_, params, check_] = expectation_ss_old_steadystate(ys_orig_, exo_,M_,options_)
    ys_=zeros(6,1);
    params=NaN(size(M_.params));
    params=M_.params;
    ys_(4)=0;
    ys_(6)=0;
    ys_(5)=0.3333333333333333;
    ys_(3)=((1/M_.params(1)-(1-M_.params(4)))/(M_.params(3)*ys_(5)^(1-M_.params(3))))^(1/(M_.params(3)-1));
    ys_(1)=ys_(5)^(1-M_.params(3))*ys_(3)^M_.params(3);
    ys_(2)=ys_(1)-M_.params(4)*ys_(3);
    params(5)=(1-M_.params(3))*ys_(1)/(ys_(2)*ys_(5)^(1+M_.params(6)));
    check_=0;
end
