@#define NumberOfCountries = 2

var log_kappa;
@#for Country in 1:NumberOfCountries
    var mu@{Country}, logit_l@{Country}, k@{Country}, a@{Country};
@#endfor

parameters alpha beta varrho delta theta rhoA sigmaA;

alpha = 0.3;
beta = 0.99;
varrho = 1.5;
delta = 0.025;
theta = 0.9;
rhoA = 0.95;
sigmaA = 0.05;

@#for Country in 1:NumberOfCountries
    varexo epsilonA@{Country};
@#endfor

model;

    #kappa = exp( log_kappa );
    #LEAD_kappa = exp( log_kappa(+1) );

    #min_A = theta ^ ( 1/alpha );
    #mean_a = log( 1 - min_A );

    @#for Country in 1:NumberOfCountries

        #K@{Country}        = exp( k@{Country} );
        #LAG_K@{Country}    = exp( k@{Country}(-1) );
        #LEAD_A@{Country}   = min_A + exp( a@{Country}(+1) );
        #A@{Country}        = min_A + exp( a@{Country} );
        #LAG_A@{Country}    = min_A + exp( a@{Country}(-1) );
        #L@{Country}        = 1 / ( 1 + exp( -logit_l@{Country} ) );
        #LEAD_L@{Country}   = 1 / ( 1 + exp( -logit_l@{Country}(+1) ) );

        #C@{Country} = kappa ^ ( -1/varrho ) - theta / ( 1-alpha ) * ( A@{Country}*(1-L@{Country}) ) ^ ( 1-alpha );
        #LEAD_phi@{Country} = ( 1 - theta * ( LEAD_A@{Country}*(1-LEAD_L@{Country}) ) ^ ( -alpha ) ) * LEAD_kappa;
        #I@{Country} = K@{Country} - ( 1-delta ) * LAG_K@{Country};

    @#endfor

    @#for Country in 1:NumberOfCountries

        ( a@{Country} - mean_a ) = rhoA * ( a@{Country}(-1) - mean_a ) - sigmaA * epsilonA@{Country};
        L@{Country} = min( LAG_K@{Country} / A@{Country}, 1 - theta ^ ( 1/alpha ) / A@{Country} );
        kappa - mu@{Country} = beta * ( ( 1-delta ) * ( LEAD_kappa - mu@{Country}(+1) ) + LEAD_phi@{Country} );
        kappa = max( beta * ( ( 1-delta ) * ( LEAD_kappa - mu@{Country}(+1) ) + LEAD_phi@{Country} ), ( theta / ( 1-alpha ) * ( A@{Country}*(1-L@{Country}) ) ^ ( 1-alpha )
        @#for OtherCountry in 1:NumberOfCountries
            + A@{OtherCountry} * L@{OtherCountry}
            @#if OtherCountry != Country
                - C@{OtherCountry} - I@{OtherCountry}
            @#endif
        @#endfor
        ) ^ ( -varrho ) );

    @#endfor

    0 = 0
    @#for OtherCountry in 1:NumberOfCountries
        + A@{OtherCountry} * L@{OtherCountry}
        - C@{OtherCountry} - I@{OtherCountry}
    @#endfor
    ;

end;

steady_state_model;

    a = log( 1 - theta ^ ( 1/alpha ) );
    mu = 0;
    L = 1 - ( 1/theta * ( 2 - 1/beta - delta ) ) ^ ( -1/alpha );
    K = L;
    I = delta * K;
    C = ( 1 - delta ) * L;

    kappa_ = ( C + theta / ( 1-alpha ) * ( 1-L ) ^ ( 1-alpha ) ) ^ ( -varrho );

    log_kappa = log( kappa_ );

    @#for Country in 1:NumberOfCountries

        mu@{Country} = mu;
        logit_l@{Country} = log( L / ( 1 - L ) );
        k@{Country} = log( K );
        a@{Country} = a;

    @#endfor

end;

steady;
check;

shocks;

var epsilonA1; periods 1; values 2;
@#for Country in 2:NumberOfCountries
    var epsilonA@{Country}; periods 1; values 0;    
@#endfor

end;

perfect_foresight_setup(periods=400);
perfect_foresight_solver(robust_lin_solve);

if ~oo_.deterministic_simulation.status
    error('Model did not solve')
end
