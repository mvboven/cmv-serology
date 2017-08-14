data {
    int<lower=0> N;                                     // number of subjects
    int<lower=1> DeltaA;                                // length of age intervals
    int<lower=1> A;                                     // number of age intervals
    int<lower=1> A_rho;                                 // (reduced) number of age intervals for reactivation
    int<lower=1, upper=A_rho> RhoClasses[A];
    matrix<lower=0>[A, A] Contact_MM;                   // gender- and age-specific contact matrices
    matrix<lower=0>[A, A] Contact_FM;
    matrix<lower=0>[A, A] Contact_MF;
    matrix<lower=0>[A, A] Contact_FF;
    int<lower=0, upper=DeltaA*A> Ages[N];               // ages of the subjects
    real Titers[N];                                     // antibody titers
    int Censor[N];                                      // 0 = normal, 1 = censored, 2 = spike
    real RightCensor;                                   // titers above this value are censored
    real MuS;                                           // means of classification mixture
    real MuL;
    real MuB;
    real<lower=0> SigmaS;                               // standard deviations of the classification mixture
    real<lower=0> SigmaL;
    real<lower=0> SigmaB;
    real<lower=0> Penalty;                              // accuracy for estimating the forces of infection
    real<lower=0, upper=1> S0;                          // fraction susceptibles at age 0
    int<lower=0, upper=1> Gender[N];                    // 0 = female, 1 = male
    int<lower=0, upper=1> mode;                         // 0 = regular sampling, 1 = sampling to compute WBIC
}

transformed data {
    real<lower=0, upper=1> watanabe_beta;               // WBIC parameter
    // some auxiliary "data" used in vector manipulations
    vector[1] Zero;
    vector[DeltaA*A] LongOnes;
    vector[DeltaA*A+1] LlongOnes;
    // indices for transforming lambda into longLambda (and back)
    int<lower=1, upper=A> ExtendIdxs[DeltaA*A];
    int<lower=1, upper=DeltaA*A+1> ReduceIdxs[A];
    int<lower=1, upper=DeltaA*A+1> ReduceIdxsRightShift[A];

    // use mode to determine the "temperature"
    if ( mode == 0 ) {
        watanabe_beta = 1.0;
    }
    else { // mode == 1
        watanabe_beta = 1.0/log(N);
    }

    Zero = rep_vector(0.0, 1);                          // a vector-version of 0 (for cumulative sums starting at 0)
    LongOnes = rep_vector(1.0, DeltaA*A);               // DeltaA*A ones (same length as longLambda and longPi)
    LlongOnes = rep_vector(1.0, DeltaA*A+1);            // DeltaA*A+1 ones (same length as S, L, and B)

    for ( aa in 1:DeltaA*A ) {                          // aa-1 = 0,1,2,3,4, 5,6,7,8,9, 10,...
        int j;
        j = (aa-1)/DeltaA + 1;                          // j = 1,1,1,1,1, 2,2,2,2,2, 3,...
        ExtendIdxs[aa] = j;
    }
    for ( j in 1:A ) {                                  // j = 1,2,3,...
        ReduceIdxs[j] = 1 + (j-1)*DeltaA;               // x-1 = 0,5,10,...,DeltaA*(A-1)
    }
    for ( j in 1:A ) {                                  // j = 1,2,3,...
        ReduceIdxsRightShift[j] = 1 + j*DeltaA;         // x-1 = 5,10,15,...,DeltaA*A
    }
}

parameters {
    vector<lower=0>[A] lambda_f;                        // force of infection for the age-intervals (female)
    vector<lower=0>[A] lambda_m;                        // force of infection for the age-intervals (male)
    real<lower=0> beta1;                                // infection rate parameter after primary infection
    real<lower=0> beta2;                                // infection rate parameter after reactivation/re-infection
    real<lower=0, upper=1> z;                           // reduction in susceptibility for reinfection
	vector<lower=0>[A_rho] shortRho_f;                  // reactivation rate (female)
    vector<lower=0>[A_rho] shortRho_m;                  // reactivation rate (male)
	real<lower=0> muRho;                                // hyperprior for rho
    real<lower=0> sigmaRho;
}

transformed parameters {
    vector<lower=0>[A] rho_f;
    vector<lower=0>[A] rho_m;
    vector<lower=0>[A] pi_f;                            // shorthand: pi = rho + z * lambda
    vector<lower=0>[A] pi_m;
    vector<lower=0, upper=1>[DeltaA*A+1] S_f;           // susceptible (female)
    vector<lower=0, upper=1>[DeltaA*A+1] S_m;           // susceptible (male)
    vector<lower=0, upper=1>[DeltaA*A+1] L_f;           // latently infected (female)
    vector<lower=0, upper=1>[DeltaA*A+1] L_m;           // latently infected (male)
    vector<lower=0, upper=1>[DeltaA*A+1] B_f;           // infected with boosted titers (female)
    vector<lower=0, upper=1>[DeltaA*A+1] B_m;           // infected with boosted titers (male)
    // auxiliary vectors for calculating S, L, and B
    vector<lower=0, upper=1>[DeltaA*A+1] X_f;           // latently infected at birth (female)
    vector<lower=0, upper=1>[DeltaA*A+1] X_m;           // latently infected at birth (male)
    vector<lower=0>[DeltaA*A+1] Y_f;                    // (L_f-X_f)/X_f
    vector<lower=0>[DeltaA*A+1] Y_m;                    // (L_m-X_m)/X_m
    // lambda hat should be very similar to lambda
    vector<lower=0>[A] lambda_hat_f;
    vector<lower=0>[A] lambda_hat_m;
    // we need long versions of lambda and pi
    vector<lower=0>[DeltaA*A] longLambda_f;
    vector<lower=0>[DeltaA*A] longLambda_m;
    vector<lower=0>[DeltaA*A] longPi_f;
    vector<lower=0>[DeltaA*A] longPi_m;
    
    // UNCOMMENT FOR THE MODEL WITH REACTIVATION FROM THE B-COMPARTMENT
    /*
    vector<lower=0>[A] aggr_B_rho_f;
    vector<lower=0>[A] aggr_B_rho_m;
    */
    // for notational convenience, transform shortRho into rho 
    rho_f = shortRho_f[RhoClasses];
    rho_m = shortRho_m[RhoClasses];
        
    // for notational convenience, define pi = rho + z * lambda
    pi_f = rho_f + z * lambda_f;
    pi_m = rho_m + z * lambda_m;

    // make long versions of lambda and pi
    longLambda_f = lambda_f[ExtendIdxs];
    longLambda_m = lambda_m[ExtendIdxs];
    longPi_f = pi_f[ExtendIdxs];
    longPi_m = pi_m[ExtendIdxs];
    
    // define S, L, B, X, and Y
    S_f = S0 * exp(-cumulative_sum(append_row(Zero, longLambda_f)));
    S_m = S0 * exp(-cumulative_sum(append_row(Zero, longLambda_m)));

    X_f = (1.0 - S0) * exp(-cumulative_sum(append_row(Zero, longPi_f)));
    X_m = (1.0 - S0) * exp(-cumulative_sum(append_row(Zero, longPi_m)));
    
    Y_f = cumulative_sum(append_row(Zero, longLambda_f .* (S_f[:DeltaA*A] ./ X_f[:DeltaA*A]) .* (LongOnes - exp(-(longLambda_f - longPi_f))) ./ (longLambda_f - longPi_f)));
    Y_m = cumulative_sum(append_row(Zero, longLambda_m .* (S_m[:DeltaA*A] ./ X_m[:DeltaA*A]) .* (LongOnes - exp(-(longLambda_m - longPi_m))) ./ (longLambda_m - longPi_m)));

    L_f = X_f .* (Y_f + LlongOnes);
    L_m = X_m .* (Y_m + LlongOnes);

    B_f = LlongOnes - S_f - L_f;
    B_m = LlongOnes - S_m - L_m;
    
    // define lambda_hat for full model
	
    lambda_hat_f = Contact_FF * (beta1 * (S_f[ReduceIdxs] - S_f[ReduceIdxsRightShift]) + beta2 * (B_f[ReduceIdxsRightShift] - B_f[ReduceIdxs]))
                 + Contact_FM * (beta1 * (S_m[ReduceIdxs] - S_m[ReduceIdxsRightShift]) + beta2 * (B_m[ReduceIdxsRightShift] - B_m[ReduceIdxs]));
    lambda_hat_m = Contact_MM * (beta1 * (S_m[ReduceIdxs] - S_m[ReduceIdxsRightShift]) + beta2 * (B_m[ReduceIdxsRightShift] - B_m[ReduceIdxs]))
                 + Contact_MF * (beta1 * (S_f[ReduceIdxs] - S_f[ReduceIdxsRightShift]) + beta2 * (B_f[ReduceIdxsRightShift] - B_f[ReduceIdxs]));
	
    // UNCOMMENT FOR THE MODEL WITH REACTIVATION FROM THE B-COMPARTMENT
    /*
    for ( a in 1:A ) { // TODO: find a quick vectorized method...
        aggr_B_rho_f[a] = rho_f[a] * sum(B_f[1+DeltaA*(a-1):DeltaA*a]);
        aggr_B_rho_m[a] = rho_m[a] * sum(B_m[1+DeltaA*(a-1):DeltaA*a]);
    }
                 
    // define lambda_hat for model with reactivation in B compartment. TODO: this is an approximation: \int B(a) rho(a) da \approx \sum_{a} B(a) rho(a)
    lambda_hat_f = Contact_FF * (beta1 * (S_f[ReduceIdxs] - S_f[ReduceIdxsRightShift]) + beta2 * (B_f[ReduceIdxsRightShift] - B_f[ReduceIdxs] + aggr_B_rho_f))
                 + Contact_FM * (beta1 * (S_m[ReduceIdxs] - S_m[ReduceIdxsRightShift]) + beta2 * (B_m[ReduceIdxsRightShift] - B_m[ReduceIdxs] + aggr_B_rho_m));
    lambda_hat_m = Contact_MM * (beta1 * (S_m[ReduceIdxs] - S_m[ReduceIdxsRightShift]) + beta2 * (B_m[ReduceIdxsRightShift] - B_m[ReduceIdxs] + aggr_B_rho_m))
                 + Contact_MF * (beta1 * (S_f[ReduceIdxs] - S_f[ReduceIdxsRightShift]) + beta2 * (B_f[ReduceIdxsRightShift] - B_f[ReduceIdxs] + aggr_B_rho_f));
    */
}

model {    
    // priors on the (hyper)parameters
    beta1 ~ normal(0, 0.1);
    beta2 ~ normal(0, 0.1);
    z ~ uniform(0, 1);
    
    muRho ~ normal(0, 0.1);
    sigmaRho ~ normal(0, 0.1);
    
    for ( a in 1:A_rho ) {
        shortRho_f[a] ~ normal(muRho, sigmaRho);
        shortRho_m[a] ~ normal(muRho, sigmaRho);
	}
  	
    // penalise the difference between lambda and lambda_hat
    // TODO: lambda_f - lambda_hat_f ~ normal(0,1/Penalty) gives a warning: do target += log(det(Jacobian))
    lambda_f ~ normal(lambda_hat_f, 1/Penalty);
    lambda_m ~ normal(lambda_hat_m, 1/Penalty);
        
    // likelihood of the data
    for ( i in 1:N ) {
        int aa;
        real pS; real pL; real pB;
        
        // make the code a little bit more readable
		aa = Ages[i] + 1; // the index for S, L and B
        
        // compute the compartment-probabilities given the subjects' age
        if ( Gender[i] == 0 ) { // 0 = female
           pS = S_f[aa];
           pL = L_f[aa];
           pB = B_f[aa];
        }
        else { // 1 = male
           pS = S_m[aa];
           pL = L_m[aa];
           pB = B_m[aa];
        }
        
        // define the likelihood
        if ( Censor[i] == 0 ) { // normal data
            target += watanabe_beta * log( pS * exp(normal_lpdf(Titers[i] | MuS, SigmaS)) +
                                           pL * exp(normal_lpdf(Titers[i] | MuL, SigmaL)) +
                                           pB * exp(normal_lpdf(Titers[i] | MuB, SigmaB)) );
        }
        else if ( Censor[i] == 1 ) { // right censored
            target += watanabe_beta * log( pS * exp(normal_lccdf(RightCensor | MuS, SigmaS)) +
                                           pL * exp(normal_lccdf(RightCensor | MuL, SigmaL)) +
                                           pB * exp(normal_lccdf(RightCensor | MuB, SigmaB)) );
        }
        else if ( Censor[i] == 2 ) { // spiked
            target += watanabe_beta * log(pS);
        }
    }
}

generated quantities {
    // save the likelihood 
    vector[N] log_lik;
    real log_like;

    // likelihood of the data
    for ( i in 1:N ) {
		int aa;
        real pS; real pL; real pB;

        // make the code a little bit more readable
		aa = Ages[i] + 1; // the index for S, L and B
        
        // compute the compartment-probabilities given the subjects' age
        if ( Gender[i] == 0 ) { // 0 = female
            pS = S_f[aa];
            pL = L_f[aa];
            pB = B_f[aa];
        }
        else{ // 1 = male
            pS = S_m[aa];
            pL = L_m[aa];
            pB = B_m[aa];
        }
        
        // define the likelihood
        if ( Censor[i] == 0 ) { // normal data
            log_lik[i] = log( pS * exp(normal_lpdf(Titers[i] | MuS, SigmaS)) +
                              pL * exp(normal_lpdf(Titers[i] | MuL, SigmaL)) +
                              pB * exp(normal_lpdf(Titers[i] | MuB, SigmaB)) );
        }
        else if ( Censor[i] == 1 ) { // right censored
            log_lik[i] = log( pS * exp(normal_lccdf(RightCensor | MuS, SigmaS)) +
                              pL * exp(normal_lccdf(RightCensor | MuL, SigmaL)) +
                              pB * exp(normal_lccdf(RightCensor | MuB, SigmaB)) );
        }
        else if ( Censor[i] == 2 ) { // spiked
            log_lik[i] = log(pS);
        }
    }
    log_like = sum(log_lik);
}

