data{
  int N;
  int OUTCOME[N];
  real TREAT[N];
  real AGED[N];
  int STUDY[N];
  real TREAT_X_AGED[N];
  int bin_total[N];
  int N_STUDY;
}

parameters{
  real Intercept;
  real beta_TREAT;
  real beta_AGED;
  real beta_TREAT_X_AGED;
  real vary_STUDY[N_STUDY];
  real<lower=0> sigma_STUDY;
}

model{
  real vary[N];
  real glm[N];
  // Priors
  Intercept ~ normal( 0 , 100 );
  beta_TREAT ~ normal( 0 , 100 );
  beta_AGED ~ normal( 0 , 100 );
  beta_TREAT_X_AGED ~ normal( 0 , 100 );
  sigma_STUDY ~ uniform( 0 , 100 );
  // Varying effects
  for ( j in 1:N_STUDY ) vary_STUDY[j] ~ normal( 0 , sigma_STUDY );
  // Fixed effects
  for ( i in 1:N ) {
    vary[i] <- vary_STUDY[STUDY[i]];
    glm[i] <- vary[i] + Intercept
    + beta_TREAT * TREAT[i]
    + beta_AGED * AGED[i]
    + beta_TREAT_X_AGED * TREAT_X_AGED[i];
    glm[i] <- inv_logit( glm[i] );
  }
  OUTCOME ~ binomial( bin_total , glm );
}

generated quantities{
  real dev;
  real vary[N];
  real glm[N];
  dev <- 0;
  for ( i in 1:N ) {
    vary[i] <- vary_STUDY[STUDY[i]];
    glm[i] <- vary[i] + Intercept
    + beta_TREAT * TREAT[i]
    + beta_AGED * AGED[i]
    + beta_TREAT_X_AGED * TREAT_X_AGED[i];
    dev <- dev + (-2) * binomial_log( OUTCOME[i] , bin_total[i] , inv_logit(glm[i]) );
  }
}
