/*
Code for estimating the original Katz and King
model for compositional data using a inverse softmax 
transformation and a multivariate t distribution.
The code extends the original model by including
a simple random intercept to allow for partial pooling.
*/
data {
  int<lower=1> N; //Number of obs.
  int<lower=1> L; //Number of congresspeople
  int<lower=1> P; //Nr. of predictors
  int<lower=2> T_full; //Nr. of cosp. portfolio dimensions 
  int<lower=2> T; //Nr. of cosp. portfolio dimensions - 1 
  matrix[N,T_full] Y; //Outcome: position on the cosp. portfolio
  matrix[N,P] X; //Design matrix
  int<lower=1> N_pred; //Cases for prediction
  int leg[N]; //congressperson indicator
  matrix[P,N_pred] X_pred; //Predictors for predicted scenarios
  int<lower=0,upper=1> fresh;//indicator for whether model is estimated only on freshmen
}

transformed data {
  vector[T] Y_raw[N]; //Projected cosp. portfolio
  for(i in 1:N)
    for(t in 1:T)
      Y_raw[i,t] = log(Y[i,(t+1)]/Y[i,1]);
}

parameters {
  matrix[L,T] alpha_l; //random intercepts by MC by dimension
  vector<lower=0>[T] sigma_alphal; //standard dev. of MC RI by dimension
  matrix[P,T] z; //std. normal location of beta coefficients
  real<lower=2> nu; //degrees of freedom scalar for MVT
  cholesky_factor_corr[T] L_Omega; //correlation cholesky factor for MVT
  cholesky_factor_corr[P] L_Beta; //correlation cholesky factor for coefficients
  vector<lower=0,upper=pi()/2>[T] tau_unif;
  vector<lower=0,upper=pi()/2>[P] tau_b_unif;
}

transformed parameters {
  matrix[P,T] beta;
  vector<lower=0>[T] tau; //standard deviations of marginals in MVT
  vector<lower=0>[P] tau_b; //standard deviations of marginals in coeff matrix
  matrix[T,T] Sigma;
  for (t in 1:T) tau[t] = 2.5 * tan(tau_unif[t]);
  for (p in 1:P) tau_b[p] = 2.5 * tan(tau_b_unif[p]);
  beta = (diag_pre_multiply(tau_b,L_Beta) * z);
  Sigma = diag_pre_multiply(tau, L_Omega) * diag_pre_multiply(tau, L_Omega)';  
}

model {
  //Likelihood
  row_vector[T] mu[N];
  for(i in 1:N)
  {
    mu[i] =  X[i,] * beta;
    if(fresh==0)
      mu[i] = mu[i] + alpha_l[leg[i],];
  }  
  
  Y_raw ~ multi_student_t(nu, mu, Sigma);

  //Random Intercepts
  for(t in 1:T)
    alpha_l[,t] ~ normal(0,sigma_alphal[t]);
  
  //Priors
  L_Omega ~ lkj_corr_cholesky(3);
  L_Beta ~ lkj_corr_cholesky(3);
  nu ~ gamma(2,0.1);
  sigma_alphal ~ cauchy(0,2.5);
  to_vector(z) ~ normal(0,1);
}

generated quantities {
  matrix[T_full,N_pred] Port_pred;
  {
    vector[T_full] Y_temp;
    matrix[T,N_pred] mu_pred; //Inner scope, so it's not stored
    mu_pred = beta' * X_pred; //No RIs, so average.
    for(i in 1:N_pred)
    {
      Y_temp[1] = 0;
      Y_temp[2:T_full] = multi_student_t_rng(nu, mu_pred[,i], Sigma); 
      Port_pred[,i] = softmax(Y_temp);
    }
  }

}
