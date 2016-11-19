data {
  int<lower=1> D;
  int<lower=1> N;
  vector[D] x[N];
  int<lower=0> y[N];
  int<lower=1> K;          // number of hyperparameters
  vector[K] phi;           // mu, log_tau
}
transformed data {
  real<lower=0> delta;
  // matrix[N,D] xm;
  // for (n in 1:N) 
  //   xm[n] = to_row_vector(x[n]);
  
  delta = 0.0001;
}
parameters {
  vector[N] eta;
  // real a;
  // real b;
}
transformed parameters {
  real<lower=0> length;
  real<lower=0> alpha;
  vector[N] f;
  length = exp(phi[1]);
  alpha = exp(phi[2]);
  {
    matrix[N, N] L_Sigma;
    matrix[N, N] Sigma;
    Sigma = cov_exp_quad(x, alpha, length);
    for (n in 1:N)
      Sigma[n, n] = Sigma[n,n] + delta^2;
    L_Sigma = cholesky_decompose(Sigma);
    f = L_Sigma * eta;
  }
}  
model {
  target += log(length);
  target += log(alpha);
  length ~ lognormal(0, 1);
  alpha ~ lognormal(0, 1);

  eta ~ normal(0, 1);
  y ~ bernoulli_logit(f);
}
// generated quantities {
//   vector[Nt] mu_gen;
//   // vector[Nt] y_gen;
//   {
//     matrix[N,N] Sigma;
//     matrix[N,Nt] K;
//     matrix[Nt,Nt] Omega;
//     matrix[Nt,N] K_transpose_div_Sigma;
//     matrix[Nt,Nt] Tau;
//     matrix[Nt,Nt] L2;
     
//     Sigma = cov_exp_quad(x, alpha, length);
//     for (n in 1:N)
//       Sigma[n, n] = Sigma[n, n] + delta^2;
//     K = cov_exp_quad(x, xt, alpha, length);
//     Omega = cov_exp_quad(xt, alpha, length);
//     for (n in 1:Nt)
//       Omega[n, n] = Omega[n, n] + delta^2;
    
//     K_transpose_div_Sigma = K' / Sigma;
//     mu_gen = K_transpose_div_Sigma * f;
    
//     Tau = Omega - K_transpose_div_Sigma * K;

//     for (i in 1:(Nt-1)){
//       for (j in (i+1):Nt){
//         Tau[i,j] = Tau[j,i];
//       }
//     }
//     // L2 = cholesky_decompose(Tau);
//   }
  
//   // y_gen = multi_normal_cholesky_rng(mu_gen, L2);
// }
