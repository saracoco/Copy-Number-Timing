data {
  // input the karyotyoe to specify the formula for theta
  int S; // number of segments
  int K; // clocks
  int N; // total number of mutations
  array[S] int karyotype; // list of karyotype associated wit the segments
  array[N] int seg_assignment; // segment_id assignment to the mutations
  array[S,2] real peaks;  //for each segment, S vectors of dim 2
  array[N] int NV;     // for all the segments
  array[N] int DP;
}

parameters {
  array[S] simplex[K] w;  //mixing proportions for each segment group;
  vector<lower=0,upper=1>[K] tau;     //clocks   ordered[K] tau  vector<lower=0,upper=1>[K]
}


transformed parameters{

  array[S,K,2] real<lower=0,upper=1> theta; //binomial mixing proportions // array[S,K] simplex[2] theta;
  
  for (s in 1:S){
   for (k in 1:K){
      if (karyotype[s] == 1) {
        theta[s,k,1] = (3 - 2*tau[k]) / (3 - tau[k]);      // 2:1
        theta[s,k,2] = tau[k] / (3 - tau[k]);
      }
      else {
        theta[s,k,1] = (2 - 2*tau[k]) / (2 - tau[k]);      // 2:0 - 2:2
        theta[s,k,2] = tau[k]/(2 - tau[k]);
      }
    }
  }

  
}



model {
  vector[K*2] contributions;
  // priors
  for (s in 1:S){
      w[s] ~ dirichlet(rep_vector(0.009, K));
  }
  tau ~ beta(2,2);
  //tau ~ uniform(0,1);  //as it is more difficult to separete events in the late/early period we could try to favour the priors on late or early periods ? 


  //likelihood
      for (i in 1:N) {
        int c = 1;
        for (k in 1:K) {
         for (j in 1:2) {
          contributions[c] = log(w[seg_assignment[i],k]) + log(theta[seg_assignment[i],k,j]) + binomial_lpmf(NV[i] | DP[i], peaks[seg_assignment[i],j]);
          c += 1;
        }
      }
      target += log_sum_exp(contributions);
    }
}


generated quantities {
  // we will generate a component identifier "comp"
  // then draw "NV_pred" using correct component theta and p

  array[N] int comp_tau;    //vector<lower=1,upper=K>[N]
  array[N] int comp_binomial; //vector<lower=1,upper=2>[N]

  vector[N] NV_pred;
  vector[2] theta_vector;

   for(i in 1:N){
     theta_vector = rep_vector(0,2);
     comp_tau[i] = categorical_rng(w[seg_assignment[i]]);
     theta_vector = to_vector(theta[seg_assignment[i], comp_tau[i]]);
     comp_binomial[i] = categorical_rng(theta_vector);
     NV_pred[i] = binomial_rng(DP[i], peaks[seg_assignment[i], comp_binomial[i]]);

   }
}

