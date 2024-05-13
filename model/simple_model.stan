data {
  int<lower=0> n;                 
  int<lower=0> num_clones;        
  matrix[n, num_clones] clone_cn_profiles; 
  vector[n] ctdna;               
  real<lower=0> scale;            
}

parameters {
  simplex[num_clones] rho;        
}

transformed parameters {
  vector[n] mu;                   

  vector[n] total_sum;  // hold sums for each clone profile row

  for (i in 1:n) {
    // convert row to vector and perform element-wise multiplication, then sum
    mu[i] = log(sum(to_vector(clone_cn_profiles[i]) .* rho));
    
    total_sum[i] = sum(to_vector(clone_cn_profiles[i]) .* rho);
  }

  // compute the mean of the total sums ,then update mu
  real mean_total_sum = mean(total_sum);
  for (i in 1:n) {
    mu[i] = mu[i] - log(mean_total_sum);
  }

}

    

model {
  rho ~ dirichlet(rep_vector(1.0, num_clones)); // prior ~ Dir(1)

  // data modeled as student's t distribution
  ctdna ~ student_t(2, mu, scale); 
}
