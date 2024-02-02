//All years estimated independently for State of Birds model

data {
  int<lower=0> n_years; // number of years of full time-series
  int<lower=0> n_species; // number of species in the group
  array[n_years,n_species] real ln_index; // log scale annual indices of abundance or annual population size
  array[n_years,n_species] real ln_index_sd; // SD of the log scale annual indices

}


parameters {
  //real<lower=0> sigma;
  array[n_years] real MU;
  array[n_years,n_species] real ln_index_true; // log scale annual indices of abundance or annual population size
  array[n_years] real<lower=0> sigma;    // sd each yearly summary

}

model {
  MU ~ student_t(3,0,1);

  for(i in 1:n_years){
  ln_index[i,] ~ normal(ln_index_true[i,], ln_index_sd[i,]);
  ln_index_true[i,] ~ normal(MU[i],sigma[i])
  }

  sigma ~ student_t(3,0,1);
}

// generated quantities {
// vector[n_years] annual_diffs;
// vector[n_years] scaled_status;
// vector[n_years] scaled_log_status;
//
// annual_diffs[1] = 0;
// scaled_status[1] = 0;
// scaled_log_status[1] = 1;
//
//
// for(y in 2:n_years){
//   scaled_status[y] = mu[y]-mu[1];
//   scaled_log_status[y] = exp(mu[y])/exp(mu[1]);
//   annual_diffs[y] = mu[y]-mu[y-1];
//
// }
//
//
// }
