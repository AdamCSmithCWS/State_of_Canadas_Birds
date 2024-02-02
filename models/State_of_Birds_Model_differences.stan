//All years estimated independently for State of Birds model

data {
  int<lower=0> n_years; // number of years of full time-series
  int<lower=0> n_species; // number of species in the group
  array[n_years] int n_species_year; // vector of number of species to include in each year
  array[n_years,n_species] int species; // matrix of species indicators included in each year allows for missing species
  array[n_years,n_species] real ln_index; // log scale annual indices of abundance or annual population size
  array[n_years,n_species] real ln_index_sd; // SD of the log scale annual indices

}


parameters {
  //real<lower=0> sigma;
  array[n_years] real MU_ann;
  array[n_years,n_species] real ln_index_true; // inclusion indicator for each species and year (allows missing data)
  array[n_years] real<lower=0> sigma;    // sd each yearly summary

}

model {
  MU_ann ~ student_t(3,0,1);

  for(i in 1:n_years){
    for(s in species[i,1:n_species_year]){ // stepping through species with data
  ln_index[i,s] ~ normal(ln_index_true[i,s], ln_index_sd[i,s]);
  ln_index_true[i,s] ~ normal(MU_ann[i],sigma[i])
    }
    for(s in species[i,n_species_year:n_species]){ //stepping through species that are missing
    ln_index_true[i,s] ~ normal(0,1)
    }
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
