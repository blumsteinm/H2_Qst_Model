data {

	// Number of Values 
   int<lower=0> Ni; // Number of level-1 observations (Replicates)
   int<lower=0> Ng; // Number of level-2 clusters (Genotypes)   
 
   // IDs
   int<lower=1> genotype_ids[Ni];
         
   // Predictors and Response
   real<lower=0> Y_ig[Ni]; // Continuous outcome - Raw Storage Data
 }
parameters {
	
   // Intercepts
   real             beta_0;        // Mean of the data
   real<lower=0.01> beta_0g[Ng]; // Genotype Intercepts 
 
   // Variances
   real<lower=0> phi_ig;     // Level-1, Variance of data within genotypes
   real<lower=0> phi_g;        // Level-2, Variance of genotypes within regions

 }
transformed parameters  {
 
 	// Define transformed parameters
 	real<lower=0> mu[Ni];         // Individual tree means
	
	real<lower=0> gamAlpha_ig[Ni]; // Replicate Gamma Coefficients
	real<lower=0> gamBeta_ig[Ni];  // Replicate Gamma Coefficients
		
	real<lower=0> gamAlpha_g[Ng]; // Genotype Gamma Coefficients
	real<lower=0> gamBeta_g[Ng];  // Genotype Gamma Coefficients
	
     // Level-2: Genotype Gamma Relationships
	for(g in 1:Ng){ // Level - 2	   
	  gamAlpha_g[g] = (beta_0 * beta_0) / phi_g; 
	  gamBeta_g[g] = beta_0 / phi_g; 
  	}
   
   // Level-1: Equation for each individual tree's storage value = genotype mean + coppiced effect
   for (i in 1:Ni) {
	 mu[i] = beta_0g[genotype_ids[i]];
	 gamAlpha_ig[i] = ( mu[i] * mu[i] ) / phi_ig; 
	 gamBeta_ig[i] = mu[i] / phi_ig;
   }
   
}
model {
	 
   // Prior part of Bayesian inference
   // Flat pior for mu (no need to specify if non-informative)
 
     for(g in 1:Ng){ // Level - 2
		 beta_0g[g] ~ gamma( gamAlpha_g[g] , gamBeta_g[g] );
	  }
 
   // Likelihood
   for (i in 1:Ni) {
     Y_ig[i] ~ gamma( gamAlpha_ig[i], gamBeta_ig[i] );
   }
}
generated quantities{
	
	real y_rep[Ni];
	
	for(i in 1:Ni){
		y_rep[i] = gamma_rng( gamAlpha_ig[i], gamBeta_ig[i] );
	}
}