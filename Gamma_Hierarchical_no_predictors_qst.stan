data {

	// Number of Values 
   int<lower=0> Ni; // Number of level-1 observations (Replicates)
   int<lower=0> Ng; // Number of level-2 clusters (Genotypes)   
   int<lower=0> Np; // Number of level-3 clusters (Populations)
 
   // IDs
   int<lower=1> genotype_ids[Ni];
   int<lower=1> population_ids[Ni];
   
   // Translate IDs across levels for loops
   int<lower=1> PopsToGens[Ng]; // Maps populations to genotypes 
         
   // Predictors and Response
   real<lower=0> Y_igp[Ni]; // Continuous outcome - Raw Storage Data
 }
parameters {
	
   // Intercepts
   real             beta_0;        // Mean of the data
   real<lower=0.01> beta_0gp[Ng]; // Genotype Intercepts 
   real<lower=0.01> beta_0p[Np];  // Population intercepts
   
   // Slopes
   real beta_1; // Slope of pedictor (coppiced or uncoppiced)
 
   // Variances
   real<lower=0> phi_igp;     // Level-1, SD of data within genotypes
   real<lower=0> phi_gp;        // Level-2, Variance of genotypes within regions
   real<lower=0> phi_p;         // Level-3, Variance of populations within regions 

 }
transformed parameters  {
 
 	// Define transformed parameters
 	real<lower=0> mu[Ni];         // Individual tree means
	
	real<lower=0> gamAlpha_igp[Ni]; // Genotype Gamma Coefficients
	real<lower=0> gamBeta_igp[Ni];  // Genotype Gamma Coefficients
		
	real<lower=0> gamAlpha_gp[Ng]; // Genotype Gamma Coefficients
	real<lower=0> gamBeta_gp[Ng];  // Genotype Gamma Coefficients
	
	real<lower=0> gamAlpha_p[Np]; // Population Gamma Coefficients
	real<lower=0> gamBeta_p[Np];  // Population Gamma Coefficients
	
     // Level-3: Population Gamma Relationships
	for(p in 1:Np){ // Level - 2	   
	  gamAlpha_p[p] = (beta_0 * beta_0) / phi_p; 
	  gamBeta_p[p] = beta_0 / phi_p; 
  	}

     // Level-2: Genotype Gamma Relationships
	for(g in 1:Ng){ // Level - 2	   
	  gamAlpha_gp[g] = (beta_0p[PopsToGens[g]] * beta_0p[PopsToGens[g]]) / phi_gp; 
	  gamBeta_gp[g] = beta_0p[PopsToGens[g]] / phi_gp; 
  	}
   
   // Level-1: Equation for each individual tree's storage value = genotype mean + coppiced effect
   for (i in 1:Ni) {
	 mu[i] = beta_0gp[genotype_ids[i]];
	 gamAlpha_igp[i] = ( mu[i] * mu[i] ) / phi_igp; 
	 gamBeta_igp[i] = mu[i] / phi_igp;
   }
   
}
model {
	 
   // Prior part of Bayesian inference
   // Flat pior for mu (no need to specify if non-informative)
 
   // Random effects distribution
	 for(p in 1:Np){  // Level-3
		 beta_0p[p] ~ gamma( gamAlpha_p[p], gamBeta_p[p] );
	 }
     for(g in 1:Ng){ // Level - 2
		 beta_0gp[g] ~ gamma( gamAlpha_gp[g] , gamBeta_gp[g] );
	  }
 
   // Likelihood
   for (i in 1:Ni) {
     Y_igp[i] ~ gamma( gamAlpha_igp[i], gamBeta_igp[i] );
   }
}
generated quantities{
	
	real y_rep[Ni];
	
	for(i in 1:Ni){
		y_rep[i] = gamma_rng( gamAlpha_igp[i], gamBeta_igp[i] );
	}
}