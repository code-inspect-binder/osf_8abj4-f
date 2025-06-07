data {
	int<lower=1> Nobs;           //Total no. of obs
	int<lower=1> Nind;           //Total no. of indvds
	int<lower=1> Nstim;          //No. of stim sizes
	real bias[Nobs];             //Bias scores 
	real stim[Nobs];             //Square size on each trial
	real RM[Nobs];               //Running mean on each trial
	real AS[Nobs];               //Anchor stim on each trial
	real B1[Nobs];               //Previous stim on each trial
	int<lower=1> id[Nobs];       //Participant IDs
	real mid;                    //Value of midpoint stim
	int stim_sizes[Nstim];       //Set of stim values
}
parameters {
	vector[Nind] a_raw;          //Raw individual-level parameters
	vector[Nind] b_raw;
	vector[Nind] wRM_raw;            
	vector[Nind] wAS_raw;            
	vector[Nind] wB1_raw;       
	real<lower=0> sig;           
	real a_mu;                   //Group means 
	real b_mu;
	real wRM_mu;                 
	real wAS_mu;                 
	real wB1_mu;                 
	real<lower=0> a_sig;         //Group standard deviations
	real<lower=0> b_sig;          
	real<lower=0> wRM_sig;       
	real<lower=0> wAS_sig;       
	real<lower=0> wB1_sig;       
}
transformed parameters {
	vector[Nind] a = a_mu + a_sig*a_raw;           //Individual-level parameters after shifting and rescaling
	vector[Nind] b = b_mu + b_sig*b_raw;           
	vector[Nind] wRM = wRM_mu + wRM_sig*wRM_raw;   
	vector[Nind] wAS = wAS_mu + wAS_sig*wAS_raw;   
	vector[Nind] wB1 = wB1_mu + wB1_sig*wB1_raw;   
}
model {
	real lambda;
	
	a_mu ~ normal(1, 5);           //Normal priors for group means
	b_mu ~ normal(0, 1);     
	wRM_mu ~ normal(0, 20);
	wAS_mu ~ normal(0, 20);
	wB1_mu ~ normal(0, 20);
	a_sig ~ normal(0, 10);         //Half-normal priors for group SDs 
	b_sig ~ normal(0, 5);
	wRM_sig ~ normal(0, 20);
	wAS_sig ~ normal(0, 20);
	wB1_sig ~ normal(0, 20);
	
	a_raw ~ std_normal();          //Implies a ~ normal(a_mu, a_sig)
	b_raw ~ std_normal();
	wRM_raw ~ std_normal();
	wAS_raw ~ std_normal();
	wB1_raw ~ std_normal();
	
	sig ~ normal(0, 20);                           
	
	for (i in 1:Nobs) {
		lambda = inv_logit(a[id[i]] + b[id[i]]*((stim[i]-mid)^2));
		bias[i] ~ normal((lambda*stim[i] + (1-lambda)*(wRM[id[i]]*RM[i] + wAS[id[i]]*AS[i] + wB1[id[i]]*B1[i]))-stim[i], sig);
	}
}
generated quantities {
	real pred_bias[Nobs];
	real pred_lambda[Nstim];
	real lambda;
	
	for (i in 1:Nobs) {
		lambda = inv_logit(a[id[i]] + b[id[i]]*((stim[i]-mid)^2));
		pred_bias[i] = normal_rng((lambda*stim[i] + (1-lambda)*(wRM[id[i]]*RM[i] + wAS[id[i]]*AS[i] + wB1[id[i]]*B1[i]))-stim[i], sig);     //Predicted y given draw from posterior 
	}

	for (s in 1:Nstim) 
		pred_lambda[s] = inv_logit(a_mu + b_mu*((stim_sizes[s]-mid)^2));      //Posteriors for lambda at various stim sizes 
}


