data {
	int<lower=1> Nobs;           //Total no. of obs
	int<lower=1> Nind;           //Total no. of indvds
	real bias[Nobs];             //Bias scores
	real stim[Nobs];             //Square sizes
	real distmean[Nobs];         //Distribution mean (constant)
	real anchor[Nobs];           //Anchor stim on each trial
	real prev[Nobs];             //Previous stim on each trial
	real freq[Nobs];             //Frequency values
	int<lower=1> id[Nobs];       //Subject identifier
	real smin;                   //Min stim magnitude
	real smax;                   //Max stim magnitude 
}
parameters {
	vector[Nind] c_raw;          //Raw individual-level parameters
	vector[Nind] wM_raw;
	vector[Nind] wA_raw;
	vector[Nind] wP_raw; 
	vector[Nind] b_raw;
	vector[Nind] w_raw;      
	real<lower=0> sig;           
	real c_mu;                   //Group means 
	real wM_mu;
	real wA_mu;
	real wP_mu;
	real b_mu;
	real w_raw_mu;               
	real<lower=0> c_sig;         //Group standard deviations 
	real<lower=0> wM_sig;
	real<lower=0> wA_sig;  
	real<lower=0> wP_sig;
	real<lower=0> b_sig;       
	real<lower=0> w_raw_sig;       
}
transformed parameters {
	vector[Nind] c = c_mu + c_sig*c_raw;           //Individual-level parameters after shifting and rescaling        
	vector[Nind] wM = wM_mu + wM_sig*wM_raw;   
	vector[Nind] wA = wA_mu + wA_sig*wA_raw;   
	vector[Nind] wP = wP_mu + wP_sig*wP_raw;
	vector[Nind] b = b_mu + b_sig*b_raw;       
	vector[Nind] w = w_raw_mu + w_raw_sig*w_raw;       
	vector[Nind] wR = inv_logit(w);                //Range-weighting on [0,1] scale
	real wR_mu = inv_logit(w_raw_mu);              //Group mean range-weighting
	real wF_mu = 1-wR_mu;                          //Group mean frequency-weighting
}
model {
	c_mu ~ normal(4, 20);       //Priors for group means
	wM_mu ~ normal(0, 20);
	wA_mu ~ normal(0, 20);
	wP_mu ~ normal(0, 20);
	b_mu ~ normal(12, 20);
	w_raw_mu ~ normal(0, 10);     
	c_sig ~ normal(0, 20);      //Priors for group SDs 
	wM_sig ~ normal(0, 20);
	wA_sig ~ normal(0, 20);
	wP_sig ~ normal(0, 20); 
	b_sig ~ normal(0, 20);
	w_raw_sig ~ normal(0, 10);
	
	c_raw ~ std_normal();       //Implies: c ~ normal(c_mu, c_sig)
	wM_raw ~ std_normal();
	wA_raw ~ std_normal();
	wP_raw ~ std_normal(); 
	b_raw ~ std_normal();
	w_raw ~ std_normal();
	
	sig ~ normal(0, 20);                           
	
	for (i in 1:Nobs) 
		bias[i] ~ normal((c[id[i]] + wM[id[i]]*distmean[i] + wA[id[i]]*anchor[i] + wP[id[i]]*prev[i] + b[id[i]]*(wR[id[i]]*((stim[i]-smin)/(smax-smin)) + (1-wR[id[i]])*freq[i]))-stim[i], sig);
	
}
generated quantities {
	real pred_bias[Nobs];

	for (i in 1:Nobs) 
		pred_bias[i] = normal_rng((c[id[i]] + wM[id[i]]*distmean[i] + wA[id[i]]*anchor[i] + wP[id[i]]*prev[i] + b[id[i]]*(wR[id[i]]*((stim[i]-smin)/(smax-smin)) + (1-wR[id[i]])*freq[i]))-stim[i], sig);
}


