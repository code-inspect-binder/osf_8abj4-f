data {
	int<lower=1> Nobs;           //Total no. of obs
	int<lower=1> Nind;           //Total no. of indvds
	real rating[Nobs];           //Magnitude ratings 
	real stim[Nobs];             //Square sizes
	real running[Nobs];          //Running mean
	real prev[Nobs];             //Previous stim on each trial
	real freq[Nobs];             //Frequency values
	real smin[Nobs];             //Running scale minimum
	real smax[Nobs];             //Running scale maximum
	int<lower=1> id[Nobs];       //Subject identifier
}
parameters {
	vector[Nind] c_raw;          //Raw individual-level parameters
	vector[Nind] wM_raw;
	vector[Nind] wP_raw; 
	vector[Nind] b_raw;
	vector[Nind] w_raw;      
	real<lower=0> sig;           
	real c_mu;                   //Group means 
	real wM_mu;
	real wP_mu;
	real b_mu;
	real w_raw_mu;               
	real<lower=0> c_sig;         //Group standard deviations 
	real<lower=0> wM_sig; 
	real<lower=0> wP_sig;
	real<lower=0> b_sig;       
	real<lower=0> w_raw_sig;       
}
transformed parameters {
	vector[Nind] c = c_mu + c_sig*c_raw;           //Individual-level parameters after shifting and rescaling        
	vector[Nind] wM = wM_mu + wM_sig*wM_raw;     
	vector[Nind] wP = wP_mu + wP_sig*wP_raw;
	vector[Nind] b = b_mu + b_sig*b_raw;       
	vector[Nind] w = w_raw_mu + w_raw_sig*w_raw;       
	vector[Nind] wR = inv_logit(w);                //Range-weighting on [0,1] scale
	real wR_mu = inv_logit(w_raw_mu);              //Group mean range-weighting
	real wF_mu = 1-wR_mu;                          //Group mean frequency-weighting
}
model {
	real smax_smin;                   //Magnitude range (SMax - SMin)
	
	c_mu ~ normal(0, 20);             //Priors for group means
	wM_mu ~ normal(0, 20);
	wP_mu ~ normal(0, 20);
	b_mu ~ normal(12, 20);            //12 pt difference between min size (2) and max size (14)
	w_raw_mu ~ normal(0, sqrt(3));    //Induces approximately uniform prior after sigmoid transform
	
	c_sig ~ normal(0, 20);            //Priors for group SDs 
	wM_sig ~ normal(0, 20);
	wP_sig ~ normal(0, 20); 
	b_sig ~ normal(0, 20);
	w_raw_sig ~ normal(0, 1);
	
	c_raw ~ std_normal();       //Implies: c ~ normal(c_mu, c_sig)
	wM_raw ~ std_normal();
	wP_raw ~ std_normal(); 
	b_raw ~ std_normal();
	w_raw ~ std_normal();
	
	sig ~ normal(0, 20);                           
	
	for (i in 1:Nobs) {
		if (smax[i]==smin[i]) {
			smax_smin = 1;
		}
		else {
			smax_smin = smax[i] - smin[i];
		}
		rating[i] ~ normal(c[id[i]] + wM[id[i]]*running[i] + wP[id[i]]*prev[i] + b[id[i]]*(wR[id[i]]*((stim[i]-smin[i])/(smax_smin)) + (1-wR[id[i]])*freq[i]), sig);
	}
	
}
generated quantities {
	real pred_rating[Nobs];
	real smax_smin;

	for (i in 1:Nobs) {
		if (smax[i]==smin[i]) {
			smax_smin = 1;
		}
		else {
			smax_smin = smax[i] - smin[i];
		}
		pred_rating[i] = normal_rng(c[id[i]] + wM[id[i]]*running[i] + wP[id[i]]*prev[i] + b[id[i]]*(wR[id[i]]*((stim[i]-smin[i])/(smax_smin)) + (1-wR[id[i]])*freq[i]), sig);
	}	
}




