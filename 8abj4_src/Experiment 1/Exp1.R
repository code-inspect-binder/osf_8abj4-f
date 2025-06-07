#==== EXPERIMENT 1 ANALYSES ====

library(BayesFactor)
library(rstan)
library(parallel)
library(ggplot2)
library(gridExtra)

dat <- read.table('Exp1_data.csv', header=T, sep=',')
dat$id <- rep(1:48, each=66)
names(dat)[2:4] <- c('AOrder','TaskOrd','Skew')

dat$Skew_F <- factor(dat$Skew, levels=c('P','N'), labels=c('Positive','Negative'))
dat$Start_F <- factor(dat$Start, levels=c(1,2), labels=c('Small','Large'))
dat$Stim_F <- factor(dat$Size)

dat$Bias <- dat$Est - dat$Size

#Remove outlier trials 
dat <- subset(dat, abs(scale(Bias)) < 3)

#==== Bayesian ANOVA for Reproduction Task ====

dat1 <- aggregate(Bias ~ Skew_F*Start_F*Stim_F*subject, data=dat, FUN=mean)
dat1$subject <- factor(dat1$subject)

set.seed(234732)

#Skew x Anchor x Stimulus ANOVA 
modBF <- anovaBF(Bias ~ subject + Skew_F*Start_F*Stim_F,  #Random effect for each participant
                 data=dat1,
                 whichModels = 'withmain',
                 whichRandom = 'subject',
                 iterations = 50000)
modBF <- sort(modBF, decreasing = T)
modBF

#BF for Anchor x Stimulus interaction: 
modBF[1]/modBF[5]

#BF for Stimulus effect: 
modBF[1]/modBF[12] 

#BF forAnchor effect:
modBF[1]/modBF[15] 

#BF in favor of the Skew effect:
modBF[1]/modBF[2]  

#Extract posterior distribution
bestmod <- modBF[1]
chains <- posterior(bestmod, iterations=50000, columnFilter="^subject$")
colnames(chains)

#Polynomial contrasts for each anchor
poly <- contr.poly(7)

sm20 <- apply(chains[,c(1,4,6,13)], 1, sum)
sm40 <- apply(chains[,c(1,4,7,14)], 1, sum)
sm60 <- apply(chains[,c(1,4,8,15)], 1, sum)
sm80 <- apply(chains[,c(1,4,9,16)], 1, sum)
sm100 <- apply(chains[,c(1,4,10,17)], 1, sum)
sm120 <- apply(chains[,c(1,4,11,18)], 1, sum)
sm140 <- apply(chains[,c(1,4,12,19)], 1, sum)
sm_mat <- cbind(sm20,sm40,sm60,sm80,sm100,sm120,sm140)
quantile(sm_mat %*% poly[,'.L'], prob=c(.025,.975))  
quantile(sm_mat %*% poly[,'.Q'], prob=c(.025,.975))  
quantile(sm_mat %*% poly[,'.C'], prob=c(.025,.975))  
quantile(sm_mat %*% poly[,'^4'], prob=c(.025,.975))  
quantile(sm_mat %*% poly[,'^5'], prob=c(.025,.975))  
quantile(sm_mat %*% poly[,'^6'], prob=c(.025,.975))  

lg20 <- apply(chains[,c(1,5,6,20)], 1, sum)
lg40 <- apply(chains[,c(1,5,7,21)], 1, sum)
lg60 <- apply(chains[,c(1,5,8,22)], 1, sum)
lg80 <- apply(chains[,c(1,5,9,23)], 1, sum)
lg100 <- apply(chains[,c(1,5,10,24)], 1, sum)
lg120 <- apply(chains[,c(1,5,11,25)], 1, sum)
lg140 <- apply(chains[,c(1,5,12,26)], 1, sum)
lg_mat <- cbind(lg20,lg40,lg60,lg80,lg100,lg120,lg140)
quantile(lg_mat %*% poly[,'.L'], prob=c(.025,.975))  
quantile(lg_mat %*% poly[,'.Q'], prob=c(.025,.975))  
quantile(lg_mat %*% poly[,'.C'], prob=c(.025,.975))  
quantile(lg_mat %*% poly[,'^4'], prob=c(.025,.975))  
quantile(lg_mat %*% poly[,'^5'], prob=c(.025,.975))  
quantile(lg_mat %*% poly[,'^6'], prob=c(.025,.975)) 

#Anchor x Stim interaction: Differences b/w polynomial contrasts for each anchor
quantile((lg_mat %*% poly[,'.L'])-(sm_mat %*% poly[,'.L']), prob=c(.025,.975))
quantile((lg_mat %*% poly[,'.Q'])-(sm_mat %*% poly[,'.Q']), prob=c(.025,.975)) 
quantile((lg_mat %*% poly[,'.C'])-(sm_mat %*% poly[,'.C']), prob=c(.025,.975))
quantile((lg_mat %*% poly[,'^4'])-(sm_mat %*% poly[,'^4']), prob=c(.025,.975))
quantile((lg_mat %*% poly[,'^5'])-(sm_mat %*% poly[,'^5']), prob=c(.025,.975))
quantile((lg_mat %*% poly[,'^6'])-(sm_mat %*% poly[,'^6']), prob=c(.025,.975))

#Skew effect
neg <- apply(chains[,c(1,3)], 1, sum)
pos <- apply(chains[,c(1,2)], 1, sum)
quantile(neg-pos, prob=c(.025, .975))

#==== Modified CA model: Hierarchical Bayes ====

dat$ASize <- numeric(nrow(dat))
dat$ASize[dat$Start==1] <- 10
dat$ASize[dat$Start==2] <- 150

#Remove 1st trial for each participant 
dat <- subset(dat, Trial != 1)

#Inputs to Stan model
stim <- dat$Size/10
bias <- dat$Bias/10
RM <- dat$Running/10
AS <- dat$ASize/10
B1 <- dat$Back1/10
id <- dat$id
Nobs <- 3080
Nind <- 48
mid <- 8
stim_sizes <- seq(20,140,20)/10
Nstim <- length(stim_sizes)

no_cores <- detectCores()-1

#Fit hierarchical Bayesian model
fit.HBM <- stan(file='CAM_full_1.stan',
                data=c('Nobs','Nind','stim','bias','RM','AS','B1','id','mid','stim_sizes','Nstim'),
                chains=4,
                iter=2500,
                cores=no_cores,
                control=list(adapt_delta=.90,
                             max_treedepth=15))
fit.HBM
sum(summary(fit.HBM)$summary[,'Rhat'] > 1.01)
pars <- extract(fit.HBM)

#Credible intervals for group means
quantile(pars$wRM_mu, prob=c(.025, .975))
quantile(pars$wAS_mu, prob=c(.025, .975))
quantile(pars$wB1_mu, prob=c(.025, .975))

#Get model predictions on each trial
pred_mat <- pars$pred_bias          #Posterior predictive distribution on each trial
preds <- apply(pred_mat, 2, mean)   #Mean of posterior predictive dist. on each trial
cor.test(bias, preds)
dat$pred_bias <- preds
dat$Bias10 <- dat$Bias/10

#Plotting individual fits
dat1 <- aggregate(cbind(Bias10, pred_bias) ~ Skew_F*Start_F*Size*subject, data=dat, FUN=mean)
dat1a <- aggregate(Bias10 ~ Skew_F*Start_F*Size*subject, data=dat, 
                   FUN=function(x)sd(x)/sqrt(length(x)))
dat1$ll <- dat1$Bias10 - dat1a$Bias10
dat1$ul <- dat1$Bias10 + dat1a$Bias10

dat1$Size10 <- dat1$Size/10
p1 <- ggplot(data=subset(dat1, Skew_F=='Positive'), aes(x=Size10)) +
  geom_pointrange(aes(y=Bias10, ymin=ll, ymax=ul, col=Start_F), fatten=2) +
  geom_line(aes(y=pred_bias, col=Start_F), lwd=1.2, alpha=.4) +
  facet_wrap(~subject, nrow=6) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        strip.text=element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.text=element_text(color='white'),
        legend.title=element_text(color='white'),
        legend.key=element_rect(fill='white')) +
  geom_abline(aes(slope=0, intercept=0), lty=3, alpha=.5) +
  xlab('Square Size') +
  ylab('Reproduction Bias') +
  ggtitle('Positive Skew') +
  scale_color_discrete(guide = guide_legend(override.aes = list(color = "white")))
p2 <- ggplot(data=subset(dat1, Skew_F=='Negative'), aes(x=Size10)) +
  geom_pointrange(aes(y=Bias10, ymin=ll, ymax=ul, col=Start_F), fatten=2) +
  geom_line(aes(y=pred_bias, col=Start_F), lwd=1.2, alpha=.4) +
  facet_wrap(~subject, nrow=6) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        strip.text=element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  geom_abline(aes(slope=0, intercept=0), lty=3, alpha=.5) +
  xlab('Square Size') +
  ylab('Reproduction Bias') +
  guides(col=guide_legend(title='Anchor')) +
  ggtitle('Negative Skew')
grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()

#==== Bayesian ANOVA for Ratings ====

dat <- read.table('Exp1_data.csv', header=T, sep=',')
dat$id <- rep(1:48, each=66)
names(dat)[2:4] <- c('AOrder','TaskOrd','Skew')

dat$Skew_F <- factor(dat$Skew, levels=c('P','N'), labels=c('Positive','Negative'))
dat$Start_F <- factor(dat$Start, levels=c(1,2), labels=c('Small','Large'))
dat$Stim_F <- factor(dat$Size)

dat1 <- aggregate(Rating ~ Skew_F*Stim_F*subject, data=dat, FUN=mean)
dat1$subject <- factor(dat1$subject)

set.seed(239813)
modBF <- anovaBF(Rating ~ subject + Skew_F*Stim_F,
                 data=dat1,
                 whichModels = 'withmain',
                 whichRandom = 'subject',
                 iterations = 50000)
modBF <- sort(modBF, decreasing = T)
modBF

#BF for Stimulus x Skew interaction
modBF[1]/modBF[2]

#BFfor Distribution effect
modBF[1]/modBF[3]

#BF for Stimulus effect
modBF[1]/modBF[4]

bestmod <- modBF[1]
chains <- posterior(bestmod, iterations=50000, columnFilter="^subject$")
colnames(chains)

#Polynomial contrasts for each skew condition
poly <- contr.poly(7)

pos20 <- apply(chains[,c(1,2,4,11)], 1, sum)
pos40 <- apply(chains[,c(1,2,5,12)], 1, sum)
pos60 <- apply(chains[,c(1,2,6,13)], 1, sum)
pos80 <- apply(chains[,c(1,2,7,14)], 1, sum)
pos100 <- apply(chains[,c(1,2,8,15)], 1, sum)
pos120 <- apply(chains[,c(1,2,9,16)], 1, sum)
pos140 <- apply(chains[,c(1,2,10,17)], 1, sum)
pos_mat <- cbind(pos20,pos40,pos60,pos80,pos100,pos120,pos140)
quantile(pos_mat %*% poly[,'.L'], prob=c(.025,.975))  
quantile(pos_mat %*% poly[,'.Q'], prob=c(.025,.975))  
quantile(pos_mat %*% poly[,'.C'], prob=c(.025,.975))  
quantile(pos_mat %*% poly[,'^4'], prob=c(.025,.975))  
quantile(pos_mat %*% poly[,'^5'], prob=c(.025,.975))  
quantile(pos_mat %*% poly[,'^6'], prob=c(.025,.975)) 

neg20 <- apply(chains[,c(1,3,4,18)], 1, sum)
neg40 <- apply(chains[,c(1,3,5,19)], 1, sum)
neg60 <- apply(chains[,c(1,3,6,20)], 1, sum)
neg80 <- apply(chains[,c(1,3,7,21)], 1, sum)
neg100 <- apply(chains[,c(1,3,8,22)], 1, sum)
neg120 <- apply(chains[,c(1,3,9,23)], 1, sum)
neg140 <- apply(chains[,c(1,3,10,24)], 1, sum)
neg_mat <- cbind(neg20,neg40,neg60,neg80,neg100,neg120,neg140)
quantile(neg_mat %*% poly[,'.L'], prob=c(.025,.975))  
quantile(neg_mat %*% poly[,'.Q'], prob=c(.025,.975))  
quantile(neg_mat %*% poly[,'.C'], prob=c(.025,.975))  
quantile(neg_mat %*% poly[,'^4'], prob=c(.025,.975))  
quantile(neg_mat %*% poly[,'^5'], prob=c(.025,.975))  
quantile(neg_mat %*% poly[,'^6'], prob=c(.025,.975))

#Skew effect
neg <- apply(chains[,c(1,3)], 1, sum)
pos <- apply(chains[,c(1,2)], 1, sum)
quantile(pos-neg, prob=c(.025, .975))

#Skew x Stimulus interaction: Difference b/w polynomial contrasts for each skew condition
quantile((pos_mat %*% poly[,'.L'])-(neg_mat %*% poly[,'.L']), prob=c(.025,.975))
quantile((pos_mat %*% poly[,'.Q'])-(neg_mat %*% poly[,'.Q']), prob=c(.025,.975))
quantile((neg_mat %*% poly[,'.C'])-(pos_mat %*% poly[,'.C']), prob=c(.025,.975))
quantile((pos_mat %*% poly[,'^4'])-(neg_mat %*% poly[,'^4']), prob=c(.025,.975))
quantile((pos_mat %*% poly[,'^5'])-(neg_mat %*% poly[,'^5']), prob=c(.025,.975))
quantile((pos_mat %*% poly[,'^6'])-(neg_mat %*% poly[,'^6']), prob=c(.025,.975))

#==== Range-frequency Model for Ratings ====

dat <- read.table('Exp1_data.csv', header=T, sep=',')
dat$id <- rep(1:48, each=66)
names(dat)[2:4] <- c('AOrder','TaskOrd','Skew')

#Compute running frequency values, SMin, and SMax
dat$freq <- numeric(nrow(dat))
dat$Smin <- numeric(nrow(dat))
dat$Smax <- numeric(nrow(dat))
ids <- unique(dat$subject)
for(id in ids) {
  for(i in 1:nrow(dat[dat$subject==id,])) {
    if(i == 1) {dat[dat$subject==id,]$freq[i] <- 1}
    else {
      dat[dat$subject==id,]$freq[i] <- (rank(dat[dat$subject==id,]$Size[1:i])[i]-1)/
                                       (length(dat[dat$subject==id,]$Size[1:i])-1)
    }
    dat[dat$subject==id,]$Smin[i] <- min(dat[dat$subject==id,]$Size[1:i])
    dat[dat$subject==id,]$Smax[i] <- max(dat[dat$subject==id,]$Size[1:i])
  }
}

stim <- dat$Size/10
rating <- dat$Rating
running <- dat$Running/10
prev <- dat$Back1/10
freq <- dat$freq
smin = dat$Smin/10
smax = dat$Smax/10
id <- dat$id
Nobs <- nrow(dat)
Nind <- length(unique(id))

no_cores <- detectCores()-1

fit.HBM <- stan(file='RF_ratings.stan',
                data=c('Nobs','Nind','stim','rating','running','prev','freq','id','smin','smax'),
                pars=c('smax_smin'),
                include=FALSE,
                chains=4,
                iter=2500,
                cores=no_cores,
                control=list(adapt_delta=.90,
                             max_treedepth=15))
fit.HBM
sum(summary(fit.HBM)$summary[,'Rhat'] > 1.01)
pars <- extract(fit.HBM)

#Posteriors for group mean parameters of interest
quantile(pars$wF_mu, prob=c(.025, .975))
quantile(pars$b_mu, prob=c(.025, .975))
quantile(pars$wM_mu, prob=c(.025, .975))
quantile(pars$wP_mu, prob=c(.025, .975))

pred_mat <- pars$pred_rating
preds <- apply(pred_mat, 2, mean)
dat$pred_rating <- preds


