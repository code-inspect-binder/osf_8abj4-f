#==== EXPERIMENT 2 ANALYSES ====

library(BayesFactor)
library(rstan)
library(parallel)
library(ggplot2)
library(gridExtra)

#==== Learning Phase Analyses ====

accdat <- read.table('Learning_data.csv', header=T, sep=',')

#Generate Table 1: Identification accuracies by skew condition 
tab1 <- round(table(subset(accdat, Skew=='Pos')$CorrResp,
                    subset(accdat, Skew=='Pos')$Resp)/224, 2)
tab1
apply(tab1, 1, sum)

tab2 <- round(table(subset(accdat, Skew=='Neg')$CorrResp,
                    subset(accdat, Skew=='Neg')$Resp)/216, 2)
tab2
apply(tab1, 1, sum)

acc_by_blk <- aggregate(Correct ~ Block*Subject, data=accdat, FUN=mean)

#Bayesian ANOVA for accuracy across blocks
acc_by_blk$Block <- factor(acc_by_blk$Block)
acc_by_blk$Subject <- factor(acc_by_blk$Subject)

set.seed(543298)
modBF <- anovaBF(Correct ~ Subject + Block,
                 data=acc_by_blk,
                 whichModels = 'withmain',
                 whichRandom = 'Subject',
                 iterations = 10000)
modBF <- sort(modBF, decreasing = T)
modBF

bestmod <- modBF[1]
chains <- posterior(bestmod, iterations=10000)
plot(chains)
dim(chains)
colnames(chains)

#Polynomial contrasts for each anchor
poly <- contr.poly(8)

blk1 <- apply(chains[,c(1,2)], 1, sum)
blk2 <- apply(chains[,c(1,3)], 1, sum)
blk3 <- apply(chains[,c(1,4)], 1, sum)
blk4 <- apply(chains[,c(1,5)], 1, sum)
blk5 <- apply(chains[,c(1,6)], 1, sum)
blk6 <- apply(chains[,c(1,7)], 1, sum)
blk7 <- apply(chains[,c(1,8)], 1, sum)
blk8 <- apply(chains[,c(1,9)], 1, sum)
blk_mat <- cbind(blk1,blk2,blk3,blk4,blk5,blk6,blk7,blk8)
quantile(blk_mat %*% poly[,'.L'], prob=c(.025,.975))  
quantile(blk_mat %*% poly[,'.Q'], prob=c(.025,.975))  
quantile(blk_mat %*% poly[,'.C'], prob=c(.025,.975))  
quantile(blk_mat %*% poly[,'^4'], prob=c(.025,.975))  
quantile(blk_mat %*% poly[,'^5'], prob=c(.025,.975))  
quantile(blk_mat %*% poly[,'^6'], prob=c(.025,.975))  
quantile(blk_mat %*% poly[,'^7'], prob=c(.025,.975))  

aggregate(Correct ~ Block, data=acc_by_blk, FUN=mean)

#Get participant accuracies in Blocks 5-8
acc_2ndhalf <- aggregate(Correct ~ Subject, data=subset(acc_by_blk, Block %in% c(5:8)), FUN=mean)
acc_2ndhalf$Z <- scale(acc_2ndhalf$Correct)
acc_2ndhalf$low <- ifelse(acc_2ndhalf$Correct < .50, 1, 0) 

#==== Named-Based Reproduction Task (LTM) ====

dat <- read.table('LTM_data.csv', header=T, sep=',')
dat$id <- rep(1:55, each=27)

dat$Actual <- numeric(nrow(dat))
dat$Actual[dat$Skew=='Pos'] <- dat$PSize[dat$Skew=='Pos'] 
dat$Actual[dat$Skew=='Neg'] <- dat$NSize[dat$Skew=='Neg']
dat$Bias <- dat$Est - dat$Actual

dat$Asize <- numeric(nrow(dat))
dat$Asize[dat$Cond <= 4] <- 20
dat$Asize[dat$Cond > 4] <- 180
dat$Anchor <- character(nrow(dat))
dat$Anchor[dat$Asize==20] <- 'Small Anchor'
dat$Anchor[dat$Asize==180] <- 'Large Anchor'

dat1 <- subset(dat, abs(scale(Bias)) < 3)  

#==== Bayesian ANOVA for Name-Based Reproduction ====

dat2 <- aggregate(Bias ~ Skew*Anchor*Actual*Subject, data=dat1, FUN=mean)

set.seed(3458934)

dat2$Skew <- factor(dat2$Skew)
dat2$Actual <- factor(dat2$Actual)
dat2$Anchor <- factor(dat2$Anchor)
dat2$Subject <- factor(dat2$Subject)

#Only analyze trials with common square sizes
dat3 <- subset(dat2, Actual %in% c(40,70,100,130,160))

modBF <- anovaBF(Bias ~ Subject + Skew*Anchor*Actual,
                 data=dat3,
                 whichModels = 'withmain',
                 whichRandom = 'Subject',
                 iterations = 50000)
modBF <- sort(modBF, decreasing = T)
modBF

modBF[1]/modBF[2]  #BF for Skew x Anchor x Stim (and Skew x Anchor) interaction
modBF[1]/modBF[4]  #BF for Anchor x Stim
modBF[1]/modBF[7]  #BF for Skew x Stim
modBF[1]/modBF[12] #BF for Stimulus effect
modBF[1]/modBF[6]  #BF for Anchor effect
modBF[1]/modBF[11] #BF for Skew effect

#Extract posteriors of effects from Bayesian ANOVA
bestmod <- modBF[1]
chains <- posterior(bestmod, iterations=50000, columnFilter="^Subject$")
colnames(chains)

#Small anchor/Neg skew means across square sizes
smN40 <- apply(chains[,c(1,2,5,6,12,15,30,40)], 1, sum)
smN70 <- apply(chains[,c(1,2,5,7,12,16,31,41)], 1, sum)
smN100 <- apply(chains[,c(1,2,5,8,12,17,32,42)], 1, sum)
smN130 <- apply(chains[,c(1,2,5,9,12,18,33,43)], 1, sum)
smN160 <- apply(chains[,c(1,2,5,10,12,19,34,44)], 1, sum)

#Small anchor/Pos skew means across square sizes
smP40 <- apply(chains[,c(1,3,5,6,14,20,30,50)], 1, sum)
smP70 <- apply(chains[,c(1,3,5,7,14,21,31,51)], 1, sum)
smP100 <- apply(chains[,c(1,3,5,8,14,22,32,52)], 1, sum)
smP130 <- apply(chains[,c(1,3,5,9,14,23,33,53)], 1, sum)
smP160 <- apply(chains[,c(1,3,5,10,14,24,34,54)], 1, sum)

#Large anchor/Neg skew means across square sizes
lgN40 <- apply(chains[,c(1,2,4,6,11,15,25,35)], 1, sum)
lgN70 <- apply(chains[,c(1,2,4,7,11,16,26,36)], 1, sum)
lgN100 <- apply(chains[,c(1,2,4,8,11,17,27,37)], 1, sum)
lgN130 <- apply(chains[,c(1,2,4,9,11,18,28,38)], 1, sum)
lgN160 <- apply(chains[,c(1,2,4,10,11,19,29,39)], 1, sum)

#Large anchor/Pos skew means across square sizes
lgP40 <- apply(chains[,c(1,3,4,6,13,20,25,45)], 1, sum)
lgP70 <- apply(chains[,c(1,3,4,7,13,21,26,46)], 1, sum)
lgP100 <- apply(chains[,c(1,3,4,8,13,22,27,47)], 1, sum)
lgP130 <- apply(chains[,c(1,3,4,9,13,23,28,48)], 1, sum)
lgP160 <- apply(chains[,c(1,3,4,10,13,24,29,49)], 1, sum)

poly <- contr.poly(5)

#Polynomial contrasts for Small/Neg
sn_mat <- cbind(smN40,smN70,smN100,smN130,smN160)
quantile(sn_mat %*% poly[,'.L'], prob=c(.025,.975))  
quantile(sn_mat %*% poly[,'.Q'], prob=c(.025,.975))  
quantile(sn_mat %*% poly[,'.C'], prob=c(.025,.975))  
quantile(sn_mat %*% poly[,'^4'], prob=c(.025,.975))  

#Polynomial contrasts for Small/Pos
sp_mat <- cbind(smP40,smP70,smP100,smP130,smP160)
quantile(sp_mat %*% poly[,'.L'], prob=c(.025,.975))  
quantile(sp_mat %*% poly[,'.Q'], prob=c(.025,.975))  
quantile(sp_mat %*% poly[,'.C'], prob=c(.025,.975))  
quantile(sp_mat %*% poly[,'^4'], prob=c(.025,.975))  

#Polynomial contrasts for Large/Neg
ln_mat <- cbind(lgN40,lgN70,lgN100,lgN130,lgN160)
quantile(ln_mat %*% poly[,'.L'], prob=c(.025,.975))  
quantile(ln_mat %*% poly[,'.Q'], prob=c(.025,.975))   
quantile(ln_mat %*% poly[,'.C'], prob=c(.025,.975))  
quantile(ln_mat %*% poly[,'^4'], prob=c(.025,.975))  

#Polynomial contrasts for Large/Pos
lp_mat <- cbind(lgP40,lgP70,lgP100,lgP130,lgP160)
quantile(lp_mat %*% poly[,'.L'], prob=c(.025,.975))  
quantile(lp_mat %*% poly[,'.Q'], prob=c(.025,.975))  
quantile(lp_mat %*% poly[,'.C'], prob=c(.025,.975))  
quantile(lp_mat %*% poly[,'^4'], prob=c(.025,.975))  

#Difference in quadratic trends for each skew condition (small anchor) 
smQuadDiff <- (sn_mat %*% poly[,'.Q']) - (sp_mat %*% poly[,'.Q'])

#Difference in quadratic trends for each skew condition (large anchor)
lgQuadDiff <- (ln_mat %*% poly[,'.Q']) - (lp_mat %*% poly[,'.Q'])

#Difference in quadratic trends reliably larger for large anchor than for small anchor
quantile(lgQuadDiff - smQuadDiff, prob=c(.025, .975))

#==== RF Model: Hierarchical Bayesian estimation ====

dat1$distmean <- numeric(nrow(dat1))
ids <- as.numeric(levels(as.factor(dat1$Subject)))

pos_sizes <- c(40,50,60,70,80,90,100,130,160)
neg_sizes <- c(40,70,100,110,120,130,140,150,160)
pos_ranks <- matrix(c(pos_sizes, (rank(pos_sizes)-1)/(length(pos_sizes)-1)), ncol=2)
neg_ranks <- matrix(c(neg_sizes, (rank(neg_sizes)-1)/(length(neg_sizes)-1)), ncol=2)

dat1$distmean[dat1$Skew=='Neg'] <- mean(neg_sizes)
dat1$distmean[dat1$Skew=='Pos'] <- mean(pos_sizes)

dat1$freqval <- numeric(nrow(dat1))
for(i in 1:nrow(dat1)){
  if(dat1$Skew[i]=='Neg') dat1$freqval[i] <- neg_ranks[which(neg_ranks[,1]==dat1$Actual[i]), 2]
  if(dat1$Skew[i]=='Pos') dat1$freqval[i] <- pos_ranks[which(pos_ranks[,1]==dat1$Actual[i]), 2]
}

dat1$prev <- numeric(nrow(dat1))
for(i in ids){
  dat1$prev[dat1$Subject==i] <- c(dat1$Est[dat1$Subject==i][1], 
                                  dat1$Est[dat1$Subject==i][1:(length(dat1$Est[dat1$Subject==i])-1)])
}

dat2 <- subset(dat1, Trial != 1)

stim <- dat2$Actual/10
bias <- dat2$Bias/10
distmean <- dat2$distmean/10
anchor <- dat2$Asize/10
prev <- dat2$prev/10
freq <- dat2$freqval
id <- dat2$id
Nobs <- 1419
Nind <- 55
smin = 4
smax = 16

no_cores <- detectCores()-1

fit.HBM <- stan(file='RF_full.stan',
                data=c('Nobs','Nind','stim','bias','distmean','anchor','prev','freq','id','smin','smax'),
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
quantile(pars$wA_mu, prob=c(.025, .975))
quantile(pars$wP_mu, prob=c(.025, .975))

pred_mat <- pars$pred_bias
preds <- apply(pred_mat, 2, mean)
plot(bias, preds)
abline(0,1)
cor.test(bias, preds)
dat2$pred_bias <- preds
dat2$Bias10 <- dat2$Bias/10

#Plotting individual fits
dat3 <- aggregate(cbind(Bias10, pred_bias) ~ Skew*Anchor*Actual*Subject, data=dat2, FUN=mean)
dat3a <- aggregate(Bias10 ~ Skew*Anchor*Actual*Subject, data=dat2, 
                   FUN=function(x)sd(x)/sqrt(length(x)))
dat3$ll <- dat3$Bias10 - dat3a$Bias10
dat3$ul <- dat3$Bias10 + dat3a$Bias10

dat3$Size10 <- dat3$Actual/10
dat3$Skew <- factor(dat3$Skew, levels=c('Pos','Neg'), labels=c('Positive','Negative'))
ggplot(data=dat3, aes(x=Size10)) +
  geom_pointrange(aes(y=Bias10, ymin=ll, ymax=ul, pch=Skew), fatten=2) +
  geom_line(aes(y=pred_bias, col=Anchor), lwd=1.2, alpha=.6) +
  facet_wrap(~Subject, nrow=5) +
  scale_shape_manual(values=c(1,2)) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        strip.text=element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  geom_abline(aes(slope=0, intercept=0), lty=3, alpha=.5) +
  xlab('Square Size') +
  ylab('Reproduction Bias') +
  guides(col=guide_legend(override.aes = list(alpha=1)))

#Comparing posterior means for c, b, and wR parameters 
c_mat <- pars$c
b_mat <- pars$b
wR_mat <- pars$wR

inddat <- aggregate(id ~ Skew*Anchor*Subject, data=dat2, FUN=mean)
inddat$c <- apply(c_mat, 2, mean)
inddat$b <- apply(b_mat, 2, mean)
inddat$wR <- apply(wR_mat, 2, mean)

inddat$Skew <- factor(inddat$Skew)
inddat$Anchor <- factor(inddat$Anchor)
inddat$Subject <- factor(inddat$Subject)

set.seed(2347239)
mod_c <- anovaBF(c ~ Skew*Anchor,
                 data=inddat,
                 whichModels = 'withmain',
                 iterations = 50000)
mod_c <- sort(mod_c, decreasing = T)
mod_c

mod_b <- anovaBF(b ~ Skew*Anchor,
                 data=inddat,
                 whichModels = 'withmain',
                 iterations = 50000)
mod_b <- sort(mod_b, decreasing = T)
mod_b

mod_b[1]/mod_b[4]

#The b parameter tends to be larger in the large anchor condition
pars.mod_b <- posterior(mod_b[1], iterations=50000)
colnames(pars.mod_b)
quantile(pars.mod_b[,'Skew-Pos']-pars.mod_b[,'Skew-Neg'], prob=c(.025, .975))
quantile(pars.mod_b[,'Anchor-Large Anchor']-pars.mod_b[,'Anchor-Small Anchor'], prob=c(.025, .975))

mod_wR <- anovaBF(wR ~ Skew*Anchor,
                  data=inddat,
                  whichModels = 'withmain',
                  iterations = 50000)
mod_wR <- sort(mod_wR, decreasing = T)
mod_wR

pars.mod_wR <- posterior(mod_wR[1], iterations=50000)
quantile(pars.mod_wR[,'Skew-Pos']-pars.mod_wR[,'Skew-Neg'], prob=c(.025, .975))

#==== Immediate Reproduction Task (STM) ====

dat <- read.table('STM_data.csv', header=T, sep=',')
dat$id <- rep(1:55, each=27)

dat$Actual <- numeric(nrow(dat))
dat$Actual[dat$Skew=='Pos'] <- dat$Psize[dat$Skew=='Pos'] 
dat$Actual[dat$Skew=='Neg'] <- dat$Nsize[dat$Skew=='Neg']
dat$Bias <- dat$Est - dat$Actual

dat$Asize <- numeric(nrow(dat))
dat$Asize[dat$Cond <= 4] <- 20
dat$Asize[dat$Cond > 4] <- 180
dat$Anchor <- character(nrow(dat))
dat$Anchor[dat$Asize==20] <- 'small'
dat$Anchor[dat$Asize==180] <- 'large'

dat1 <- subset(dat, abs(scale(Bias)) < 3)  

#==== Bayesian ANOVA for Immediate Reproduction ====

dat2 <- aggregate(Bias ~ Skew*Anchor*Actual*Subject, data=dat1, FUN=mean)

set.seed(346534)

dat2$Skew <- factor(dat2$Skew)
dat2$Actual <- factor(dat2$Actual)
dat2$Anchor <- factor(dat2$Anchor)
dat2$Subject <- factor(dat2$Subject)

#Only analyze trials with common square sizes
dat3 <- subset(dat2, Actual %in% c(40,70,100,130,160))

modBF <- anovaBF(Bias ~ Subject + Skew*Anchor*Actual,
                 data=dat3,
                 whichModels = 'withmain',
                 whichRandom = 'Subject',
                 iterations = 50000)
modBF <- sort(modBF, decreasing = T)
modBF

modBF[1]/modBF[15]  #BF for Stim effect

bestmod <- modBF[1]
chains <- posterior(bestmod, iterations=50000, columnFilter="^Subject$")
colnames(chains)

#Polynomial contrasts for Stim effect
poly <- contr.poly(5)

stim40 <- apply(chains[,c(1,2)], 1, sum)
stim70 <- apply(chains[,c(1,3)], 1, sum)
stim100 <- apply(chains[,c(1,4)], 1, sum)
stim130 <- apply(chains[,c(1,5)], 1, sum)
stim160 <- apply(chains[,c(1,6)], 1, sum)
stim_mat <- cbind(stim40,stim70,stim100,stim130,stim160)
quantile(stim_mat %*% poly[,'.L'], prob=c(.025,.975))  
quantile(stim_mat %*% poly[,'.Q'], prob=c(.025,.975))  
quantile(stim_mat %*% poly[,'.C'], prob=c(.025,.975))  
quantile(stim_mat %*% poly[,'^4'], prob=c(.025,.975))  

#==== CA Model for Immediate Reproduction Bias =====

dat1$Running <- numeric(nrow(dat1))
ids <- as.numeric(levels(as.factor(dat1$Subject)))

for(i in ids){
  dat1$Running[dat1$Subject==i] <- cumsum(dat1$Actual[dat1$Subject==i])/(1:length(dat1$Actual[dat1$Subject==i]))
}

dat1$Back1 <- numeric(nrow(dat1))
for(i in ids){
  dat1$Back1[dat1$Subject==i] <- c(dat1$Actual[dat1$Subject==i][1], 
                                   dat1$Actual[dat1$Subject==i][1:(length(dat1$Actual[dat1$Subject==i])-1)])
}


dat2 <- subset(dat1, Trial != 1)

stim <- dat2$Actual/10
bias <- dat2$Bias/10
RM <- dat2$Running/10
AS <- dat2$Asize/10
B1 <- dat2$Back1/10
id <- dat2$id
Nobs <- 1419
Nind <- 55
mid <- 10
stim_sizes <- seq(40,160,10)/10
Nstim <- length(stim_sizes)

library(parallel)
no_cores <- detectCores()-1

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

#Posteriors for group mean parameters of interest
quantile(pars$wRM_mu, prob=c(.025, .975))
quantile(pars$wAS_mu, prob=c(.025, .975))
quantile(pars$wB1_mu, prob=c(.025, .975))

#Get trial-by-trial predictions
pred_mat <- pars$pred_bias
preds <- apply(pred_mat, 2, mean)
plot(bias, preds)
abline(0,1)
cor.test(bias, preds)
dat2$pred_bias <- preds
dat2$Bias10 <- dat2$Bias/10

#Plotting individual fits
dat3 <- aggregate(cbind(Bias10, pred_bias) ~ Skew*Anchor*Actual*Subject, data=dat2, FUN=mean)
dat3a <- aggregate(Bias10 ~ Skew*Anchor*Actual*Subject, data=dat2, 
                   FUN=function(x)sd(x)/sqrt(length(x)))
dat3$ll <- dat3$Bias10 - dat3a$Bias10
dat3$ul <- dat3$Bias10 + dat3a$Bias10

dat3$Size10 <- dat3$Actual/10
dat3$Skew <- factor(dat3$Skew, levels=c('Pos','Neg'), labels=c('Positive','Negative'))
ggplot(data=dat3, aes(x=Size10)) +
  geom_pointrange(aes(y=Bias10, ymin=ll, ymax=ul, pch=Skew), fatten=2) +
  geom_line(aes(y=pred_bias, col=Anchor), lwd=1.2, alpha=.6) +
  facet_wrap(~Subject, nrow=5) +
  scale_shape_manual(values=c(1,2)) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        strip.text=element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  geom_abline(aes(slope=0, intercept=0), lty=3, alpha=.5) +
  xlab('Square Size') +
  ylab('Reproduction Bias') +
  guides(col=guide_legend(override.aes = list(alpha=1)))
