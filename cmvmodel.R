##### Stan code for CMV transmission model #####

setwd("mypath")
library(readr)
library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

### import CMV data: Titers, Sexes, Ages ###
cmvdata <- read_csv("data/cmvdata.csv")

# filter non-western subjects
cmvdata <- cmvdata[cmvdata$nl == '1',]

# filter ages <6 months
cmvdata <- cmvdata[cmvdata$lftinmnd2 > 6,]

Ages <- cmvdata$LFTB
Titers <- cmvdata$boxcoxvalues
N <- length(Ages)

### deal with censoring and spike ###
RightCensor <- 3.41
Spike <- -2.91
Censor <- rep(0,N)
for (i in 1:N){
  if (Titers[i] > RightCensor){
    Censor[i] <- 1
  }
  else if (Titers[i] <= Spike){
    Censor[i] <- 2
  }
}

### distinguish gender ###
Gender <- rep(0,N)
Gender[cmvdata$gesl2 == 'male'] <- 1 # female = 0

### import Contact Matrix ###
ContactData <- read_csv("data/contact_intensities_aggregated.csv")
MM <- ContactData$mMM # length 256 = 16*16
FM <- ContactData$mFM
MF <- ContactData$mMF
FF <- ContactData$mFF
Contact_MM <- matrix(MM,nrow=16,ncol=16)
Contact_FM <- matrix(FM,nrow=16,ncol=16)
Contact_MF <- matrix(MF,nrow=16,ncol=16)
Contact_FF <- matrix(FF,nrow=16,ncol=16)

# some values
DeltaA <- 5
A <- 16
A_rho <- 3
S0 <- 0.81
RhoClasses <- c(1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3)
mode <- 0 ### 0=normal sampling, 1=WBIC sampling

# prepare a data dictionary an initial values for Stan
DataDict <- list('N'= N,'A'= A,'A_rho'= A_rho,'DeltaA'= DeltaA,'Titers'= Titers,
  'Ages'= Ages,'Censor'= Censor,'RhoClasses'= RhoClasses,'RightCensor'= RightCensor,
  'Contact_MM'= Contact_MM,'Contact_FM'= Contact_FM,'Contact_MF'= Contact_MF,
  'Contact_FF'= Contact_FF,'MuS'= -1.68,'MuL'= 1.03,'MuB'= 2.6,'SigmaS'= 1/6.27, 
  'SigmaL'= 1/1.52,'S0'= S0,'SigmaB'= 1/2.03,'Penalty'= 1e4,'Gender'= Gender,
  'mode'= mode) 

# initial values for parameters to be estimated
initials = function(){
  return(list(beta1=0.0024,beta2=0.073, z=0.41,
             lambda_f=c(0.015,0.019,0.021,0.016,0.011,0.012,0.013,0.013,0.013,0.012,0.011,0.009,0.0081,0.00083,0.0071,0.0058),
             lambda_m=c(0.014,0.019,0.017,0.013,0.0091,0.01,0.012,0.014,0.013,0.011,0.0094,0.0083,0.009,0.0076,0.0066,0.0065),
             shortRho_f=c(0.015,0.023,0.023),
             shortRho_m=c(0.0051,0.011,0.012),
             muRho=0.015,sigmaRho=0.0093
            )
        )
}

# run the Stan model
fitW <- stan(file = 'cmvmodel (15062017).stan', data = DataDict, init = initials, iter = 3000, 
                   warmup = 500, thin=5, chains = 10, control = list(adapt_delta = 0.99))

# some parameter outputs 
print(fit, pars=c("beta1", "beta2", "z",
  "shortRho_m[1]", "shortRho_m[2]", "shortRho_m[3]",
  "shortRho_f[1]", "shortRho_f[2]", "shortRho_f[3]", 
  "muRho", "sigmaRho"), digits=5)

# calculate WAIC
LL <- extract_log_lik(fit, parameter_name = 'log_lik')
waic(LL)

# write selected output
output <- as.data.frame(fit, pars=c("beta1", "beta2", "z",
  "shortRho_m[1]", "shortRho_m[2]", "shortRho_m[3]",
  "shortRho_f[1]", "shortRho_f[2]", "shortRho_f[3]", 
  "muRho", "sigmaRho", "lambda_m", "lambda_f",
  "S_m", "S_f", "L_m", "L_f", "B_m", "B_f")) 
write.csv(output, file = "results/myoutput.csv", row.names = FALSE)

# calculate WBIC; note set mode=1 in separate run
print(fit, pars="log_like", digits=5) 
