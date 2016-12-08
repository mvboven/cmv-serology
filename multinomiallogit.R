#
# Init
#

# Load packages
library(sas7bdat)
library(splines)
library(runjags)
library(coda)

# Set working directory
folder <- "myfolder"
switch(Sys.info()["nodename"],
  "node1" = setwd(file.path("path1", folder)),
  "node2" = setwd(file.path("path2", folder)))

#
# Read and modify CMV data
#

# Read data
cmv.data <- #import your data

# Some modifications
cmv.data <- within(cmv.data, {
  # Make nl and gesl2 categorical variables with informative labels
  nl <- factor(nl, levels = 1:2, labels = c("yes", "no"))
  sex <- factor(gesl2, levels = 1:2, labels = c("male", "female"))
  # Rename average_IU_ml to cmv
  cmv <- average_IU_ml
  rm(average_IU_ml, gesl2)
  # Calcluate age in decimal years from lftinmnd2
  age <- lftinmnd2/12
})

# Remove all records with age <= 6 months & not NL
cmv.data <- subset(cmv.data, subset = lftinmnd2 > 6 & nl == "yes")

# Histogram of log(cmv)
with(cmv.data, hist(log(cmv), breaks = 100))

# We see a spike on the left
# There are right censored observations
# The log-transformation results in a left skewed distribution -> Log transformation is too strong
# -> Apply the Box-Cox transformation instead of a log-transformation
# -> Model the spike separately

# Spike is where cmv = min(cmv), that is where cmv = 1e-04
# Right censored observations are where cmv > max(cmv[cmv < 12])
min.cmv <- with(cmv.data, min(cmv)) # 1e-4
max.cmv <- 10
abline(v = log(c(min.cmv, max.cmv)), col = 2)

# Add variables spike and cens for spike and right censored data
cmv.data <- within(cmv.data, {
  spike <- as.numeric(cmv == min.cmv)
  cens <- as.numeric(cmv > max.cmv)
})

# Check
head(subset(cmv.data, subset = spike == 0 & cens == 0)) # Observed
head(subset(cmv.data, subset = spike == 1 & cens == 0)) # Spike
head(subset(cmv.data, subset = spike == 0 & cens == 1)) # Right censored
head(subset(cmv.data, subset = spike == 1 & cens == 1)) # Empty set

# Apply Box-Cox transformation
lambda <- 0.3
cmv.data <- within(cmv.data, {
  cmv.trans <- (cmv^lambda-1)/lambda
})

# Histogram of box-cox transformed cmv
with(cmv.data, hist(cmv.trans, breaks = 100))

# Histograms by age category
cmv.data.split <- with(cmv.data, split(x = cmv.data, f = cut(age, breaks = seq(0, 80, length = 9+1), include.lowest = TRUE)))
par(mar = c(2.1, 2.1, 0.1, 0.1), mfrow = c(3, 3))
lapply(cmv.data.split, function(x) {
  with(x, hist(cmv.trans, 
    prob = TRUE,
    col = "grey", breaks = 40,
    xlim = c((min.cmv^lambda-1)/lambda, (max.cmv^lambda-1)/lambda),
    ylim = c(0, 1), main = NULL))
  box()
})

# Inform components
cmv.data <- within(cmv.data, {
  comp <- ifelse(test = cmv.trans < -1, yes = NA,
    no = ifelse(test = cmv.trans > 0 & cmv.trans < 3 & age < 10, yes = NA,
      no = ifelse(test = cens == 1, yes = NA, no = NA)))
})

# Scatterplot
with(cmv.data, {
  plot(age, cmv.trans, col = 8)
  points(age, cmv.trans, col = comp)
})

#
# Analysis
#

# JAGS model
model.string <- "model {
  
  # Likelihood
  for (i in 1:n) {

    # Observations are a mixture of three normals and a spike
    y[i] ~ dnorm(mu.y[i], tau.y[i])
    mu.y[i]  <- min.y*(comp[i] == 1 && spike[i] == 1) +  mu.neg*(comp[i] == 1 && spike[i] == 0) +  mu.inf*(comp[i] == 2) +  mu.boo*(comp[i] == 3)
    tau.y[i] <-   1e6*(comp[i] == 1 && spike[i] == 1) + tau.neg*(comp[i] == 1 && spike[i] == 0) + tau.inf*(comp[i] == 2) + tau.boo*(comp[i] == 3)
    # Some are censored
    cens[i] ~ dinterval(y[i], max.y)
    
    # Component number has a categorical distribution
    comp[i] ~ dcat(p[i, 1:3])
    p[i, 1] <- p.neg[i]; p[i, 2] <- p.inf[i]; p[i, 3] <- p.boo[i]

    # Multinomial logit with splines
    p.neg[i] <-              1 / (1 + exp.eta.inf[i] + exp.eta.boo[i])
    p.inf[i] <- exp.eta.inf[i] / (1 + exp.eta.inf[i] + exp.eta.boo[i])
    p.boo[i] <- exp.eta.boo[i] / (1 + exp.eta.inf[i] + exp.eta.boo[i])
    log(exp.eta.inf[i]) <- sum(Q[i, 1:n.gamma]*gamma.inf[1:n.gamma])
    log(exp.eta.boo[i]) <- sum(Q[i, 1:n.gamma]*gamma.boo[1:n.gamma])

  }

  # Priors for mixture components
  mu.neg <- mu[1]; mu.inf <- mu[2]; mu.boo <- mu[3]
  mu[1:3] <- sort(mu0[1:3])
  for (j in 1:3) {
    mu0[j] ~ dnorm(0, 0.001)
  }
  tau.neg ~ dgamma(0.5, 0.005)
  tau.inf ~ dgamma(0.5, 0.005)
  tau.boo ~ dgamma(0.5, 0.005)

  # Priors for spline in multinomial logit model
  for (k in 1:n.gamma) {
    gamma.inf[k] ~ dnorm(0, 0.001)
    gamma.boo[k] ~ dnorm(0, 0.001)
  }

}"

# Create model matrix
# We have age from -1 to 81 in order to calculate the derivative at 0 to 80.
cmv.pred <- expand.grid(age = -1:81, sex = c("male", "female"))
X      <- with(cmv.data, model.matrix(~ ns(age, knots = c(20, 40, 60), Boundary.knots = c(0, 80)) + sex))
X.pred <- with(cmv.pred, model.matrix(~ ns(age, knots = c(20, 40, 60), Boundary.knots = c(0, 80)) + sex))

# For better mixing properties, we calculate the QR decomposition of X
# see http://sites.stat.psu.edu/~jls/stat511/lectures/lec17.pdf
# eta = X%*%beta = Q%*%gamma
# gamma = R%*%beta
qr.X <- qr(X); Q <- qr.Q(qr.X); R <- qr.R(qr.X)

# Data list
data.list <- with(cmv.data, list(
  n = nrow(cmv.data),
  y = ifelse(cens == 0, cmv.trans, NA),
  min.y = (min.cmv^lambda-1)/lambda, max.y = (max.cmv^lambda-1)/lambda,
  spike = spike, cens = cens, comp = comp,
  Q = Q, n.gamma = ncol(Q)))

# Inits function
inits.fun <- with(data.list, function() list(
  mu0 = c(-1.5, 1, 3), tau.neg = 1, tau.inf = 1, tau.boo = 1,
  y = ifelse(cens == 0, NA, max.y+0.01),
  gamma.inf = rep(0, n.gamma), gamma.boo = rep(0, n.gamma),
  .RNG.name="base::Mersenne-Twister", .RNG.seed = sample(1:10000, 1)))

# Run Model
n.cores <- 10
post.runjags <- run.jags(model = model.string, data = data.list, inits = inits.fun,
  n.chains = n.cores, adapt = 500, burnin = 1000, sample = round(10000/n.cores), thin = 1, method = "parallel", modules = "glm",
  monitor = c("mu.neg", "mu.inf", "mu.boo", "tau.neg", "tau.inf", "tau.boo", "gamma.inf", "gamma.boo"))

# Traceplots
plot(post.runjags, vars = c("mu.neg", "mu.inf", "mu.boo", "tau.neg", "tau.inf", "tau.boo", "gamma.boo", "gamma.inf"),
  plot.type = "trace", layout = c(4, 4))

#
# Post processing
#

# Extract parameters
post.mat <- as.matrix(as.mcmc(post.runjags))
n.iter <- nrow(post.mat)
gamma.boo <- post.mat[, grep("gamma.boo", colnames(post.mat))]
gamma.inf <- post.mat[, grep("gamma.inf", colnames(post.mat))]

# Calculate beta's
beta.boo <- solve(R)%*%t(gamma.boo)
beta.inf <- solve(R)%*%t(gamma.inf)

# Make predictions
exp.eta.inf <- exp(t(X.pred%*%beta.inf))
exp.eta.boo <- exp(t(X.pred%*%beta.boo))
p.neg.MF <-           1 / (1 + exp.eta.inf + exp.eta.boo)
p.inf.MF <- exp.eta.inf / (1 + exp.eta.inf + exp.eta.boo)
p.boo.MF <- exp.eta.boo / (1 + exp.eta.inf + exp.eta.boo)

# Combine males and females
age <- -1:81; n <- length(age)
p.neg.M <- p.neg.MF[, 1:n]; p.neg.F <- p.neg.MF[, (n+1):(2*n)]; p.neg <- (p.neg.M + p.neg.F)/2
p.inf.M <- p.inf.MF[, 1:n]; p.inf.F <- p.inf.MF[, (n+1):(2*n)]; p.inf <- (p.inf.M + p.inf.F)/2
p.boo.M <- p.boo.MF[, 1:n]; p.boo.F <- p.boo.MF[, (n+1):(2*n)]; p.boo <- (p.boo.M + p.boo.F)/2

#
# Plots
#

# Function to calculate s
polygon.ci <- function(x, boundaries, ...) polygon(c(x, rev(x)), c(boundaries[, 1], rev(boundaries[, 2])), border = NA, ...)

# Plot prevalences
cols <- RColorBrewer::brewer.pal(n = 3, name = "Set1")
cols.trans <- adjustcolor(cols, alpha = 0.5)
# Open figure
pdf(file = "results/prevalences.pdf", width = 29.7/2.54, height = 21/2.54, paper = "a4r")
par(mar = c(4.5, 4.5, 1.5, 0.1), mfrow = c(1, 3))
# Males
plot.new(); plot.window(xlim = c(0, 80), ylim = c(0, 1))
grid(); axis(1); axis(2); box()
title(main = "Males", xlab = "Age", ylab = "Prevalence")
polygon.ci(age, t(apply(p.neg.M, 2, quantile, prob = c(0.025, 0.975))), col = cols.trans[1])
polygon.ci(age, t(apply(p.inf.M, 2, quantile, prob = c(0.025, 0.975))), col = cols.trans[2])
polygon.ci(age, t(apply(p.boo.M, 2, quantile, prob = c(0.025, 0.975))), col = cols.trans[3])
lines(age, apply(p.neg.M, 2, median), col = cols[1])
lines(age, apply(p.inf.M, 2, median), col = cols[2])
lines(age, apply(p.boo.M, 2, median), col = cols[3])
# Males
plot.new(); plot.window(xlim = c(0, 80), ylim = c(0, 1))
grid(); axis(1); axis(2); box()
title(main = "Females", xlab = "Age", ylab = "Prevalence")
polygon.ci(age, t(apply(p.neg.F, 2, quantile, prob = c(0.025, 0.975))), col = cols.trans[1])
polygon.ci(age, t(apply(p.inf.F, 2, quantile, prob = c(0.025, 0.975))), col = cols.trans[2])
polygon.ci(age, t(apply(p.boo.F, 2, quantile, prob = c(0.025, 0.975))), col = cols.trans[3])
lines(age, apply(p.neg.F, 2, median), col = cols[1])
lines(age, apply(p.inf.F, 2, median), col = cols[2])
lines(age, apply(p.boo.F, 2, median), col = cols[3])
# Males
plot.new(); plot.window(xlim = c(0, 80), ylim = c(0, 1))
grid(); axis(1); axis(2); box()
title(main = "Males and Females", xlab = "Age", ylab = "Prevalence")
legend("topright", legend = c("Negative", "Infected", "Boosted"), fill = cols, bty = "n")
polygon.ci(age, t(apply(p.neg, 2, quantile, prob = c(0.025, 0.975))), col = cols.trans[1])
polygon.ci(age, t(apply(p.inf, 2, quantile, prob = c(0.025, 0.975))), col = cols.trans[2])
polygon.ci(age, t(apply(p.boo, 2, quantile, prob = c(0.025, 0.975))), col = cols.trans[3])
lines(age, apply(p.neg, 2, median), col = cols[1])
lines(age, apply(p.inf, 2, median), col = cols[2])
lines(age, apply(p.boo, 2, median), col = cols[3])
# Close figure
dev.off()

# Histograms by age category
post.vec <- apply(post.mat, 2, median)
mu.neg <- post.vec[grep("mu.neg", names(post.vec))]
tau.neg <- post.vec[grep("tau.neg", names(post.vec))]
mu.inf <- post.vec[grep("mu.inf", names(post.vec))]
tau.inf <- post.vec[grep("tau.inf", names(post.vec))]
mu.boo <- post.vec[grep("mu.boo", names(post.vec))]
tau.boo <- post.vec[grep("tau.boo", names(post.vec))]

median[mu.neg]

x <- seq(-3.5, 3.5, length = 101)
par(mar = c(2.1, 2.1, 0.1, 0.1), mfrow = c(3, 3))
lapply(cmv.data.split, function(data) {
  # Histogram
  with(data, hist(ifelse(cens == 0, cmv.trans, NA),
    col = "grey", breaks = 30, prob = TRUE,
    xlim = c((min.cmv^lambda-1)/lambda, (max.cmv^lambda-1)/lambda), ylim = c(0, 0.85),
    main = NULL))
  box()
  # Density
  index <- round(mean(data$age))+1
  d <-
    median(p.neg[, index])*dnorm(x, mu.neg, sqrt(1/tau.neg)) +
    median(p.inf[, index])*dnorm(x, mu.inf, sqrt(1/tau.inf)) +
    median(p.boo[, index])*dnorm(x, mu.boo, sqrt(1/tau.boo))
  lines(x, d, col = 2, lwd = 2)
})

#
# Export
#

# Save matrices
write.csv(x = round(p.neg,   6), file = "results/p_neg.csv", row.names = FALSE)
write.csv(x = round(p.inf,   6), file = "results/p_inf.csv", row.names = FALSE)
write.csv(x = round(p.boo,   6), file = "results/p_boo.csv", row.names = FALSE)
write.csv(x = round(p.neg.M, 6), file = "results/p_neg_M.csv", row.names = FALSE)
write.csv(x = round(p.inf.M, 6), file = "results/p_inf_M.csv", row.names = FALSE)
write.csv(x = round(p.boo.M, 6), file = "results/p_boo_M.csv", row.names = FALSE)
write.csv(x = round(p.neg.F, 6), file = "results/p_neg_F.csv", row.names = FALSE)
write.csv(x = round(p.inf.F, 6), file = "results/p_inf_F.csv", row.names = FALSE)
write.csv(x = round(p.boo.F, 6), file = "results/p_boo_F.csv", row.names = FALSE)
