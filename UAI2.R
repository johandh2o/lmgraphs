###############################################################################
# Libraries

rm(list = ls())    # Clear environment
library(MASS)      # Miscellaneous
library(mice)      # For imputations
library(twang)     # For AIPW
library(lmtest)    # For confidence intervals
library(margins)   # For marginal effects
library(ggplot2)   # For plots
library(tidyr)     # For pivoting tables
library(dplyr)     # For transforming data

################################################################################
# Auxiliary functions

# Function to compute logit for a given probability p
logit <- function(p) {
  return(log(p / (1 - p)))
}

###############################################################################
# Vessels

# Oracle vessels
oracleFATE = c() # Oracle sample FATE
oracleNATE = c() # Oracle sample NATE

# Estimators
imputation = c() # Imputation estimatot
zeroed = c()     # Missing indicator (MIM) estimator
mineFATE = c()   # Computed FATE
mineNATE = c()   # Computed NATE

###############################################################################
# Causal mechanisms

# Exposure mechanism
generate_A <- function(W, X, R, param_set1, param_set2) {
  logit_p <- R * (param_set1[1] + param_set1[2] * W + param_set1[3] * X) + 
    (1-R) * (param_set2[1] + param_set2[2] * W)
  p <- 1 / (1 + exp(-logit_p))
  return(rbinom(length(W), 1, p))
}

# Outcome mechanism
generate_Y <- function(W, X, R, A) {
  return( R*(3 + (7 * W - 8 * A - 6 * X - 3 * A * W + 4 * A * X )/4) +
            (1-R)*(4 + (6 * W + 8 * A - 8 * W * A)))
}

# Function to generate R using a logit model
generate_R <- function(W, param_set) {
  logit_p <- param_set[1] + param_set[2] * W 
  p <- 1 / (1 + exp(-logit_p))
  return(rbinom(length(W), 1, p))
}

################################################################################

# Manually adjust parameters to achieve desired proportions for R
param_sets_R <- list(
  c(-4.5, 3)/2.5,   # 20%
  c(-2, 2)/2.5,     # 35%
  c(-1, 4)/2.5,     # 40%
  c(0.5, 3)/2.5     # 50%
)

################################################################################
# Data generation

# Parameters
N <- 5e3  # Number of observations
M <- 200  # Number of replications

# Data generation in loop
# Generate M number or replications for the computation of the estimators
# For each of the K scenarios
for(K in 4:1){
for(J in 1:M){
  
  # Set seed  
  set.seed(1000*K+J)
  
  # Generate W from Uniform(0, 1)
  W <- rnorm(N, mean = 0, sd=1)
  
  # Generate X as Gaussian with mean as a quadratic function of W 
  mean_X <- (1 * W^2 - 2 * W - 3)       # Quadratic function of W
  X <- rnorm(N, mean = mean_X, sd = 7)  # Gaussian noise
  
  # Generate R with the new parameter sets
  MissingV <- sapply(param_sets_R, function(param_set) generate_R(W, param_set))
  
  # Check the proportions of R for each case
  colMeans(MissingV)
  
  ##############################################################################
  # First scenario
  
  R = 1-MissingV[,K]
  #summary(lm(R ~ W)) 
  #mean(R)
  
  param_set_A_1 <- c(-0, 2, 0.5)/2
  param_set_A_2 <- c(-1, 3)/2
  A_R <- generate_A(W, X, R, param_set_A_1, param_set_A_2)
  
  # Generate noise uT[R == 1] for R = 1
  uT <- rnorm(N, mean = 0, sd = 7)
  
  # Generate counterfactuals Y1.0 and Y1.1 for R = 1
  YF_0 <- generate_Y(W, X, 1, 0) + uT
  YF_1 <- generate_Y(W, X, 1, 1) + uT
  
  YN_0 <- generate_Y(W, X, R, 0) + uT
  YN_1 <- generate_Y(W, X, R, 1) + uT
  
  # Dataframe
  dataF = data.frame(W, X, R=R, A=A_R,
                     Y=YN_1*A_R + (1-A_R)*YN_0, 
                     dF=YF_1-YF_0, dN=YN_1-YN_0)
  
  #mean(A_R[R==0])
  #mean(A_R[R==1])
  
  #summary(lm(A ~  W + X, data = dataF, subset = (R == 0))) 
  #summary(lm(Y ~ W + A + X + A*W + A*X, data = dataF, subset = (R == 0)))
  
  #summary(lm(A ~  W, data = dataF, subset = (R == 0))) 
  #summary(lm(Y ~ W + A + A:W, data = dataF, subset = (R == 0)))
  
  oracleFATE = c(oracleFATE , mean(dataF$dF))
  oracleNATE = c(oracleNATE , mean(dataF$dN))
  
  #############################################################################
  # estimators: Imputation
  
  impudata = dataF
  impudata$X[impudata$R==0] <- NA
  impudata = impudata[,c(1,2,3,4,5)]
  
  Im = 10
  imp_data <- mice(impudata, m = Im, method = 'pmm',
                   maxit = 5, seed = 123, printFlag = FALSE)
  
  cumsuma = rep(0,6)
  
  # Generate multiple imputations
  for(k in 1:Im){
    completed_data <- complete(imp_data, k)  # Using the first imputed dataset
    
    # Run regression with interactions between A, W, and X
    model1 <- lm(Y ~ W + A + X + A*W + A*X, data = completed_data)
    
    # Compute the marginal effects of A on Y
    marginal1 <- margins(model1, variables = "A")
    
    cumsuma = cumsuma + summary(marginal1)[2:6]
  }
  finsuma = cumsuma / Im
  finsuma
  
  imputation = c(imputation, as.numeric(finsuma$AME))
  
  ##############################################################################
  # estimators: Missingness Indicator Method (MIM)
  
  data_zeroed <- impudata
  data_zeroed$X[is.na(data_zeroed$X)] <- 0
  
  # Run regression with interactions between A, W, and X (with X=0 where missing)
  model2 <- lm(Y ~ W + A + X + A:W + A:X + A:R + A:W, data = data_zeroed)
  marginal2 <- margins(model2, variables = "A")
  #summary(marginal2)
  zeroed = c(zeroed, as.numeric(summary(marginal2)$A))
  
  #############################################################################
  # estimators: AIPW FATE
  
  dataR1 = subset(dataF, R == 1)
  
  # Step 1: Fit separate outcome models
  outcome_model_R1 <- lm(Y ~ W + A + X + A:W + A:X, data = dataF, subset = (R == 1))
  
  # Step 2: Fit separate propensity score models
  ps_model_R1 <- glm(A ~ X + W, data = dataF, subset = (R == 1), family = binomial())
  
  # Step 3: Predict treatment probabilities
  ps_A_R1 <- predict(ps_model_R1, type = "response")
  
  # Step 4: Compute AIPW estimates for both R=1 and R=0 cases
  dummyData = dataR1
  
  dummyData$A = 1
  mu_hat_Y1.R1 <- predict(outcome_model_R1, newdata = dummyData)  
  
  dummyData$A = 0
  mu_hat_Y0.R1 <- predict(outcome_model_R1, newdata = dummyData)
  
  aipw_R1 = (mu_hat_Y1.R1 - mu_hat_Y0.R1) +
    (dataR1$A / ps_A_R1  - (1 - dataR1$A) / (1 - ps_A_R1)) * (dataR1$Y - (dataR1$A * mu_hat_Y1.R1 + (1-dataR1$A) * mu_hat_Y0.R1))
  
  dummyData = dataF
  dummyData$aipwPseudo = NA
  dummyData$aipwPseudo[R==1] = aipw_R1
  
  # Step 5: Pseudo-outcome regression
  outcome_extra <- lm(aipwPseudo ~ W, data = dummyData, subset = (R == 1))
  ps_model_extra <- glm(R ~ W, data = dataF, family = binomial())
  
  out.pre = predict(outcome_extra, newdata=dataF)
  ps.pre = predict(ps_model_extra, type = "response")
  
  aipw_Pre = mean(out.pre) +
    mean((dataF$R / ps.pre) * (dummyData$aipwPseudo - out.pre), na.rm = TRUE)
  
  mineFATE = c( mineFATE, aipw_Pre)

  #############################################################################
  # estimators: AIPW NATE
  
  dataR1 = subset(dataF, R == 1)
  dataR0 = subset(dataF, R == 0)
  
  # Step 1: Fit separate outcome models
  outcome_model_R1 <- lm(Y ~ W + A + X + A:W + A:X, data = dataF, subset = (R == 1))
  outcome_model_R0 <- lm(Y ~ W + A + A:W, data = dataF, subset = (R == 0)) 
  
  # Step 2: Fit separate propensity score models
  ps_model_R1 <- glm(A ~ X + W, data = dataF, subset = (R == 1), family = binomial())
  ps_model_R0 <- glm(A ~ W, data = dataF, subset = (R == 0), family = binomial())
  
  # Step 3: Predict treatment probabilities
  ps_A_R1 <- predict(ps_model_R1, type = "response")
  ps_A_R0 <- predict(ps_model_R0, type = "response")
  
  # Step 4: Compute AIPW estimates for both R=1 and R=0 cases
  
  # Stratum R=1
  dummyData = dataR1
  
  dummyData$A = 1
  mu_hat_Y1.R1 <- predict(outcome_model_R1, newdata = dummyData)  
  
  dummyData$A = 0
  mu_hat_Y0.R1 <- predict(outcome_model_R1, newdata = dummyData)
  
  # Stratum R=0
  dummyData = dataR0
  
  dummyData$A = 1
  mu_hat_Y1.R0 <- predict(outcome_model_R0, newdata = dummyData)
  
  dummyData$A = 0
  mu_hat_Y0.R0 <- predict(outcome_model_R0, newdata = dummyData)
  
  # Compute stratum-specific estimates
  aipw_R1 = (mu_hat_Y1.R1 - mu_hat_Y0.R1) +
    (dataR1$A / ps_A_R1  - (1 - dataR1$A) / (1 - ps_A_R1)) * 
    (dataR1$Y - (dataR1$A * mu_hat_Y1.R1 + (1-dataR1$A) * mu_hat_Y0.R1))
  
  aipw_R0 = (mu_hat_Y1.R0 - mu_hat_Y0.R0) +
    (dataR0$A / ps_A_R0  - (1 - dataR0$A) / (1 - ps_A_R0)) * 
    (dataR0$Y - (dataR0$A * mu_hat_Y1.R0 + (1-dataR0$A) * mu_hat_Y0.R0))
  
  # Step 5: Marginalize R
  model4_effect <- mean(R) * mean(aipw_R1) + (1-mean(R)) * mean(aipw_R0)
  mineNATE = c(mineNATE, model4_effect)
  
  print(1000*K+J)
}}

# Final tables
final.table = data.frame(
  Imput = imputation,      # Multiple imputation estimator
  MIM  = zeroed,           # Missing indicator method
  FATE = mineFATE,         # Computed FATE
  NATE = mineNATE,         # Computed NATE
  K = c(rep(4,M),rep(3,M),rep(2,M),rep(1,M))
)

pivot.table <- pivot_longer(
  final.table,
  cols = c(Imput, MIM, FATE, NATE),
  names_to = "Estimator",
  values_to = "Value"
)

# Plot
mr <- c(20, 35, 40, 50)
label_map <- setNames(
  paste0("Missing rate ", mr, "%"),
  c(1, 2, 3, 4)
)

ggplot(pivot.table, aes(x = Estimator, y = Value, fill = Estimator)) +
  geom_boxplot() +
  facet_wrap(~ K, nrow = 2, labeller = labeller(K = label_map)) +
  theme_minimal() +
  labs(
    title = "Distribution of estimators by scenario",
    x = "Estimator",
    y = "Value"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )