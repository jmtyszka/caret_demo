# Analyse headcase motion statistics for ASD patients
#
# AUTHOR : Mike Tyszka, Ph.D.
# PLACE  : Caltech
# DATES  : 2017-07-31 JMT From scratch

# Import packages
library(splitstackshape)
library(plyr)
library(doBy)
library(gdata)
library(ggplot2)
library(Metrics)

# Flags
do_scatterplot = FALSE
do_preproc_checks = FALSE

# Regressors to tune
do_mvtb   = FALSE  # Multivariate tree boosting
do_gbm    = TRUE   # Gradient boosted trees
do_lrback = TRUE   # Linear regression with backwards selection

# User selected file from Qualtrics
data_fname <- "Results.txt"

# Read TSV data
X <- read.table(data_fname, sep="\t", header=TRUE, fill=TRUE)

# The variables of interest are:
# (Q3 (Age))
# NIHTB_Emotional_Support
# NIHTB_Friendship
# NIHTB_Loneliness
# NIHTB_Instrumetnal_Support
# NIHTB_Perceived_Hostility
# NIHTB_Perceived_Rejection
# SNI_Ex_People_Count
# SNI_Ex_Embed_Networks
# SNI_Ex_Network_Diversity

myvars <- c('NIHTB_Emotional_Support',
            'NIHTB_Friendship',
            'NIHTB_Loneliness',
            'NIHTB_Instrumental_Support',
            'NIHTB_Perceived_Hostility',
            'NIHTB_Perceived_Rejection',
            'SNI_Ex_People_Count',
            'SNI_Ex_Embed_Networks',
            'SNI_Ex_Network_Diversity')

#-------------------
# Exclusion Criteria
#-------------------

# Exclude initial debugging runs (Excluded = "SAMPLE_NONE")
X_good <- X[X$Excluded != "SAMPLE_NONE", myvars]

# Exclude any rows containing NA in any column (incomplete cases)
X_good <- X_good[complete.cases(X_good),]

if (do_scatterplot) {
  
  require(car)

  # Style constants
  plot.width = 10
  plot.height = 8
  
  # Pairwise scatterplots
  scatterplotMatrix(X_good)

}

# Extract included NIHTB and SNI measures
NTB <- X_good[, colnames(X_good) %like% "NIHTB"] * 1.0
SNI <- X_good[, colnames(X_good) %like% "SNI"] * 1.0

# Standardize variables [mean,sd] = [0,1]
NTBs <- scale(NTB)
SNIs <- scale(SNI)

#----------------------------------------------
# Multivariate gradient boosted tree regression
# - predicting SNI measures from NIHTB Social Measures
#----------------------------------------------

if (do_mvtb) {
  
  library(mvtboost)
  
  # Run MV Tree Regression with cross validation
  mvtb.model <- mvtb(Y = SNIs, X = NTBs,
                     n.trees = 250, shrinkage = 0.05,
                     interaction.depth = 1, cv.folds = 10)
  
  cat('\nMVTB Summary\n')
  summary(mvtb.model)
  
  # Test SNI predictions from NTB
  SNIs.hat <- predict(mvtb.model, newdata=NTBs)
  
  cat('\nSNI Prediction R2\n')
  SNIs.R2 <- diag(var(SNIs.hat)/var(SNIs))
  print(SNIs.R2)

  cat('\nSNI Prediction RMSE\n')
  SNIs.RMSE <- rmse(SNIs.hat, SNIs)
  print(SNIs.RMSE)
  
  # Relative influences of NIHTB Social on SNI
  cat('\nRelative Influences of NTB Social on SNI\n')
  par(mar=c(15,12,2,2))
  relinfl <- t(mvtb.ri(mvtb.model, relative="tot"))
  mvtb.heat(relinfl, clust.method="ward.D")
   
  # Outcome covariance explained
  cat('\nOutcome Covariance Explained')
  covex <- mvtb.covex(mvtb.model, Y=SNIs, X=NTBs)
  par(mar=c(12,25,2,2))
  mvtb.heat(covex, clust.method="ward.D")
  
  # Training and CV errors
  Trees <- 1:length(mvtb.model$trainerr)
  errs <- data.frame(Trees, mvtb.model$trainerr, mvtb.model$cv.err)
  errs <- melt(errs,id="Trees")
  vp <- ggplot(errs, aes(x=Trees, y=value, color=variable))
  vp <- vp + geom_line()
  print(vp)

}

# Compare a variety of regression models using Caret
library(caret)
library(doMC)
registerDoMC(cores = 8)

#
# Collect data
#
social <- merge(NTBs, SNIs)

#
# Preprocessing checks
#

if (do_preproc_checks) {

  cat('\nNear-zero Variance Predictors\n')
  print(nearZeroVar(NTBs, saveMetrics = TRUE))
  
  cat('\nSummary of NTB correlations\n')
  corNTB <- cor(NTBs)
  print(summary(corNTB[upper.tri(corNTB)]))
  
  cat('\nNumber of highly correlated predictors\n')
  print(findCorrelation(corNTB, cutoff = .75))
  
  cat('\nLinear Combinations within Predictors\n')
  print(findLinearCombos(NTBs))

}


#------------------------------------------
# Setup 10-fold cross-validation
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

# Univariate regression formulae
f1 = SNI_Ex_People_Count ~ . - SNI_Ex_Embed_Networks - SNI_Ex_Network_Diversity
f2 = SNI_Ex_Embed_Networks ~ . - SNI_Ex_People_Count - SNI_Ex_Network_Diversity
f3 = SNI_Ex_Network_Diversity ~ . - SNI_Ex_People_Count - SNI_Ex_Embed_Networks

#------------------------------------------
# Gradient Boosted Tree (GBM)

if (do_gbm) {
  
  cat('\nGradient Boosted Trees Optimization\n')

  # Customize the training grid
  gbmGrid <-  expand.grid(interaction.depth = 1:5, 
                          n.trees = 50 * (1:10), 
                          shrinkage = c(0.05, 0.01, 0.005),
                          n.minobsinnode = 10)
  
  set.seed(1966)

  cat('\nGBM : SNI People Count\n')   
  gbmFit.SNIcount <- train(f1, data = social,
                           method = "gbm", metric = "Rsquared",
                           trControl = fitControl, tuneGrid = gbmGrid,
                           verbose = FALSE)
  
  print(summary(gbmFit.SNIcount))
  print(plot(gbmFit.SNIcount))

  cat('\nGBM : SNI Embedded Networks\n') 
  gbmFit.SNIemb <- train(f2, data = social,
                         method = "gbm", metric = "Rsquared",
                         trControl = fitControl, tuneGrid = gbmGrid,
                         verbose = FALSE)
  
  print(summary(gbmFit.SNIemb))
  print(plot(gbmFit.SNIemb))
  
  cat('\nGBM : SNI Network Diversity\n')
  gbmFit.SNIdiv <- train(f3, data = social,
                         method = "gbm", metric = "Rsquared",
                         trControl = fitControl, tuneGrid = gbmGrid,
                         verbose = FALSE)
  
  print(summary(gbmFit.SNIdiv))
  print(plot(gbmFit.SNIdiv))
  
}

#------------------------------------------
# Linear Regression with Backwards Selection

if (do_lrback) {
  
  cat('\nLinear Regression with Backwards Selection\n')
  
  lrbackGrid <-  expand.grid(nvmax = 1:6)

  set.seed(1966)
  
  cat('\nLR Backwards : SNI People Count\n')
  lrbackFit.SNIcount <- train(f1, data = social,
                              method = "leapBackward", metric = "Rsquared",
                              trControl = fitControl, tuneGrid = lrbackGrid)
  
  print(lrbackFit.SNIcount)
  print(plot(lrbackFit.SNIcount))

  cat('\nLR Backwards : SNI Embedded Networks\n')
  lrbackFit.SNIemb <- train(f2, data = social, 
                            method = "leapBackward", metric = "Rsquared",
                            trControl = fitControl, tuneGrid = lrbackGrid)
  
  print(lrbackFit.SNIemb)
  print(plot(lrbackFit.SNIemb))
  
  cat('\nLR Backwards : SNI Network Diversity\n')
  lrbackFit.SNIdiv <- train(f3, data = social,
                            method = "leapBackward", metric = "Rsquared",
                            trControl = fitControl, tuneGrid = lrbackGrid)
  
  print(lrbackFit.SNIdiv)
  print(plot(lrbackFit.SNIdiv))
  
}


