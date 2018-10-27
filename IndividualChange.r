## Eduardo Estrada. Psychology UC Davis
## eestradaa@ucdavis.edu  -  eduardo.estrada.rs@gmail.com

# R code for the analysis described in:
# Estrada, Ferrer, Shaywitz, Holahan, & Shaywitz (2018).
# Identifying atypical change at the individual level from childhood to adolescence.
# Developmental Psychology 54(11), 2193-2206.
# http://dx.doi.org/10.1037/dev0000583


d <- d.calib       # A data frame composed of a set of observations coming from 2 known populations/groups
d.test <- d.test   # A data frame composed of new observations. They come from either one of the 2 groups, 


vars <- c("gviq1","gviq3", "gviq5","gviq7","gviq9") # Repeated measures of the variable at grades 1,3,5,7,9

d.typ <- d[d$group=="typical", vars]
d.atp <- d[d$group=="atypical", vars]


## Identifying typical and atypical trajectories (>2 time points) ------------------

## Estimating LCS model ##
if ("package:lavaan" %in% search() == FALSE) {require(lavaan)}

LCSm = "
# Declaring latent level
viq1 =~ 1* gviq1
viq3 =~ 1* gviq3
viq5 =~ 1* gviq5
viq7 =~ 1* gviq7
viq9 =~ 1* gviq9

# Auto-regression
viq3 ~ 1* viq1
viq5 ~ 1* viq3
viq7 ~ 1* viq5
viq9 ~ 1* viq7

# Define latent change
d3 =~ 1* viq3
d5 =~ 1* viq5
d7 =~ 1* viq7
d9 =~ 1* viq9

# Auto-proportions
d3 ~ b * viq1
d5 ~ b * viq3
d7 ~ b * viq5
d9 ~ b * viq7

# Latent intercept and slope
viq0 =~ 1* viq1
viqS =~ 1*d3 + 1*d5 + 1*d7 + 1*d9
viq0 ~ InMn * 1
viqS ~ SlMn * 1

# Latent variances and covariances
viq0 ~~ InV * viq0
viqS ~~ SlV * viqS
viq0 ~~ InSlCv * viqS # Covariance

# Measurement errors
gviq1 ~~ er * gviq1
gviq3 ~~ er * gviq3
gviq5 ~~ er * gviq5
gviq7 ~~ er * gviq7
gviq9 ~~ er * gviq9
"

# Fit the model for each of the two groups
LCS.fit.tp <- lavaan(model=LCSm, d.typ, missing = "fiml")
LCS.fit.at <- lavaan(model=LCSm, d.atp, missing = "fiml")

# Extract model info from both groups
modelInfo <- list()
modelInfo$sampstat$typical <- inspect(LCS.fit.tp, what="sampstat.h1", missing = "fiml")
modelInfo$mod_cov$typical <- inspect(LCS.fit.tp, what="cov.ov")
modelInfo$mod_mn$typical <- inspect(LCS.fit.tp, what="mean.ov")
modelInfo$InV_t <- coef(LCS.fit.tp)["InV"]

modelInfo$sampstat$atypical <- inspect(LCS.fit.at, what="sampstat.h1", missing = "fiml")
modelInfo$mod_cov$atypical <- inspect(LCS.fit.at, what="cov.ov")
modelInfo$mod_mn$atypical <- inspect(LCS.fit.at, what="mean.ov")



# Function for computing individual likelihoods
ILi <- function(data, varnames, referenceModel = modelInfo,
                classes = c("typical", "atypical"),
                FullInfo = TRUE) {
  Yi <- data[,names(data) %in% varnames]
  Occ <- which(colnames(referenceModel$sampstat$typical$cov) %in% varnames)
  
  if ((FullInfo == TRUE) & (sum(is.na(Yi))>0)) { # Full Information computations
    r <- 1 * is.na(Yi) # Create matrix of missing data points
    mdp <- (r %*% (2^((seq_len(ncol(Yi))) - 1))) + 1 # Assign one number to each missing data pattern
    
    Yi <- cbind(Yi, mdp)
    
    indstats <- sapply(unique(mdp), function(x) {
      ids <- data$id[Yi$mdp==x]
      cYi <- Yi[Yi$mdp==x,varnames]
      cr <- matrix(r[Yi$mdp==x,], ncol=length(varnames))
      if (nrow(cr)>1) {cr <- unique(cr)}
      occasions <- which(as.logical(1-cr))
      
      # Person Maximum Likelihood for the saturated model
      ILsat_tp <-
        as.vector((-1/2) * ( length(referenceModel$sampstat$typical$mean[occasions]) * log(2*pi) +
                               log(det(referenceModel$sampstat$typical$cov[occasions,occasions])) +
                               mahalanobis(cYi[,occasions],
                                           center=referenceModel$sampstat$typical$mean[occasions],
                                           cov=referenceModel$sampstat$typical$cov[occasions,occasions]) ) )
      ILsat_at <-
        as.vector((-1/2) * ( length(referenceModel$sampstat$atypical$mean[occasions]) * log(2*pi) +
                               log(det(referenceModel$sampstat$atypical$cov[occasions,occasions])) +
                               mahalanobis(cYi[,occasions],
                                           center=referenceModel$sampstat$atypical$mean[occasions],
                                           cov=referenceModel$sampstat$atypical$cov[occasions,occasions]) ) )
      
      # Person Maximum Likelihood for the test model
      IL_tp <-
        as.vector((-1/2) * ( length(occasions) * log(2*pi) +
                               log(det(referenceModel$mod_cov$typical[occasions,occasions])) +
                               mahalanobis(cYi[,occasions],
                                           center=referenceModel$mod_mn$typical[occasions],
                                           cov=referenceModel$mod_cov$typical[occasions,occasions]) ) )
      
      IL_at <-
        as.vector((-1/2) * ( length(occasions) * log(2*pi) +
                               log(det(referenceModel$mod_cov$atypical[occasions,occasions])) +
                               mahalanobis(cYi[,occasions],
                                           center=referenceModel$mod_mn$atypical[occasions],
                                           cov=referenceModel$mod_cov$atypical[occasions,occasions]) ) )
      
      return(cbind(ids,
                   as.integer(row.names(cYi[,occasions])),
                   ILsat_tp, ILsat_at, IL_tp, IL_at))
    } ) 
    
    indstats <- as.data.frame(do.call(rbind,indstats))
    colnames(indstats)[2] <- "id"
    
    data = merge(data, indstats, by = "id")
    
  } else { # Listwise deletion if missing data
    occasions <- Occ
    
    # Person Maximum Likelihood for the saturated model
    data$ILsat_tp <-
      as.vector((-1/2) * ( length(referenceModel$sampstat$typical$mean[occasions]) * log(2*pi) +
                             log(det(referenceModel$sampstat$typical$cov[occasions,occasions])) +
                             mahalanobis(Yi,
                                         center=referenceModel$sampstat$typical$mean[occasions],
                                         cov=referenceModel$sampstat$typical$cov[occasions,occasions]) ) )
    data$ILsat_at <-
      as.vector((-1/2) * ( length(referenceModel$sampstat$atypical$mean[occasions]) * log(2*pi) +
                             log(det(referenceModel$sampstat$atypical$cov[occasions,occasions])) +
                             mahalanobis(Yi,
                                         center=referenceModel$sampstat$atypical$mean[occasions],
                                         cov=referenceModel$sampstat$atypical$cov[occasions,occasions]) ) )
    
    # Person Maximum Likelihood for the test model
    data$IL_tp <-
      as.vector((-1/2) * ( length(occasions) * log(2*pi) +
                             log(det(referenceModel$mod_cov$typical[occasions,occasions])) +
                             mahalanobis(Yi,
                                         center=referenceModel$mod_mn$typical[occasions],
                                         cov=referenceModel$mod_cov$typical[occasions,occasions]) ) )
    data$IL_at <-
      as.vector((-1/2) * ( length(occasions) * log(2*pi) +
                             log(det(referenceModel$mod_cov$atypical[occasions,occasions])) +
                             mahalanobis(Yi,
                                         center=referenceModel$mod_mn$atypical[occasions],
                                         cov=referenceModel$mod_cov$atypical[occasions,occasions]) ) )
  }
  
  # Classifying
  data$ILsat_class <- factor((data$ILsat_at - data$ILsat_tp) < 0,
                             levels = c(TRUE, FALSE), labels = classes)
  data$IL_class <- factor((data$IL_at - data$IL_tp) < 0,
                          levels = c(TRUE, FALSE), labels = classes)
  
  return(data)
} # ILi <- function(...)


# Identifying typical and atypical trajectories in a new test sample - - - - - -
d.test <- ILi(d.test, varnames= vars)
with(d.test, table(group, IL_class)) # Check results




## Identifying typical and atypical change between >2 time points ------------------
vars <- c("gviq1", "gviq3") # Pre and Post measures


## Simple linear regression - - - - - -
# Estimate the model coefficients from the typical calibration sample
regModel <- lm(paste0(vars[2], " ~ ", vars[1]), data = d.typ)
regCoefs <- regModel$coefficients

# Compute residuals in the test sample
d.test[paste0(vars[2],"_preds")] <- regCoefs[[1]] + d.test[vars[1]] * regCoefs[[2]]
d.test[paste0(vars[2],"_obsRes")] <- d.test[vars[2]] - d.test[paste0(vars[2],"_preds")]



## 2 occasion LCS - - - - - -
LCSi <- function(data.calib, refGrp = "typical",
                 data.test, y1, y2,
                 var0, # We must provide either the variance of the latent intercept
                 reliab = NULL) # Or the test reliability observed in previous works
  {
  
  data.calib <- data.calib[data.calib$group==refGrp,c(y1, y2)]
  data.calib <- data.calib[complete.cases(data.calib),]
  vary1 <- var(data.calib[,y1], use = "pairwise.complete.obs")
  
  if (is.null(var0)) {
    var0 <- ""
    varEr <- vary1*(1-reliab) 
    varEr <- as.character(varEr)
    varEr <- paste0(varEr, " * ")
  } else {
    varEr <- vary1 - var0
    varEr <- as.character(varEr)
    varEr <- paste0(varEr, " * ")
    var0 <- ""
  }
  
  ly1 <- substring(y1, 2)
  ly2 <- substring(y2, 2)
  
  # Fitting model to typical calibration sample ---
  LCS2m = paste0("
                 # Declaring latent level
                 ",ly1," =~ 1* ",y1,"
                 ",ly2," =~ 1* ",y2,"
                 
                 ",ly2," ~ 1* viq1 # Auto-regression
                 
                 # Define latent change scores
                 d =~ 1* ",ly2," 
                 d ~ b * ",ly1," # Change-regression
                 
                 # Mean / intercept structure
                 ",ly1," ~ mu1 * 1
                 d ~ tauD * 1 # replace for muD when modeling latent covariance
                 
                 # Latent variances and covariances
                 ",ly1," ~~ ",var0,ly1,"  
                 d ~~ d
                 #",ly1," ~~ d
                 
                 # Measurement error
                 ",y1," ~~ ",varEr,y1,"
                 ",y2," ~~ ",varEr,y2,"
                 ")
  
  LCS2.fit <- lavaan(model=LCS2m, data.calib)
  
  # Extract parameters
  calib.pars <- coef(LCS2.fit)
  
  E_d <- calib.pars["tauD"] + calib.pars["b"] * calib.pars["mu1"]
  
  c.p.atr <- attributes(calib.pars)$names
  calib.pars <- as.character(calib.pars)
  attributes(calib.pars)$names <- c.p.atr
  
  # Fitting model to test sample ---
  LCS2new = paste0("
                   # Declaring latent level
                   ",ly1," =~ 1* ",y1,"
                   ",ly2," =~ 1* ",y2,"
                   
                   ",ly2," ~ 1* viq1 # Auto-regression
                   
                   # Define latent change scores
                   d =~ 1* ",ly2," 
                   d ~ ",calib.pars['b']," * ",ly1," # Change-regression
                   
                   # Mean / intercept structure
                   ",ly1," ~ 1
                   #d ~ ",calib.pars['muD']," * 1
                   d ~ ",calib.pars['tauD']," * 1
                   
                   # Latent variances and covariances
                   #",ly1," ~~ ",var0,ly1,"
                   ",ly1," ~~ ",ly1,"
                   d ~~ d
                   #",ly1," ~~ d
                   
                   # Measurement error
                   ",y1," ~~ ",varEr,y1,"
                   ",y2," ~~ ",varEr,y2,"
                   ")
  
  data.test.fit <- data.test[,c("id",y1, y2)]
  data.test.fit <- data.test.fit[complete.cases(data.test.fit),]
  
  LCS2new.fit <- lavaan(model=LCS2new, data.test.fit)
  fs <- as.data.frame(lavPredict(LCS2new.fit, type = "lv", method = "EBM"))
  
  fs[paste0(y2, "_Lres")] <- fs["d"] - E_d
  names(fs)[3] <- paste0("d",substring(ly1, 1,3), substring(ly1, 4), substring(ly2, 4))
  
  data.test.fit[colnames(fs)] <- fs
  data.test.fit <- data.test.fit[,-which(names(data.test.fit) %in% c(y1, y2))]
  
  data.test <- merge(data.test, data.test.fit, by="id", all=TRUE)
  
  return(data.test)
  
} # LCSi <- function(...)


## Assessing atypical pre-post change between 2 occasions in a new test sample - - - - - -
d.test <- LCSi(data.calib = d,
               var0 = modelInfo$InV_t,
               data.test = d.test,
               y1="gviq1", y2="gviq3")

with(d.test, plot(gviq3_obsRes, gviq3_Lres)) # Compare observed and latent residuals


