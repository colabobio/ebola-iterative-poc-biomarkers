# Extract single frame from MI data
# https://stackoverflow.com/a/42820029
getImpute <- function(impute, df, im) {
  cbind.data.frame(impute.transcan(x = impute, 
                                   imputation = im, 
                                   data = df, 
                                   list.out = TRUE, 
                                   pr = FALSE))
}

# Select variables using Elastic net
selectVariables <- function(rseed, nimp, fimp, train_data, inc_ct) {
  set.seed(rseed)
  
  # Generate imputed datasets
  formula <- as.formula(fimp)
  imp_data <- aregImpute(formula, data=train_data, n.impute=nimp)
  
  # Setup the LASSO
  
  # alpha=0 is Ridge Regression (L1 norm penalty only)
  # alpha=0.5 is elastic net (mixture of L1 and L2 at a 50%)
  # alpha=1 is lasso (L2 norm penalty only)
  aelast = 0.5 # actually, we are using the elastic net
  
  # Set binomial as the prediction family so we run logistic regression
  predfam = "binomial"
  
  # Iterate over imputations
  coeffs <- list()
  for (i in 1:nimp) {
    df <- getImpute(impute=imp_data, df=train_data, im=i)
    # Get outcome and predictor variables
    y <- as.matrix(df[,1])    
    x <- as.matrix(df[,2:ncol(df)])
    
    # Finds optimal lambda by cross-validation
    cv <- cv.glmnet(x, y, family=predfam, alpha=aelast, nfolds=10)
    lbest <- cv$lambda.min
    
    # Fit model
    fit <- glmnet(x, y, family=predfam, alpha=aelast, lambda=lbest)
    
    # Store coefficients  
    coeffs[[i]] <- coef(fit, s = "lambda.min")
  }
  
  ncoeff <- length(rownames(coeffs[[1]]))
  counts <- list()
  for (v in 1:ncoeff) {
    name <- rownames(coeffs[[1]])[v]
    c <- 0.0  
    for (i in 1:nimp) {
      if (0 < coeffs[[i]][v]) {
        c <- c + 1    
      }
    }
    counts[[name]] <- c/nimp  
  }
  
  # Sorting the counts https://stackoverflow.com/a/30651395
  counts <- counts[order(unlist(counts), decreasing=TRUE)]
  
  # And print out...
  res_names <- c()
  res_counts <- c()    
  for (k in names(counts)) {    
    n <- counts[[k]]
    cat(k, counts[[k]], '\n')
    res_names <- c(res_names, k)
    res_counts <- c(res_counts, n)
  }  
  
  newList <- list("names" = res_names, "counts" = res_counts)
}

# Trims leading and trailing whitespaces
# https://stackoverflow.com/a/21882152
trim <- function(x) {
  return(gsub("(^[[:space:]]+|[[:space:]]+$)", "", x))
}


saveToTXT <- function(obj, dir, fn) {
  sink(paste0(dir, "/", fn), append=FALSE, split=FALSE)
  print(obj)
  sink()    
}

saveEvaluation <- function(an, val, cal, dir) {
  print(val)
  plot(cal)
  
  saveToTXT(an, dir, "anova.txt")
  saveToPDF(an, dir, "anova.pdf")
  
  saveToTXT(val, dir, "validation.txt")
  saveToPDF(cal, dir, "calibration.pdf")
}

saveToPDF <- function(obj, dir, fn) {
  pdf(paste0(dir, "/", fn), useDingbats=FALSE)
  plot(obj)
  dev.off()
}

saveModelToCSV <- function(f, vars, dir) {
  terms = names(f$coefficients)
  coeff = unname(f$coefficients)
  
  # Extract RCS knots from model specs. 
  # This code is now specific to get the knots for age and fever temperature, 
  # but it could be generalized easily by providing a list with all the 
  # variables modeled as RCS.
  
  spec = specs(f)
  # print(str(spec)) # This is useful to understand the structure of the object
  
  ageIdx <- which(vars == "PatientAge")
  ctIdx <- which(vars == "CT")
  
  ageKnots <- spec$how.modeled[ageIdx, 2]
  ctKnots <- spec$how.modeled[ctIdx, 2]
  ageKnots <- trim(ageKnots)
  ctKnots <- trim(ctKnots)
  
  types <- rep.int("linear", length(terms))
  knots = rep.int("none", length(terms))
  
  age0 <- which(terms == "PatientAge")
  age1 <- which(terms == "PatientAge'")
  
  ct0 <- which(terms == "CT")
  ct1 <- which(terms == "CT'")
  
  types[age0] <- "RCS0"
  types[age1] <- "RCS1"
  knots[age1] <- ageKnots
  
  types[ct0] <- "RCS0"
  types[ct1] <- "RCS1"
  knots[ct1] <- ctKnots
  
  model <- data.frame("Term" = terms, "Coefficient" = coeff, "Type" = types, "Knots" = knots)
  write.table(model, file = paste0(dir, "/model.csv"), sep = ",", row.names=FALSE, qmethod = "double")    
}

saveDescription <- function(f, vars, dir) {
  print(f)
  
  saveToTXT(f, dir, "model.txt")
  saveToTXT(specs(f), dir, "specs.txt")
  
  saveModelToCSV(f, vars, dir)    
}

generateModel <- function(rseed, nimp, nboot, imp_formula, model_formula, train_data, vars, dir) {
  set.seed(rseed)
  
  # Impute data and fit pooled model
  imp_data <- aregImpute(as.formula(imp_formula), data=train_data, n.impute=nimp)
  model <- fit.mult.impute(as.formula(model_formula), lrm, imp_data, data=train_data)
  
  save(imp_formula, model_formula, train_data, imp_data, model, rseed, nimp, nboot, file =  paste0(dir, "/model.RData"))
  
  # Calculate ANOVA and validation/calibration
  mdl_anv <- anova(model)
  mdl_upd <- update(model, x=TRUE, y=TRUE)
  mdl_val <- validate(mdl_upd, B=nboot)
  mdl_cal <- calibrate(mdl_upd, B=nboot)
  
  # Calculate distribution summaries for potential predictor variables
  saveDescription(f=model, vars=vars, dir=dir)
  saveEvaluation(an=mdl_anv, val=mdl_val, cal=mdl_cal, dir=dir)
  
  return(model)
}