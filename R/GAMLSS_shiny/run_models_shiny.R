run_models <- function(data, covariates = NULL, criterion = "BIC", updateProgress = NULL, lambda = 1, family = "NO"){
  
  # Progress bar function
  update = function(n){
    if(is.function(updateProgress)) {
      updateProgress(value = n)
    }
  }
  
  # Using metaprogramming expressions to define models
  ML_expr = expr(pb(xvar, method = "ML"))
  AIC_expr = expr(pb(xvar, method = "GAIC", k = 2))
  BIC_expr = expr(pb(xvar, method = "GAIC", k = log(length(data$xvar))))
  # BIC_expr = expr(fp(xvar, npoly = 3))
  GAIC3_expr = expr(pb(xvar, method = "GAIC", k = 3))
  stopifnot(length(c(covariates)) == length(colnames(data)[grepl("covariate",colnames(data))]))

  # Adding covariates to model definitions
  for (col in colnames(data)[grepl("covariate",colnames(data))]) {
    col = sym(col)
    ML_expr = expr(!!ML_expr + pb(!!col, method = "ML"))
    AIC_expr = expr(!!AIC_expr + pb(!!col, method = "GAIC", k = 2))
    BIC_expr = expr(!!BIC_expr + pb(!!col, method = "GAIC", k = log(length(data$xvar))))
    GAIC3_expr = expr(!!GAIC3_expr + pb(!!col, method = "GAIC", k = 3))
  }
  
  # Adding random subject effects
  ML_expr_sub = expr(!!ML_expr + random(subject_id, lambda = lambda))
  AIC_expr_sub = expr(!!AIC_expr + random(subject_id, lambda = lambda))
  BIC_expr_sub = expr(!!BIC_expr + random(subject_id, lambda = lambda))
  GAIC3_expr_sub = expr(!!GAIC3_expr + random(subject_id, lambda = lambda))
  
  # Running models
  tryCatch(
    {
      if(criterion == "AIC"){
        model = gamlss(eval(expr(yvar ~ !!AIC_expr_sub)), 
                       sigma.fo = eval(expr(~ !!AIC_expr)), 
                       tau.formula = eval(expr(~ !!AIC_expr)), 
                       nu.formula = eval(expr(~ !!AIC_expr)), 
                       family = family,
                       method = mixed(10, 50),
                       data = as.data.frame(data))
      } else if(criterion == "ML"){
        model = gamlss(eval(expr(yvar ~ !!ML_expr_sub)), 
                             sigma.fo = eval(expr(~ !!ML_expr)), 
                             tau.formula = eval(expr(~ !!ML_expr)), 
                             nu.formula = eval(expr(~ !!ML_expr)),
                             family = family,
                             method = mixed(10,50),
                             data = as.data.frame(data))
      } else if(criterion == "BIC"){
        model = gamlss(eval(expr(yvar ~ !!BIC_expr_sub)), 
                              sigma.fo = eval(expr(~ !!BIC_expr)), 
                              tau.formula = eval(expr(~ !!BIC_expr)), 
                              nu.formula = eval(expr(~ !!BIC_expr)), 
                              family = family,
                              method = mixed(10, 50), 
                              data = as.data.frame(data))
      } else if(criterion == "GAIC3"){
        model = gamlss(eval(expr(yvar ~ !!GAIC3_expr_sub)), 
                                sigma.fo = eval(expr(~ !!GAIC3_expr)), 
                                tau.formula = eval(expr(~ !!GAIC3_expr)), 
                                nu.formula = eval(expr(~ !!GAIC3_expr)), 
                                family = family,
                                method = mixed(10, 50),
                                data = as.data.frame(data))
        
        model$formula = as.character(eval(expr(yvar ~ !!GAIC3_expr)))
      }
    },
    error = function(e){
      message("An error occurred - most likely the model failed to converge. Try a different model or criterion. Here's the error message:\n",e)
      model = NULL
    },
    finally = {
      update(4)
      return(model)
    }
  )
}

