run_models <- function(data, covariates = NULL, criterion = "BIC", updateProgress = NULL, lambda = 1, family = "NO", random = FALSE){
  
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
  
  ML_expr_sub = ML_expr
  AIC_expr_sub = AIC_expr
  BIC_expr_sub = BIC_expr
  GAIC3_expr_sub = GAIC3_expr
  
  # Adding random subject effects
  if(random == TRUE){
    ML_expr_sub = expr(!!ML_expr + random(subject_id, lambda = lambda))
    AIC_expr_sub = expr(!!AIC_expr + random(subject_id, lambda = lambda))
    BIC_expr_sub = expr(!!BIC_expr + random(subject_id, lambda = lambda))
    GAIC3_expr_sub = expr(!!GAIC3_expr + random(subject_id, lambda = lambda))
  }
  
  # Running models
  tryCatch(
    {
      if(criterion == "AIC"){
        model_call = call("gamlss", expr(yvar ~ !!AIC_expr_sub),
                          sigma.fo = expr(~ !!AIC_expr), 
                          tau.formula = expr(~ !!AIC_expr),
                          nu.formula = expr(~ !!AIC_expr),
                          family = family,
                          method = expr(mixed(10, 50)),
                          data = data)
        model = eval(model_call)
      } else if(criterion == "ML"){
        model_call = call("gamlss", expr(yvar ~ !!ML_expr_sub),
                          sigma.fo = expr(~ !!ML_expr), 
                          tau.formula = expr(~ !!ML_expr),
                          nu.formula = expr(~ !!ML_expr),
                          family = family,
                          method = expr(mixed(10, 50)),
                          data = data)
        model = eval(model_call)
      } else if(criterion == "BIC"){
        model_call = call("gamlss", expr(yvar ~ !!BIC_expr_sub),
                        sigma.fo = expr(~ !!BIC_expr), 
                        tau.formula = expr(~ !!BIC_expr),
                        nu.formula = expr(~ !!BIC_expr),
                        family = family,
                        method = expr(mixed(10, 50)),
                        data = data)
        model = eval(model_call)
        # model = gamlss(eval(expr(yvar ~ !!BIC_expr_sub)), 
        #                       sigma.fo = eval(expr(~ !!BIC_expr)), 
        #                       tau.formula = eval(expr(~ !!BIC_expr)), 
        #                       nu.formula = eval(expr(~ !!BIC_expr)), 
        #                       family = family,
        #                       method = mixed(10, 50), 
        #                       data = data)
      } else if(criterion == "GAIC3"){
        model_call = call("gamlss", expr(yvar ~ !!GAIC3_expr_sub),
                          sigma.fo = expr(~ !!GAIC3_expr), 
                          tau.formula = expr(~ !!GAIC3_expr),
                          nu.formula = expr(~ !!GAIC3_expr),
                          family = family,
                          method = expr(mixed(10, 50)),
                          data = data)
        model = eval(model_call)
      }
    },
    error = function(e){
      message("An error occurred in run_models_shiny.R - most likely the model failed to converge. Try a different model or criterion. Here's the error message:\n",e)
    },
    finally = {
      update(4)
      return(model)
    }
  )
}

