library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)


getFile_debug = function(datapath = "C:/Users/chris/OneDrive - Yale University/Projects/GAMLSS/data/FINALIZEDVEP_Without6MSleep_VEP_Longitudinal_Modeling_Data_11-15-23.csv",
                         xvar = "vep_age_days", yvar = "peak_latency_P1", covariates = NULL){
  read_csv(datapath) %>% 
    dplyr::mutate(subject_id = as.factor(subject_id)) %>%
    dplyr::select(subject_id:peak_latency_N2)
    # dplyr::filter(!if_any(where(is.numeric), ~ . == 99999))
}



modifyDf_debug = function(data, xvar = "vep_age_days", yvar = "peak_latency_P1", covariates = NULL){
  data = data %>% dplyr::filter(!if_any(all_of(c(xvar, yvar, covariates)), ~ . == 99999))
  if(!is.null(covariates)){
    df = data[,c("subject_id", xvar, yvar,covariates)]
    covnames = sprintf("covariate%d", 1:(length(covariates)))
    colnames(df) = c("subject_id","xvar","yvar",covnames)
  } else{
    df = data[,c("subject_id", xvar, yvar)]
    colnames(df) = c("subject_id","xvar","yvar")
  }
  df
}

facet_plot = function(data, lines = F){
  data_long = data %>% pivot_longer(cols = grep("peak",colnames(data))) %>% 
    filter(child_sex_0female_1male_sa != 99999)
  if(lines){
    ggplot(data_long, aes(vep_age_days, value, 
                          color = child_sex_0female_1male_sa, 
                          group = child_sex_0female_1male_sa)) + 
      geom_point() + geom_line(alpha = .5, aes(group = subject_id)) + 
      geom_smooth(color = "red", method = "loess", aes(group = NULL)) + 
      facet_wrap(vars(name), scales = "free_y")
  } else
    ggplot(data_long, aes(vep_age_days, value, 
                          color = child_sex_0female_1male_sa, 
                          group = child_sex_0female_1male_sa)) + 
      geom_point() + 
      geom_smooth(color = "red", method = "loess", aes(group = NULL)) + 
      facet_wrap(vars(name), scales = "free_y")
}

fit_model = function(model, df){
  ypred = predictAll(model, data = df)$mu
  # sigmapred = predictAll(model, data = df)$sigma[order(df$xvar)]
  sigmapred = fitted(model, "sigma")
  
  ggdata = cbind(ypred, sigmapred, df[,c("xvar", "yvar")])
  colnames(ggdata) = c("ypred","sigmapred","xvar","yvar")
  ggplot(ggdata, aes(x = xvar)) + 
    geom_line(linewidth = 1, aes(y = ypred, color = 'Predicted Mean')) + 
    geom_line(aes(y = sigmapred)) + 
    geom_point(color = "#87A5C0", aes(y = yvar)) + 
    geom_ribbon(fill = "#6C0E23", alpha = .3, aes(ymin = ypred - sigmapred, 
                                                  ymax = ypred + sigmapred, color = 'Predicted SD')) + scale_color_manual(name = "",
                                                             breaks = c('Predicted Mean', 'Predicted SD'), 
                                                             values = c('Predicted Mean' = "darkblue",'Predicted SD' = "#6C0E23")) + 
    theme_bw() + ggtitle('Mu and Sigma Model Fit')
}

centile_plot = function(model, df)({
  centiles(model, xvar = df$xvar, save = T)
})

gamlss_debug = function(datapath = "C:/Users/chris/OneDrive - Yale University/Projects/GAMLSS/data/FINALIZEDVEP_Without6MSleep_VEP_Longitudinal_Modeling_Data_11-15-23.csv",
                        xvar = "vep_age_days", yvar = "peak_latency_P1", covariates = NULL, lines = F, criterion = "BIC", lambda = 1, family = "NO",
                        random = F) {
  data = getFile_debug(datapath, xvar, yvar, covariates)
  df = modifyDf_debug(data, xvar, yvar, covariates)
  model = run_models(data = df, covariates = covariates, criterion = criterion, lambda = lambda, family = family, random = random)
  fp = facet_plot(data, lines)
  fit = fit_model(model, df)
  cp = centile_plot(model, df)
  centile_zscores = NULL
  centile_curves = NULL
  try({
      centile_zscores = centiles.pred(model,type = "z-scores", xname = "xvar", 
                                      xvalues = df$xvar, yval = df$yvar, data = df, plot = F)
      centile_curves = centiles.pred(model, xname = "xvar", 
                                      xvalues = df$xvar, data = df, plot = F)
    }, silent = T)
    return(list(data = data, df = df, model = model, 
              facet_plot = fp, centile_plot = cp, fit_plot = fit, 
              centile_zscores = centile_zscores,
              centile_curves = centile_curves, datapath = datapath, 
              xvar = xvar, yvar = yvar, covariates = covariates, 
              lines = lines, criterion = criterion, lambda = lambda, 
              family = family, random = random))
}
