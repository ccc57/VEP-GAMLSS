# GAMLSS for VEP centile curves
# Chris C. Camp
# 12/2023
# 
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#


rm(list = ls())
library(gamlss)
library(readr)
library(dplyr)
library(ggplot2)
library(shiny)


source("run_models_shiny.R") # Contains model functions
# source("debugging_functions.R") To be removed

ui = fluidPage(
  titlePanel("GAMLSS for Brain Centile Curves"),
  sidebarLayout(
    sidebarPanel(
      fileInput("dataFile", "Data"),
      p("Data file should be in .csv format. Check column names in the getFile function of the app.R file."),
      uiOutput("xvarSelect"),
      uiOutput("yvarSelect"),
      uiOutput("randomSelect"),
      uiOutput("covSelect"),
      uiOutput("critSelect"),
      uiOutput("lambdaSelect"),
      uiOutput("familySelect"),
      uiOutput("run"),
      conditionalPanel(condition = "input.run > 0 && !$('html').hasClass('shiny-busy')",
        h5("Model Summary"),
        verbatimTextOutput("call"),
        h5("Effective degrees of freedom for each model parameter:"),
        verbatimTextOutput("edf"),
      ),
    ),
    conditionalPanel(condition = "input.run > 0 && !$('html').hasClass('shiny-busy')",
      mainPanel(
        h4("Model prediction mapped onto real data"),
        plotOutput("line"),
        htmlOutput("fit"),
        br(),
        h4("Centile Plots"),
        p("Adjust 'splits' parameter to separate plots by age ranges"),
        plotOutput("centilePlot"),
        plotOutput("centileFanPlot"),
        br(),
        h4("Coefficients for each distribution parameter"),
        plotOutput("parameters"),
        br(),
        downloadButton("downloadModel", "Save Model"),
        downloadButton("downloadScores", "Save Centile Z Scores"),
        downloadButton("downloadCurves", "Save Centile Data"),
        br(),
        h3("Model Diagnostics"),
        h4("Worm Plot"),
        plotOutput("wp"),
        p("See ", a("Buuren & Fredriks 2001", href = "https://onlinelibrary.wiley.com/doi/10.1002/sim.746"),
          "Table II for interpretation:"),
        plotOutput("wpInterpretation"),
        p("Excluded points in each plot - more exclusions indicates greater # of outliers and worse model fit."),
        verbatimTextOutput("wpText"),
        br(),
        h4('Another Plot of Residuals'),
        p("Residuals should be normally distributed"),
        plotOutput("residuals"),
        br(),
        h4('Q-Stats'),
        p("Significant Q values for each column (Z1, Z2, Z3, and Z4) suggest possible model 
        inadequacy for mu, sigma, nu, and tau, respectively.  
        Large Z values for a given age range identify the areas of worst fit for that 
        parameter. GAMLSS authors suggest |Z| > 2 as a guideline."),
        tableOutput("QStats"),
        plotOutput("QStatsPlot")
      )
    )
  )
)


server = function(input, output, session) {
  
  # Load and process data
  # You may need to make changes here if your variable names are different
  getFile = reactive({
    file = input$dataFile
    req(file)
    
    read_csv(file$datapath) %>% 
      dplyr::mutate(subject_id = as.factor(subject_id)) %>%
      dplyr::select(subject_id:peak_latency_N2)
  })
  
  # Runs when "Run GAMLSS" is pressed
  # Filters out missing data points for selected variables and strips variable names from data
  modifyDf = eventReactive(input$run, {
    df = getFile()
    df = df %>% dplyr::filter(!if_any(c(input$xvar, input$yvar, all_of(input$covariates)), ~ . == 99999))
    if(!is.null(input$covariates)){
      df = df[,c("subject_id", input$xvar, input$yvar, input$covariates)]
      covnames = sprintf("covariate%d", 1:(length(input$covariates)))
      colnames(df) = c("subject_id","xvar","yvar",covnames)
    } else{
      df = df[,c("subject_id", input$xvar, input$yvar)]
      colnames(df) = c("subject_id","xvar","yvar")
    }
    df
  })
  
  #User inputs
  output$xvarSelect = renderUI({
    df = getFile()
    selectInput("xvar", "Select explanatory variable", choices = colnames(df), selected = "vep_age_days")
  })
  output$yvarSelect = renderUI({
    df = getFile()
    selectInput("yvar", "Select predicted variable", choices = colnames(df), selected = "peak_latency_P1")
  })
  output$randomSelect = renderUI({
    df = getFile()
    checkboxInput("random",label = "Random effect of subject? (Centile predictions will be unavailable)", value = FALSE)
  })
  output$covSelect = renderUI({
    df = getFile()
    checkboxGroupInput("covariates", "Choose covariates if desired",choices = colnames(df))
  })
  output$critSelect = renderUI({
    df = getFile()
    selectInput("crit","Choose criterion to minimize:", choices = c("AIC","BIC","GAIC3","ML"), selected = "BIC")
  })
  output$lambdaSelect = renderUI({
    df = getFile()
    numericInput("lambda","Choose shrinkage parameter for random effect of subject", min = 0, max = 30, value = 1)
  })
  output$familySelect = renderUI({
    df = getFile()
    textInput("family","Choose distribution family", value = "NO")
  })
  output$run = renderUI({
    req(getFile())
    actionButton("run","Run GAMLSS")
  })
  
  # Generate models using run_models function defined in run_models_shiny.R
  model = eventReactive(input$run, {
    progress = shiny::Progress$new()
    progress$set(message = "Running models", value = 0)
    on.exit(progress$close())
    updateProgress = function(value) {
      progress$set(value = value)
    }
    run_models(data = modifyDf(), 
               covariates = input$covariates, 
               criterion = input$crit, 
               updateProgress = updateProgress,
               lambda = input$lambda,
               family = input$family, 
               random = input$random)
  })
  
  # Outputs
  
  # Model information
  output$call = renderPrint({
    summary(model())
    print(model()$mu.formula)
  })
  output$edf = renderPrint({
    edfAll(model())
  })
  
  # Plot of model fit
  output$line = renderPlot({
    # ypred = predict(model())
    # sigmapred = predict(model(), what = "sigma")
    ypred = predictAll(model(), data = modifyDf())$mu
    sigmapred = predictAll(model(), data = modifyDf())$sigma
    
    ggdata = cbind(ypred, sigmapred, modifyDf()[,c("xvar", "yvar")])
    ggplot(ggdata, aes(x = xvar)) + 
      geom_line(linewidth = 1, aes(y = ypred, color = 'Predicted Mean')) + 
      geom_point(color = "#87A5C0", aes(y = yvar)) + 
      geom_ribbon(fill = "#6C0E23", alpha = .3, aes(ymin = ypred - sigmapred, 
                                                    ymax = ypred + sigmapred, color = 'Predicted SD')) + 
      ylab(input$yvar) + xlab(input$xvar) + scale_color_manual(name = "",
                                                   breaks = c('Predicted Mean', 'Predicted SD'), 
                                                   values = c('Predicted Mean' = "darkblue",'Predicted SD' = "#6C0E23")) + 
      theme_bw() + ggtitle('Mu and Sigma Model Fit')
  })
  output$fit = renderText({
    paste("<b>Model fit:</b>",
          "<br>","r = ", round(cor(predict(model()),modifyDf()$yvar),2), 
           "<br>","R² = ", round(Rsq(model()), digits = 2),
           "<br>","AIC: ", round(GAIC(model())), 
           "<br>","BIC: ", round(GAIC(model(), k = log(length(modifyDf())))))
  })
  
  # Centile plots
  output$centilePlot = renderPlot({
    centiles(model(), xvar = modifyDf()$xvar, 
             xlab = input$xvar, ylab = input$yvar)
  })
  output$centileFanPlot = renderPlot({
    centiles.fan(model(), xvar = modifyDf()$xvar, xlab = input$xvar,ylab = input$yvar,
                 points = T, colors = "cm")
  })
  
  # Plot of model coefficients
  output$parameters = renderPlot({
    fittedPlot(model(), x = modifyDf()$xvar, xlab = input$xvar)
  })
  
  # Worm plots
  output$wp = renderPlot({
    wp(xvar = modifyDf()$xvar, resid = resid(model()), n.inter = 12, mar = c(1, 1, 1, 1), col = "#6C0E23", bg = "#87A5C0", bar.bg = "#6C0E23")
  })
  output$wpText = renderPrint({
    wp(xvar = modifyDf()$xvar, resid = resid(model()), n.inter = 12, mar = c(1, 1, 1, 1), col = "#6C0E23", bg = "#87A5C0", bar.bg = "#6C0E23")
  })
  output$wpInterpretation = renderImage({
    req(input$run)
    width  = session$clientData$output_wpInterpretation_width
    height = width * 0.2942686 # Aspect ratio of image
    
    fname = "TableII.png"
    list(src = fname, width = width, height = height)
    }, deleteFile = F)
  output$residuals = renderPlot({
    plot(model(), xvar = modifyDf()$xvar, parameters = par(mfrow = c(2, 2), mar = par("mar") + c(0, 1, 0,0),
                                                 col.axis = "blue4", col = "blue4", col.main = "blue4",
                                                 col.lab = "blue4", pch = 20, cex = 0.45, cex.lab = 1.2,
                                                 cex.axis = 1, cex.main = 1.2))
  })
  
  # Q stats
  output$QStats = renderTable({
    Q.stats(model(), xvar = modifyDf()$xvar, n.inter = 12)
  })
  output$QStatsPlot = renderPlot({
    Q.stats(model(), xvar = modifyDf()$xvar, n.inter = 12)
  })
  
  output$downloadModel = downloadHandler(
    filename = function() {
      # Use the selected dataset as the suggested file name
      paste0("GAMLSS_model_",input$yvar,".RDS")
    },
    content = function(file) {
      # Write the dataset to the `file` that will be downloaded
      saveRDS(model(), file)
    }
  )
  
  output$downloadScores = downloadHandler(
    filename = function() {
      # Use the selected dataset as the suggested file name
      paste0("centile_zscores_",input$yvar,".csv")
    },
    content = function(file) {
      # Write the dataset to the `file` that will be downloaded
      write.csv(centiles.pred(model(),type = "z-scores", xname = "xvar", 
                              xvalues = modifyDf()$xvar, yval = modifyDf()$yvar, data = modifyDf(), plot = F), file)
    }
  )
  
  output$downloadCurves = downloadHandler(
    filename = function() {
      # Use the selected dataset as the suggested file name
      paste0("centile_data_",input$yvar,".csv")
    },
    content = function(file) {
      # Write the dataset to the `file` that will be downloaded
      write.csv(centiles.pred(model(), xname = "xvar", 
                              xvalues = modifyDf()$xvar, data = modifyDf(), plot = F), file)
    }
  )
}





# Run the application 
shinyApp(ui = ui, server = server)
