# Functions used in the analysis of GLAD task force, namely 
# GLAD_Analysis_RegressionModels_Mood_on_EnergyAdjustedResidualDietaryIntake.R,
# GLAD_Analysis_Model_Disgnosis.R,
# and Analysis_GLAD2025.R.
# Written by Ehsan Motazedi, Amsterdam UMC, 12-01-2024
# Last modified: 27-09-2025
library(car)
library(dplyr)
library(sandwich)
library(lmtest)
library(ggplot2)
library(grid)
library(gtable)
# 1) Function to perform the planned regression analyses of the outcome
# on dietary variables. Interactions with the dietary variables can also 
# be tested, as well as other interaction terms in the model.
# if energy_intake_var is given, energy adjusted residuals
# of the dietary intake variables will be used according to 
# Willett, Howe & Kushi (1997) instead of the raw values.
run_analysis<-function(outcomes = c(), dietary_var_names = c(),
                       energy_intake_var = c(),
                       covariates = c(), 
                       interaction_terms = c(),
                       interaction_with_dietary_var = c(),
                       interaction_contrasts = c(),
                       data = NULL,
                       level_for_no_outcome = c('No', '0'),
                       type_of_analysis = c('logistic', 'poisson'), 
                       maxit = 100, tol = 1e-08, ...){
  if(is.null(data) || length(outcomes)==0 || length(dietary_var_names)==0){
    warning('No analysis was performed as part of the model or data was missing.')
    return(NULL)
  }
  if(length(interaction_terms)>0){
    covariates<-unique(c(covariates, # add missing main effects of the interaction terms
                         unlist(strsplit(split=':|\\*|\\.', interaction_terms))))
    interaction_terms<-gsub('\\*|\\.', ':', interaction_terms)
  }
  if(length(interaction_with_dietary_var)>0){
    covariates<-unique(c(covariates, # add missing main effects of the interaction terms
                         interaction_with_dietary_var))
  }
  if(!all(c(outcomes, energy_intake_var, dietary_var_names, covariates)%in%colnames(data))){
    stop('Some model variables were missing in the given data.')
  }
  type_of_analysis<-unique(type_of_analysis)
  if(!all(type_of_analysis%in%c('logistic', 'poisson'))){
    stop('Analysis in type_of_analysis must be either "logistic" or "poisson".')
  }
  if(length(energy_intake_var)==0){
    warning('No variable has been specified for the total energy intake. The raw dietary intake values will be used.')
  } else if(length(energy_intake_var)>1){
    stop('Only one variable is allowed to be given for the total energy intake.')
  }
  if(length(c(covariates, interaction_terms))>0){
    confounding_adjustment<-c(covariates, interaction_terms)|>paste(sep='', collapse='+')
  } else {
    confounding_adjustment<-NULL
  }
  fitted_models<-list()
  data<-data[, c(outcomes, dietary_intake_variables, energy_intake_var, covariates)]
  for(outcome in outcomes){
    fitted_models[[outcome]]<-list()
    for(analysis in type_of_analysis){
      fitted_models[[outcome]][[analysis]]<-list()
      for(predictor in dietary_var_names){
        analysis_data<-data[, c(outcome, energy_intake_var, covariates), drop = FALSE]
        if(length(energy_intake_var)>0){
          analysis_data[, predictor]<-tryCatch({
            mod<-lm(as.formula(paste(predictor, energy_intake_var, sep='~')), data=data)
            newdata<-data[1, ]
            newdata[, predictor]<-NA
            newdata[, energy_intake_var]<-mean(data[, energy_intake_var], na.rm=TRUE)
            resid(mod)+predict(mod, newdata=newdata, type='response') # return 'standardized' residuals by adding the predicted dietary value at the mean energy intake          
          }, error = function(err){
            write(sprintf('An error occurred in getting energy adjusted residuals, with %s predictor and %s for the energy intake:',
                          predictor, energy_intake_var), stderr())
            message(err)
            write('\nRaw dietary intake values will be used!\n--------\n', stderr())
            return(data[, predictor])
          }
          )
        } else {
          analysis_data[, predictor]<-data[, predictor]
        }
        analysis_data <- analysis_data[complete.cases(analysis_data), ]
        if(length(interaction_with_dietary_var)==0){
          interaction_terms_with_dietary<-c()
        } else {
          interaction_terms_with_dietary<-paste(interaction_with_dietary_var,
                                                rep(predictor, length(interaction_with_dietary_var)), sep = ':')
          i<-0 # Use the global contrasts set for the current R session if no interaction contrasts is given, otherwise make use of the given contrast function. The interaction_contrasts argument can be a vector giving contrasts for some or all of the factors with interaction.
          while(i<length(interaction_contrasts) && i<length(interaction_with_dietary_var)){ # Set the specific given contrasts for the interactions with the dietary variable in the regression model (for example, orthonormal contrasts)
            i<-i+1
            ff<-interaction_with_dietary_var[i] # Match the given contrasts one-by-one to the factors with interaction, until the given contrasts are exhausted or all the factors are matched to a given contrast
            analysis_data[, ff]<-droplevels(analysis_data[, ff]) # drop levels non-existing in the analysis dataset from the interaction factors
            ff_interaction_contrasts<-interaction_contrasts[[i]]
            contrasts(analysis_data[, ff])<-round(do.call(ff_interaction_contrasts, list(levels(analysis_data[, ff]))), 3) # for the factors that are going to be used in interactions
            colnames(attr(analysis_data[, ff], 'contrasts'))<-paste0('contrast_level_', 1:dim(attr(analysis_data[, ff], 'contrasts'))[2])
          }
        }
        fitted_models[[outcome]][[analysis]][[predictor]]<-tryCatch({
          if(analysis=='logistic'){
            glm(as.formula(paste(outcome,
                                 paste(c(energy_intake_var, confounding_adjustment, predictor, interaction_terms_with_dietary),
                                       sep='', collapse='+'), sep='~')),
                data = analysis_data,
                family = binomial(link = 'logit'),
                control = list(maxit = maxit, epsilon = tol), ...)
          } else {
            do_poisson_glm_for_binary_outcome(as.formula(paste(outcome,
                                                               paste(c(energy_intake_var, confounding_adjustment, predictor, interaction_terms_with_dietary),
                                                                     sep='', collapse='+'), sep='~')),
                                              data = analysis_data, outcome = outcome, level_for_no_outcome = level_for_no_outcome,
                                              maxit = maxit, tol = tol, ...)
          }
        }, error = function(err){
          write(sprintf('An error occurred in %s regression, with %s outcome and %s predictor:', 
                        analysis, outcome, predictor), stderr())
          message(err)
          write('\n--------\n', stderr())
          return(NULL)
        })
      }
    }
  }
  return(fitted_models)
}

# 2) Get robust variance estimation and confidence intervals
my_summary <- function(model, robust = TRUE, type='HC1', alpha = 0.05, decimals = 5, take_exp = FALSE){
  df <- model$df.residual
  summ <- summary(model)
  if(robust){
    coeftest <- lmtest::coeftest(model, vcov = vcovHC(model, type=type),
                                 df = df)
    summ$coefficients[,2] <- coeftest[,2]
    summ$coefficients[,3] <- coeftest[,3]
    summ$coefficients[,4] <- coeftest[,4]
  }
  if(is.null(df)){
    qci = qnorm(1-alpha/2) 
  } else {
    qci = qt(1-alpha/2, df=df)
  }
  summ$coefficients<-as.data.frame(summ$coefficients)
  summ$coefficients[, sprintf('Low_%.2fCI', 1-alpha)]<-summ$coefficients[, 'Estimate'] - qci*summ$coefficients[, 'Std. Error']
  summ$coefficients[, sprintf('Up_%.2fCI', 1-alpha)]<-summ$coefficients[, 'Estimate'] + qci*summ$coefficients[, 'Std. Error']
  if(!is.null(df)){
    colnames(summ$coefficients)[colnames(summ$coefficients)=='z value']<-'t-value'
    colnames(summ$coefficients)[colnames(summ$coefficients)=='Pr(>|z|)']<-'Pr(>|t|)'
  }
  if(take_exp){
    summ$exp_coefficients<-data.frame(exp(summ$coefficients[,1]),
                                      exp(summ$coefficients[, sprintf('Low_%.2fCI', 1-alpha)]),
                                      exp(summ$coefficients[, sprintf('Up_%.2fCI', 1-alpha)]))
    colnames(summ$exp_coefficients)<-c(paste0('Exp_', colnames(summ$coefficients)[1]),
                                       sprintf('Exp_Low_%.2fCI', 1-alpha),
                                       sprintf('Exp_Up_%.2fCI', 1-alpha))
    rownames(summ$exp_coefficients)<-rownames(summ$coefficients)
    summ$exp_coefficients[summ$exp_coefficients > 1000]<- Inf
    summ$coefficients<-round(cbind(summ$coefficients, summ$exp_coefficients), decimals)
    summ$exp_coefficients<-NULL
  }
  return(summ)
}

# 3) Write a summary of the fitted regression model, i.e. the results, to the specified output file
#    or to the standard output. The my_summary function (defined function #2 in this file) will be used
#    with the provided parameters. If do.anova is set to true (default false), an anova analysis will follow
#    with the type specified by type.anova (default to type III).
write_results<-function(results, outfile=NULL, analysis_name = 'performed analysis',
                        sepp_adjusted_unadjusted_analysis = paste(rep('#', getOption('width')), collapse = '', sep=''),
                        sepp_dietary_analysis = paste(rep('-', floor(getOption('width')/2)), collapse = '', sep=''),
                        type_of_analysis = 'logistic',
                        alpha = 0.05,
                        decimals = 5,
                        robust = FALSE,
                        take_exp = TRUE,
                        append = FALSE,
                        do.anova = FALSE,
                        type.anova = 'III',
                        print.digits = 10){
  if(do.anova){
    if(length(type.anova)==0){
      type.anova = 'III'
    } else {
      for(i in 1:length(type.anova)){
        if(!(type.anova[i]%in%c('I', 'II', 'III') || type.anova[i]%in%1:3)){
          warning(sprintf("Type of anova must be 1, 2, 3 or 'I', 'II', 'III'. The given type %s is not recognized and type I will therefore be used instead.", type.anova[i]))
          type.anova[i] = 1
        }
      }
    }
    type.anova = as.character(type.anova)
  }
  tryCatch({
    if(!is.null(outfile)){
      sink(outfile, append = append)
    }
    write(sprintf('\n%s\n%s\nResults of the %s:\n', 
                  sepp_adjusted_unadjusted_analysis, sepp_adjusted_unadjusted_analysis, analysis_name), stdout())
    for(name_of_outcome in names(results)){
      write(sprintf('\n#### Results of the outcome: %s\n', name_of_outcome), stdout())
      for(predictor in names(results[[name_of_outcome]][[type_of_analysis]])){
        write(sprintf('\n%s\n>>>Dietary intake: %s\nOutcome is %s\n', sepp_dietary_analysis, predictor, name_of_outcome), stdout())
        print(my_summary(results[[name_of_outcome]][[type_of_analysis]][[predictor]],
                         robust = robust, decimals = decimals, alpha = alpha, take_exp = take_exp),
              digits = print.digits)
        if(do.anova){
          for(i in 1:length(type.anova)){
            anova_type = type.anova[i]
            write(sprintf('\n-----\nResults of type %s ANOVA analysis (analysis of deviance residuals):\n', anova_type), stdout())
            ff<-results[[name_of_outcome]][[type_of_analysis]][[predictor]]$family
            if(anova_type%in%c('1', 'I')){
              if(!is.null(ff) && ff$family%in%c('binomial', 'poisson')){
                print(anova(results[[name_of_outcome]][[type_of_analysis]][[predictor]], test = 'Chisq'), digits = print.digits)
              } else {
                print(anova(results[[name_of_outcome]][[type_of_analysis]][[predictor]], test = 'F'), digits = print.digits)
              }  
            }  else {
              if(!is.null(ff) && ff$family%in%c('binomial', 'poisson')){
                print(car::Anova(results[[name_of_outcome]][[type_of_analysis]][[predictor]], type = anova_type, test.statistic = 'LR'), digits = print.digits)
              } else {
                print(car::Anova(results[[name_of_outcome]][[type_of_analysis]][[predictor]], type = anova_type, test.statistic = 'F'), digits = print.digits)
              } 
            }
          }
        }
      }
    }
  },
  error = function(err){
    message(err)
  }, 
  finally = {
    if(!is.null(outfile)){
      sink()  
    }
    return(invisible(NULL))
  }
  )
}
# 4) Create scatterplots for the fitted dependent variable values 
#    on the glm link scale or the original outcome scale, vs. the 
#    values of a predictor in the fitted glm model.
plot.scatter<-function(model, predictor,
                       type = c('link', 'response'),
                       xlab = NULL, ylab=NULL,
                       plot = FALSE){
  if(!all(type%in%c('link', 'response'))){
    stop('The glm estimate type must be "link" or "response".')
  }
  type = type[1]
  if(is.null(xlab)){
    xlab = predictor
  }
  if(is.null(ylab)){
    ylab = sprintf('Predicted values (%s scale)', type)
  }
  newdata<-model$data
  newdata<-newdata[complete.cases(newdata), ]
  plotdata<-data.frame(x = newdata[, predictor], 
                       y = predict(model, newdata=newdata,  type = type))
  p<-ggplot(plotdata, aes(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = 'loess') +
    xlab(xlab) + ylab(ylab)
  if(plot){
    print(p)
  }
  return(p)
}
# 5) Plot scatterplots for the fitted dependent variable values
#    on the glm link scale or the original outcome scale, vs. the 
#    values of "all" predictors in the fitted glm model, and plot 
#    them to file or standard plotting device.
plot.scatter.all_model_predictors<-function(model,
                                            type = c('link', 'response'), ylab = NULL, 
                                            plot.per.page = 3, outfile = NULL){
  tryCatch({
    outcome<-as.character(model$formula)[2]
    if(is.na(outcome) || length(outcome)==0){
      stop('The given model is not a valid glm object.')
    }
    predictors<-names(model$model)
    predictors<-predictors[!predictors==outcome]
    pp<-lapply(predictors, function(x) ggplotGrob(plot.scatter(model = model, predictor = x,
                                                               type = type, ylab = ylab,
                                                               xlab = x, plot = FALSE)))
    
    if(!is.null(outfile)){
      create_device(outfile)
    }
    depicted_plots<-0
    plot.per.page<-min(plot.per.page, length(predictors))
    while(depicted_plots<length(predictors)){
      gg<-do.call(rbind, c(pp[(depicted_plots+1):min(depicted_plots+plot.per.page, length(predictors))], size='first')) 
      gg$widths<-do.call(unit.pmax, lapply(pp[(depicted_plots+1):min(depicted_plots+plot.per.page, length(predictors))], function(g) g$widths))
      grid.newpage()
      grid.draw(gg)
      depicted_plots<-depicted_plots+plot.per.page
    }
    return(pp)
  }, 
  error = function(e){
    message(e)
    return(NULL)
  }, finally = {
    if(!is.null(attributes(.Device)$filepath) && attributes(.Device)$filepath==outfile){
      dev.off()
    }
  }
  )
}

# 6) Function to create a device based on file extension
create_device <- function(file_name) {
  extension <- tools::file_ext(file_name)
  if (extension == 'pdf') {
    pdf(file_name)
  } else if (extension == 'png') {
    png(file_name)
  } else if (extension == 'jpeg' || extension == 'jpg') {
    jpeg(file_name)
  } else {
    stop('Error: Unsupported file extension.')
  }
}

# 7) Calculate the Cook's distance, standardized residuals and differences in the fit Beta (DFBETA)
#    for a given model. Report outliers based on the given thresholds for each measure of individual
#    influence. If thresholds are not given, 4/sqrt(n-p) and 2/sqrt(n) will be used for
#    the Cook's distance and for the DFBETA, respectively, and 3 will be used for the standardized residuals
#    to report the outliers in the data used to fit the model. These thresholds can of course be set to Inf to prevent
#    any outlier report based on one or more criteria. If 'model_update' argument is set to TRUE (default),
#    the given model will also be update by removing the detected outliers, and the updated model is added to the 
#    output. This option is useful to perform a sensitivity analysis by excluding outliers.
#    For dfbeta, one can specify the main predictor(s) in the model whose Beta is of ineterest. If it is not given,
#    one will consider the maximum reported dfbeta for all parameters.
get_outliers<-function(model, mainPredictor = NULL, threshold.Cooks = NULL, threshold.stdres = NULL, threshold.dfbeta = NULL,
                       update_model = TRUE, type.resid = 'deviance'){
  if(!any(class(model)%in%c('glm', 'lm'))){
    stop('Input model must be of "glm" or "lm" class!')
  }
  MyData<-model$data
  n<-dim(MyData)[[1]] # Number of included individuals
  predictorSet<-strsplit(split='\\s+\\+\\s+', as.character(model$formula)[3])[[1]] # Set of all predictors in the model
  p<-length(predictorSet) # number of parameters in the model, number of variables minus 1 for the output
  if(length(mainPredictor)==0){
    mainPredictor<-predictorSet
  } else if(any(!mainPredictor%in%predictorSet)){
    stop('Given predictors are not included in the given model.')
  }
  if(length(threshold.Cooks) == 0){
    threshold.Cooks = 4/sqrt(n-p)
  }
  if(length(threshold.dfbeta) == 0){
    threshold.dfbeta = 2/sqrt(n)
  }
  if(length(threshold.stdres) == 0){
    threshold.stdres = 3
  }
  output<-list(threshold.Cooks = cooks.distance(model),
               threshold.dfbeta = dfbeta(model), threshold.stdres = rstandard(model, type = type.resid)) # Get the influence measures
  output[['threshold.Cooks']]<-which(output[['threshold.Cooks']]>threshold.Cooks) # Get outliers on Cook's distance
  output[['threshold.stdres']]<-which(abs(output[['threshold.stdres']])>threshold.stdres*sd(output[['threshold.stdres']])) # Get outliers on standardized residuals
  output[['threshold.dfbeta']]<-output[['threshold.dfbeta']][, colnames(output[['threshold.dfbeta']])[unlist(lapply(colnames(output[['threshold.dfbeta']]),
                                                                                                                    function(beta) any(unlist(lapply(mainPredictor, function(x) grepl(x, beta))))))], drop = FALSE] # Just look into the given main predictors and their interactions
  output[['threshold.dfbeta']]<-which(apply(output[['threshold.dfbeta']], 1, function(df_beta_i) max(abs(df_beta_i), na.rm = TRUE)>threshold.dfbeta)) # Get outliers on dfbetas
  if(update_model){
    filtered_subjects <- setdiff(1:dim(MyData)[[1]], do.call(c, output))
    output$model.update <- update(model, formula = model$formula,
                                  data = MyData[filtered_subjects,,drop = FALSE], control = model$control) # Return modified model having removed the outliers from the data
  }
  return(output)
}

# 8) Generate residual plots against predictor values for the given glm model.
#    If no predictor is given, using the argument predictors, residuals will be plotted against all main predictors 
#    present in the model. If plot.against.fitted is TRUE, then residuals will also be plotted
#    against the fitted linear predictor values.
#    Optionally, results can be exported to the given outfile.
plot.residual_predictors<-function(model, predictors = NULL,
                                   plot.against.fitted = TRUE,
                                   residual.type = c('deviance', 'predictive', 'pearson'), 
                                   plot.per.page = 3, outfile = NULL){
  tryCatch({
    outcome <- as.character(model$formula)[2]
    if(is.na(outcome) || length(outcome)==0){
      stop('The given model is not a valid glm object.')
    }
    all.predictors <- names(model$model)
    all.predictors <- all.predictors[!all.predictors==outcome]
    if(length(predictors)>0){
      predictors <- unique(do.call(c, strsplit(split=':\\|*', predictors))) # Get rid of interaction terms, only keep main effect predictors
      if(any(!predictors%in%all.predictors)){
        stop('The given predictor(s) do not exist in the given model.')
      }
    } else {
      predictors <- all.predictors
    }
    plot.data <- model$data[, predictors, drop = FALSE]
    plot.data[, 'Residuals'] <- rstandard(model, type = residual.type) 
    if(plot.against.fitted){
      found <- FALSE
      x<-0
      LP.name<-'LP'
      PV.name<-'PV'
      while(!found){
        if(LP.name %in% predictors || PV.name %in% predictors){
          x<-x+1
          LP.name<-paste0('LP', x)
          PV.name<-paste0('PV', x)
        } else {
          found <- TRUE
          rm(x)
        }
      }
      plot.data[, LP.name] <- predict(model, type = 'link')
      predictors<-c(predictors, LP.name)
      names(predictors)<-c(predictors[1:(length(predictors)-1)], 'Linear predictor values')
      plot.data[, PV.name] <- predict(model, type = 'response')
      predictors<-c(predictors, PV.name)
      names(predictors)<-c(predictors[1:(length(predictors)-1)], 'Predicted values')      
    }
    plot.N <- length(predictors)
    pp<-lapply(1:plot.N, function(i) {
      dd <- data.frame(x = plot.data[, predictors[i]], y = plot.data[, 'Residuals'])
      return(ggplotGrob(ggplot(data = dd, aes(x = x, y = y))+geom_point(color='black')+
                          geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') + geom_smooth(method = 'loess') +
                          labs(title = 'Standardized residual plot', x = names(predictors)[i], 
                               y = sprintf('Std residuals (%s)', residual.type))))
    })
    if(!is.null(outfile)){
      create_device(outfile)
    }
    depicted_plots<-0
    plot.per.page<-min(plot.per.page, plot.N)
    while(depicted_plots<plot.N){
      gg<-do.call(rbind, c(pp[(depicted_plots+1):min(depicted_plots+plot.per.page, length(predictors))], size='first')) 
      gg$widths<-do.call(unit.pmax, lapply(pp[(depicted_plots+1):min(depicted_plots+plot.per.page, length(predictors))], function(g) g$widths))
      grid.newpage()
      grid.draw(gg)
      depicted_plots<-depicted_plots+plot.per.page
    }
    return(pp)
  }, 
  error = function(e){
    message(e)
    return(NULL)
  }, finally = {
    if(!is.null(attributes(.Device)$filepath) && attributes(.Device)$filepath==outfile){
      dev.off()
    }
  }
  )
}

# 9) Function to perform Poisson regression for binary outcomes
# This function is needed to convert the dichotomous regression outcome 
# to a numeric one, assigning zero to the level representing the absence
# and assigning 1 otherwise. In case no outcome is specified in the function input,
# Poisson regression will be performed juts using the given formula.
do_poisson_glm_for_binary_outcome<-function(formula, data,
                                            outcome = NULL, level_for_no_outcome = c('No', '0'),
                                            tol = 1e-08, maxit = 100, ...){
  if(!is.null(outcome)){
    if(is.numeric(data[, outcome])){
      warning(sprintf('The given outcome "%s" is already numeric, no change will be made to it!', outcome))
    } else {
      binary_outcome<-as.factor(data[, outcome])
      if(length(intersect(toupper(levels(binary_outcome)), toupper(as.character(level_for_no_outcome))))==0){
        stop(sprintf('Error: Invalid or no outcome level(s) given for the absence of the outcome %s!', outcome))
      }
      data[, outcome]<-ifelse(toupper(levels(binary_outcome)[binary_outcome])%in%toupper(as.character(level_for_no_outcome)), 0, 1)
    }
  }
  return(glm(formula = formula, data = data, family = poisson(link = 'log'), control = list(maxit = maxit, epsilon = tol), ...))
}


# 10) Function to extract the p-values of type I and type III ANOVA analyses, 
# effect estimates, standard errors, confidence intervals (log or logit scale
# as well as the response scale) for each dietary intake predictor, parsing
# the text files that include the regression results.

collect_info_of_dietary_vars_from_txt<-function(file, analysis.title, outcome.title, predictor.title){
  out.tables<-list() # An empty list to store the predictor effects (main & interactive) and p-values
  text <- readLines(file) #Read the text file
  text.analyses <- split(text, cumsum(grepl(analysis.title, text)))[-1]  # Split the text into blocks of each analysis setting
  names(text.analyses)<-unlist(lapply(text.analyses, function(x) sub(':$', '', x[1])))
  for(analysis in names(text.analyses)){
    text.outcomes<-split(text.analyses[[analysis]], cumsum(grepl(outcome.title, text.analyses[[analysis]])))[-1]   # Split the text of each analysis setting into blocks for each outcome
    names(text.outcomes)<-lapply(text.outcomes, function(x) sub('^#+\\s+', '', sub(paste0(outcome.title,':\\s+'), '', x[1])))
    for(outcome in names(text.outcomes)){
      text.predictors<-split(text.outcomes[[outcome]], cumsum(grepl(predictor.title, text.outcomes[[outcome]])))[-1] # Split the text of each outcome into blocks for the models of the outcome on each predictor
      names(text.predictors)<-lapply(text.predictors, function(x) sub(paste0(predictor.title,':\\s+'), '', x[1]))
      for(predictor in names(text.predictors)){
        blank_lines <- grep('^[[:blank:]]*$', text.predictors[[predictor]])
        coefficients_start <- grep('^Coefficients|^coefficients', text.predictors[[predictor]])
        coefficients_end <- blank_lines[blank_lines>coefficients_start][1]
        coefficients_table <- as.data.frame(do.call(rbind, lapply((coefficients_start+2):(coefficients_end-1), function(r, table.txt){
          return(strsplit(split='\\s+', table.txt[r])[[1]])
        }, table.txt = text.predictors[[predictor]])))
        rownames(coefficients_table)<-coefficients_table[,1]
        coefficients_table<-coefficients_table[, -c(1), drop=FALSE]
        colnames(coefficients_table)<-sub('Error', 'Std.Error', strsplit(split= '\\s+', text.predictors[[predictor]][coefficients_start+1])[[1]][-c(1,3)])
        coefficients_table<-cbind(coefficients_table, # Add ANOVA p-values to the table
                                  data.frame(TypeI_Pvalue = rep(NA, dim(coefficients_table)[[1]]),
                                             TypeIII_Pvalue = rep(NA, dim(coefficients_table)[[1]])))
        anovaI_start <- grep("Results of type 1 ANOVA analysis", text.predictors[[predictor]])
        anovaIII_start <- grep("Results of type 3 ANOVA analysis", text.predictors[[predictor]])
        ff <- grep(predictor, text.predictors[[predictor]])
        if(length(ff)>0){
          if(length(anovaI_start)>0){
            rr <- ff[ff>anovaI_start & ff<ifelse(length(anovaIII_start)>0, anovaIII_start, Inf)]
            for(r in rr){
              values <- strsplit(split="\\s+", text.predictors[[predictor]][r])[[1]]
              row_indication <- paste(strsplit(split=':|\\*', values[1])[[1]], sep='', collapse='*.*')
              coefficients_table[grepl(row_indication, rownames(coefficients_table)), 'TypeI_Pvalue']<-ifelse(!grepl('^\\.$|^\\*+$', values[length(values)]), values[length(values)], paste(values[length(values)-1]))#, values[length(values)]))
            }
          }
          if(length(anovaIII_start)>0){
            rr <- ff[ff>anovaIII_start]
            for(r in rr){
              values<-strsplit(split="\\s+", text.predictors[[predictor]][r])[[1]]
              row_indication <- paste(strsplit(split=':|\\*', values[1])[[1]], sep='', collapse='*.*')
              coefficients_table[grepl(row_indication, rownames(coefficients_table)), 'TypeIII_Pvalue']<-ifelse(!grepl('^\\.$|^\\*+$', values[length(values)]), values[length(values)], paste(values[length(values)-1]))#, values[length(values)])) 
            }
          }
        }
        predictor.table<-coefficients_table[which(grepl(predictor, rownames(coefficients_table))),,drop=FALSE]
        predictor.table<-cbind(data.frame(Analysis = rep(analysis, dim(predictor.table)[[1]]),
                                          Outcome = rep(outcome, dim(predictor.table)[[1]]),
                                          Predictor = rownames(predictor.table)), predictor.table)
        rownames(predictor.table)<-NULL
        out.tables[[paste(analysis, outcome, predictor, sep = '_')]]<-predictor.table
      }
    }
  }
  return(do.call(rbind, out.tables))
}

# 11) Function to make the first letter of a string a capital letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

# 12) Function to perform Poisson regression for imputed data
poisson_robust_imputed <- function(data, outcome, predictor,
                                   covariates = c(), 
                                   outcome_no_levels = c('No', 0), vcov.type = 'HC1') {
  dd = data
  if(!is.numeric(dd[, outcome])){
    dd[, outcome] <- ifelse(dd[, outcome]%in%outcome_no_levels, 0, 1)
  }
  Imps<-unique(dd$.imp)
  if(length(Imps)==0){ # in case no imputation is done, make the function still work
    dd$.imp = 1
    Imps = 1
  }
  output<-list()
  for(ii in sort(Imps)){
    model <- glm(as.formula(paste(outcome,
                                  paste(c(predictor, covariates), sep='', collapse= '+'), sep='~')),
                 family = poisson(link = 'log'), data = subset(dd, .imp==ii))
    robust_se <- sqrt(diag(vcovHC(model, type = vcov.type)))  
    output[[paste0('.imp_', ii)]] <- data.frame(Estimate = coef(model), 
                                                Robust_SE = robust_se)
  }
  return(output)
}

# 13) Function to pool the results of Poisson regression for multiple imputed
# data using Rubin's rules. Works on the output of Function 12.
pool_results_poisson <- function(result_imp, ci = 0.95){
  if(length(result_imp)>1){
    coef_se_matrix <- do.call(cbind, lapply(1:length(result_imp), function(cc){
      dd <- result_imp[[names(result_imp)[cc]]]
      colnames(dd) <- paste(colnames(dd),
                            names(result_imp)[cc], sep='_')
      return(dd)
    }))
  } else {
    coef_se_matrix <- result_imp[[names(result_imp)[1]]] 
  }
  pooled_coefs <- coef_se_matrix[, !grepl('se', colnames(coef_se_matrix), ignore.case = TRUE), drop=FALSE]
  pooled_var <- coef_se_matrix[, grepl('se', colnames(coef_se_matrix), ignore.case = TRUE), drop=FALSE]^2
  pooled_coefs <- rowMeans(pooled_coefs)
  pooled_var <- rowMeans(pooled_var)
  if(length(result_imp)>1){
    Lambda <- pooled_var
    pooled_var <- pooled_var + (1+1/length(result_imp))*apply(coef_se_matrix[, !grepl('se', colnames(coef_se_matrix), ignore.case = TRUE)], 1, var)
    Lambda <- 1-Lambda/pooled_var
  } else {
    Lambda <- NA
  }
  pooled_se <- sqrt(pooled_var)
  pooled_results <- data.frame(
    Estimate = pooled_coefs,
    Robust_SE = pooled_se,
    Lower_CI = pooled_coefs - qnorm((1-ci)/2, lower.tail = FALSE) * pooled_se,
    Upper_CI = pooled_coefs + qnorm((1-ci)/2, lower.tail = FALSE) * pooled_se,
    pval = pchisq((pooled_coefs/pooled_se)^2, df = 1, lower.tail = FALSE),
    Lambda = Lambda
  )
  return(pooled_results)
}

# 14) Provide a summary of missing in the input dataframe,
# giving N(%) for each variable.
missing_summary <- function(df) {
  missing_counts <- colSums(is.na(df))
  missing_percent <- (missing_counts / dim(df)[[1]]) * 100
  result <- data.frame(
    Variable = names(df),
    Missing = paste0(missing_counts, " (", round(missing_percent, 3), "%)")
  )
  colnames(result)[which(colnames(result)=='Missing')]<-'Missing N(%)'
  return(result)
}
# 15) Add stars (* for <0.05, ** for <0.01 and *** for <0.001) to p-values.
# Round the p-values by d meaningful digits and show all pval<0.001 as '<0.001'.
add_star_to_p <- function(x, d=3) {
  if(is.null(x) || is.na(x)){
    return(x)
  }
  if(!is.numeric(x)){
    stop('Non-numerical p-value is not accepted!')
  }
  if(x<0.001){
    return("<0.001***")
  } else  if(x<0.01){
    return(paste0(round(x, d), '**'))
  } else if(x<0.05){
    return(paste0(round(x, d), '*'))
  } else {
    return(as.character(round(x, d)))
  }
}

# 16) Summarize the dataset in a table. Works only on selected variables and 
# reports for each variable: Total N. Number missing, % missing, 
# If numeric: Mean (SD), If categorical: Counts per level (%).
summarize_vars <- function(data, vars, decimal = 2) {
  results <- lapply(vars, function(v) {
    x <- data[[v]]
    total_n <- length(x)
    missing_n <- sum(is.na(x))
    missing_pct <- round(100 * missing_n / total_n, decimal - 1)
    
    if (is.numeric(x)) {
      summary_val <- paste0(
        round(mean(x, na.rm = TRUE), decimal), 
        " (", round(sd(x, na.rm = TRUE), decimal - 1), ")"
      )
      
      tibble(
        variable = v,
        total_n = total_n,
        missing_n = missing_n,
        missing_pct = missing_pct,
        summary = summary_val
      )
      
    } else {
      counts <- table(x, useNA = "no")
      props <- round(100 * prop.table(counts), decimal - 1)
      cat_summary <- paste0(
        names(counts), ": ", counts, " (", props, "%)",
        collapse = "; "
      )
      
      tibble(
        variable = v,
        total_n = total_n,
        missing_n = missing_n,
        missing_pct = missing_pct,
        summary = cat_summary
      )
    }
  })
  
  bind_rows(results)
}
