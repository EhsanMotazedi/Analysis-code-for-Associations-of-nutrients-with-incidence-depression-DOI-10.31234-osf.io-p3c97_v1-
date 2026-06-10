# Look into the effect of Mediterranean Diet Score (MDS) measured at the
# baseline in the HELIUS study, and later in LASA study, on major depressive disorder (MDD).
# MDD is diagnosed using H1_PHQ9_deprsymp and H2_PHQ9_deprsymp scores at H1 and H2 waves, respectively.
# This analysis is an update of the GLAD_2024 analysis, using the latest update
# of the GLAD consortium guidelines.Specifically, We filter individuals at H1, so
# that only those healthy at the baseline are included. Besides, no clinical diagnosis
# of MDD (H1_MDD and H2_MDD in 2024 analysis) is used at output. Also, marital status
# has now only 4 levels instead of the original 5. The dietary intake values are first
# analysis as such, and later the energy adjusted residuals of the dietary variables are used
# in a sensitivity analysis. Moreover, only Poisson regression analysis is performed
# to obtain risk ratios as the population summary measure of the dietary effects.
# More details about the required models and analyses are documented in:
# "GLAD_Taskforce_Data_Analysis_Plan_24.06.2024.pdf" and
# "confounders and model 2025_EM.docx".
# Previous analyses (GLAD 2024) are done according to:
# "GLAD_Analysis_RegressionModels_Mood_on_EnergyAdjustedResidualDietaryIntake.R"
# Written by Ehsan Motazedi, Amsterdam UMC, 08-03-2025
# Last modified: 27-09-2025
library(foreign)
library(dplyr)
library(mice)
rm(list=ls())
# Specify which parts of the code need to run, e.g. all, only imputation, only writing summary, etc.
do.impute <-  FALSE # Do we need to do imputation or it has been done once before?
run_analyses <- FALSE # Set to TRUE to run the regression analyses, otherwise just load the already fitted models
write.summary <- FALSE # Set to TRUE to write the summary of regression results to the text files given in outfiles
## Change the source path to your own when using
source('C:/Users/P077588/OneDrive - Amsterdam UMC/Documenten/GLAD_HELIUS_Mary_Nicolaou/GLAD_Functions.R')
## Set the parameters for the analyses
alpha = 0.05 # Significance level
decimals = 5 # The number of decimal digits to print in the regression effect estimates
ow = getOption('width')
scipen = 999
# Imputation parameters
Imp_m = 50 # Number of imputations
Imp_max_iter = 300 # Max number of iterations for MI
Imp_burnin = 100 # Burn-in period for MI 
Imp_seed =  710651 # Seed for MI
writefolder <- 'C:\\Users\\P077588\\OneDrive - Amsterdam UMC\\Documenten\\GLAD_HELIUS_Mary_Nicolaou\\New Analysis_March 2025'
datafolder <- 'G:\\divjk\\sg\\GLAD\\Update 2025'
datafilename <-'GLAD_analyses_2025.sav'
AS <- list(MI = 'Imp_Data', CC = 'MyData') # Analysis dataset names for the MI and completers analyses, respectively
studyName <- 'HELIUS' # set to your study name
outfiles<-list(CC = file.path(writefolder, paste0('Poisson_regression_results_', studyName, '_GLAD2025_Completers')),
               MI = file.path(writefolder, paste0('Poisson_regression_results_', studyName, '_GLAD2025_SensitivityMI')))
outcomes <- c(AVALCAT1 = 'H1_PHQ9_deprsymp', # Categorical outcome variable at baseline
              AVALCAT2 = 'H2_PHQ9_deprsymp', # Categorical outcome variable at H2
              AVAL1 = 'H1_PHQ9_sumscore', # Survey sum score (used for imputation) at H1
              AVAL2 = 'H2_PHQ9_sumscore' # Survey sum score (used for imputation) at H2
)
demographics<-c(Age = 'H1_lft', # Name of the baseline covariates used in the analyses
                Sex = 'H1_geslacht',
                `Highest educational level` = 'H1_Opleid',
                `Marital status` = 'H1_Marital_Status',
                `Physical activity` = 'H1_Squash_rlbew',
                `Diabetes` = 'H1_Diabetes_Self',
                `Cerebrovascular accident (CVA)` = 'H1_CVA_Self',
                `Heart attack` = 'H1_Infarct',
                `Alcohol consumption` = 'H1_AlcoholConsumption', 
                Smoking = 'H1_Roken',
                `BMI` = 'H1_LO_BMI')
energy_intake <- c('ENKcal_Sum') # Name of the variable containing energy intake values
# Name of the variable for ethnicity. Set to NULL if ethnicity is not used in the cohort
ethnicity_variable <- 'H1_etniciteit'
# Names of the dietary intake variables in the dataset
dietary_intake_variables <- c('Fruit', 'Vegetables','Legumes', 'Whole_grains','Nuts_seeds', 'Milk','Red_meat', 'Processed_meat',
                              'Sugar_sweetened_beverages','Fibre_Sum', 'Ca_Sum', 'omega_3', # Ca_Sum: Calcium
                              'PUFA_Sum', 'Trans_Fatty_Acids', 'Fish', 'Total_MDS') #PUFA: Polysaturated fatty acids, Total_MDS: Mediterrainian diet score
MyData <- read.spss(file=file.path(datafolder, datafilename), # Read the HELIUS data
                    use.value.labels = FALSE, to.data.frame = TRUE, 
                    use.missings = TRUE, trim_values = FALSE)
for(ff in outcomes){ # Convert the dichotomous outcomes to factors
  aa<-attributes(MyData[, ff])$value.labels
  if(length(aa)>0){
    MyData[, ff]<-factor(MyData[, ff], 
                         levels = unname(aa),
                         labels = names(aa))
    if(length(intersect(toupper(levels(MyData[, ff])), c('NO', '0')))>0){
      MyData[, ff]<-relevel(MyData[, ff], ref = levels(MyData[, ff])[toupper(levels(MyData[, ff]))%in%c('NO', '0')])
    }
  }
}
MyData[, demographics['Highest educational level']] <- factor(MyData[, demographics['Highest educational level']],
                                                              levels = unname(attributes(MyData[, demographics['Highest educational level']])$value.labels),
                                                              labels = names(attributes(MyData[, demographics['Highest educational level']])$value.labels))

MyData[, demographics['Physical activity']] <- factor(MyData[, demographics['Physical activity']],
                                                      levels = unname(attributes(MyData[, demographics['Physical activity']])$value.labels),
                                                      labels = names(attributes(MyData[, demographics['Physical activity']])$value.labels))
MyData[, demographics['Sex']]<-factor(MyData[, demographics['Sex']],
                                      levels = unname(attributes(MyData[, demographics['Sex']])$value.labels),
                                      labels = names(attributes(MyData[, demographics['Sex']])$value.labels))
MyData[, demographics['Sex']] <- relevel(MyData[, demographics['Sex']], ref = 'man')

MyData[, demographics['Marital status']] <- factor(MyData[, demographics['Marital status']],
                                                   levels = unname(attributes(MyData[, demographics['Marital status']])$value.labels),
                                                   labels = names(attributes(MyData[, demographics['Marital status']])$value.labels))
MyData[, demographics['Marital status']] <- relevel(MyData[, demographics['Marital status']],
                                                    ref = 'married/registered partnership/Living together')
MyData[, demographics['Smoking']] <- factor(MyData[, demographics['Smoking']],
                                            levels = unname(attributes(MyData[, demographics['Smoking']])$value.labels),
                                            labels = gsub('\\s+', '_', names(attributes(MyData[, demographics['Smoking']])$value.labels)))
MyData[, demographics['Smoking']] <- relevel(MyData[, demographics['Smoking']], ref = 'Nee,_ik_heb_nooit_gerookt')

for(disease in c('Cerebrovascular accident (CVA)', 'Diabetes', 'Heart attack')){
  MyData[, demographics[disease]] <-  factor(MyData[, demographics[disease]],
                                             levels = unname(attributes(MyData[, demographics[disease]])$value.labels),
                                             labels = names(attributes(MyData[, demographics[disease]])$value.labels))
  MyData[, demographics[disease]] <- relevel(MyData[, demographics[disease]], ref = 'Nee')
}
MyData[, demographics['Alcohol consumption']] <- factor(MyData[, demographics['Alcohol consumption']],
                                                        levels = unname(attributes(MyData[, demographics['Alcohol consumption']])$value.labels),
                                                        labels = gsub(pattern = ' \\(.*', replacement = '', names(attributes(MyData[, demographics['Alcohol consumption']])$value.labels)))
MyData[, demographics['Alcohol consumption']] <- relevel(MyData[, demographics['Alcohol consumption']],
                                                         ref = levels(MyData[, demographics['Alcohol consumption']])[which(grepl('low', levels(MyData[, demographics['Alcohol consumption']])))[1]])

if(length(ethnicity_variable)>0){
  ethnicities<-names(attributes(MyData[, ethnicity_variable])$value.labels)
  # Set non-specific ethnicity to NA
  ethnicities<-unique(c('Nederlands', ethnicities[!ethnicities%in%c('Anders', 'Onbekend', 'Other', 'Unknown')]))
  MyData[, ethnicity_variable]<-factor(MyData[, ethnicity_variable],
                             levels = unname(attributes(MyData[, ethnicity_variable])$value.labels[ethnicities]),
                             labels = ethnicities)
}


# As described in the methods section: 
# implausible energy intake is defined as <500 kcal/day or >3,500 kcal/day for women,
# and <800 kcal/day or >4,000 kcal/day for men. We exclude this from LASA and HELIUS.
MyData[which(MyData[, demographics['Sex']]=='man' & 
              (MyData[, energy_intake]<800 | MyData[, energy_intake]>4000)), energy_intake]<-NA #182 individuals in HELIUS
MyData[which(MyData[, demographics['Sex']]=='vrouw' & 
              (MyData[, energy_intake]<500 | MyData[, energy_intake]>3500)), energy_intake]<-NA #152 individuals in HELIUS
MyData <- MyData[which(!is.na(MyData[, energy_intake])),] # Exclude implausible energy intake levels
MyData <- MyData[which(!is.na(MyData[, outcomes['AVALCAT1']])), ] # Exclude those with unknown depression at the baseline
print(sum(MyData[, outcomes['AVALCAT1']]=='No')) # 4109 with negative depression at baseline
# 4750 individuals remain
# Calculate the energy adjusted residuals for the dietary variables
# using Willett, Howe and Kushi (1997) method
residual_dietary_intake <- paste('resid', dietary_intake_variables, sep = '_')
for(i in 1:length(dietary_intake_variables)){
  MyData[, residual_dietary_intake[i]]<-tryCatch({
    mod <- lm(as.formula(paste(dietary_intake_variables[i], energy_intake, sep='~')), data=MyData)
    newdata<-MyData[1, ] # This newdata is only used in the following predict function to get the lsmeans average of the dietary intake variable
    newdata[, dietary_intake_variables[i]]<-NA # Set it to NA to later predict it
    newdata[, energy_intake]<-mean(MyData[, energy_intake], na.rm=TRUE) # Set the predictor value to the population mean
    resid(mod)+predict(mod, newdata=newdata, type='response') # return 'standardized' residuals by adding the predicted dietary value at the mean energy intake          
  }, error = function(err){
    write(sprintf('An error occurred in getting energy adjusted residuals, with %s predictor and %s for the energy intake:',
                  dietary_intake_variables[i], energy_intake), stderr())
    message(err)
    write('\nRaw dietary intake values will be used!\n--------\n', stderr())
    return(MyData[, dietary_intake_variables[i]])
  }
  )
}

#head(MyData[, c(demographics, dietary_intake_variables, energy_intake, residual_dietary_intake)])
print(missing_summary(MyData[, c(demographics, 
                                 ethnicity_variable,
                                 energy_intake,
                                 dietary_intake_variables,
                                 outcomes)]), row.names = FALSE)
# Full analysis set of multiple imputation (FAS-MI)
for(diet_var in dietary_intake_variables){
  MyData <- MyData[which(!is.na(MyData[, diet_var])), ]
}
print(dim(MyData))
print(sum(!is.na(MyData[, outcomes['AVALCAT2']])))
if(do.impute){## Perform multiple imputation of the data using the full analysis set
  missingness <- apply(MyData[, c(demographics, dietary_intake_variables, energy_intake)],
                       2, FUN = function(x,N) sum(is.na(x))/N*100, N=dim(MyData)[[1]])
  Imputed_vars <- names(missingness[missingness>0])
  predictors <- c(Imputed_vars, c(demographics[c('Age', 'Sex')], 
                                  dietary_intake_variables, energy_intake))
  Imp_matrix <- make.predictorMatrix(MyData)
  Imp_matrix[!rownames(Imp_matrix)%in%predictors, ] <- 0
  Imp_matrix[, !colnames(Imp_matrix)%in%Imputed_vars] <- 0
  Imp_matrix[c(predictors, outcomes['AVAL1']), outcomes['AVAL1']] <- 1
  Imp_matrix[c(predictors, outcomes['AVAL1'], outcomes['AVAL2']), outcomes['AVAL2']] <- 1
  Imp_matrix[c(predictors, outcomes['AVAL1']), outcomes['AVALCAT1']] <- 1 # Impute outcome at the baseline using PHQ9 score at baseline
  Imp_matrix[c(predictors, outcomes['AVAL1'], outcomes['AVAL2']), outcomes['AVALCAT2']] <- 1 # Impute outcome at H2 using PHQ9 score at baseline and at H2
  Imp_methods <- rep('', dim(Imp_matrix)[[2]])
  names(Imp_methods) <- colnames(Imp_matrix)
  Imputed_vars_nlevel <-  unlist(lapply(Imputed_vars, FUN = function(vv) ifelse(is.factor(MyData[, vv]), length(levels(MyData[, vv])), -1)))
  Imp_methods[Imputed_vars] <- ifelse(Imputed_vars_nlevel ==-1, 'pmm', # Set the method to Predicted mean matching for cintinuous vars, polyreg for multinomial and logistic reg for binary vars
                                      ifelse(Imputed_vars_nlevel>2, 'polyreg', 'logreg'))
  Imp_methods[outcomes] <- 'logreg'
  Imp_methods[c('H1_PHQ9_sumscore', 'H2_PHQ9_sumscore')] <- 'pmm'
  Imp_Data <- mice(MyData, 
                   m = Imp_m,                   # Number of imputed datasets
                   method = Imp_methods,    # Imputation methods for each variable
                   predictorMatrix = Imp_matrix, # Predictor matrix to specify relationships
                   maxit = Imp_max_iter,        # Maximum number of iterations
                   burnin = Imp_burnin,         # Burn-in period to discard initial iterations
                   seed = Imp_seed)              # Seed for reproducibility
  saveRDS(Imp_Data, file.path(writefolder, 'Imp_Data.rds'))
} else {
  Imp_Data <- readRDS(file.path(writefolder, 'Imp_Data.rds'))
}
# Imp_Data$imp$H1_Marital_Status # Look at the imputed missing values for marital status
# head(complete(Imp_Data, action='long')) # Look at the aggregated multiple imputed data
Imp_Data <- complete(Imp_Data, action='long')
Imp_Data[which(MyData[, demographics['Sex']]=='man' & MyData[, energy_intake]<800), energy_intake]<-800
Imp_Data[which(MyData[, demographics['Sex']]=='man' & MyData[, energy_intake]>4000), energy_intake]<-4000
Imp_Data[which(MyData[, demographics['Sex']]=='vrouw' & MyData[, energy_intake]<500), energy_intake]<-500
Imp_Data[which(MyData[, demographics['Sex']]=='vrouw' & MyData[, energy_intake]>3500), energy_intake]<-3500

# AFTER the imputation, filter out individuals with depressive symptoms at the baseline
# table(MyData$H1_PHQ9_deprsymp, useNA = "always")
# No  Yes <NA>
# 4109  616 0
# Full analysis set (FAS) for the dietary analysis, augmented by MI
print(dim(MyData[which(MyData[, outcomes['AVALCAT1']]=='No'),]))
Imp_Data <- Imp_Data[with(Imp_Data, ID%in%MyData$ID[which(MyData[, outcomes['AVALCAT1']]=='No')]), ]
# Completers' data: modified full analysis set (mFAS)
MyData <- MyData[which(MyData[, outcomes['AVALCAT1']]=='No' & !is.na(MyData[, outcomes['AVALCAT2']])), ]
w=options('width')$width
options(width=500)
print(summarize_vars(MyData, 
                     unname(c(outcomes, 
                         energy_intake, dietary_intake_variables, 
                         demographics))), width = 100000, n = 100) %>%
        capture.output() %>%
        writeLines(file.path(writefolder, "summary_HELIUS_all.txt"))
print(summarize_vars(MyData[which(MyData[,demographics["Sex"]]=='vrouw'),], 
                     unname(c(outcomes, 
                         energy_intake, dietary_intake_variables, 
                         demographics))), width = 100000, n = 100) %>%
        capture.output() %>%
        writeLines(file.path(writefolder, "summary_HELIUS_women.txt"))
print(summarize_vars(MyData[which(MyData[,demographics["Sex"]]=='man'),], 
                     unname(c(outcomes, 
                         energy_intake, dietary_intake_variables, 
                         demographics))), width = 100000, n = 100) %>%
        capture.output() %>%
        writeLines(file.path(writefolder, "summary_HELIUS_men.txt"))
options(width=w)
### Start the analyses here
### The main analysis must be using the multiplly imputed data (MI)
### Also a sensitivity analysis will be run using only the complete-case data set
analysis_sets <- c('all', levels(MyData[, demographics['Sex']])) # Perform the analyses for the full analysis set and also per gender subgroup 
#analysis_sets <- c(levels(MyData[, demographics['Sex']])) # Perform the analyses only per gender subgroup
for(sex in analysis_sets){
  if(run_analyses){
    Results_1_raw_unadjusted_H2_PHQ9_original <- Results_1_raw_unadjusted_H2_PHQ9_pooled <- lapply(dietary_intake_variables, function(x) NULL)
    names(Results_1_raw_unadjusted_H2_PHQ9_original) <- names(Results_1_raw_unadjusted_H2_PHQ9_pooled) <- dietary_intake_variables
    Results_2_raw_minimallyadjusted_H2_PHQ9_original <- Results_2_raw_minimallyadjusted_H2_PHQ9_pooled <- lapply(dietary_intake_variables, function(x) NULL)
    names(Results_2_raw_minimallyadjusted_H2_PHQ9_original) <- names(Results_2_raw_minimallyadjusted_H2_PHQ9_pooled) <- dietary_intake_variables
    Results_3_energyadjusted_minimallyadjusted_H2_PHQ9_original <- Results_3_energyadjusted_minimallyadjusted_H2_PHQ9_pooled <- lapply(dietary_intake_variables, function(x) NULL)
    names(Results_3_energyadjusted_minimallyadjusted_H2_PHQ9_original) <- names(Results_3_energyadjusted_minimallyadjusted_H2_PHQ9_pooled) <- paste('resid', dietary_intake_variables, sep = '_')  
    Results_4_energyadjusted_adjusted_H2_PHQ9_original <- Results_4_energyadjusted_adjusted_H2_PHQ9_pooled <- lapply(dietary_intake_variables, function(x) NULL)
    names(Results_4_energyadjusted_adjusted_H2_PHQ9_original) <- names(Results_4_energyadjusted_adjusted_H2_PHQ9_pooled) <- paste('resid', dietary_intake_variables, sep = '_')
    Results_5_energyadjusted_adjustedwithBMI_H2_PHQ9_original <- Results_5_energyadjusted_adjustedwithBMI_H2_PHQ9_pooled <- lapply(dietary_intake_variables, function(x) NULL)
    names(Results_5_energyadjusted_adjustedwithBMI_H2_PHQ9_original) <- names(Results_5_energyadjusted_adjustedwithBMI_H2_PHQ9_pooled) <- paste('resid', dietary_intake_variables, sep = '_')
    output_results <- lapply(names(AS), function(x) return(list()))
    names(output_results) <- names(AS)
    for(aset in names(AS)){
      output_results[[aset]] <- list(Results_1_raw_unadjusted_H2_PHQ9_original = Results_1_raw_unadjusted_H2_PHQ9_original,
                                   Results_1_raw_unadjusted_H2_PHQ9_pooled = Results_1_raw_unadjusted_H2_PHQ9_pooled,
                                   Results_2_raw_minimallyadjusted_H2_PHQ9_original = Results_2_raw_minimallyadjusted_H2_PHQ9_original,
                                   Results_2_raw_minimallyadjusted_H2_PHQ9_pooled = Results_2_raw_minimallyadjusted_H2_PHQ9_pooled,
                                   Results_3_energyadjusted_minimallyadjusted_H2_PHQ9_original = Results_3_energyadjusted_minimallyadjusted_H2_PHQ9_original,
                                   Results_3_energyadjusted_minimallyadjusted_H2_PHQ9_pooled = Results_3_energyadjusted_minimallyadjusted_H2_PHQ9_pooled,
                                   Results_4_energyadjusted_adjusted_H2_PHQ9_original = Results_4_energyadjusted_adjusted_H2_PHQ9_original,
                                   Results_4_energyadjusted_adjusted_H2_PHQ9_pooled = Results_4_energyadjusted_adjusted_H2_PHQ9_pooled,
                                   Results_5_energyadjusted_adjustedwithBMI_H2_PHQ9_original = Results_5_energyadjusted_adjustedwithBMI_H2_PHQ9_original,
                                   Results_5_energyadjusted_adjustedwithBMI_H2_PHQ9_pooled = Results_5_energyadjusted_adjustedwithBMI_H2_PHQ9_pooled)
      AnalysisSet = get(AS[[aset]], envir = .GlobalEnv) # Get the data and filter for gender if necessary
      dietary_intake_variables_per_analysis <- dietary_intake_variables # The set of all dietary intake variables will be used for the full analysis set 
      if(sex%in%levels(AnalysisSet[, demographics['Sex']])){
        AnalysisSet <- AnalysisSet[AnalysisSet[, demographics['Sex']]==sex, ]
        #dietary_intake_variables_per_analysis <-  c('Total_MDS') # Perform per gender analysis only for Total_MDS
      }
      print(sprintf('Analysis for: %s-%s', sex, aset))
      for(crude_diet_var in dietary_intake_variables_per_analysis){
        ## a) Unadjusted analysis of raw dietary intake variables
        print(sprintf('Analysis a for: %s', crude_diet_var))
        output_results[[aset]][['Results_1_raw_unadjusted_H2_PHQ9_original']][[crude_diet_var]] <- poisson_robust_imputed(outcome = outcomes['AVALCAT2'],
                                                                                                         predictor = crude_diet_var,
                                                                                                         data = AnalysisSet)
        output_results[[aset]][['Results_1_raw_unadjusted_H2_PHQ9_pooled']][[crude_diet_var]] <- pool_results_poisson(output_results[[aset]][['Results_1_raw_unadjusted_H2_PHQ9_original']][[crude_diet_var]])
        ## b) Minimally adjusted analyses on the raw dietary intake variable (include age, sex, education level, ethnicity as covariates)
        print(sprintf('Analysis b for: %s', crude_diet_var))
        if(sex=='all'){
          output_results[[aset]][['Results_2_raw_minimallyadjusted_H2_PHQ9_original']][[crude_diet_var]] <- poisson_robust_imputed(outcome = outcomes['AVALCAT2'],
                                                                                predictor = crude_diet_var,
                                                                                covariates = c(unname(demographics[c('Age', 'Sex', 'Highest educational level')]), ethnicity_variable),
                                                                                data = AnalysisSet)
        } else {
        output_results[[aset]][['Results_2_raw_minimallyadjusted_H2_PHQ9_original']][[crude_diet_var]] <- poisson_robust_imputed(outcome = outcomes['AVALCAT2'],
                                                                                predictor = crude_diet_var,
                                                                                covariates = c(unname(demographics[c('Age', 'Highest educational level')]), ethnicity_variable),
                                                                                data = AnalysisSet)    
        }                                        
        output_results[[aset]][['Results_2_raw_minimallyadjusted_H2_PHQ9_pooled']][[crude_diet_var]] <- pool_results_poisson(output_results[[aset]][['Results_2_raw_minimallyadjusted_H2_PHQ9_original']][[crude_diet_var]])
      }
      for(resid_diet_var in paste('resid', dietary_intake_variables_per_analysis, sep = '_')){
        ## c) Minimally adjusted analyses on energy adjusted dietary intakes, covariates as in Model (b) as well as total energy intake
        print(sprintf('Analysis c for: %s', resid_diet_var))
        if(sex=='all'){
        output_results[[aset]][['Results_3_energyadjusted_minimallyadjusted_H2_PHQ9_original']][[resid_diet_var]] <- poisson_robust_imputed(outcome = outcomes['AVALCAT2'],
                                                                                                                                          predictor = resid_diet_var,
                                                                                                                                          covariates = c(unname(demographics[c('Age', 'Sex', 'Highest educational level')]), ethnicity_variable, energy_intake),
                                                                                                                                          data = AnalysisSet)
        } else {
          output_results[[aset]][['Results_3_energyadjusted_minimallyadjusted_H2_PHQ9_original']][[resid_diet_var]] <- poisson_robust_imputed(outcome = outcomes['AVALCAT2'],
                                                                                                                                          predictor = resid_diet_var,
                                                                                                                                          covariates = c(unname(demographics[c('Age', 'Highest educational level')]), ethnicity_variable, energy_intake),
                                                                                                                                          data = AnalysisSet)
        }
        output_results[[aset]][['Results_3_energyadjusted_minimallyadjusted_H2_PHQ9_pooled']][[resid_diet_var]] <- pool_results_poisson(output_results[[aset]][['Results_3_energyadjusted_minimallyadjusted_H2_PHQ9_original']][[resid_diet_var]])
        ## d) Adjusted analyses on energy adjusted dietary intakes, covariates as in Model (c),
        ## as well as marital status, PA, chronic disease (diabetes, Heart attack and CVA as three components in model),
        ## alcohol, smoking
        print(sprintf('Analysis d for: %s', resid_diet_var))
        if(sex=='all'){
          output_results[[aset]][['Results_4_energyadjusted_adjusted_H2_PHQ9_original']][[resid_diet_var]] <- poisson_robust_imputed(outcome = outcomes['AVALCAT2'],
                                                                                                                                predictor = resid_diet_var,
                                                                                                                                covariates = c(unname(demographics[c('Age', 'Sex', 'Highest educational level',
                                                                                                                                                         'Marital status', 'Physical activity', 'Diabetes', 
                                                                                                                                                          'Cerebrovascular accident (CVA)', 'Heart attack', 
                                                                                                                                                          'Alcohol consumption', 'Smoking')]), ethnicity_variable, energy_intake),
                                                                                                                                data = AnalysisSet)
        } else {
          output_results[[aset]][['Results_4_energyadjusted_adjusted_H2_PHQ9_original']][[resid_diet_var]] <- poisson_robust_imputed(outcome = outcomes['AVALCAT2'],
                                                                                                                                predictor = resid_diet_var,
                                                                                                                                covariates = c(unname(demographics[c('Age', 'Highest educational level',
                                                                                                                                                  'Marital status', 'Physical activity', 'Diabetes', 
                                                                                                                                                  'Cerebrovascular accident (CVA)', 'Heart attack', 
                                                                                                                                                  'Alcohol consumption', 'Smoking')]), ethnicity_variable, energy_intake),
                                                                                                                                data = AnalysisSet)
          
        }
        output_results[[aset]][['Results_4_energyadjusted_adjusted_H2_PHQ9_pooled']][[resid_diet_var]] <- pool_results_poisson(output_results[[aset]][['Results_4_energyadjusted_adjusted_H2_PHQ9_original']][[resid_diet_var]])
        ## e) Adjusted analyses on energy adjusted dietary intakes, covariates as in Model (d), as well as BMI
        print(sprintf('Analysis e for: %s', resid_diet_var))
        if(sex=='all'){
          output_results[[aset]][['Results_5_energyadjusted_adjustedwithBMI_H2_PHQ9_original']][[resid_diet_var]] <- poisson_robust_imputed(outcome = outcomes['AVALCAT2'],
                                                                                                                                            predictor = resid_diet_var,
                                                                                                                                            covariates = c(unname(demographics),
                                                                                                                                                          ethnicity_variable, energy_intake),
                                                                                                                                            data = AnalysisSet)
                                                                                                                                            
        } else {
          output_results[[aset]][['Results_5_energyadjusted_adjustedwithBMI_H2_PHQ9_original']][[resid_diet_var]] <- poisson_robust_imputed(outcome = outcomes['AVALCAT2'],
                                                                                                                                            predictor = resid_diet_var,
                                                                                                                                            covariates = c(unname(demographics[which(names(demographics)!='Sex')]),
                                                                                                                                                           ethnicity_variable, energy_intake),
                                                                                                                                            data = AnalysisSet)
        }
        output_results[[aset]][['Results_5_energyadjusted_adjustedwithBMI_H2_PHQ9_pooled']][[resid_diet_var]] <- pool_results_poisson(output_results[[aset]][['Results_5_energyadjusted_adjustedwithBMI_H2_PHQ9_original']][[resid_diet_var]])
      }
    }
    saveRDS(output_results, file = file.path(writefolder, paste0('GLAD_Regression_Results_', sex, '.rds')))
    rm(list="output_results")
  }
  output_results <- readRDS(file = file.path(writefolder, paste0('GLAD_Regression_Results_', sex, '.rds'))) # Load the regression models  
  if(write.summary){
    ow = getOption('width')
    options(width = 3*ow)
    for(aset in names(AS)){
      PooledResults <-  names(output_results[[aset]])[grepl('Results_.*_pooled', names(output_results[[aset]]))]
      tables_adjusted <- list()
      for(pooledres in PooledResults){
        tmp <- get(pooledres, envir = as.environment(output_results[[aset]]))
        tmp <- tmp[!unlist(lapply(tmp, is.null))]
        tables_adjusted[[pooledres]] <- do.call(rbind, lapply(names(tmp), function(xx) return(tmp[[xx]] %>% 
                                                                                              filter(grepl(xx, rownames(tmp[[xx]]))))))
        tables_adjusted[[pooledres]][, 'Dietary_var'] <- rownames(tables_adjusted[[pooledres]])
        tables_adjusted[[pooledres]][, 'adjusted_pval_BH'] <- round(p.adjust(tables_adjusted[[pooledres]][, 'pval'],
                                                                            method = 'BH'), 5) 
        #tables_adjusted[[pooledres]] <- tables_adjusted[[pooledres]][order(tables_adjusted[[pooledres]][, 'adjusted_pval_BH']), ]
        tables_adjusted[[pooledres]][, 'adjusted_pval_BH'] <- unlist(lapply(tables_adjusted[[pooledres]][, 'adjusted_pval_BH'], add_star_to_p))
        tables_adjusted[[pooledres]][, 'pval'] <- unlist(lapply(tables_adjusted[[pooledres]][, 'pval'], add_star_to_p))
        tables_adjusted[[pooledres]] <- tables_adjusted[[pooledres]] %>% select(Dietary_var, Estimate, 
                            Robust_SE, Lower_CI, Upper_CI, pval, adjusted_pval_BH, Lambda)
      }
      wfile = paste0(outfiles[[aset]], '_', sex, '.txt')
      sink(wfile, append = FALSE)
      sink()
      nn <- 0 # analysis number
      for(pooledres in PooledResults){
        nn <- nn + 1
        write(sprintf('%s) %s:', letters[nn], pooledres), file = wfile, append = TRUE)
        suppressWarnings(write.table(tables_adjusted[[pooledres]], row.names = FALSE, sep = '\t',
                                   file = wfile, append = TRUE, quote = FALSE))
        write(paste0(rep('*', ow), collapse = ''), file = wfile, append = TRUE)
      }
    }
    rm(list=c("output_results", "tables_adjusted"))
    options(width = ow)
  }
}
