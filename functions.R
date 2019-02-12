# General functions -------------------------------------------------------
on_Midway <- function(){
  if (Sys.getenv('USER')=='pilosofs'){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

detect_locale <- function(){
  require(stringr)
  if (str_detect(Sys.info()[4],'midway2')){return('Midway')}
  if (str_detect(Sys.info()[4],'Shais-')){return('Mac')}
  if (str_detect(Sys.info()[4],'ee-pascual')){return('Lab')}
}

prep.packages <- function(package.list, verbose=T) {
  loaded = package.list %in% .packages()
  if ( all(loaded) ) return(invisible())
  
  package.list = package.list[!loaded]
  installed = package.list %in% .packages(TRUE)
  if ( !all(installed) ) install.packages(package.list[!installed], repos="http://cran.rstudio.com/")
  for ( p in package.list ){
    print(paste("Loading package:",p))
    if (verbose){
      library(p,character.only=TRUE)
    } else {
      suppressMessages(library(p,character.only=TRUE))  
    }
  }
}

collapse_with_commas <- function(x, use_brackets=T){
  if (use_brackets){
    paste('[',paste(x,collapse=','),']',sep='')
  } else {
    paste(x,collapse=',')
  }
}

sink.reset <- function(){
  for(i in seq_len(sink.number())){
    sink(NULL)
  }
}

chunk2 <- function(x,n) {split(x, cut(seq_along(x), n, labels = FALSE))}

## A function used to get the repertoire names from unique repertoire copy
## names.
splitText <- function(str,after=T,splitchar='\\.'){
  if (after){
    return(sapply(strsplit(str, split=splitchar), tail, 1))
  }
  if (after==F){
    return(sapply(strsplit(str, split=splitchar), head, 1))
  }
}

# A function to remove sqlite and parameter files, or any other file for
# specific combinations of parameter space, scenario and experiment. Will remove
# across all runs.
clear_previous_files <- function(parameter_space=NULL, scenario=NULL, experiment=NULL, run = NULL, exclude_sqlite=T, exclude_CP=T, exclude_control=T, test=F){
  files <- list.files(path = '~/Documents/malaria_interventions_data/', full.names = T)
  if (!is.null(parameter_space)){
    files <- files[str_detect(files,paste('PS',parameter_space,sep=''))]
  }
  if (!is.null(scenario)){
    files <- files[str_detect(files,paste('_',scenario,'_',sep='')) | str_detect(files,paste(scenario,'E',sep=''))]
  }
  if (!is.null(experiment)){
    files <- files[str_detect(files,paste('E',experiment,sep=''))]
  }
  if (!is.null(run_range)){
    files <- files[str_detect(files,paste('_R',run,sep=''))]
  }
  if(exclude_CP){
    files <- files[!str_detect(files,'E000')]
  }
  if(exclude_control){
    files <- files[!str_detect(files,'E001')]
  }
  if(exclude_sqlite){
    files <- files[!str_detect(files,'\\.sqlite')]
  }
  if (test){
    print('test mode, not actually removing')
    print(files)
  } else {
    print(files)
    file.remove(files)
  }
}

# Requires package googlesheets
loadExperiments_GoogleSheets <- function(local=F, workBookName='malaria_interventions_design',sheetID=4){
  if (local){
    col_types <- read.csv('~/Documents/malaria_interventions/malaria_interventions_design_col_types.csv',header = T,stringsAsFactors = F)
    col_t <- unname(as.list(col_types[1,]))
    experiments <- read_csv(paste('~/Documents/malaria_interventions/',workBookName,sep=''), col_types = as.list(col_types))
  }
  if (!local){
    GS <- gs_title(workBookName)
    col_types <- GS %>% gs_read(ws=1, col_names=T)
    col_t <- unname(as.list(col_types[1,]))
    experiments <- GS %>% gs_read(ws=sheetID, col_names=T, col_types=col_t)
  }
  print(experiments)
  return(experiments)
}


# Manage parameters in parameter files ------------------------------------

# Function to get reference parameter file and output a data frame of parameters
get_parameter_reference <- function(parameter_file_ref='parameter_file_ref.py'){
  reference <- readLines(parameter_file_ref)
  for (l in 1:length(reference)){
    if (str_detect(reference[l],'#')){
      reference[l] <- str_sub(reference[l],1,str_locate(reference[l],'#')[1]-1)
    }
  }
  
  lines <- 1:length(reference)
  parameters <- map(lines, function(l){
    tmp <- str_sub(reference[l], 1, str_locate(reference[l], '=')[1]-1)
    if(!is.na(tmp)){return(tmp)}
  }) %>% dplyr::combine()
  
  param_values <- map(lines, function(l){
    tmp <- str_trim(str_sub(reference[l], str_locate(reference[l], '=')[1]+1, 10^6), side = 'both')
    if(!is.na(tmp)){return(tmp)}
  }) %>% dplyr::combine()
  
  return(data.frame(param=str_trim(unlist(parameters)), value=unlist(param_values), stringsAsFactors = F))
}

# Function to set a parameter in the parameters data frame
set_parameter <- function(param_data, parameter, value){
  param_data$value[param_data$param==parameter] <- value
  return(subset(param_data, param==parameter))
}

get_random_seed <- function(PS, scenario, experiment='000', run_range, folder){
  seeds <- c()
  for (run in run_range){
    x <- try(readLines(paste(folder,'PS',PS,'_',scenario,'_','E',experiment,'_R',run,'.py',sep='')))
    if (inherits(x, "try-error")){
      print('Cannot find parameter file to extract random number. Most probably a generalized-immunity experiment/run that does not exist.')
      seeds <- c(seeds, NA)
    } else {
      seeds <- c(seeds, parse_number(x[1]))  
    }
  }
  return(seeds)
}

# This function extracts the biting rates from the parameter file, which are in
# daily resolution and returns a vector of the averaged biting rates for a given
# period. (e.g. 30 would be biting rates averaged over a month and 1 will return
# the daily biting rates).
get_biting_rate <- function(parameter_file, sampling_period=30){
  x <- readLines(parameter_file)
  y <- x[grep('BITING_RATE_MEAN',x)[1]]
  BITING_RATE_MEAN <- parse_number(y)
  y <- x[grep('DAILY_BITING_RATE_DISTRIBUTION',x)[1]]
  DAILY_BITING_RATE_DISTRIBUTION <- eval(parse(text=paste('c(',(str_sub(y, str_locate(y, '\\[(.*?)\\]')[1]+1, str_locate(y, '\\[(.*?)\\]')[2]-1)),')',sep='')))
  BITING_RATE <- BITING_RATE_MEAN*DAILY_BITING_RATE_DISTRIBUTION
  BITING_RATE <- chunk2(BITING_RATE, 360/sampling_period)
  sapply(BITING_RATE, mean)
}

# This function creates a data frame with an experimental IRS design. This is a
# convenience function to avoid putting the desing in the Google sheet.
create_intervention_scheme_IRS <- function(PS_benchmark, scenario_benchmark, IRS_START_TIMES=NULL, immigration_range, length_range, coverage_range, poolsize_bounce='False', write_to_file=T, design_ref=NULL){
  # This fixes the 1. in the file name produced in Mathematica
  files <- list.files(path = '~/Documents/malaria_interventions_data/', full.names = T, pattern = '1\\._')
  if(length(files)>0){sapply(files, function(x){file.rename(from = x, to = str_replace(x, '\\.', ''))})}
  
  # PS_benchmark is the parameter space for which intervention is tested. it should already have the checkpoint (000) and control (001) experiments in the google sheet.
  # IRS_START_TIMES can also be something like: '29000,35000'
  if(is.null(design_ref)){
    design_ref <- loadExperiments_GoogleSheets() # Get data design
  }
  design_ref <- subset(design_ref, PS==PS_benchmark & scenario==scenario_benchmark)
  reference_row <- which(grepl(PS_benchmark, design_ref$PS))[2] # reference_row is the row in the design data set which is the reference for this set of IRS exepriments. Should be the line of the control experiment (not the checkpoint)
  
  scheme <- expand.grid(PS=design_ref$PS[reference_row],
                        scenario=design_ref$scenario[reference_row],
                        BITING_RATE_MEAN=design_ref$BITING_RATE_MEAN[reference_row],
                        DAILY_BITING_RATE_DISTRIBUTION=design_ref$DAILY_BITING_RATE_DISTRIBUTION[reference_row],
                        N_GENES_INITIAL=design_ref$N_GENES_INITIAL[reference_row],
                        N_LOCI=design_ref$N_LOCI[reference_row],
                        N_ALLELES_INITIAL=design_ref$N_ALLELES_INITIAL[reference_row],
                        N_POPULATIONS=design_ref$N_POPULATIONS[reference_row],
                        populations_dist=design_ref$populations_dist[reference_row],
                        T_BURNIN=design_ref$T_BURNIN[reference_row],
                        T_END=design_ref$T_END[reference_row],
                        IRS_START_TIMES=IRS_START_TIMES,
                        IRS_IMMIGRATION=immigration_range,
                        IRS_length=length_range,
                        IRS_coverage=coverage_range,
                        MDA_START=NA, # This so to not have an MDA (for function generate_files)
                        POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION=poolsize_bounce,
                        wall_time=design_ref$wall_time[reference_row],
                        mem_per_cpu=design_ref$mem_per_cpu[reference_row],
                        stringsAsFactors = F)
  scheme$exp <- sprintf('%0.3d', 2:(nrow(scheme)+1)) # Start with 002 because 000 and 001 are checkpoint and control
  # Important!!! The files for the IRS in the next line are generated separately in Mathematica, and should already exist.
  scheme$IRS_input <- paste('IRS_120_300_',scheme$IRS_coverage,'_',scheme$IRS_length,'.csv',sep='')
  if (write_to_file){
    write_csv(scheme, paste('PS',scheme$PS[1],'_',scheme$scenario[1],'_IRS_scheme.csv',sep=''))
  }
  return(scheme)
}

# Manage parameters in neutral models -------------------------------------

# This function obtains the duration of infection from the selection mode
# counterpart simulations to calculate the var switching rate in the neutral
# scenario. it makes more sense to take the duration from the control experiment
# (001), otherwise interventions or transience can affect it, so I set
# exepriment=001 as default. sqlite_path_global is a global parameter with the
# path to the sqlite files. It is needed because different projects may have
# different paths to keep sqlite files.
set_transition_rate_neutral <- function(parameter_space, run, N_GENES_PER_STRAIN = 60){
  if (on_Midway()){
    sqlite_file <- list.files(path = 'sqlite/', pattern=paste('PS',parameter_space,'_S_E001_R',run,'.sqlite',sep=''), full.names = T)
  } else {
    sqlite_file <- list.files(path = sqlite_path_global, pattern=paste('PS',parameter_space,'_S_E001_R',run,'.sqlite',sep=''), full.names = T)
  }
  db <- dbConnect(SQLite(), dbname = sqlite_file)
  sampled_duration <- dbGetQuery(db, 'SELECT duration FROM sampled_duration')
  mean_doi <- mean(sampled_duration$duration-14)
  TRANSITION_RATE_NOT_IMMUNE <- 1/(mean_doi/N_GENES_PER_STRAIN)
  # setwd('~/Documents/malaria_interventions/')
  dbDisconnect(db)
  return(TRANSITION_RATE_NOT_IMMUNE)
}

# A function to obtain the transition rate from an existing checkpoint parameter
# file.
get_transition_rate_neutral <- function(PS, run){
  transition_rate <- c()
  x <- try(readLines(paste(parameter_files_path_global,'/','PS',PS,'_N_','E000_R',run,'.py',sep='')))
  if (inherits(x, "try-error")){
    print(paste('Cannot find parameter file to extract transition rate! (','PS',PS,'_N_','E000_R',run,'.py)',sep=''))
    transition_rate <- NA
  } else {
    transition_rate <- parse_number(x[10])
  }
  return(transition_rate)
}

#A function to set the parameters for generazlied immunity scenario. Unlike the
#set_transition_rate_neutral, which returns a single value, this function
#returns a list which describes a fit of a curve of a particular run. It makes
#most sense to fit the curve to an experiment without interevntions and in
#stable state, so 001 is by default.
# REQUIRES rPython
set_generalized_immunity <- function(parameter_space, run){
  # Prepare the python file for the curve-fittign code
  if (on_Midway()){
    sqlite_file <- paste('sqlite/','PS',parameter_space,'_S_E001_R',run,'.sqlite',sep='')
    pyFile <- readLines('generalized_immunity_fitting.py')
  } else {
    sqlite_file <- paste(sqlite_path_global,'/','PS',parameter_space,'_S_E001_R',run,'.sqlite',sep='')
    pyFile <- readLines('/home/shai/Documents/malaria_interventions/generalized_immunity_fitting.py')
  }
  pyFile[11] <- paste('path=','"',sqlite_file,'"',sep='')
  writeLines(pyFile,'generalized_immunity_fitting_experiment.py')
  
  # Try to run the code. Sometimes the curve does not fit
  res <- try(python.load('generalized_immunity_fitting_experiment.py'))  
  if(inherits(res, "try-error")){
    return(NA) # If fit is not sucessful return a NA
  } else {
    N_INFECTIONS_FOR_GENERAL_IMMUNITY <- python.get('infectionTimesToImmune')
    GENERAL_IMMUNITY_PARAMS <- python.get('generalImmunityParams')
    GENERAL_IMMUNITY_PARAMS <- paste('[',paste(GENERAL_IMMUNITY_PARAMS, collapse = ','),']',sep='')
    CLEARANCE_RATE_IMMUNE <- python.get('clearanceRateConstantImmune')
    # file.remove('generalized_immunity_fitting_experiment.py')
    return(list(N_INFECTIONS_FOR_GENERAL_IMMUNITY=N_INFECTIONS_FOR_GENERAL_IMMUNITY,
                GENERAL_IMMUNITY_PARAMS=GENERAL_IMMUNITY_PARAMS,
                CLEARANCE_RATE_IMMUNE=CLEARANCE_RATE_IMMUNE))
  }
}

# A function to obtain the GI parameters from an existing checkpoint parameter
# file.
get_generalized_immunity <- function(PS, run){
  x <- try(readLines(paste(parameter_files_path_global,'/','PS',PS,'_G_','E000_R',run,'.py',sep='')))
  if (inherits(x, "try-error")){
    print(paste('Cannot find parameter file to extract parameters! (','PS',PS,'_G_','E000_R',run,'.py)',sep=''))
    params <- NA
  } else {
    params <- list(N_INFECTIONS_FOR_GENERAL_IMMUNITY=parse_number(x[7]),
                   GENERAL_IMMUNITY_PARAMS=str_sub(x[8],25,str_length(x[8])),
                   CLEARANCE_RATE_IMMUNE=parse_number(x[9]))
  }
  return(params)
}

# This function sets the parameters for a single IRS. It is possible to run it
# several times, once for each IRS scheme. It uses a parameter file which has
# already has the rest of the parameters in place, and just changes the relevant
# parameters for IRS.
set_IRS <- function(design_ID, run, IRS_START_TIME, IRS_input, IRS_IMMIGRATION_RATE_FACTORS, POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION, experimental_design){
  # Regime
  parameter_space <- experimental_design$PS[design_ID]
  scenario <- experimental_design$scenario[design_ID]
  experiment <- experimental_design$exp[design_ID] # Use 000 for checkpoints and 001 for control
  base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
  output_file <- paste(base_name,'.py',sep = '')
  param_data <- get_parameter_reference(paste(parameter_files_path_global,output_file,sep='/'))
  
  # Turn on IRS
  param_data[param_data$param=='IRS_ON',] <- set_parameter(param_data, 'IRS_ON', 'True')
  
  # Add the start time to the vector IRS_START_TIMES
  x <- substr(param_data[param_data$param=='IRS_START_TIMES',]$value, 1, nchar(param_data[param_data$param=='IRS_START_TIMES',]$value)-1)
  x <- paste(x,',',IRS_START_TIME,']',sep='')
  x <- str_replace(x, '\\[,','\\[')
  param_data[param_data$param=='IRS_START_TIMES',] <- set_parameter(param_data, 'IRS_START_TIMES', x)
  # Add the IRS scheme to the vector BITING_RATE_FACTORS
  IRS_scheme <- read_csv(IRS_input, col_names = c('day','num_mosquitos'))
  x <- substr(param_data[param_data$param=='BITING_RATE_FACTORS',]$value, 1, nchar(param_data[param_data$param=='BITING_RATE_FACTORS',]$value)-1)
  tmp1 <- paste(IRS_scheme$num_mosquitos, collapse=',')
  tmp2 <- paste('[',tmp1,']',sep='')
  x <- paste(x,',',tmp2,']',sep='')
  x <- str_replace(x, '\\[,','\\[')
  param_data[param_data$param=='BITING_RATE_FACTORS',] <- set_parameter(param_data, 'BITING_RATE_FACTORS', x)
  # Add to the IRS_IMMIGRATION_RATE_FACTORS
  x <- substr(param_data[param_data$param=='IRS_IMMIGRATION_RATE_FACTORS',]$value, 1, nchar(param_data[param_data$param=='IRS_IMMIGRATION_RATE_FACTORS',]$value)-1)
  x <- paste(x,',',IRS_IMMIGRATION_RATE_FACTORS,']',sep='')
  x <- str_replace(x, '\\[,','\\[')
  param_data[param_data$param=='IRS_IMMIGRATION_RATE_FACTORS',] <- set_parameter(param_data, 'IRS_IMMIGRATION_RATE_FACTORS', x)
  # Add POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION
  param_data[param_data$param=='POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION',] <- set_parameter(param_data, 'POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION', POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION)
  
  
  # Write parameter file
  param_data$output <- paste(param_data$param,param_data$value,sep='=')
  write_lines(param_data$output, paste(parameter_files_path_global,output_file,sep='/'))
}

# This function sets the parameters for a single MDA. It uses a parameter file which has
# already has the rest of the parameters in place, and just changes the relevant
# parameters for MDA.
set_MDA <- function(design_ID, run, experimental_design){
  # Regime
  parameter_space <- experimental_design$PS[design_ID]
  scenario <- experimental_design$scenario[design_ID]
  experiment <- experimental_design$exp[design_ID] # Use 000 for checkpoints and 001 for control
  base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
  output_file <- paste(base_name,'.py',sep = '')
  param_data <- get_parameter_reference(output_file)
  
  # Turn on MDA
  param_data[param_data$param=='MDA_ON',] <- set_parameter(param_data, 'MDA_ON', 'True')
  
  # Set parameters
  param_data[param_data$param=='MDA_START_TIMES',] <- set_parameter(param_data, 'MDA_START_TIMES', paste('[',experimental_design$MDA_START[design_ID],']',sep=''))
  param_data[param_data$param=='HOST_FAIL_RATE',] <- set_parameter(param_data, 'HOST_FAIL_RATE', paste('[',experimental_design$HOST_FAIL_RATE[design_ID],']',sep=''))
  param_data[param_data$param=='DRUG_EFF_DURATION',] <- set_parameter(param_data, 'DRUG_EFF_DURATION', paste('[',experimental_design$DRUG_EFF_DURATION[design_ID],']',sep=''))
  param_data[param_data$param=='MDA_IMMIGRATION_RATE_FACTORS',] <- set_parameter(param_data, 'MDA_IMMIGRATION_RATE_FACTORS', paste('[',experimental_design$MDA_IMMIGRATION_RATE_FACTORS[design_ID],']',sep=''))
  
  # Write parameter file
  param_data$output <- paste(param_data$param,param_data$value,sep='=')
  write_lines(param_data$output, output_file)
}



# File generation ---------------------------------------------------------

# Function to create the necesary files and pipeline for a single run of an experiment.
# Each run has its own random seed across experiments and scenarios.
create_run <- function(design_ID, run, RANDOM_SEED, experimental_design, biting_rate_mathematica=NULL, params_GI=NULL, target_folder=NULL){
  if(is.null(target_folder)){target_folder <- getwd()}
  
  # Regime
  parameter_space <- experimental_design$PS[design_ID]
  scenario <- experimental_design$scenario[design_ID]
  experiment <- experimental_design$exp[design_ID] # Use 000 for checkpoints and 001 for control
  base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
  
  if (detect_locale()=='Midway'){param_data <- get_parameter_reference('parameter_file_ref.py')}
  if (detect_locale()=='Mac'){param_data <- get_parameter_reference('~/GitHub/malaria_interventions/parameter_file_ref.py')}
  if (detect_locale()=='Lab'){param_data <- get_parameter_reference('~/Documents/malaria_interventions/parameter_file_ref.py')}

  # General parameters
  param_data[param_data$param=='RANDOM_SEED',] <- set_parameter(param_data, 'RANDOM_SEED', RANDOM_SEED)
  T_END <- experimental_design$T_END[design_ID]
  param_data[param_data$param=='T_END',] <- set_parameter(param_data, 'T_END', T_END)
  param_data[param_data$param=='VERIFICATION_ON',] <- set_parameter(param_data, 'VERIFICATION_ON', 'False')
  param_data[param_data$param=='VERIFICATION_PERIOD',] <- set_parameter(param_data, 'VERIFICATION_PERIOD', T_END)
  
  if(scenario=='S'){
  param_data[param_data$param=='TRANSITION_RATE_NOT_IMMUNE',] <- set_parameter(param_data, 'TRANSITION_RATE_NOT_IMMUNE', experimental_design$TRANSITION_RATE_NOT_IMMUNE[design_ID])
  }
  
  # Scenario
  if(scenario=='N'){
    param_data[param_data$param=='SELECTION_MODE',] <- set_parameter(param_data, 'SELECTION_MODE', "\'NEUTRALITY\'")
    # For each run, there is a need to calculate the transition rate once, for
    # the CP experiment. the rest of the experiments use the same value.
    if (experiment=='000'){ 
      print(paste('Calculating transition rate from PS',parameter_space,'_S_E001_R',run,sep = ''))
      TRANSITION_RATE_NOT_IMMUNE <- set_transition_rate_neutral(parameter_space, run)
    } else {
      print(paste('Getting transition rate from PS',parameter_space,'_N_E000_R',run,sep = ''))
      TRANSITION_RATE_NOT_IMMUNE <- get_transition_rate_neutral(parameter_space, run)
    }
    param_data[param_data$param=='TRANSITION_RATE_NOT_IMMUNE',] <- set_parameter(param_data, 'TRANSITION_RATE_NOT_IMMUNE', TRANSITION_RATE_NOT_IMMUNE)
  }
  if(scenario=='G'){
    param_data[param_data$param=='SELECTION_MODE',] <- set_parameter(param_data, 'SELECTION_MODE', "\'GENERAL_IMMUNITY\'")
    if (is.null(params_GI)){
      # For each run, there is a need to calculate the general immunity parameters once, for
      # the CP experiment. the rest of the experiments use the same values.
      if (experiment=='000'){ 
        print(paste('Calculating generalized immunity parameters from PS',parameter_space,'_S_E001_R',run,sep = ''))
        params_GI <- set_generalized_immunity(parameter_space, run)
      } else {
        print(paste('Getting generalized immunity parameters from PS',parameter_space,'_G_E000_R',run,sep = ''))
        params_GI <- get_generalized_immunity(parameter_space, run)
      }
      if (is.na(params_GI)[1]){
        print(paste('ERROR (create_run): Could not fit function for generalized immunity, skipping and NOT PRODUCING FILE ',base_name,'.py',sep = ''))
        return(NULL)
      }
      param_data[param_data$param=='GENERAL_IMMUNITY_PARAMS',] <- set_parameter(param_data, 'GENERAL_IMMUNITY_PARAMS', params_GI$GENERAL_IMMUNITY_PARAMS)
      param_data[param_data$param=='N_INFECTIONS_FOR_GENERAL_IMMUNITY',] <- set_parameter(param_data, 'N_INFECTIONS_FOR_GENERAL_IMMUNITY', params_GI$N_INFECTIONS_FOR_GENERAL_IMMUNITY)
      param_data[param_data$param=='CLEARANCE_RATE_IMMUNE',] <- set_parameter(param_data, 'CLEARANCE_RATE_IMMUNE', params_GI$CLEARANCE_RATE_IMMUNE)
    } else {
      print(paste('Using external GI parameters for file ',base_name,'.py',sep = ''))
      param_data[param_data$param=='GENERAL_IMMUNITY_PARAMS',] <- set_parameter(param_data, 'GENERAL_IMMUNITY_PARAMS', params_GI$GENERAL_IMMUNITY_PARAMS)
      param_data[param_data$param=='N_INFECTIONS_FOR_GENERAL_IMMUNITY',] <- set_parameter(param_data, 'N_INFECTIONS_FOR_GENERAL_IMMUNITY', params_GI$N_INFECTIONS_FOR_GENERAL_IMMUNITY)
      param_data[param_data$param=='CLEARANCE_RATE_IMMUNE',] <- set_parameter(param_data, 'CLEARANCE_RATE_IMMUNE', params_GI$CLEARANCE_RATE_IMMUNE)
    }
  }
  # Populations. This duplicates the following parameters as the number of populations.
  N_POPULATIONS <- experimental_design$N_POPULATIONS[design_ID]
  param_data[param_data$param=='N_POPULATIONS',] <- set_parameter(param_data, 'N_POPULATIONS', N_POPULATIONS)
  # This creates the distance matrix, assuming equal distances between all populations (parameter populations_dist)
  DISTANCE_MAT <- matrix(experimental_design$populations_dist[design_ID], ncol=N_POPULATIONS, nrow=N_POPULATIONS)
  diag(DISTANCE_MAT) <- 1
  DISTANCE_MAT_vectorized <- c()
  for (i in 1:N_POPULATIONS){
    DISTANCE_MAT_vectorized <- paste(DISTANCE_MAT_vectorized, collapse_with_commas(DISTANCE_MAT[i,]), sep=',')
  }
  DISTANCE_MAT <- paste(str_replace(DISTANCE_MAT_vectorized,',','['),']',sep='')
  param_data[param_data$param=='DISTANCE_MAT',] <- set_parameter(param_data, 'DISTANCE_MAT', DISTANCE_MAT)
  # Here are some parameters with fixed values
  param_data[param_data$param=='N_HOSTS',] <- set_parameter(param_data, 'N_HOSTS', collapse_with_commas(rep(10000,N_POPULATIONS)))
  param_data[param_data$param=='N_INITIAL_INFECTIONS',] <- set_parameter(param_data, 'N_INITIAL_INFECTIONS', collapse_with_commas(rep(20,N_POPULATIONS)))
  param_data[param_data$param=='BITING_RATE_RELATIVE_AMPLITUDE',] <- set_parameter(param_data, 'BITING_RATE_RELATIVE_AMPLITUDE', collapse_with_commas(rep(0,N_POPULATIONS)))
  param_data[param_data$param=='BITING_RATE_PEAK_PHASE',] <- set_parameter(param_data, 'BITING_RATE_PEAK_PHASE', collapse_with_commas(rep(0,N_POPULATIONS)))
  param_data[param_data$param=='IMMIGRATION_RATE',] <- set_parameter(param_data, 'IMMIGRATION_RATE', collapse_with_commas(rep(1,N_POPULATIONS)))
  
  # Biting rates
  BITING_RATE_MEAN <- experimental_design$BITING_RATE_MEAN[design_ID]
  param_data[param_data$param=='BITING_RATE_MEAN',] <- set_parameter(param_data, 'BITING_RATE_MEAN', paste('[',paste(rep(BITING_RATE_MEAN,N_POPULATIONS),collapse = ','),']',sep=''))
  # It is possible to provide a specific file for biting rates for each
  # experiments, but usually that is not necessary so no need to re-read this
  # when producing batches of files. Just provide it in the call to the
  # function.
  if (is.null(biting_rate_mathematica)){
    mathematica_file <- experimental_design$DAILY_BITING_RATE_DISTRIBUTION[design_ID]
    biting_rate_mathematica <- read_csv(mathematica_file, col_names = c('day','num_mosquitos'))
  }
  DAILY_BITING_RATE_DISTRIBUTION <- biting_rate_mathematica$num_mosquitos
  param_data[param_data$param=='DAILY_BITING_RATE_DISTRIBUTION',] <- set_parameter(param_data, 'DAILY_BITING_RATE_DISTRIBUTION', paste('[',paste(DAILY_BITING_RATE_DISTRIBUTION, collapse=','),']',sep=''))
  
  # Genetic diversity
  N_GENES_INITIAL <- experimental_design$N_GENES_INITIAL[design_ID]
  param_data[param_data$param=='N_GENES_INITIAL',] <- set_parameter(param_data, 'N_GENES_INITIAL', N_GENES_INITIAL)
  N_ALLELES_INITIAL <- experimental_design$N_ALLELES_INITIAL[design_ID]
  param_data[param_data$param=='N_ALLELES_INITIAL',] <- set_parameter(param_data, 'N_ALLELES_INITIAL', paste('N_LOCI*[',N_ALLELES_INITIAL,']',sep=''))
  
  # Checkpoints
  T_BURNIN <- experimental_design$T_BURNIN[design_ID]
  if (experiment == '000'){
    param_data[param_data$param=='SAMPLE_DB_FILENAME',] <- set_parameter(param_data, 'SAMPLE_DB_FILENAME', paste('\'\"',base_name,'.sqlite\"\'',sep='')) # The run ID will be added while running the job (in the sbatch execution).
    param_data[param_data$param=='SAVE_TO_CHECKPOINT',] <- set_parameter(param_data, 'SAVE_TO_CHECKPOINT', 'True')
    param_data[param_data$param=='CHECKPOINT_SAVE_PERIOD',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_PERIOD', T_END) # The save period should be the T_END
    param_data[param_data$param=='CHECKPOINT_SAVE_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_FILENAME', paste('\'\"',base_name,'_CP.sqlite\"\'',sep=''))
    param_data[param_data$param=='LOAD_FROM_CHECKPOINT',] <- set_parameter(param_data, 'LOAD_FROM_CHECKPOINT', 'False')
    param_data[param_data$param=='CHECKPOINT_LOAD_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_LOAD_FILENAME', '\'\"\"\'')
    param_data[param_data$param=='T_BURNIN',] <- set_parameter(param_data, 'T_BURNIN', T_BURNIN)
  }
  if (experiment != '000'){ # This section prepares a parameter file to load a checkpoint and run an experiment
    param_data[param_data$param=='SAMPLE_DB_FILENAME',] <- set_parameter(param_data, 'SAMPLE_DB_FILENAME', paste('\'\"',base_name,'.sqlite\"\'',sep='')) # The run ID will be added while running the job (in the sbatch execution).
    param_data[param_data$param=='SAVE_TO_CHECKPOINT',] <- set_parameter(param_data, 'SAVE_TO_CHECKPOINT', 'False')
    param_data[param_data$param=='CHECKPOINT_SAVE_PERIOD',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_PERIOD', 0)
    param_data[param_data$param=='CHECKPOINT_SAVE_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_FILENAME', '\'\"\"\'')
    param_data[param_data$param=='LOAD_FROM_CHECKPOINT',] <- set_parameter(param_data, 'LOAD_FROM_CHECKPOINT', 'True')
    param_data[param_data$param=='CHECKPOINT_LOAD_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_LOAD_FILENAME', paste('\'\"PS',parameter_space,'_',scenario,'_E000','_R',run,'_CP.sqlite\"\'',sep=''))
    param_data[param_data$param=='T_BURNIN',] <- set_parameter(param_data, 'T_BURNIN', T_BURNIN) # The burnin value should be the value where the checkpoint was taken
  }
  
  # Write parameter file
  output_file=paste(target_folder,'/',base_name,'.py',sep = '')
  param_data$output <- paste(param_data$param,param_data$value,sep='=')
  write_lines(param_data$output, output_file)
}

# This function generates parameter files and corresponding job files (to run on
# Midway), for several experiments and runs. It keeps the random seed for each
# RUN across experiments AND PARAMETER SPACES the same. Also possible to provide
# a random seed. row_range is the row numbers in the design data frame.
generate_files <- function(row_range, run_range, random_seed=NULL, 
                           experimental_design, 
                           biting_rate_file='mosquito_population_seasonality.csv', 
                           params_GI=NULL, 
                           target_folder=NULL,
                           only_sbatch=F){
  
  if(is.null(target_folder)){target_folder <- getwd()}
  
  if (!is.null(biting_rate_file)){
    biting_rate_mathematica <- read_csv(biting_rate_file, col_names = c('day','num_mosquitos'))
  }
  
  # Using this data structure makes it easier to control specific combinations
  # of runs in experimets, for eaxmple for parameter files that cannot be
  # produced because of ill fitting of curves in the generlized immunity.
  cases <- expand.grid(run = run_range, design_ID=row_range, file_created=T)
  
  # Random seeds are equal for each run across experiments AND PARAMETER SPACES
  # (e.g., the same seed for run 1 in experiments '000', '001', '002' in
  # parameter spaces 01 and 02).
  if (is.null(random_seed)){
    cases$seed <- rep(round(runif(n = length(run_range), min=1, max=10^7),0),length(row_range))
  } else {
    cases$seed <- rep(random_seed,length(row_range))
  }
  
  cases$param_file <- NA
  
  if (!only_sbatch){
    for (idx in 1:nrow(cases)){
      RANDOM_SEED <- cases$seed[idx]
      RUN <- cases$run[idx]
      design_ID <- cases$design_ID[idx]
      # print(paste('Run: ',RUN,' | Row: ',design_ID,sep=''))
      
      # Create parameter file with main parameters. If the scenario is generalized
      # immunity then the fit for some runs may have not convereged (functions
      # 'create_run' and 'set_generalized_immunity'). Parameter files cannot be
      # produced in these cases and so the function skips to the next case.
      try_run <- create_run(design_ID, RUN, RANDOM_SEED, experimental_design, biting_rate_mathematica, params_GI=params_GI, target_folder = target_folder)
      if (is.null(try_run)){
        cases$file_created[idx] <- F
        next
      }
      
      # Set IRS
      if (!is.na(experimental_design$IRS_START_TIMES[design_ID])){
        IRS_scheme <- data.frame(IRS_START_TIME=str_split(experimental_design$IRS_START_TIMES[design_ID], ',')[[1]],
                                 IRS_input=str_split(experimental_design$IRS_input[design_ID], ',')[[1]],
                                 IRS_IMMIGRATION_RATE_FACTORS=str_split(experimental_design$IRS_IMMIGRATION_RATE_FACTORS[design_ID], ',')[[1]],
                                 POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION=experimental_design$POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION[design_ID],
                                 stringsAsFactors = F)
        for (i in 1:nrow(IRS_scheme)){
          set_IRS(design_ID = design_ID, 
                  run = RUN, 
                  IRS_START_TIME = IRS_scheme$IRS_START_TIME[i], 
                  IRS_IMMIGRATION_RATE_FACTORS = IRS_scheme$IRS_IMMIGRATION_RATE_FACTORS[i], 
                  IRS_input = IRS_scheme$IRS_input[i], 
                  POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION=IRS_scheme$POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION[i], 
                  experimental_design)
        }
      }
      # Set MDA
      if (!is.na(experimental_design$MDA_START[design_ID])){
        set_MDA(design_ID, RUN, experimental_design)
      }
      
      cases$param_file[idx] <- paste('PS',experimental_design$PS[design_ID],'_',experimental_design$scenario[design_ID],'_E',experimental_design$exp[design_ID],'_R',RUN,sep='')
    }
  }
  
  # Generate batch file for each experiment
  for (design_ID in row_range){
    parameter_space <- experimental_design$PS[design_ID]
    scenario <- experimental_design$scenario[design_ID]
    experiment <- experimental_design$exp[design_ID] # Use 000 for checkpoints and 001 for control
    base_name <- paste('PS',parameter_space,scenario,'E',experiment,sep='')
    
    SLURM_ARRAY_RANGE <- cases$file_created[which(cases$design_ID==design_ID)]
    if (any(SLURM_ARRAY_RANGE==F)){  # If some runs are missing (e.g., runs for which the fit has not converged in generalized immunity)
      SLURM_ARRAY_RANGE <-  paste("\'",collapse_with_commas(run_range[SLURM_ARRAY_RANGE],F),"\'",sep='')
    } else {
      SLURM_ARRAY_RANGE <- paste("\'",min(run_range),'-',max(run_range),"\'",sep='')
    }
    
    # Write the job file for the exepriment
    if (detect_locale()=='Midway'){job_lines <- readLines('job_file_ref.sbatch')}
    if (detect_locale()=='Mac'){job_lines <- readLines('~/GitHub/malaria_interventions/job_file_ref.sbatch')}
    if (detect_locale()=='Lab'){job_lines <- readLines('~/Documents/malaria_interventions/job_file_ref.sbatch')}
    wall_time <-  experimental_design$wall_time[design_ID]
    mem_per_cpu <-  experimental_design$mem_per_cpu[design_ID]
    job_lines[2] <- paste('#SBATCH --job-name=',paste(parameter_space,scenario,'E',experiment,sep=''),sep='')
    job_lines[3] <- paste('#SBATCH --time=',wall_time,sep='')
    job_lines[4] <- paste('#SBATCH --output=slurm_output/',base_name,'_%A_%a.out',sep='')
    job_lines[5] <- paste('#SBATCH --error=slurm_output/',base_name,'_%A_%a.err',sep='')
    job_lines[6] <- paste('#SBATCH --array=',SLURM_ARRAY_RANGE,sep='')
    job_lines[9] <- paste('#SBATCH --mem-per-cpu=',mem_per_cpu,sep='')
    job_lines[19] <- paste("PS='",parameter_space,"'",sep='')
    job_lines[20] <- paste("scenario='",scenario,"'",sep='')
    job_lines[21] <- paste("exp='",experiment,"'",sep='')
    
    if (str_detect(target_folder,'PLOS_Biol')){
      job_lines[23] <- "base_folder='/scratch/midway2/pilosofs/PLOS_Biol/'"
    }
    
    if (experiment == '000'){
      job_lines[26] <- "CHECKPOINT='create'"
    }
    if (experiment != '000'){
      job_lines[26] <- "CHECKPOINT='load'"
    }
    output_file=paste(target_folder,'/',base_name,'.sbatch',sep = '')
    write_lines(job_lines, output_file)
  }
  
  print(cases)
  cat('Run this on Midway: \n')
  for (e in row_range){
    cat(paste('sbatch PS',experimental_design$PS[e],experimental_design$scenario[e],'E',experimental_design$exp[e],'.sbatch','\n',sep=''))
  }
}

# A function to generate sbatch files on demand
generate_sbatch <- function(ps, scen, experiment, runs, unzip_py_files){
  base_name <- paste('PS',ps,scen,'E',experiment,sep='')
  SLURM_ARRAY_RANGE <- collapse_with_commas(runs, use_brackets = F)
  
  # This is to get the correct walltime and memory for Midway
  if (experiment == '000'){
    experimental_design <- subset(design, PS==ps & scenario==scen & exp=='000')
  }
  if (experiment != '000'){
    experimental_design <- subset(design, PS==ps & scenario==scen & exp=='001')
  }  
  
  # Write the job file for the exepriment
  job_lines <- readLines('~/Documents/malaria_interventions/job_file_ref.sbatch')
  wall_time <-  experimental_design$wall_time
  mem_per_cpu <-  experimental_design$mem_per_cpu
  job_lines[2] <- paste('#SBATCH --job-name=',paste(ps,scen,'E',experiment,sep=''),sep='')
  job_lines[3] <- paste('#SBATCH --time=',wall_time,sep='')
  job_lines[4] <- paste('#SBATCH --output=slurm_output/',base_name,'_%A_%a.out',sep='')
  job_lines[5] <- paste('#SBATCH --error=slurm_output/',base_name,'_%A_%a.err',sep='')
  job_lines[6] <- paste('#SBATCH --array=',SLURM_ARRAY_RANGE,sep='')
  job_lines[9] <- paste('#SBATCH --mem-per-cpu=',mem_per_cpu,sep='')
  job_lines[19] <- paste("PS='",ps,"'",sep='')
  job_lines[20] <- paste("scenario='",scen,"'",sep='')
  job_lines[21] <- paste("exp='",experiment,"'",sep='')
  if (experiment == '000'){
    job_lines[26] <- "CHECKPOINT='create'"
  }
  if (experiment != '000'){
    job_lines[26] <- "CHECKPOINT='load'"
  }  
  
  output_file=paste(base_name,'.sbatch',sep = '')
  write_lines(job_lines, output_file)
  cat(paste('sbatch ',base_name,'.sbatch',sep=''));cat('\n')
  
  if (unzip_py_files){
    unzip(paste(ps,'_',scen,'_py.zip',sep=''), files = paste('PS',ps,'_',scen,'_','E',experiment,'_R',runs,'.py',sep=''))
  }
}


make_sbatch_get_data <- function(sbatch_arguments,
                                 make_networks=F,
                                 repertoire_persistence=F,
                                 prepare_infomap=F,
                                 run_Infomap=F,
                                 read_infomap_results=F,
                                 temporal_diversity=F,
                                 module_Fst=F,
                                 run_experiments_file='/media/Data/PLOS_Biol/parameter_files/run_experiments') {
  
  if (detect_locale()=='Lab'){
    setwd('/home/shai/Documents/malaria_interventions')
    sqlite_path_global <- '/media/Data/PLOS_Biol/sqlite'
    parameter_files_path_global <- '/media/Data/PLOS_Biol/parameter_files'
  }
  if (detect_locale()=='Mac'){
    setwd('~/GitHub/malaria_interventions')
    sqlite_path_global <- '~/GitHub/PLOS_Biol/sqlite'
    parameter_files_path_global <- '~/GitHub/PLOS_Biol/parameter_files'
  }
  
  for (i in 1:nrow(sbatch_arguments)){
    ps <- sbatch_arguments$PS[i]
    scenario <- sbatch_arguments$scen[i]
    e <- sbatch_arguments$exp[i]
    modularity_exp <- sbatch_arguments$modularity_exp[i]
    cutoff_prob <- sbatch_arguments$cutoff_prob[i]
    write_edge_weights <- sbatch_arguments$write_edge_weights[i]
    if (detect_locale()=='Lab'){
      x <- readLines('~/Documents/malaria_interventions/PLOS_Biol/get_data_midway_plosbiol.sbatch')
    }
    if (detect_locale()=='Mac'){
      x <- readLines('~/GitHub/malaria_interventions/PLOS_Biol/get_data_midway_plosbiol.sbatch')
    }
    str_sub(x[2],20,22) <- paste(ps,scenario,e,cutoff_prob,sep='')
    str_sub(x[3],16,18) <- sbatch_arguments$time[i]
    str_sub(x[4],31,33) <- paste(ps,scenario,e,cutoff_prob,sep='')
    str_sub(x[5],30,32) <- paste(ps,scenario,e,cutoff_prob,sep='')
    str_sub(x[6],17,20) <- sbatch_arguments$array[i]
    str_sub(x[9],23,25) <- sbatch_arguments$mem_per_cpu[i]
    str_sub(x[19],5,7) <- ps
    str_sub(x[20],11,13) <- scenario
    str_sub(x[21],6,8) <- e
    str_sub(x[22],9,11) <- sbatch_arguments$layers[i]
    str_sub(x[23],13,16) <- cutoff_prob
    str_sub(x[25],17,19) <- modularity_exp
    x[26] <- paste("write_edge_weights=",write_edge_weights,sep='')
    x[length(x)+1] <- ""
    if (make_networks){x[length(x)+1] <- "Rscript $prog $PS $scenario $exp $cutoff_prob $layers 'make_networks' $modularity_exp $write_edge_weights"}
    if (repertoire_persistence){x[length(x)+1] <- "Rscript $prog $PS $scenario $exp $cutoff_prob $layers 'repertoire_persistence' $modularity_exp $write_edge_weights"}
    if (prepare_infomap){x[length(x)+1] <- "Rscript $prog $PS $scenario $exp $cutoff_prob $layers 'prepare_infomap' $modularity_exp $write_edge_weights"}
    if (run_Infomap){x[length(x)+1] <- "infomap_name='PS'$PS'_'$scenario'_E'$exp'_R'$SLURM_ARRAY_TASK_ID'_'$cutoff_prob'_Infomap_multilayer'"}
    if (run_Infomap){x[length(x)+1] <- "./Infomap_v01926 $infomap_name'.txt' . -i multilayer -d -N 10 --rawdir --two-level --tree --expanded"}
    if (read_infomap_results){x[length(x)+1] <- "Rscript $prog $PS $scenario $exp $cutoff_prob $layers 'read_infomap_results' $modularity_exp $write_edge_weights"}
    if (temporal_diversity){x[length(x)+1] <- "Rscript $prog $PS $scenario $exp $cutoff_prob $layers 'temporal_diversity' $modularity_exp $write_edge_weights"}
    if (module_Fst){x[length(x)+1] <- "Rscript $prog $PS $scenario $exp $cutoff_prob $layers 'module_Fst' $modularity_exp $write_edge_weights"}
    
    writeLines(x, paste(parameter_files_path_global,'/','PS',ps,'_',scenario,'_','E',e,'_',cutoff_prob,'_',modularity_exp,'_get_data_midway.sbatch',sep=''))
  }
  # Write a file to execute all the sbatch files
  chunks <- split(1:nrow(sbatch_arguments), ceiling(seq_along(1:nrow(sbatch_arguments))/500))
  for (ch in 1:length(chunks)){
    sink(paste(run_experiments_file,'_',ch,'.sh',sep=''))    
    for (i in chunks[[ch]]){
      ps <- sbatch_arguments$PS[i]
      scenario <- sbatch_arguments$scen[i]
      cutoff_prob <- sbatch_arguments$cutoff_prob[i]
      exp <- sbatch_arguments$exp[i]
      modularity_exp <- sbatch_arguments$modularity_exp[i]
      if ('after_job'%in%names(sbatch_arguments)){
        cat('sbatch -d afterok:',sbatch_arguments$after_job[i], ' PS',ps,'_',scenario,'_','E',exp,'_',cutoff_prob,'_',modularity_exp,'_get_data_midway.sbatch',sep='');cat('\n')
      } else {
        cat('sbatch PS',ps,'_',scenario,'_','E',exp,'_',cutoff_prob,'_',modularity_exp,'_get_data_midway.sbatch',sep='');cat('\n')
      }
    }
    sink.reset()
  }
}

# Plotting ----------------------------------------------------------------
library(ggplot2)
mytheme <- theme_bw() + theme(
  legend.title  = element_text(colour = "black", size=17),
  # legend.position = "none",
  #	legend.direction = "horizontal",
  legend.key = element_blank(),
  legend.text  = element_text(colour = "black", size=17),
  panel.background = element_blank(),
  # panel.grid.major = element_blank(),
  # panel.grid.minor = element_blank(),
  axis.text = element_text(color='black', family="Helvetica", size=14),
  axis.title = element_text(color='black', family="Helvetica", size=14),
  strip.text.x = element_text(family = "Helvetica", size = 14),
  strip.text.y = element_text(family = "Helvetica", size = 14),
  panel.border = element_rect(colour = "black", size=1.3),
  axis.ticks = element_line(size = 1.3),
  strip.background = element_rect( fill = "transparent", size = 1.3, colour = "black"  ),
  strip.text = element_text(size = 19)
)

mytheme_no_legend <- theme_bw() + theme(
  legend.title  = element_text(colour = "black", size=17),
  legend.position = "none",
  #	legend.direction = "horizontal",
  legend.key = element_blank(),
  legend.text  = element_text(colour = "black", size=17),
  panel.background = element_blank(),
  # panel.grid.major = element_blank(),
  # panel.grid.minor = element_blank(),
  axis.text = element_text(color='black', family="Helvetica", size=14),
  axis.title = element_text(color='black', family="Helvetica", size=14),
  strip.text.x = element_text(family = "Helvetica", size = 14),
  strip.text.y = element_text(family = "Helvetica", size = 14),
  panel.border = element_rect(colour = "black", size=1.3),
  axis.ticks = element_line(size = 1.3),
  strip.background = element_rect( fill = "transparent", size = 1.3, colour = "black"  ),
  strip.text = element_text(size = 19)
)

manuscript_theme <- theme_bw() + theme(
  legend.title  = element_text(colour = "black", size=17),
  legend.position = "none",
  #	legend.direction = "horizontal",
  legend.key = element_blank(),
  legend.text  = element_text(colour = "black", size=17),
  panel.background = element_blank(),
  # panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text = element_text(color='black', family="Helvetica", size=14),
  axis.title = element_text(color='black', family="Helvetica", size=14),
  strip.text.x = element_blank(),
  strip.text.y = element_blank(),
  panel.border = element_rect(colour = "black", size=1.3),
  axis.ticks = element_line(size = 1.3),
  strip.background = element_rect( fill = "transparent", size = 1.3, colour = "black"  ),
  strip.text = element_text(size = 19)
)


gg_color_hue <- function(n, hue_min = 10, hue_max = 280, l = 62, c = 100) {
  hues = seq(hue_min, hue_max, length=n+1)
  hcl(h=hues, l=l, c=c)[1:n]
}

# This function uses the list obtained by get_data() and generates relevant plots
generate_plots <- function(data, time_range=NULL){
  data1 <- data[[1]]
  data2 <- data[[2]]
  if (is.null(time_range)){
    time_range <- c(min(data1$time), max(data1$time))
  }
  plot_variables <- data1 %>% 
    select(-n_infected, -year, -month, -PS, -exp, -run, -scenario) %>% 
    filter(time>time_range[1]&time<time_range[2]) %>% gather(variable, value, -time) %>% 
    ggplot(aes(time, value, color=variable))+
    geom_line()+
    facet_wrap(~variable, scales = 'free')+
    mytheme+theme(legend.position = 'none')
  
  plot_eir <- data1 %>% 
    ggplot(aes(x=month,y=EIR))+
    geom_boxplot()+
    geom_point(stat='summary', fun.y=mean, color='red')+
    stat_summary(fun.y=mean, geom="line")+mytheme
  
  # Plot the age structure of infected hosts
  plot_age_structure <- data2 %>% group_by()%>% ggplot(aes(x=host_age))+geom_histogram() + 
    labs(x='Infected host age (months)') + 
    geom_vline(xintercept = 60) +
    mytheme
  return(list(plot_variables, plot_eir, plot_age_structure))
}


plotLayer <- function(network_object, l, edge_weight_multiply=1, remove.loops=T, ver.col=NULL, coords=NULL,...) { 
  g <- network_object$temporal_network[[l]]
  if (class(g)=='matrix'){g <- graph.adjacency(g, weighted = T, mode = 'directed')}
  if(remove.loops){g <- simplify(g, remove.multiple = F, remove.loops = T)}
  # g <- delete_edges(g, which(E(g)$weight<quantile(E(g)$weight, cutoff_g))) # remove all edges smaller than the cutoff
  # layout <-layout.kamada.kawai(g)
  # plot.new()
  # par(mar=c(0,3,3,3))
  if (is.null(ver.col)){ # All nodes have the same color?
    V(g)$color <- 'deepskyblue2'
  } else {
    if (class(ver.col)=='character') {V(g)$color <- ver.col}
    if (class(ver.col)=='data.frame') {
      V(g)$color <- ver.col$color[match(V(g)$name, ver.col$node_name)]
    }
  }
  if (!is.null(coords)){
    V(g)$repertoire <- splitText(V(g)$name,splitchar = '_', after = F)
    V(g)$x <- coords$x[match(V(g)$repertoire, coords$node_name)]
    V(g)$y <- coords$y[match(V(g)$repertoire, coords$node_name)]
  }
  plot(g, 
       vertex.color=V(g)$color,
       # vertex.label=V(g)$module,
       vertex.size=3,
       vertex.label=NA,
       # edge.arrow.mode='-', 
       edge.arrow.width=0.2,
       edge.arrow.size=0.2,
       edge.curved=0.5, 
       edge.width=E(g)$weight*edge_weight_multiply,
       # asp=0, # This needed if changing margins
       ...)
}

# Data extraction/analysis ------------------------------------------------


# This function obtains data from an sqlite file and prepares them for further analysis.
# requires sqldf
get_data <- function(parameter_space, scenario, experiment, run, cutoff_prob=0.9, sampling_period=30, host_age_structure=F, use_sqlite=T, tables_to_get=c('summary_general','sampled_infections')){
  # Initialize
  if (use_sqlite){
    base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
    
    if (detect_locale()=='Midway'){
      sqlite_file <- paste('/scratch/midway2/pilosofs/PLOS_Biol/sqlite/',base_name,'.sqlite',sep='')
    }
    if (detect_locale()=='Mac'){
      sqlite_file <- paste('~/GitHub/PLOS_Biol/sqlite/',base_name,'.sqlite',sep='')
    }
    if (detect_locale()=='Lab'){
      sqlite_file <- paste('/media/Data/PLOS_Biol/sqlite/',base_name,'.sqlite',sep='')
    }
    print(sqlite_file)
    print(file.exists(sqlite_file))
    if (!file.exists(sqlite_file)) {
      print (paste(sqlite_file, ' does not exist, ignoring and returning NULL'))
      return(NULL)
    }
    # parameter_file <- paste(base_name,'.py',sep='') # This may be necessary so I leave it
    
    # Extract data from sqlite. variable names correspond to table names
    print('Connecting to sqlite file...')
    db <- dbConnect(SQLite(), dbname = sqlite_file)
    summary_general <- dbGetQuery(db, 'SELECT * FROM summary')
    summary_general$PS <- parameter_space
    summary_general$exp <- experiment
    summary_general$scenario <- scenario
    summary_general$run <- run
    summary_general$year <- rep(1:(nrow(summary_general)/12), each=12)
    summary_general$month <- gl(n = 12, k = 1, length = nrow(summary_general),labels = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'), ordered = F)
    
    summary_alleles <- dbGetQuery(db, 'SELECT * FROM summary_alleles')
    summary_alleles %<>% group_by(time) %>% summarise(n_alleles=sum(n_circulating_alleles))
    
    summary_general <- suppressMessages(left_join(summary_general, summary_alleles))
    
    if ('sampled_infections'%in%tables_to_get){
      sampled_infections <- dbGetQuery(db, 'SELECT * FROM sampled_infections')
      sampled_infections$PS <- parameter_space
      sampled_infections$exp <- experiment
      sampled_infections$scenario <- scenario
      sampled_infections$run <- run
    }
    
    # Extract statistics
    if (nrow(summary_general)%%sampling_period==1){# The beginning of the data set in E00 has a timestap of 0 and this line is unnecesary. So remove it
      summary_general <- summary_general[-1,]
    }
    
    ## Prevalence
    summary_general$prevalence <- summary_general$n_infected/10^4
    
    ## Get EIR from the table in the sqlite
    summary_general$EIR <- summary_general$n_infected_bites/10000 # 10000 is the size of the human population
    
    # ##EIR can also be calculated manually as EIR=biting_rate * prevalence
    # biting_rate <- get_biting_rate(parameter_file)
    # summary_general$b <- NA
    # summary_general$b[] <- biting_rate # The [] is for recycling the biting_rate
    # summary_general$EIR1 <- summary_general$prevalence*summary_general$b*sampling_period
    
    # # Annual EIR is given by dividing by the sampling period to get a daily EIR, then multiplying by 360
    # summary_general %>% mutate(eir_y=EIR2/sampling_period*360) %>% group_by(year) %>% summarise(EIR_Year=mean(eir_y))
    
    # MOI#
    if ('sampled_infections'%in%tables_to_get){
      meanMOI <- sampled_infections %>% group_by(time, host_id) %>% dplyr::summarise(MOI=length(strain_id)) %>% group_by(time) %>% summarise(meanMOI=mean(MOI))
      summary_general <- suppressMessages(inner_join(summary_general, meanMOI)) # Add MOI to the summary
    }
    
    # Host age structure
    if(host_age_structure) {
      if (!'sampled_infections'%in%tables_to_get) {
        print('Cannot extract host age structre without sampled_infections. Make sure sampled_infections is included in the output. Returning NULL')
        return(NULL)
      }
      hosts <- dbGetQuery(db, 'SELECT * FROM hosts')
      names(hosts)[1] <- 'host_id'
      hosts$lifespan <- round((hosts$death_time-hosts$birth_time))
      hosts <- subset(hosts, host_id%in%sampled_infections$host_id)
      sampled_infections <- suppressMessages(left_join(sampled_infections, hosts, by='host_id'))
      sampled_infections$host_age <- round((sampled_infections$time-sampled_infections$birth_time)/30)    
    }
    dbDisconnect(db)
  }
  
  if (!use_sqlite){
    file <- paste('/media/Data/PLOS_Biol/Results/',parameter_space,'_',scenario,'/PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,'_',cutoff_prob,'_summary_general.csv',sep='')
    if (!file.exists(file)){
      print(paste(file, 'does not exist, ignoring and returning NULL'))
      return(NULL)
    }
    summary_general <- fread(file)
    summary_general$exp <- sprintf('%03d',summary_general$exp)
    summary_general$PS <- sprintf('%02d',summary_general$PS)
    summary_general$year <- as.factor(summary_general$year)
    summary_general$month <- factor(summary_general$month, levels = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))
    
    
    if ('sampled_infections'%in%tables_to_get){
      file <- paste(parameter_space,'_',scenario,'/PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,'_sampled_infections.csv',sep='')
      if (!file.exists(file)){
        print(paste(file, 'does not exist, ignoring and returning NULL'))
        return(NULL)
      }
      sampled_infections <- fread(file)
      sampled_infections$exp <- sprintf('%03d',sampled_infections$exp)
      sampled_infections$PS <- sprintf('%02d',sampled_infections$PS)
    }
  }
  
  
  if ('sampled_infections'%in%tables_to_get){
    return(list(summary_general=as.tibble(summary_general), sampled_infections=as.tibble(sampled_infections)))
  } else {
    return(list(summary_general=as.tibble(summary_general)))
  }
}





# This function compares an experiment to control and calculates statistics for the intervention
post_intervention_stats <- function(PS, scenario='S', exp, run, post_intervention_lag=360, control_data=NULL, plot.it=F, design_irs, use_sqlite=F){
  
  plots_out <- list()
  
  # If a data frame with control data is not included (NULL), then load the
  # control experiment
  if (is.null(control_data)){
    control_data <- get_data(parameter_space = PS, scenario = scenario, experiment = '001', run = run, use_sqlite = use_sqlite, tables_to_get = 'summary_general')[[1]]
  }
  
  # Calculate variable means (across time) for control
  control_means <- control_data %>% select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id) %>% 
    group_by(PS, exp, variable) %>% summarise(mean_value=mean(value, na.rm = T))
  
  # Plot control and means
  if (plot.it){
    plots_out$control <- control_data %>% select(-year, -month, -n_infected) %>%
      gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id) %>%
      ggplot(aes(x=time, y=value))+
      geom_line()+
      facet_wrap(~variable,scales='free')+
      geom_hline(aes(yintercept=mean_value), data=control_means, color='blue')+
      mytheme
  }
  
  # Load experiment data
  experiment_data <- get_data(parameter_space = PS, scenario = scenario, experiment = exp, run = run, use_sqlite = use_sqlite, tables_to_get = 'summary_general')[[1]]
  if (is.null(experiment_data)){
    return(NULL)
  }
  # Plot control vs experiment
  if (plot.it){
    plots_out$control_experiment <- bind_rows(control_data, experiment_data) %>%
      select(-year, -month, -n_infected) %>% 
      gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id) %>% 
      ggplot(aes(x=time, y=value, color=exp))+
      geom_vline(xintercept = c(29160,29160+1800), linetype='dashed', color='tomato')+
      geom_line()+
      scale_color_manual(values=c('#504E4E','purple'))+
      scale_x_continuous(breaks=pretty(x=subset(d, time>time_range[1]&time<time_range[2])$time,n=5))+
      geom_hline(aes(yintercept=mean_value), data=control_means, color='tomato')+
      facet_wrap(~variable, scales = 'free')+mytheme
  }
  
  # Calculate absolute difference from control at every timepoint
  x <- control_data %>%
    select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id)
  y <- experiment_data %>%
    select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id)
  diff_control <- suppressMessages(inner_join(x,y,by=c('pop_id','PS','scenario','run','time','variable'))) %>% 
    mutate(diff=value.y-value.x, ratio=value.y/value.x) %>% 
    select(variable, time, diff, ratio)
  
  if (plot.it){
    plots_out$control_diff <- diff_control %>%
      ggplot(aes(time, diff))+
      geom_line()+
      facet_wrap(~variable,scales='free')+
      geom_hline(aes(yintercept=mean_value), data=control_means, color='blue')+
      geom_vline(xintercept = c(29160,29160+3600), color='red')+
      mytheme
  }
  
  # Maximum amplitudes after an intervention
  amplitude <- experiment_data %>% left_join(subset(design_irs, select=c(IRS_START_TIMES,IRS_length,exp)), by='exp') %>% 
    mutate(intervention_lift=as.numeric(IRS_START_TIMES)+as.numeric(IRS_length)) %>% 
    filter(time>intervention_lift) %>% 
    select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id,-IRS_START_TIMES,-IRS_length) %>% 
    group_by(variable) %>% summarise(max_value=max(value))
  
  # New mean from some lag post intervention
  new_mean <- experiment_data %>% left_join(subset(design_irs, select=c(IRS_START_TIMES,IRS_length,exp)), by='exp') %>% 
    mutate(intervention_lift=as.numeric(IRS_START_TIMES)+as.numeric(IRS_length)) %>% 
    filter(time>intervention_lift+post_intervention_lag) %>% 
    select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id,-IRS_START_TIMES,-IRS_length) %>% 
    group_by(variable) %>% summarise(mean_value=mean(value, na.rm = T))
  
  # Time when extinction occured (last time point with data)
  time_extinct <- experiment_data %>% left_join(subset(design_irs, select=c(IRS_START_TIMES,IRS_length,exp)), by='exp') %>% 
    select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id,-IRS_START_TIMES,-IRS_length) %>% 
    group_by(variable) %>% summarise(time_extinct=max(time))
  
  summary_stats <- suppressMessages(left_join(time_extinct, new_mean))
  summary_stats <- suppressMessages(left_join(summary_stats, amplitude))
  summary_stats$PS <- PS
  summary_stats$scenario <- scenario
  summary_stats$exp <- exp
  summary_stats$run <- run  
  
  diff_control$PS <- PS
  diff_control$scenario <- scenario
  diff_control$exp <- exp
  diff_control$run <- run 
  
  if (plot.it){
    return(list(plots=plots_out,stats=list(diff_control=diff_control, summary_stats=summary_stats)))
  } else {
    return(list(diff_control=diff_control, summary_stats=summary_stats))  
  }
}


get_duration_infection <- function(parameter_space, scenario, experiment, run){
  if (on_Midway()){
    sqlite_file <- list.files(path = 'sqlite/', pattern=paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,'.sqlite',sep=''), full.names = T)
  } else {
    sqlite_file <- list.files(path = paste('/media/Data/malaria_interventions_data/sqlite_',scenario,'/',sep=''), pattern=paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,'.sqlite',sep=''), full.names = T)
  }
  db <- dbConnect(SQLite(), dbname = sqlite_file)
  sampled_duration <- dbGetQuery(db, 'SELECT * FROM sampled_duration')
  dbDisconnect(db)
  return(as.tibble(sampled_duration))
}


build_calendar <- function(num_years = 25, year_to_start=10){
  # This calendar helps match the days inthe ABM to the days in the empirical
  # data. This is especially important for simulations of seasonality, to match
  # the rain cycle and the IRS rounds.
  
  # Infomration on IRS:
  #
  # Survey 1: October 2012/End of Wet
  # Survey 2: June 2013/ End of Dry
  # IRS Round 1: Between October-December 2013 (End of the wet to Beginning of dry): 80% of Compounds sprayed with Organophosphates
  # Survey 3: June 2014/End of Dry
  # IRS Round 2: Between May-July 2014 (End of Dry/Beginning of Wet): 97% of Compounds sprayed with Organophosphates
  # Survey 4: October 2014/End of Wet
  # IRS Round 3: Between December 2014-Febuary 2015 (Middle of Dry): 96% of Compounds sprayed with Organophosphates Actellic 300CS
  # Survey 5: October 2015/End of Wet
  # Survey 6: June 2016/End of Dry
  
  # IRS interventions are supposed to be timed and completed before the start of
  # the wet season so that the insecticide is on the walls of the homes before the
  # mosquito population expands.An IRS round means that every home is sprayed once
  # during the time frame indicated. Insectisides are expected to last between
  # ~2-3 months (Organophosphates) or ~4-6 months (Organophosphates newer
  # formulation: Actellic 300CS) depending on the insecticide and/or formualtion
  # used.  This impacts the effect on reducing the mosquito population as females
  # rest on walls ~2-3 days post feeding, and therefore will be exposed to the
  # insecticide leading to death and therefore restricting/preventing transmission
  # of P.fal.
  # num_years: This is the number of years from start to end, not including pre-burnin
  months_in_year <- rep(c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'), each=30)
  calendar <- data.frame(running_day=seq(from = 1,to = 360*num_years,by=1),
                         year_sim=rep(1:num_years, each=360),
                         month_sim=rep(months_in_year,num_years),
                         day_sim=rep(1:30,num_years*12),
                         layer=rep(1:300, each=30),
                         stringsAsFactors = F)
  
  calendar$survey <- NA
  calendar$IRS <- NA
  
  # Set the dates corresponding to the empirical survey dates. year_to_start is
  # the year in the simulation where simulated data is starting to
  # be matched to the empirical data.
  calendar$survey[calendar$year_sim==year_to_start & calendar$month_sim=='Oct'] <- 'S1'
  calendar$survey[calendar$year_sim==year_to_start+1 & calendar$month_sim=='Jun'] <- 'S2'
  calendar$survey[calendar$year_sim==year_to_start+2 & calendar$month_sim=='Jun'] <- 'S3'
  calendar$survey[calendar$year_sim==year_to_start+2 & calendar$month_sim=='Oct'] <- 'S4'
  calendar$survey[calendar$year_sim==year_to_start+3 & calendar$month_sim=='Oct'] <- 'S5'
  calendar$survey[calendar$year_sim==year_to_start+4 & calendar$month_sim=='Jun'] <- 'S6'
  
  IRS_design <- data.frame(IRS=paste('IRS',1:3,sep='_'),
                           start_day=c(1,1,1),
                           start_month=c('Nov','Jun','Jan'),
                           start_year=c(year_to_start+1,year_to_start+2,year_to_start+3)
  )
  IRS_design %<>% rowwise() %>% 
    mutate(running_day_start=extract_from_calendar(calendar, run_day=NULL, d=start_day, m=start_month, y=start_year)$running_day)
  IRS_design$running_day_end <- IRS_design$running_day_start+c(75,75,150)-1
  
  calendar[calendar$running_day>=IRS_design$running_day_start[1] & calendar$running_day<=IRS_design$running_day_end[1],'IRS'] <- 'IRS1'
  calendar[calendar$running_day>=IRS_design$running_day_start[2] & calendar$running_day<=IRS_design$running_day_end[2],'IRS'] <- 'IRS2'
  calendar[calendar$running_day>=IRS_design$running_day_start[3] & calendar$running_day<=IRS_design$running_day_end[3],'IRS'] <- 'IRS3'
  
  return(calendar)
}

extract_from_calendar <- function(cal, run_day=NULL, y=NULL, m=NULL, d=NULL){
  if(!is.null(run_day)){
    return(subset(cal, running_day==run_day))
  } else {
    if (is.null(d)){
      x <- subset(cal, year_sim==y & month_sim==m)
      return(x)
    } else {
      x <- subset(cal, year_sim==y & month_sim==m & day_sim==d)
      return(x)
    }  
  }
}


# Generate networks -------------------------------------------------------

# Calculate edge values in networks of repertories
overlapAlleleAdj<-function(mat){
  newmat<-tcrossprod(mat>0)
  newmat<-newmat/rowSums(mat>0)
  return(newmat)  
}


# Function to calculate survival probability of repertoires between layers. This
# is used to rescale the interlayer edges.
rescale_by_divergence <- function(ps,scenario,exp,run,layers_to_include,interlayer_edges){
  
  ## rescale function for similarity-----------
  ## simMat: vector of inter-layer edges similarities that needs to be rescaled
  ## Nt0: baseline of interval poptime we are rescaling to,
  ## for example, if baseline population size is 1000, and interval is 1 month
  ## then Nt0 = 1000*1 = 1000
  #   -- That would be the N_e (Harmonic mean) of the population before intervention 
  ## popsize: a vector of population size per month between two layers
  ## t: time intervals between two layers, if it's 4 months, then t is 4
  divRescale<-function(simMat, Nt0, popsize, t){
    Ne<-1/mean(1/popsize)
    rescale_factor<-Ne*t/Nt0
    return(1-(1-simMat)/rescale_factor)
  }
  
  # Get infection data
  base_name <- paste('PS',ps,'_',scenario,'_E',exp,'_R',run,sep='')
  if (on_Midway()){
    sqlite_file <- paste('/scratch/midway2/pilosofs/PLOS_Biol/sqlite/',base_name,'.sqlite',sep='')
  } else {
    sqlite_file <- paste('/media/Data/PLOS_Biol/sqlite/',base_name,'.sqlite',sep='')
  }
  
  db <- dbConnect(SQLite(), dbname = sqlite_file)
  print('Getting genetic data from sqlite...')
  summary_table <- as.tibble(dbGetQuery(db, 'SELECT time, n_infections FROM summary'))
  dbDisconnect(db)
  summary_table$layer <- group_indices(summary_table, time)
  
  # Make Nt0 the minimum Ne*t of the 5 layers
  Nt0 <- 10^9 # Just a very large number
  for (l in 1:(length(layers_to_include)-1)){
    popsize <- subset(summary_table, layer>=layers_to_include[l] & layer <= layers_to_include[l+1])$n_infections
    t <- layers_to_include[l+1]-layers_to_include[l]
    Ne <- 1/mean(1/popsize)
    print(Ne*t)
    Nt0 <- ifelse(Ne*t<Nt0,Ne*t,Nt0)
  }
  
  # Make Nto the minimum number of infections
  # Nt0 <- min(summary_table$n_infections)

  # Calculate survival probability
  w_rescaled <- c()
  for (l in 1:(length(layers_to_include)-1)){
    popsize <- subset(summary_table, layer>=layers_to_include[l] & layer <= layers_to_include[l+1])$n_infections
    t <- layers_to_include[l+1]-layers_to_include[l]
    simMat <- subset(interlayer_edges, layer_s==l)$w
    w_rescaled <- c(w_rescaled, divRescale(simMat, Nt0, popsize, t))
  }
  interlayer_edges$w_rescaled <- w_rescaled
  return(interlayer_edges)
}


# Function to calculate survival probability of repertoires between layers. This
# is used to rescale the interlayer edges.
calculate_rep_survival <- function(ps,
                                   scenario,
                                   exp,
                                   run,
                                   cutoff_prob,
                                   layers_to_include){
  
  Gn<-function(t, z=0, r){ # # Prob of having a given frequency at time t
    if (t<=1) {
      return(exp(r[t]*(z-1)));
    }else{
      return(exp(r[t]*(Gn(t-1,z,r)-1)))
    }
  }
  surprob<-function(popsizes) {
    r<-c()
    for (t in 2:length(popsizes)){
      r<-c(r, popsizes[t]/popsizes[t-1])
    }
    return(1-Gn(length(r),0,r))
  }
  
  # Get infection data
  base_name <- paste('PS',ps,'_',scenario,'_E',exp,'_R',run,sep='')
  if (on_Midway()){
    sqlite_file <- paste('/scratch/midway2/pilosofs/PLOS_Biol/sqlite/',base_name,'.sqlite',sep='')
  } else {
    sqlite_file <- paste('/media/Data/PLOS_Biol/sqlite/',base_name,'.sqlite',sep='')
  }
  
  db <- dbConnect(SQLite(), dbname = sqlite_file)
  print('Getting genetic data from sqlite...')
  summary_table <- as.tibble(dbGetQuery(db, 'SELECT time, n_infections FROM summary'))
  dbDisconnect(db)
  summary_table$layer <- group_indices(summary_table, time)
  
  # Calculate survival probability
  print('Calculating survival prob...')
  surv_prob <- c()
  for (l in 1:(length(layers_to_include)-1)){
    x <- subset(summary_table, layer>=layers_to_include[l] & layer <= layers_to_include[l+1])
    surv_prob <- c(surv_prob, surprob(x$n_infections))
  }
  return(surv_prob)
}

# A function to build the similarity matrix for a single layer and calculate some summary stats
build_layer <- function(infection_df, unit_for_edges, write_to_files=F, base_filename=NULL){
  infection_df %<>% group_by(strain_id) %>%
    mutate(freq = n()/120) %>% # strain frequency (the number of strain copies should be equal to the frequency)
    arrange(strain_id_unique) 
  # Calculate the edges
  if (unit_for_edges=='alleles'){
    bipartite_layer <- table(infection_df$strain_id_unique, infection_df$allele_locus)
  }
  if (unit_for_edges=='genes'){
    bipartite_layer <- table(infection_df$strain_id_unique, infection_df$gene_id)/2 # divide by two beause each hene appears twice because there are two loci  
  }
  similarity_matrix <- overlapAlleleAdj(bipartite_layer)
  # Some summary
  layer_summary <- with(infection_df, 
                        data.frame(hosts=length(unique(host_id)),
                                   repertoires_unique=length(unique(strain_id)),
                                   repertoires_total=length(unique(strain_id_unique)),
                                   num_units=ncol(as.data.frame.matrix(bipartite_layer)),
                                   unit_for_edges=unit_for_edges
                        ))
  if (write_to_files){
    write.csv(as.data.frame.matrix(bipartite_layer),paste(base_filename,'_bipartite_layer.csv',sep=''), row.names = T)
    write.csv(as.data.frame.matrix(similarity_matrix),paste(base_filename,'_similarity_layer.csv',sep=''), row.names = T)
  }
  return(list(similarity_matrix=similarity_matrix, infections=infection_df, layer_summary=layer_summary))
}


createTemporalNetwork <- function(ps, 
                                  scenario, 
                                  exp, 
                                  run, 
                                  cutoff_prob, 
                                  cutoff_value=NULL,
                                  layers_to_include=NULL, 
                                  sampled_infections=NULL, 
                                  unit_for_edges='alleles',
                                  repertoires_to_sample=NULL, 
                                  write_to_files=F){
  # Define the sqlite file to use
  base_name <- paste('PS',ps,'_',scenario,'_E',exp,'_R',run,sep='')
  if (on_Midway()){
    sqlite_file <- paste('/scratch/midway2/pilosofs/PLOS_Biol/sqlite/',base_name,'.sqlite',sep='')
  } else {
    sqlite_file <- paste('/media/Data/PLOS_Biol/sqlite/',base_name,'.sqlite',sep='')
  }
  
  # Extract data from sqlite. variable names correspond to table names
  db <- dbConnect(SQLite(), dbname = sqlite_file)
  print('Getting genetic data from sqlite...')
  sampled_strains <- as.tibble(dbGetQuery(db, 'SELECT id, gene_id FROM sampled_strains'))
  names(sampled_strains)[1] <- c('strain_id')
  sampled_alleles <- as.tibble(dbGetQuery(db, 'SELECT * FROM sampled_alleles'))
  names(sampled_alleles)[3] <- c('allele_id')
  sampled_strains <- full_join(sampled_strains, sampled_alleles)
  sampled_strains$allele_locus <- paste(sampled_strains$allele_id,sampled_strains$locus,sep='_') # each allele in a locus is unique
  dbDisconnect(db)
  ## Get infection data
  if (is.null(sampled_infections)){
    print('Getting infection data from sqlite...')
    sampled_infections <- get_data(parameter_space = ps, experiment = exp, scenario = scenario, run = run)$sampled_infections
  }
  
  # Build the data set
  print('Building data set...')
  ## Add layer numbers. Each time point is a layer
  sampled_infections$layer <- group_indices(sampled_infections, time) 
  sampled_infections %<>% arrange(layer,strain_id,host_id) %>% 
    group_by(layer,strain_id) %>% 
    mutate(strain_copy = row_number()) # add unique id for each strain copy within each layer. A strain copy is an instance of a strain in a particualr host
  sampled_infections$strain_id_unique <- paste(sampled_infections$strain_id,sampled_infections$strain_copy,sep='_') # Create the unique strains
  ## Integrate the strain composition into the infections table
  if (all(unique(sampled_strains$strain_id)%in%unique(sampled_infections$strain_id))==F || all(unique(sampled_infections$strain_id)%in%unique(sampled_strains$strain_id))==F) {
    stop('There may be a mismatch in repertoires between the sampled_strains and sampled_infections data sets. Revise!')
  }
  sampled_infections %<>% select(time, layer, host_id, strain_id, strain_copy, strain_id_unique) %>% 
    left_join(sampled_strains)
  
  # Build layers
  print('Building layers...')
  Layers <- list()
  if (is.null(layers_to_include)){
    layers_to_include <- sort(unique(sampled_infections$layer))
  }
  for (l in layers_to_include){
    cat(paste('[',Sys.time(), '] building layer ',l,'\n',sep=''))
    sampled_infections_layer <- subset(sampled_infections, layer==l) # This is MUCH faster than sampled_infections_layer <- sampled_infections %>% filter(layer==l)
    # Sub-sample repertoires (for empirical data)
    if (!is.null(repertoires_to_sample)){
      # If there are not enough repertoires to sample (because the intervention was too strong)
      if(length(unique(sampled_infections_layer$strain_id))<repertoires_to_sample[which(layers_to_include==l)]){
        sampled_repertoires <- unique(sampled_infections_layer$strain_id)
      } else {
        sampled_repertoires <- sample(unique(sampled_infections_layer$strain_id),repertoires_to_sample[which(layers_to_include==l)],F)
      }
      sampled_infections_layer <- subset(sampled_infections_layer, strain_id%in%sampled_repertoires)
    }
    # This line makes the layer
    Layers[[which(layers_to_include==l)]] <- build_layer(infection_df = sampled_infections_layer,
                                                         unit_for_edges = unit_for_edges,
                                                         write_to_files = write_to_files,
                                                         base_filename = paste(base_name,cutoff_prob,str_pad(l,3,'left','0'),sep='_'))
  }
  intralayer_matrices <- lapply(Layers, function(x) x$similarity_matrix)   # Get just the matrices
  layer_summary <- do.call(rbind, lapply(Layers, function(x) x$layer_summary)) # Get the layer summary
  layer_summary$layer <- layers_to_include
  
  # Build "interlayer networks"
  print('Building inter-layer networks...')
  interlayer_matrices <- list()
  for (current_layer in layers_to_include[-length(layers_to_include)]){
    next_layer <- layers_to_include[which(layers_to_include==current_layer)+1]
    current_layer_idx <- which(layers_to_include==current_layer)
    next_layer_idx <- which(layers_to_include==current_layer)+1
    
    strain_copies_t <- rownames(intralayer_matrices[[current_layer_idx]]) # repertoires at time t
    strain_copies_t1 <- rownames(intralayer_matrices[[next_layer_idx]]) # repertoires at time t+1
    # need minimum of 2 strains in t and t+1 to build a matrix
    if (length(strain_copies_t)<2 | length(strain_copies_t1)<2){
      print(paste('No interlayer edges between layers ',current_layer,' and ',next_layer,' because there are < 2 repertoires in one of them.',sep=''))
      return(NULL)
    } 
    sampled_infections_interlayer <- subset(sampled_infections, layer%in%c(current_layer,next_layer)) # This is MUCH faster than sampled_infections_layer <- sampled_infections %>% filter(layer==l)
    x <- build_layer(sampled_infections_interlayer,unit_for_edges)$similarity_matrix # This is the similarity matrix for all the repertoires in both layers.
    # Pull only the similarity values between the repertoires from the correct layers (i.e. create a bipartite)
    inter_layer_edges_matrix <- x[strain_copies_t,strain_copies_t1]
    interlayer_matrices[[current_layer_idx]] <- inter_layer_edges_matrix
    print(paste('Built interlayer edges for layers: ',current_layer,'-->',next_layer,sep=''))
  }
  
  # Create cutoff. Note that zeroes SHOULD BE included in the distribution of edge weights
  # Cutoff is based on all the intralayer and interlayer edges.
  if (is.null(cutoff_value)){
    # Raw values of intralayer edges, in case needed to plot edge weight distributions; and only produced when cutoff is not already defined
    intralayer_edges_no_cutoff <- tibble(layer=rep(1:length(intralayer_matrices),times=(sapply(intralayer_matrices, function(x) nrow(x)*ncol(x)))),
                                         w=unlist(sapply(intralayer_matrices, as.vector)))
    # Raw values of interlayer edges, in case needed to plot edge weight distributions; and only produced when cutoff is not already defined           
    interlayer_edges_no_cutoff <- tibble(layer=rep(1:length(interlayer_matrices),times=(sapply(interlayer_matrices, function(x) nrow(x)*ncol(x)))),
                                         w=unlist(sapply(interlayer_matrices, as.vector)))
    edges <- c(intralayer_edges_no_cutoff$w,interlayer_edges_no_cutoff$w)
    cutoff_value <- quantile(edges, probs = cutoff_prob)
    #as.tibble(edges) %>% ggplot(aes(value))+geom_density()+geom_vline(xintercept = cutoff_value)
  }
  
  # Apply cutoff to layers
  layer_summary$density <- layer_summary$density_interlayer <- NA
  print('Applying cutoff to layers...')
  for (i in 1:length(intralayer_matrices)){
    # print(i)
    x <- intralayer_matrices[[i]]
    x[x<cutoff_value] <- 0
    intralayer_matrices[[i]] <- x
    layer_summary$density[i] <- sum(x>0)/(nrow(x)*ncol(x))
  }
  # Apply cutoff to interlayer blocks
  print('Applying cutoff to interlayer blocks...')
  for (i in 1:length(interlayer_matrices)){
    # print(i)
    x <- interlayer_matrices[[i]]
    x[x<cutoff_value] <- 0
    interlayer_matrices[[i]] <- x
    layer_summary$density_interlayer[i] <- sum(x>0)/(nrow(x)*ncol(x))
  }
  
  layer_summary <- as.tibble(layer_summary)
  
  print('Done!')
  
  return(list(intralayer_matrices=intralayer_matrices, 
              interlayer_matrices=interlayer_matrices,
              intralayer_edges_no_cutoff=intralayer_edges_no_cutoff, # Raw values of intralayer edges, incase needed to plot edge weight distributions
              interlayer_edges_no_cutoff=interlayer_edges_no_cutoff, # Raw values of interlayer edges, incase needed to plot edge weight distributions
              cutoff_prob=cutoff_prob, 
              cutoff_value=cutoff_value, 
              layer_summary=layer_summary, 
              base_name = base_name, ps=ps, scenario=scenario, experiment=exp, run=run))
}


# Get network data --------------------------------------------------------


get_edge_disributions <- function(PS, scenario, exp, run, cutoff_prob, get_inter=T){
  x <- readLines(paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_network_info.csv',sep=''))
  cutoff_value <- as.numeric(x[6])
  intra <- fread(paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_intralayer_no_cutoff.csv',sep=''))  
  if (get_inter){
    inter <- fread(paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_interlayer_no_cutoff.csv',sep=''))  
    x <- rbind(intra,inter)
  } else {
    x <- intra
  }
  x$PS=PS
  x$scenario=scenario
  x$exp=exp
  x$run=run
  x$cutoff_value=cutoff_value
  x$cutoff_prob=cutoff_prob
  # p <- as.tibble(x) %>% ggplot(aes(value))+geom_density()+geom_vline(xintercept = cutoff_value, color='red')+
  #   labs(title = paste('Cut-off quantile:',cutoff_prob,'Cut-off value:',cutoff_value))
  return(x)
}


get_results_for_cutoff <- function(cutoff_prob_seq=seq(0.25,0.95,0.05), scenario='S', run_range=1:10){
  design_cutoff <- expand.grid(PS=sprintf('%0.2d', 4:6),
                               scenario=scenario, 
                               exp='001',
                               run_range=run_range,
                               cutoff_prob=cutoff_prob_seq,
                               stringsAsFactors = F)
  for (run in unique(design_cutoff$run_range)){
    results_cutoff <- c()
    design_cutoff_run <- subset(design_cutoff, run_range==run)
    print(design_cutoff_run)
    for (i in 1:nrow(design_cutoff_run)){
      PS <- design_cutoff_run[i,1]
      scenario <- design_cutoff_run[i,2]
      exp <- design_cutoff_run[i,3]
      # run <- design_cutoff_run[i,4]
      cutoff_prob <- design_cutoff_run[i,5]
      file <- paste('Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_modules.csv',sep='')
      if(file.exists(file)){
        print(paste(PS,scenario,exp,run,cutoff_prob,sep=' | '))
        x <- read_csv(file, col_types = 'iiccciccccd')  
        results_cutoff <- rbind(results_cutoff, x)
      } else {
        print(paste('File does not exist: ',file,sep=''))
      }
    }
    write_csv(results_cutoff, paste('Results/results_cutoff_',scenario,'_R',run,'.csv',sep=''))
  }
  
  results_cutoff <- c()
  for (run in unique(design_cutoff$run_range)){
    file <- paste('Results/results_cutoff_',scenario,'_R',run,'.csv',sep='')
    if(file.exists(file)){
      x <- read_csv(file, col_types = 'iicccicccid')
      results_cutoff <- rbind(results_cutoff, x)  
    } else {
      print(paste('File does not exist: ',file,sep=''))
    }
  }
  results_cutoff <- as.tibble(results_cutoff)
  return(results_cutoff)
}




# This function builds the network object from the result files produced on Midway
get_network_structure <- function(ps, scenario, exp, run, cutoff_prob, layers_to_include, parse_interlayer=T, plotit=F, folder='/media/Data/'){
  require(utils)
  basename <- paste('PS',ps,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,sep='')
  network_structure <- list()
  
  
  print('Getting network info...')
  network_structure$layer_summary <- suppressMessages(read_csv(paste(folder,'/Results/',ps,'_',scenario,'/',basename,'_layer_summary.csv',sep='')))
  info <- fread(paste(folder,'/Results/',ps,'_',scenario,'/',basename,'_network_info.csv',sep=''))
  network_structure$cutoff_prob <- as.numeric(info[5,1])
  network_structure$cutoff_value <- as.numeric(info[6,1])
  network_structure$base_name <- basename
  network_structure$source <- 'CSV'
  
  
  print('Parsing nodes and edges...')
  network_structure$temporal_network <- list()
  intralayer <- fread(paste(folder,'/Results/',ps,'_',scenario,'/',basename,'_intralayer.csv',sep=''))
  nodelist <- fread(paste(folder,'/Results/',ps,'_',scenario,'/',basename,'_node_list.csv',sep=''))
  print('Getting layers...')
  pb <- txtProgressBar(min = 1, max=max(layers_to_include))
  for (l in layers_to_include){
    setTxtProgressBar(pb, l)
    tmp <- subset(intralayer, layer_s==l)
    if (nrow(tmp)==0){
      network_structure$temporal_network[[l]] <- matrix(0,0,0)
    } else {
      names(tmp)[5] <- 'weight'
      tmp_nodelist <- subset(nodelist, nodeID%in%union(tmp$node_s,tmp$node_t))
      g <- graph.data.frame(tmp[,c(2,4,5)], vertices = tmp_nodelist)
      V(g)$name <- V(g)$nodeLabel
      network_structure$temporal_network[[l]] <- g
    }
  }
  close(pb)
  
  
  if(parse_interlayer){
    print('Parsing interlayer edges...')
    network_structure$interlayer <- list()
    interlayer <- fread(paste(folder,'/malaria_interventions_data/Results/',ps,'_',scenario,'/',basename,'_interlayer.csv',sep=''))
    print('Getting interlayer edge structures as bipartite graphs...')
    pb <- txtProgressBar(min = 1, max=max(layers_to_include)-1)
    for (l in head(layers_to_include,-1)){
      setTxtProgressBar(pb, l)
      tmp <- subset(interlayer, layer_s==l)
      if (nrow(tmp)==0){
        network_structure$interlayer[[l]] <- matrix(0,0,0)
      } else {
        names(tmp)[5] <- 'weight'
        tmp_nodelist <- subset(nodelist, nodeID%in%union(tmp$node_s,tmp$node_t))
        # print('--> Applying a cutoff to the interlayer edges...')
        tmp <- subset(tmp, weight>=network_structure$cutoff_value)
        # print('--> Adjusting state nodes...')
        # Need to define state nodes
        tmp$from <- paste(tmp$layer_s,tmp$node_s,sep = '_')
        tmp$to <- paste(tmp$layer_t,tmp$node_t,sep = '_')
        # Define state nodes in the node list data frame
        tmp_nodelist_from <- tmp_nodelist_to <- tmp_nodelist
        tmp_nodelist_from$nodeID <- paste(l,tmp_nodelist_from$nodeID,sep='_')
        tmp_nodelist_from <- subset(tmp_nodelist_from, nodeID%in%tmp$from)
        tmp_nodelist_from$layer <- l
        tmp_nodelist_to$nodeID <- paste(l+1,tmp_nodelist_to$nodeID,sep='_')
        tmp_nodelist_to <- subset(tmp_nodelist_to, nodeID%in%tmp$to)
        tmp_nodelist_to$layer <- l+1
        tmp_nodelist <- rbind(tmp_nodelist_from,tmp_nodelist_to)
        
        # Create the graph
        # print('--> parsing interlayer edges as a bipartite graph...')
        g <- graph.data.frame(tmp[,c(6,7,5)], vertices = tmp_nodelist, directed = T)
        tmp <- as.data.frame(tmp)
        V(g)$type <- V(g)$name %in% tmp[,6] # make it bipartite
        network_structure$interlayer[[l]] <- g
      }
    }
    close(pb)
  }
  
  if (plotit){
    x <- as.tibble(network_structure$layer_summary)
    print(x %>% gather(variable, value, -layer) %>% 
            ggplot(aes(layer, value, color=variable))+
            geom_line()+
            scale_x_continuous(breaks = seq(1,max(x$layer)+5,5))+
            mytheme+
            theme(panel.grid.minor = element_blank())
    )
  } 
  print('Done')
  return(network_structure)
}

# This function gets a network object and calculates network properties poer layer
get_network_properties <- function(network_object, layers_to_include, num_properties=25){
  network_properties <- matrix(NA,ncol=num_properties, nrow=length(layers_to_include))
  rownames(network_properties) <- layers_to_include
  print('Calculating network properties for layers...')
  pb <- txtProgressBar(min = 1, max=max(layers_to_include))
  for (l in layers_to_include){
    setTxtProgressBar(pb, l)
    # print(l)
    if (class(network_object$temporal_network[[l]])=='matrix'){
      if (nrow(network_object$temporal_network[[l]])==0){
        print(paste('Layer ',l,' does not exist, putting NA',sep=''))
        network_properties[as.character(l),] <- NA
        next
      }
    }
    tmp <- calculateFeatures(network_object, l)
    network_properties[as.character(l),] <- tmp
  }
  close(pb)
  colnames(network_properties) <- names(tmp)
  network_properties <- as.tibble(network_properties)
  network_properties$layer <- layers_to_include
  return(network_properties)
}

get_interlayer_properties <- function(network_object, layers_to_include, num_properties=7){
  interlayer_properties <- matrix(NA,ncol=num_properties, nrow=length(layers_to_include)-1)
  rownames(interlayer_properties) <- head(layers_to_include, -1)
  print('Calculating network properties for interlayer edges...')
  pb <- txtProgressBar(min = 1, max=max(layers_to_include)-1)
  for (l in head(layers_to_include, -1)){
    setTxtProgressBar(pb, l)
    g <- network_object$interlayer[[l]]
    # print(l)
    if (class(g)=='matrix'){
      if (nrow(g)==0){
        print(paste('Layer ',l,' does not exist, putting NA',sep=''))
        interlayer_properties[as.character(l),] <- NA
        next
      }
    }
    close(pb)
    interlayer_properties[as.character(l),1] <- sum(V(g)$type==T)
    interlayer_properties[as.character(l),2] <- sum(V(g)$type==F)
    interlayer_properties[as.character(l),3] <- graph.density(g)
    interlayer_properties[as.character(l),4] <- ecount(g)
    interlayer_properties[as.character(l),5] <- mean(E(g)$weight)
    interlayer_properties[as.character(l),6] <- mean(strength(g, vids = V(g)[V(g)$type==T]))
    interlayer_properties[as.character(l),7] <- mean(strength(g, vids = V(g)[V(g)$type==F]))
  }
  colnames(interlayer_properties) <- c('num_nodes_s','num_nodes_t','density_il','num_edges_il','mean_edge_weight_il','mean_strength_s','mean_strength_t')
  interlayer_properties <- as.tibble(interlayer_properties)
  interlayer_properties$layer <- head(layers_to_include, -1)
  return(interlayer_properties)
}


# This is a wrapper function
analyze_network <- function(ps, scenario, exp, run, layers_to_include){
  network_object <- get_network_structure(ps,scenario,exp,run, layers_to_include)
  network_properties <- get_network_properties(network_object, layers_to_include)
  network_properties$exp <- exp
  interlayer_properties <- get_interlayer_properties(network_object, layers_to_include)
  interlayer_properties$exp <- exp
  return(list(network_object=network_object,network_properties=network_properties,interlayer_properties=interlayer_properties))
}

# This is a wrapper function
analyze_networks_multiple <- function(ps, scenario, experiments=c('001','002','003','004'), runs, layers_to_include, parse_interlayer=T, folder='/media/Data/'){
  results <- c()
  for (exp in experiments){
    for (run in runs){
      network_object <- get_network_structure(ps,scenario,exp,run, layers_to_include, folder = folder, parse_interlayer = parse_interlayer)
      
      network_properties <- get_network_properties(network_object, layers_to_include)
      network_properties$PS <- ps
      network_properties$scenario <- scenario
      network_properties$exp <- exp
      network_properties$run <- run
      
      if(parse_interlayer){
        interlayer_properties <- get_interlayer_properties(network_object, layers_to_include)
        interlayer_properties$PS <- ps
        interlayer_properties$scenario <- scenario
        interlayer_properties$exp <- exp
        interlayer_properties$run <- run
        x <- full_join(network_properties, interlayer_properties, by=c('PS','scenario','exp','run','layer'))
      } else {
        x <- network_properties
      }
      results <- rbind(results, x)
     }
  }
  return(results)
}

# Network properties ------------------------------------------------------
density_bipartite <- function(x){
  return(
    sum(x>0)/(nrow(x)*ncol(x))
  )  
}

f_01_averageLocalClusteringCoeff <- function(g,GC=F){
  if (GC) {
    gc <- giant.component(g)
    return(mean(transitivity(gc, type = 'local'), na.rm = T))
  } else {
    return(mean(transitivity(g, type = 'local'), na.rm = T))
  }
}

f_02_averageLocalClusteringCoeffWeighted <- function(g,GC=F){
  if (GC) {
    gc <- giant.component(g)
    return(mean(transitivity(gc, type = 'barrat'), na.rm = T))
  } else {
    return(mean(transitivity(g, type = 'barrat'), na.rm = T))
  }
}


f_03_globalClusteringCoeff <- function(g, GC=F){
  if (GC) {
    gc <- giant.component(g)
    return(transitivity(gc, type = 'global'))
  } else {
    return(transitivity(g, type = 'global'))
  }
}

f_04_gdensity <- function(g){
  return(ecount(g)/(vcount(g)*(vcount(g)-1)))
}

f_05_proportionSingletons <- function(g){
  sum(degree(g)==0)/length(V(g))
}

f_06_proportionEndpoints <- function(g){
  sum(degree(g)==1)/length(V(g))
}

f_07_meanDegree <- function(g){
  mean(degree(g))
}

f_08_meanDegreeNotSingletons <- function(g){
  mean(degree(g)[degree(g)!=0])
}

f_09_meanStrength <- function(g){
  mean(strength(g))
}

f_10_meanStrengthNotSingletons <- function(g){
  mean(strength(g)[degree(g)!=0])
}

f_11_entropyDegreeDistribution <-  function(g,verbose=F){
  y=degree(g)
  freq=prop.table(table(y))
  if (verbose){print(freq)}
  -sum(freq * log(freq, base = 2))
}

f_12_ratioComponents <- function(g) { 
  cl <- clusters(g) 
  cl$no/length(V(g))
}

giant.component <- function(g) { 
  cl <- clusters(g) 
  induced.subgraph(g, which(cl$membership == which.max(cl$csize)))
}

f_13_averageComponentSize <- function(g) { 
  cl <- clusters(g) 
  mean(cl$csize)
}

f_14_entropyComponentSize <-  function(g,verbose=F){
  cl <- clusters(g)
  y=cl$csize
  freq=prop.table(table(y))
  if (verbose){print(freq)}
  -sum(freq * log(freq, base = 2))
}

f_15_giantConnectedRatio <- function(g){
  cl <- clusters(g)
  max(cl$csize)/length(V(g))
}

f_16_meanEccentricity <- function(g){
  mean(eccentricity(g))
}

f_17_gdiameter <- function(g){
  return(diameter(g))
}

f_18_meanDiameterComponents <- function(g){
  cl <- clusters(g)
  d <- 0
  for (m in 1:(cl$no)){
    d <- d + diameter(induced.subgraph(g, which(cl$membership == m)))
  }
  return(d/cl$no)
}

f_19_globalEfficiency <- function(g){
  d_ij <- shortest.paths(g)
  d_ij <- d_ij[lower.tri(d_ij)]
  d_ij <- d_ij[!is.infinite(d_ij)]
  N=length(V(g))
  1/(N*(N-1))*sum(1/d_ij)
}

f_20_averageClosenessCentrality <- function(g){
  mean(closeness(g, weights = NULL))
}

f_21_averageEdgeWeight <- function(g){
  mean(E(g)$weight)
}

graphDensityDirected <- function(m){ # density of directed graphs
  E = sum(m!=0)-length(diag(m)) # in our networks the diagonal is 1 so this should be considered.
  N=nrow(m)
  return(E/(N*(N-1)))
}




motifsProportion <- function(g){ # Calculate the proportion of each of the 16 motifs out of the total motifs found
  
  # This code shows the motifs by order
  # par(mfrow=c(4,4), mar=c(.75,.75,.75,.75))
  # for (i in 0:15){ # note that counting starts at 0
  #   g <- graph_from_isomorphism_class(size = 3, number = i)
  #   V(g)$name <- c('B','A','C')
  #   plot(g,
  #        edge.arrow.size = 0.4,
  #        edge.color='black',
  #        main = i + 1)
  #   box(col='red')
  # }
  
  motifs <- graph.motifs(g, size = 3)
  motifs.prop <- motifs/sum(motifs, na.rm = T)
  names(motifs.prop) <- c('M01--A_B_C', #1
                          'M02--A->B_C', #2
                          'M03--A->B<-C', #3
                          'M04--A<->B_C', #4
                          'M05--A->B->C', #5
                          'M06--A<->B<-C',#6
                          'M07--A<-B->C', #7
                          'M08--A->B<-C_A->C', #8
                          'M09--A<-B->C_A<->C', #9
                          'M10--A<->B->C', #10
                          'M11--A<->B<->C', #11
                          'M12--A<-B<-C_A->C', #12
                          'M13--A->B->C_A<->C', #13
                          'M14--A->B<-C_A<->C', #14
                          'M15--A->B<->C_A<->C', #15
                          'M16--A<->B<->C_A<->C') #16
  return(motifs.prop)
}


calculateFeatures <- function(network_object, l, remove.loops=F){
  g <- network_object$temporal_network[[l]]
  if (class(g)!='matrix' & class(g)!='igraph'){
    stop('Network is not a matrix neither an igraph object')
  }
  if (class(g)=='matrix'){
    g <- graph.adjacency(g, weighted = T, mode = 'directed')
  }
  if(remove.loops){g <- simplify(g, remove.multiple = F, remove.loops = T)}
  
  featureVector <- NULL
  featureVector <- c(featureVector, vcount(g))
  featureVector <- c(featureVector, ecount(g))
  featureVector <- c(featureVector, f_03_globalClusteringCoeff(g))
  featureVector <- c(featureVector, f_04_gdensity(g))
  featureVector <- c(featureVector, f_07_meanDegree(g))
  featureVector <- c(featureVector, f_12_ratioComponents(g))
  featureVector <- c(featureVector, f_15_giantConnectedRatio(g))
  featureVector <- c(featureVector, f_17_gdiameter(g))
  featureVector <- c(featureVector, f_21_averageEdgeWeight(g))
  featureVector <- c(featureVector, motifsProportion(g))
  names(featureVector)[1:9] <- c('Num_nodes','Num_edges','GCC','density','mean_degree','ratio_comp','nodes_in_giant_comp','diameter','mean_edge_weight')
  
  # featureVector <- vector(length=32)
  # Diagnostics of transitivity
  # featureVector[1] <- f_01_averageLocalClusteringCoeff(g,F)    # Clustering coefficient averaged across all nodes
  # featureVector[2] <- f_02_averageLocalClusteringCoeffWeighted(g,F)    # Barrat's clustering coefficient averaged across all nodes
  # featureVector[3] <- f_03_globalClusteringCoeff(g,F)
  # Diagnostics of degree/sterngth
  # featureVector[4] <- f_04_gdensity(g)                    # Graph density
  # featureVector[5] <- f_05_proportionSingletons(g)             # Proportion of nodes with degree 0 of all the nodes
  # featureVector[6] <- f_06_proportionEndpoints(g)              # Proportion of nodes with degree 1 of all the nodes
  # featureVector[7] <- f_07_meanDegree(g)                       # Average degree 
  # featureVector[8] <- f_08_meanDegreeNotSingletons(g)
  # featureVector[9] <- f_09_meanStrength(g)                     # Average strength
  # featureVector[10] <- f_10_meanStrengthNotSingletons(g)
  # featureVector[11] <- f_11_entropyDegreeDistribution(g)       # Average measurement of the heterogeneity of the network. See eq. 14 in: da F. Costa, L., Rodrigues, F. A., Travieso, G. & Boas, P. R. V. Characterization of complex networks: A survey of measurements. Adv. Phys. 56, 167–242 (2007).
  # featureVector[12] <- f_12_ratioComponents(g)                 # Number of components relative to networks size
  # featureVector[13] <- f_13_averageComponentSize(g)
  # featureVector[14] <- f_14_entropyComponentSize(g)
  # featureVector[15] <- f_15_giantConnectedRatio(g)             # Proportion of nodes in the giant component
  # #Diagnostics of shortest-paths
  # featureVector[16] <- f_16_meanEccentricity(g)                 # Eccentricity is the maximum shortest distance from a node to all other nodes. This is averaged across all nodes
  # featureVector[17] <- f_17_gdiameter(g)                        # length of the longest geodesic
  # featureVector[18] <- f_18_meanDiameterComponents(g)
  # featureVector[19] <- f_19_globalEfficiency(g)                # See eq. 14 in: da F. Costa, L., Rodrigues, F. A., Travieso, G. & Boas, P. R. V. Characterization of complex networks: A survey of measurements. Adv. Phys. 56, 167–242 (2007).
  # featureVector[20] <- f_20_averageClosenessCentrality(g)      # 
  # # Diagnostics of motifs
  # featureVector[21:32] <- motifsProportion(g)
  return(featureVector)
}


# Infomap -----------------------------------------------------------------


##  A function that gets the layer as a matrix and writes it for infomap as an edge list
# network_object is a list of matrices, each element in the list is a layer.
# requires igraph
matrix_to_infomap_intralayer <- function(l, nodeList, network_object, remove_self_links=T){
  current_layer <- network_object$intralayer_matrices[[l]]
  if(nrow(current_layer)<2){
    print(paste('Less than 2 repertoires in layer',l,'!!! skipping layer.'))
    return(NULL)
  }
  if (all(current_layer==0)){
    print(paste('All edges in layer',l,' are 0!!! skipping layer.'))
    return(NULL)
  }
  g <- graph.adjacency(current_layer, mode = 'directed', weighted = T)
  current_layer_el <- as.tibble(as_data_frame(g, what = 'edges'))
  names(current_layer_el) <- c('node_s','node_t','w')
  current_layer_el$layer_s <- l
  current_layer_el$layer_t <- l
  current_layer_el$node_s <- nodeList$nodeID[match(current_layer_el$node_s,nodeList$nodeLabel)]
  current_layer_el$node_t <- nodeList$nodeID[match(current_layer_el$node_t,nodeList$nodeLabel)]
  current_layer_el %<>% select(layer_s, node_s, layer_t, node_t, w) # Re-arrange columns for the infomap input order
  if(remove_self_links){
    current_layer_el %<>% filter(node_s != node_t)
  }
  print(paste('[',Sys.time(), '] Created edge list of layer ',l,' for Infomap | ', nrow(current_layer_el),' edges',sep=''))
  return(current_layer_el)
}

##  A function that gets the layer as a matrix and writes it for infomap as an edge list
# network_object is a list of matrices, each element in the list is a layer.
# requires igraph
matrix_to_infomap_interlayer <- function(l, nodeList, network_object){
  interlayer_block <- network_object$interlayer_matrices[[l]]
  if (all(interlayer_block==0)){
    print(paste('There are no interlayer edges in block',l,'-->',l+1,'!!! skipping layer.'))
    return(NULL)
  }
  g <- graph.incidence(incidence = interlayer_block, mode = 'out', weighted = T)
  current_layer_el <- as.tibble(as_data_frame(g, what = 'edges'))
  names(current_layer_el) <- c('node_s','node_t','w')
  current_layer_el$layer_s <- l
  current_layer_el$layer_t <- l+1
  current_layer_el$node_s <- nodeList$nodeID[match(current_layer_el$node_s,nodeList$nodeLabel)]
  current_layer_el$node_t <- nodeList$nodeID[match(current_layer_el$node_t,nodeList$nodeLabel)]
  current_layer_el %<>% select(layer_s, node_s, layer_t, node_t, w) # Re-arrange columns for the infomap input order
  print(paste('[',Sys.time(), '] Created interlayer edge list for layers ',l,'-->',l+1,' for Infomap | ', nrow(current_layer_el),' edges',sep=''))
  return(current_layer_el)
}


# This function takes a list of temporal matrices and returns an intralayer and
# interlayer edge lists in a format: [layer_source, node_source, layer_target node_target, weight].
# It also returns the list of node names.
# requires igraph
build_infomap_objects <- function(network_object, 
                                  write_to_infomap_file=T, 
                                  infomap_file_name, 
                                  return_objects=T, 
                                  repertoire_survival_prob=NULL,
                                  rescale_by_survival_prob=F){
  require(data.table)
  intralayer_matrices <- network_object$intralayer_matrices
  base_name <- network_object$base_name

  # Get the node list
  nodeLabel <- sort(unique(unlist(lapply(intralayer_matrices,rownames))))
  nodeList <- data.frame(nodeID=1:length(nodeLabel), nodeLabel)
  
  # Create intralayer edge lists.
  print('Creating intralayer edge lists')
  layers <- 1:length(intralayer_matrices)
  infomap_intralayer <- lapply(layers, function (x) matrix_to_infomap_intralayer(x, nodeList = nodeList, network_object = network_object))
  print(head(infomap_intralayer[[1]]))
  print('Creating a DF of intralayer edges')
  infomap_intralayer <- do.call("rbind", infomap_intralayer)
  print(head(infomap_intralayer))
  
  # Create interlayer edge lists.
  print('Creating interlayer edge lists')
  layers <- layers[-length(layers)]
  infomap_interlayer <- lapply(layers, function (x) matrix_to_infomap_interlayer(x, nodeList = nodeList, network_object = network_object))
  
  if (rescale_by_survival_prob){
    print('Rescaling interlayer edges by survival probabilities...')
    # Rescale interlayer edges using survival probability
    surv_prob <- calculate_rep_survival(network_object$ps, 
                                        network_object$scenario, 
                                        network_object$experiment, 
                                        network_object$run, 
                                        network_object$cutoff_prob,
                                        network_object$layer_summary$layer)
    rescaling_factors <- tibble(layer_s=1:length(surv_prob),surv_prob=surv_prob)
    print(rescaling_factors)
    for (l in 1:length(surv_prob)){
      # Divide edge weights by the probability of persistence. those that persisted
      # for longer will have stronger values
      x <- infomap_interlayer[[l]]
      suppressMessages(x %<>% left_join(rescaling_factors) %>% mutate(w_rescaled=w/surv_prob))
      infomap_interlayer[[l]] <- x
    }
  }
  
  # Rescale edges using probabilities recorded in a file, calculated from the neutral scenario
  if (!is.null(repertoire_survival_prob)){
    print('Rescaling interlayer edges...')
    # 8 months between S1 and S2
    # 12 months between S2 and S3
    # 4 months between S3 and S4
    # 12 months between S4 and S5
    # 8 months between S5 and S6
    rescaling_factors <- tibble(layer_s=1:length(repertoire_survival_prob),surv_prob=repertoire_survival_prob)
    print(rescaling_factors)
    for (l in 1:length(repertoire_survival_prob)){
      # Divide edge weights by the probability of persistence. those that persisted
      # for longer will have stronger values
      x <- infomap_interlayer[[l]]
      suppressMessages(x %<>% left_join(rescaling_factors) %>% mutate(w_rescaled=w/surv_prob))
      infomap_interlayer[[l]] <- x
    }
  }

  
  print('Creating a DF of interlayer edges')
  infomap_interlayer <- do.call("rbind", infomap_interlayer)
  print(head(infomap_interlayer))
  

  
  if (write_to_infomap_file){
    ## Write file for infomap
    print('Writing Infomap files')
    print(paste('Infomap file:',infomap_file_name))
    if (file.exists(infomap_file_name)){unlink(infomap_file_name)}
    
    if (rescale_by_survival_prob | !is.null(repertoire_survival_prob)){
      x <- infomap_interlayer %>% select(layer_s,node_s,layer_t,node_t,w_rescaled) %>% rename(w=w_rescaled)
      edges_to_write <- infomap_intralayer %>% bind_rows(x)
    } else {
      edges_to_write <- infomap_intralayer %>% bind_rows(infomap_interlayer)
    }
    fwrite(edges_to_write, infomap_file_name, sep=' ', col.names = F)
  }
  
  if (return_objects){
    return(list(infomap_intralayer=infomap_intralayer, infomap_interlayer=infomap_interlayer, nodeList=nodeList))
  }
}



infomap_readTreeFile <- function(PS, scenario, exp, run, cutoff_prob, folder='/media/Data/PLOS_Biol/Results/cutoff_to_use/', infomap_file=NULL){
  print ('Reading infomap file...')
  if (is.null(infomap_file)){
    if (on_Midway()){
      infomap_file <-       paste('PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_Infomap_multilayer_expanded.tree',sep='')
      node_list <- read_csv(paste('PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_node_list.csv',sep=''), col_types = list(col_character(),col_character()))
    } else {
      infomap_file <-       paste(folder,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_Infomap_multilayer_expanded.tree',sep='')
      node_list <- read_csv(paste(folder,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_node_list.csv',sep=''), col_types = list(col_character(),col_character()))
    }
  }
  lines <- readLines(infomap_file)
  cat(lines[1]);cat('\n')
  # x <- fread(infomap_file, skip = 2, stringsAsFactors = F) # Read results of Infomap
  x <- read_delim(infomap_file, 
                  delim = ' ',
                  col_types = list(col_character(), col_double(), col_character(), col_integer(), col_integer()), 
                  col_names = c('path', 'flow', 'name', 'layer', 'node'),
                  skip = 2) # Read results of Infomap
  print(x)
  
  # Create a data frame to store results
  modules <- tibble(module=rep(NA,nrow(x)),
                    nodeID=rep(NA,nrow(x)),
                    layer=rep(NA,nrow(x)),
                    path=x$path)
  
  print('Creating module data frame...')
  modules$module <- as.numeric(str_split(string = modules$path, pattern = ':', simplify = T)[,1])
  modules$nodeID <- str_trim(str_split(string = x$name, pattern = '\\|', simplify = T)[,1])
  modules$layer <- as.numeric(str_trim(str_split(string = x$name, pattern = '\\|', simplify = T)[,2])) # can also use x$layer
  
  
  # modules %>% filter(module==1) %>%  ggplot(aes(x=layer,y=module))+geom_point()

  # Rename modules because Infomap gives random names
  print('Adding information on strains...')
  modules2 <- modules %>% 
    distinct(module,layer) %>% 
    arrange(module,layer)
  x <- c(1,table(modules2$module))
  module_birth_layers <- modules2 %>% slice(cumsum(x)) %>% arrange(layer,module)
  module_renaming <- data.frame(module=module_birth_layers$module, module_renamed = 1:max(module_birth_layers$module)) 
  modules2 %<>% left_join(module_renaming)
  modules2 %<>% full_join(modules) 
  modules2 %<>% select(-module, -path)
  names(modules2)[2] <- 'module'
  modules2 %<>% arrange(module, layer, nodeID)
  
  # modules2 %>% ggplot(aes(x=layer,y=module))+geom_point()
  
  # Change node IDs to repertoire names
  modules2$nodeID <- as.integer(modules2$nodeID)
  node_list$nodeID <- as.integer(node_list$nodeID)
  print(paste('Same strains in module and the node_list strain data frames?',setequal(modules2$nodeID,node_list$nodeID)))
  modules2 %<>% left_join(node_list) %>% 
    rename(strain_unique=nodeLabel) %>% 
    mutate(strain_id=str_split(strain_unique,'_', simplify = T)[,1])
  
  # Add information
  modules2$PS <- PS
  modules2$scenario <- scenario
  modules2$exp <- exp
  modules2$run <- run
  modules2$cutoff_prob <- cutoff_prob
  
  print('Getting information on genes and alleles...')
  # Get the gene list for the repertoires
  if (on_Midway()){
    db <- dbConnect(SQLite(), dbname = paste('/scratch/midway2/pilosofs/PLOS_Biol/sqlite/','PS',PS,'_',scenario,'_E',exp,'_R',run,'.sqlite',sep=''))
  } else {
    db <- dbConnect(SQLite(), dbname = paste('/media/Data/PLOS_Biol/sqlite/','PS',PS,'_',scenario,'_E',exp,'_R',run,'.sqlite',sep=''))
  }
  
  sampled_strains <- as.tibble(dbGetQuery(db, 'SELECT id,gene_id FROM sampled_strains'))
  names(sampled_strains)[1] <- 'strain_id'
  sampled_alleles <- as.tibble(dbGetQuery(db, 'SELECT * FROM sampled_alleles'))
  dbDisconnect(db)
  
  print('Done!')
  return(list(modules=modules2,sampled_strains=sampled_strains,sampled_alleles=sampled_alleles))
}



# Analysis of modules------------------------------------------------------
calculate_module_diversity <- function(PS, scenario, exp, run, cutoff_prob){
  if (on_Midway()){
    file_modules <- paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_modules.csv',sep='')
    file_strains <- paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_sampled_strains.csv',sep='')
  } else {
    file_modules <- paste('Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_modules.csv',sep='')
    file_strains <- paste('Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_sampled_strains.csv',sep='')  
  }
  # file_modules <- paste('PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_modules.csv',sep='')
  # file_strains <- paste('PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_sampled_strains.csv',sep='')  
  
  if(file.exists(file_modules) & file.exists(file_strains)){
    print(paste(PS,scenario,exp,run,cutoff_prob,sep=' | '))
    modules <- read_csv(file_modules, col_types = 'iiccciccccd')
    
    max_layer_to_persist <- max(modules$layer)-min(modules$layer)+1 # This is important for analyses that do not have sequential number of layers (1:300), like in the empirical IRS.
      
    module_persistence <- modules %>% 
      select(scenario, PS, run, cutoff_prob, layer, module) %>% 
      group_by(scenario, PS,run,cutoff_prob,module) %>% 
      summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
      mutate(relative_persistence=persistence/(max_layer_to_persist-birth_layer+1))
    
    sampled_strains <- read_csv(file_strains, col_types = 'ccc')
    sampled_strains <-  sampled_strains[,-3]
    suppressMessages(modules %<>% select(scenario, PS, scenario, exp, run, cutoff_prob, module, strain_id) %>% left_join(sampled_strains))
    allele_freq <- xtabs(~module+allele_locus, modules)
    module_diversity <- vegan::diversity(allele_freq)
    module_diversity_normalized <- vegan::diversity(allele_freq)/log(ncol(allele_freq))
    
    module_persistence$D <- module_diversity
    module_persistence$D_normalized <- module_diversity_normalized
    module_persistence$statistic <- module_diversity_normalized/module_persistence$relative_persistence
    return(module_persistence)
  } else {
    print(paste('One file does not exist:',file_modules,file_strains))
  }
}

pairFST<-function(totalSt, totalAllele, pop1 = x, pop2 = y){
  
  subLoc<-totalSt%>%filter(module %in% c(pop1, pop2))
  subLoc<-t(subLoc[,-1])
  subAllele<-totalAllele%>%filter(module %in% c(pop1, pop2))
  
  subAllele<-t(subAllele[,-1])
  
  
  #get total number of copy of genes per locus, a vector of dimension, number of locus
  indPLocus<-rowSums(subLoc)
  #get number of populations per locus, a vector of dimension, number of locus
  popPLocus<-rep(2, length(indPLocus))
  nc<- (indPLocus - apply(subLoc^2, 1, sum, na.rm = TRUE)/indPLocus)/(popPLocus - 1)
  #expand it to each allele, replicate the values
  indPAllele <- rep(indPLocus, 2)
  ncPAllele<-rep(nc, 2)
  indPAllelePop<-rbind(subLoc,subLoc)
  #calculate allele frequency per pop
  mhom<-subAllele
  p<-mhom/indPAllelePop
  #calculate allele frequency across pop
  pb<-rowSums(mhom)/rowSums(indPAllelePop)
  
  SSG<-rep(0,length(indPLocus))
  
  dum<-indPAllelePop * (p - 2 * p^2) + mhom
  SSi <- rowSums(dum, na.rm = TRUE) #per allele, sum across location level
  dum1 <- indPAllelePop * (sweep(p, 1, pb))^2
  SSP <- 2 * rowSums(dum1, na.rm = TRUE)
  #how many populations per allele
  popPAllele <- rep(2, length(indPAllele))
  
  MSG <- SSG/indPAllele
  MSP <- SSP/(popPAllele - 1)#population effect
  MSI <- SSi/(indPAllele - popPAllele)#individual effect
  sigw <- MSG
  sigb <- 0.5 * (MSI - MSG)
  siga <- 1/2/ncPAllele * (MSP - MSI)#across each allele
  Fxy <- function(x) x[1]/sum(x, na.rm = TRUE)
  
  tsiga <- sum(siga, na.rm = TRUE)
  tsigb <- sum(sigb, na.rm = TRUE)
  tsigw <- sum(sigw, na.rm = TRUE)
  tFST <- Fxy(c(tsiga, tsigb, tsigw))
  return(tFST)
}

pairFSTMat<-function(dat,maxModule = 100){
  print("starting FST calculation")
  dat<-arrange(dat, module)
  totalSt<- dat %>%group_by(module)%>%
    summarise_all(length)
  al1<-dat %>%group_by(module)%>%
    summarise_all(sum)
  al0<-totalSt-al1
  totalAllele<-cbind(al1,al0[,-1])
  colnames(totalAllele)[-1]<-paste('al', c(1:(ncol(totalAllele)-1)),sep="")
  
  moduleList<-unique(dat$module)
  moduleNumber<-length(moduleList)
  print(paste("total number of modules are", moduleNumber))
  if (moduleNumber<2){
    print("no sub communities detected")
    return(matrix(NA,comNumber,comNumber))
  }
  moduleSize<-as.vector(table(dat$module))
  moduleList<-moduleList[moduleSize>=5]
  moduleSize<-moduleSize[moduleSize>=5]
  print(paste("total number of modules having larger than 5 strains are", length(moduleList)))
  if (length(moduleSize)>maxModule){
    moduleList <- sort(sample(moduleList,maxModule))
    print(moduleList)
  }
  moduleNumber<-length(moduleList)
  
  outMat<-matrix(NA,moduleNumber,moduleNumber)
  for (i in 1:(moduleNumber-1)){
    for (j in (i+1):moduleNumber){
      outMat[i,j]<-outMat[j,i]<-pairFST(totalSt, totalAllele, pop1 = moduleList[i], pop2 = moduleList[j])
    }
  }
  return(outMat)
}

calculate_mFst <- function(PS, scenario, exp, run, cutoff_prob, maxModule=100){
  if (on_Midway()){
    file_modules <- paste('PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_modules.csv',sep='')
    file_strains <- paste('PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_sampled_strains.csv',sep='')
  } else {
    file_modules <- paste('Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_modules.csv',sep='')
    file_strains <- paste('Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_sampled_strains.csv',sep='')  
  }
  if(file.exists(file_modules) & file.exists(file_strains)){
    print(paste(PS,scenario,exp,run,cutoff_prob,sep=' | '))
    modules <- read_csv(file_modules, col_types = 'iiccciccccd')
    
    # modules_in_last_layer <- modules %>% filter(layer>290) %>% select(module)
    # modules_in_last_layer <- unique(modules_in_last_layer$module)
    # 
    module_list <- modules %>% 
      # filter(module%in%modules_in_last_layer) %>% 
      select(module, strain_id)
    
    sampled_strains <- read_csv(file_strains, col_types = 'ccc')
    sampled_strains <-  sampled_strains[,-3]
    sampled_strains %<>%
      filter(strain_id %in% module_list$strain_id)
    
    # print(setequal(sampled_strains$strain_id, module_list$strain_id))
    
    allele_freq <- xtabs(~strain_id+allele_locus, sampled_strains)
    
    # print(setequal(rownames(allele_freq), module_list$strain_id))
    
    # write.csv(allele_freq, paste('~/Dropbox/Qixin_Shai_Malaria/PLOS_Biol/mFst/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_mFst_strains_alleles_matrix.csv',sep=''))
    # write.csv(sampled_strains, paste('~/Dropbox/Qixin_Shai_Malaria/PLOS_Biol/mFst/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_mFst_strains_alleles_list.csv',sep=''))
    # write.csv(module_list, paste('~/Dropbox/Qixin_Shai_Malaria/PLOS_Biol/mFst/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_mFst_strains_modules.csv',sep=''), row.names=F)
    
    allele_freq_mat <- matrix(allele_freq, nrow=nrow(allele_freq), ncol=ncol(allele_freq), dimnames = list(rownames(allele_freq), colnames(allele_freq)))
    allele_freq_df <- as.data.frame(allele_freq_mat)
    allele_freq_df$strain_id <- rownames(allele_freq_mat)
    
    # Calculate the mFst
    module_list = suppressMessages(left_join(module_list,allele_freq_df))
    
    
    FstMat<-pairFSTMat(module_list[,-2], maxModule = maxModule)
    
    
    # Results
    x <- FstMat[upper.tri(FstMat)]
    result <- tibble(PS, scenario, exp, run, cutoff_prob, mFst=x)
    return(result)
  } else {
    print(paste('One file does not exist:',file_modules,file_strains,'--returning NULL.'))
    return(NULL)
  }
}


get_modularity_results <- function(PS,scenario,exp,run,cutoff_prob,folder='/media/Data/PLOS_Biol/Results/cutoff_to_use/'){
  file <- paste(folder,'PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_modules.csv',sep='')
  if(file.exists(file)){
    print(paste(PS,scenario,exp,run,cutoff_prob,sep=' | '))
    # x <- read_csv(file, col_types = 'iicccicccid')  
    x <- fread(file, colClasses=c('integer','integer','character','character','character','integer','character','character','character','integer','double'))
    return(x)
  } else {
    print(paste('File does not exist: ',file,sep=''))
    return()
  }
}

get_temporal_diversity <- function(PS,scenario,exp,run,cutoff_prob,folder='/media/Data/PLOS_Biol/Results/cutoff_to_use/'){
  file <- paste(folder,'PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_temporal_diversity.csv',sep='')
  if(file.exists(file)){
    print(paste(PS,scenario,exp,run,cutoff_prob,sep=' | '))
    # x <- fread(file, colClasses=c('character','character','integer','double','integer','integer','integer','integer','double','double','double'))
    x <- fread(file)
    return(x)
  } else {
    print(paste('File does not exist: ',file,sep=''))
    return()
  }
}

get_mFst <- function(PS,scenario,exp,run,cutoff_prob,folder='/media/Data/PLOS_Biol/Results/cutoff_to_use/'){
  file <- paste(folder,'PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_mFst.csv',sep='')
  if(file.exists(file)){
    print(paste(PS,scenario,exp,run,cutoff_prob,sep=' | '))
    x <- fread(file, colClasses=c('character','character','character','integer','double','double'))
    return(x)
  } else {
    print(paste('File does not exist: ',file,sep=''))
    return()
  }
}

