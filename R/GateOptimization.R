#COMPLEX - USER-EXPOSED ---------------------
#' Run Gate Optimization
#'
#' Executes a gated model-based optimization starting from a pre-selected treatment parameter set.
#' The continuous-point objective value is used as a gate reference to identify optimized parameter
#' sets with penalty values that are better than or equal to the starting-point penalty.
#'
#' @import dplyr
#' @import reshape2
#' @importFrom pso psoptim
#' @importFrom optimx optimx
#' @importFrom GA ga
#' @importFrom optimization optim_sa
#' @importFrom DEoptim DEoptim
#' @importFrom ABCoptim abc_optim
#'
#' @param opt_object An optimization task object as created by \code{initiateOptimizationTask}.
#' Must include a valid population sample, search space definition, and PKPD model.
#' @param starting_point Named numeric vector defining the initial treatment parameter values
#' used to initialize and gate the optimization search.
#' @param optimization_method Character. Name of the optimization algorithm to use.
#' Supported methods are the same as for \code{runOptimization}, with \code{"PSO"}
#' recommended for gated optimization.
#' @param optimization_control Named list of control parameters specific to the selected
#' optimization method.
#' @param return_optimization Logical. If TRUE, the optimization result is returned together
#' with the starting point. If FALSE, the updated \code{opt_object} is returned. Default is FALSE.
#' @param verbose_level Integer (0–3). Controls verbosity during optimization.
#' Higher values provide more detailed console output.
#' @param output_filename Character or NULL. Optional filename to which optimization history
#' or results may be written.
#' @param save_optimization Logical. If TRUE, the optimization results and updated object
#' are saved to disk. Default is TRUE.
#' @param multiple Logical. If TRUE, suppresses console messages for batch execution of
#' multiple gate optimizations. Default is FALSE.
#'
#' @return If \code{return_optimization = FALSE}, the updated \code{opt_object} is returned.
#' If TRUE, a named vector combining the starting point and optimization result is returned.
#'
#' @seealso \code{\link{runOptimization}}
#' @export
runGateOptimization <- function(opt_object, starting_point, 
                                optimization_method = 'PSO',
                                optimization_control=list(), return_optimization = F,
                                verbose_level = 2, output_filename = NULL, save_optimization=T,
                                multiple = F){
  # get continuous point penalty value
  pointPenalty_continuous <- objFn(as.numeric(starting_point), names(starting_point),
                                   opt_object[['simulationSettings']],
                                   model=opt_object[['pkpdModel']], 
                                   evaluation_matrix=opt_object[['treatmentObjectives']],
                                   pkpd_pars=opt_object[['modelParameters']], 
                                   fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']],
                                   population_samples=opt_object[['populationSample']], 
                                   n_sub=opt_object[['simulationSettings']]['n_subjects'],
                                   fun_postprocessing=opt_object[['fun_postProcessing']],
                                   verbose_level=0, write_file=NULL,
                                   dose_scaling = opt_object[['fun_doseScaling']],
                                   eval_format = 'continuous')
  
  # adjust search space starting point
  opt_object[['searchSpace_gated']] <- opt_object[['searchSpace']]
  opt_object[['searchSpace_gated']]$start <- sapply(c(1:length(starting_point)), function(i){
    starting_point[names(starting_point)==rownames(opt_object[['searchSpace_gated']])[i]]
  })
  
  # replace value with 0.5 of 1 failed attainment if initial penalty was exactly zero
  if(any(pointPenalty_continuous==0)){
    pointPenalty_continuous[pointPenalty_continuous==0] <- 0.5*(1/opt_object[['simulationSettings']]['n_subjects'])
  }
  
  # run optimization
  if(!multiple){print("Starting gate optimization . . . ")}
  opt_result <- coreOptimization_gate(opt_object, pointPenalty_continuous,
                                      method_name = optimization_method, 
                                      kurtosis = 50, 
                                      verbose = verbose_level, opt_control=optimization_control, 
                                      save_object = save_optimization,
                                      output_file = output_filename,
                                      return_optimization_result=return_optimization)
  if(!multiple){print("Gate optimization complete!")}
  
  # return value
  if(!return_optimization){
    return(opt_object)
  }else{
    names(starting_point) <- paste0("startingPoint_", names(starting_point))
    return_value <- c(starting_point, opt_result)
    return(return_value)
  }
}

#' Run Gate Optimization for Multiple Starting Points
#'
#' Executes gate optimization for a set of candidate starting points (separate runs).
#' Each row of dataframe \code{start_points} is treated as a distinct starting-point vector and is passed
#' to \code{runGateOptimization}. Runs can be executed sequentially or in parallel (via
#' \code{parallel::mclapply}) depending on the number of available cores.
#'
#' The function creates an output directory (if needed), writes the starting points to file,
#' and assigns a dedicated output filename to each optimization run.
#'
#' @importFrom parallel mclapply
#'
#' @param opt_object An optimization task object as created by \code{initiateOptimizationTask}.
#' Must include a valid population sample, search space definition, and PKPD model.
#' @param start_points Data frame of starting points, where each row corresponds to a named
#' treatment parameter vector used as \code{starting_point} input to \code{runGateOptimization}.
#' Columns are matched and reordered to \code{rownames(opt_object[['searchSpace']])}.
#' @param n_cores Integer. Number of CPU cores to use for parallel execution. If greater than 1,
#' runs are distributed using \code{parallel::mclapply}. Default is 1.
#' @param optimization_method Character. Name of the optimization algorithm to use for each run
#' (passed to \code{runGateOptimization}). Default is \code{"PSO-LBFGSB"}.
#' @param output_directory_name Character or NULL. Name of the subdirectory within
#' \code{opt_object[['WorkDirectory']]} where outputs will be written. If provided, the directory
#' is created if it does not exist.
#' @param optimization_control Named list of control parameters for the selected optimization
#' method (passed to \code{runGateOptimization}).
#'
#' @return The (unchanged) \code{opt_object}. Results are produced as side effects (files written
#' to disk); each run returns its optimization output internally but is not propagated by this
#' wrapper function.
#'
#' @seealso \code{\link{runGateOptimization}}
#' @export
runGateOptimization__Multiple <- function(opt_object, start_points, n_cores=1,
                                          optimization_method='PSO-LBFGSB',
                                          output_directory_name=NULL, optimization_control=list()){
  # order starting point names
  start_points <- start_points[,rownames(opt_object[['searchSpace']])]
  
  # prepare output directory and filenames
  output_dir <- paste0(opt_object[['WorkDirectory']], output_directory_name, "/")
  if(!file.exists(output_dir)){dir.create(output_dir)}
  
  # print out starting points
  rownames(start_points) <- c(1:nrow(start_points))
  write.csv(start_points, paste0(output_dir, "GateOptimization_StartingPoints.csv"), row.names=T)
  
  # prepare output filenames
  output_filenames <- paste0(output_dir, "GateOptimizationRun_", c(1:nrow(start_points)))
  
  # run optimization (parallel)
  print("Starting gate optimization (multiple) . . . ")
  n_cores <- min(n_cores, nrow(start_points))
  if(n_cores>1){
    gate_optimization_results <- mclapply(c(1:nrow(start_points)), function(xi){
      runGateOptimization(opt_object, start_points[xi,], 
                          optimization_method = optimization_method,
                          optimization_control=optimization_control,
                          output_filename = output_filenames[xi], 
                          return_optimization = T, save_optimization = F, multiple=T)}, mc.cores=n_cores)
  }else{
    gate_optimization_results <- lapply(c(1:nrow(start_points)), function(xi){
      runGateOptimization(opt_object, start_points[xi,], 
                          optimization_method = optimization_method,
                          optimization_control=optimization_control,
                          output_filename = output_filenames[xi], 
                          return_optimization = T, save_optimization = F, multiple=T)})
  }
  print("Gate optimization complete (multiple)!")
  
  # return object
  return(opt_object)
}

#NON-USER ---------------------
# Optimization - units
activationFunction_gateOptimization <- function(x, x_threshold, kurtosis = 50){
  0.5*(1 + exp(-kurtosis * (x_threshold/x - 1)))/(x_threshold/x) # 0.5 multiplier so that value is 1 at x=x_threshold
}

# Optimization - Complex
#' @import dplyr
#' @import tidyr
#' @import rxode2
objFn_gate <- function(pars, par_names,
                       simulation_settings,
                       model, evaluation_matrix,
                       pkpd_pars, fun_generalEventMatrix, 
                       population_samples, n_sub, 
                       gate_penalties,
                       gamma = 50,
                       fun_postprocessing=NULL,
                       verbose_level=2,
                       write_file=NULL,
                       dose_scaling = NULL){
  
  # generate event matrix
  event_matrix <- fun_generalEventMatrix(pars, par_names, simulation_settings) %>%
    eventMatrix_appendEtaDemographicsSample(population_samples, n_subjects = n_sub,
                                            dose_scaling_function = dose_scaling)
  
  # run simulation
  sim_res <- model$solve(pkpd_pars, event_matrix) %>% data.frame()
  
  # apply post-processing function if given
  if(!is.null(fun_postprocessing)){
    sim_res <- fun_postprocessing(sim_res)
  }
  
  # run evaluation
  penalty_value_continuous <- evaluation_allObjectives(sim_res, evaluation_matrix,
                                                       gamma_kurtosis = gamma,
                                                       event_matrix_eval = event_matrix,
                                                       evaluation_format = 'continuous') %>% unlist()
  
  # calculate gate penalties
  gate_penalty <- activationFunction_revSwish(penalty_value_continuous, gate_penalties,
                                              kurtosis=gamma) %>% sum()
  
  # calculate sum penalty
  penalty_sum <- as.numeric(sum(penalty_value_continuous) + gate_penalty)
  
  # reporting
  if(verbose_level>=2){
    # visual report
    names(penalty_value_continuous) <- paste(evaluation_matrix$evaluation_column, evaluation_matrix$type, sep="_")
    print(round(penalty_value_continuous, 3))
    
    # write to history file at verbose level = 3
    if(verbose_level == 3){
      write.table(data.frame(t(round(penalty_value_continuous, 3))), write_file,
                  sep=",", row.names=F, col.names=F, append=T)
    }
    
  }else if(verbose_level==1){
    print(penalty_sum)
  }
  
  # return sum of weighted penalty values
  if(verbose_level > 0){
    return(penalty_sum)
  }else{
    return(penalty_value_continuous)
  }
  
  return(penalty_sum)
}

#' @import dplyr
#' @import reshape2
#' @importFrom pso psoptim
#' @importFrom optimx optimx
#' @importFrom GA ga
#' @importFrom optimization optim_sa
#' @importFrom DEoptim DEoptim
#' @importFrom ABCoptim abc_optim
coreOptimization_gate <- function(opt_object, gate_penalty,
                                  method_name = 'PSO-LBFGSB', 
                                  kurtosis = 50, 
                                  verbose = 2, opt_control=list(), 
                                  save_object=T, output_file=NULL,
                                  return_optimization_result=F){
  # check control stream
  if(grepl("PSO", method_name)){
    if(!'pso_maxit' %in% names(opt_control)){
      print(" !!! PSO number of iterations not defined. Default value used (100)")
      opt_control[['pso_maxit']] <- 100
    }
  }
  if(method_name=="PSO-LBFGSB"){
    if(!'use_previous_PSO' %in% names(opt_control) | !'use_PSO' %in% names(opt_control)){
      print(" !!! PSO prefit options incomplete. Default options used (previous PSO used if available)")
      opt_control[['use_PSO']] <- T
      opt_control[['use_previous_PSO']] <- T
    }
  }else if(method_name=="GA"){
    if(!'GA_maxit' %in% names(opt_control)){
      print(" !!! GA number of iterations not defined. Default value used (100)")
      opt_control[['GA_maxit']] <- 100
    }
    if(!'GA_parallel' %in% names(opt_control)){
      print(" !!! GA parallel processing value not defined. Single-core precess selected by default.")
      opt_control[['GA_parallel']] <- F
    }
  }else if(method_name=='ABC'){
    if(!'ABC_maxit' %in% names(opt_control)){
      print(" !!! ABC number of iterations not defined. Default value used (500).")
      opt_control[['ABC_maxit']] <- 500
    }
    if(!'ABC_nFoodSource' %in% names(opt_control)){
      print(" !!! ABC number of food source undefined. Default value used (20).")
      opt_control[['ABC_nFoodSource']] <- 20
    }
    if(!'ABC_foodLimit' %in% names(opt_control)){
      print(" !!! ABC limit number of food undefined. Default value used (50).")
      opt_control[['ABC_foodLimit']] <- 50
    }
  }
  
  # process
  if(is.null(opt_object[['populationSample']])){
    print("Population sample required!")
  }else{
    # setup loading
    sim_settings <- opt_object[['simulationSettings']]
    
    # directory setup
    if(is.null(output_file)){output_file <- paste0(opt_object[['WorkDirectory']], "GateOptimization")}
    
    # load previous fit if prompted and if file existed
    if(method_name=='PSO-LBFGSB'){
      if(opt_control[['use_previous_PSO']]){
        if(file.exists(paste0(output_file, "_PSOprefit.Rdata"))){
          # if file exist, load previous prefit and use as init for BFGS
          load(paste0(output_file, "_PSOprefit.Rdata"))
        }else{
          # if file does not exist; override 'use_previous_prefit' request
          opt_control[['use_previous_PSO']]=F
        }
      }
    }
    
    # initiate output history track if verbose=3
    if(verbose==3){
      empty_df <- data.frame(matrix(ncol = nrow(opt_object[['treatmentObjectives']]), nrow = 0))
      colnames(empty_df) <- opt_object[['treatmentObjectives']]$evaluation_column
      write_file_name <- paste0(output_file, method_name, "_OptimizationHistory.csv")
      write.csv(empty_df, write_file_name, row.names = FALSE)
    }else{
      write_file_name <- NULL
    }
    
    # method selection
    if(method_name=='PSO-LBFGSB'){
      # run PSO pre-fit
      # library(pso)
      print('### PSO ###')
      if(!opt_control[['use_previous_PSO']] & opt_control[['use_PSO']]){
        pso_prefit <- psoptim(as.numeric(opt_object[['searchSpace_gated']]$start), fn=objFn_gate,
                              par_names = rownames(opt_object[['searchSpace']]),
                              simulation_settings=sim_settings,
                              model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
                              pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']], 
                              population_samples=opt_object[['populationSample']], n_sub=sim_settings['n_subjects'],
                              gate_penalties = gate_penalty,
                              gamma = kurtosis,
                              verbose_level = verbose,
                              write_file = write_file_name,
                              dose_scaling = opt_object[['fun_doseScaling']],
                              fun_postprocessing=opt_object[['fun_postProcessing']],
                              lower=opt_object[['searchSpace']]$lb,
                              upper=opt_object[['searchSpace']]$ub,
                              control=list(maxit=opt_control[['pso_maxit']]))
        pso_prefitPars <- pso_prefit$par
        save(pso_prefit, file=paste0(output_file, "_PSOprefit.Rdata"))
      }else{
        if(!opt_control[['use_PSO']]){
          # get default starting value if use_PSO is not prompted
          pso_prefitPars <- opt_object[['searchSpace']]$start
        }else{
          # get previously fitted PSO parameter
          pso_prefitPars <- pso_prefit$par
        }
      }
      
      # run L-BFGS-B
      # library(optimx)
      print('### L-BFGS-B ###')
      optRes <- optimx(pso_prefitPars, fn=objFn_gate,
                       par_names = rownames(opt_object[['searchSpace']]),
                       simulation_settings=sim_settings,
                       model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
                       pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']], 
                       population_samples=opt_object[['populationSample']], n_sub=sim_settings['n_subjects'],
                       gate_penalties = gate_penalty,
                       gamma = kurtosis,
                       verbose_level = verbose,
                       write_file = write_file_name,
                       dose_scaling = opt_object[['fun_doseScaling']],
                       fun_postprocessing=opt_object[['fun_postProcessing']],
                       lower=opt_object[['searchSpace']]$lb,
                       upper=opt_object[['searchSpace']]$ub,
                       method = 'L-BFGS-B')
      save(optRes, file=paste0(output_file, "_PSO_LBFGSB_mainfit.Rdata"))
    }else if(method_name=='PSO'){
      # library(pso)
      print('### PSO ###')
      optRes <- psoptim(as.numeric(opt_object[['searchSpace_gated']]$start), fn=objFn_gate,
                        par_names = rownames(opt_object[['searchSpace']]),
                        simulation_settings=sim_settings,
                        model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
                        pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']], 
                        population_samples=opt_object[['populationSample']], n_sub=sim_settings['n_subjects'],
                        gate_penalties = gate_penalty,
                        gamma = kurtosis,
                        verbose_level = verbose,
                        write_file = write_file_name,
                        dose_scaling = opt_object[['fun_doseScaling']],
                        fun_postprocessing=opt_object[['fun_postProcessing']],
                        lower=opt_object[['searchSpace']]$lb,
                        upper=opt_object[['searchSpace']]$ub,
                        control=list(maxit=opt_control[['pso_maxit']]))
      save(optRes, file=paste0(output_file, method_name, "_mainfit.Rdata"))
    }else if(method_name=='GA'){
      # library(GA)
      print('### GA ###')
      optRes <- ga(type='real-valued', 
                   fitness = objFn_gate,
                   par_names = rownames(opt_object[['searchSpace']]),
                   simulation_settings=sim_settings,
                   model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
                   pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']], 
                   population_samples=opt_object[['populationSample']], n_sub=sim_settings['n_subjects'],
                   gate_penalties = gate_penalty,
                   gamma = kurtosis,
                   verbose_level = verbose,
                   write_file = write_file_name,
                   dose_scaling = opt_object[['fun_doseScaling']],
                   fun_postprocessing=opt_object[['fun_postProcessing']],
                   maxiter=opt_control[['GA_maxit']],
                   lower=opt_object[['searchSpace']]$lb,
                   upper=opt_object[['searchSpace']]$ub,
                   parallel=opt_control[['GA_parallel']])
      save(optRes, file=paste0(output_file, method_name, "_mainfit.Rdata"))
    }else if(method_name %in% c('SA', 'DE')){
      # methods requiring wrapper function
      # create wrapper function
      wrapperFunction <- function(pars){
        objFn_gate(pars, 
                   par_names = rownames(opt_object[['searchSpace']]),
                   simulation_settings=sim_settings,
                   model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
                   pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']], 
                   population_samples=opt_object[['populationSample']], n_sub=sim_settings['n_subjects'],
                   gate_penalties = gate_penalty,
                   gamma = kurtosis,
                   verbose_level = verbose,
                   write_file = write_file_name,
                   dose_scaling = opt_object[['fun_doseScaling']],
                   fun_postprocessing=opt_object[['fun_postProcessing']])
      }
      
      # run SA
      if(method_name=='SA'){
        # library(optimization)
        print('### SA ###')
        optRes <- optim_sa(fun=wrapperFunction,
                           start=as.numeric(opt_object[['searchSpace_gated']]$start),
                           maximization = F, trace = F,
                           lower=opt_object[['searchSpace']]$lb,
                           upper=opt_object[['searchSpace']]$ub)
      }else{
        # library(DEoptim)
        print("### DE ###")
        optRes <- DEoptim(wrapperFunction,
                          lower=opt_object[['searchSpace']]$lb,
                          upper=opt_object[['searchSpace']]$ub) # currently no control function (add later!)
      }
    }else if(method_name=='ABC'){
      # library(ABCoptim)
      print("### ABC ###")
      optRes <- abc_optim(as.numeric(opt_object[['searchSpace_gated']]$start), fn=objFn_gate,
                          par_names = rownames(opt_object[['searchSpace']]),
                          simulation_settings=sim_settings,
                          model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
                          pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']], 
                          population_samples=opt_object[['populationSample']], n_sub=sim_settings['n_subjects'],
                          gate_penalties = gate_penalty,
                          gamma = kurtosis,
                          verbose_level = verbose,
                          write_file = write_file_name,
                          dose_scaling = opt_object[['fun_doseScaling']],
                          fun_postprocessing=opt_object[['fun_postProcessing']],
                          lb=opt_object[['searchSpace']]$lb,
                          ub=opt_object[['searchSpace']]$ub,
                          FoodNumber = opt_control[['ABC_nFoodSource']],
                          maxCycle = opt_control[['ABC_maxit']],
                          limit = opt_control[['ABC_foodLimit']])
    }else{
      # library(optimx)
      print(paste('###', method_name, '###', sep=" "))
      
      # rename LBFGSB
      use_method_name <- ifelse(method_name=='LBFGSB', "L-BFGS-B", method_name)
      
      # other methods by optimx - remove boundary when method is BFGS
      if(method_name != 'BFGS'){
        lower_bound <- opt_object[['searchSpace']]$lb
        upper_bound <- opt_object[['searchSpace']]$ub
      }else{
        lower_bound <- NULL
        upper_bound <- NULL
      }
      
      # run optimization
      optRes <- optimx(as.numeric(opt_object[['searchSpace_gated']]$start), fn=objFn_gate,
                       par_names = rownames(opt_object[['searchSpace']]),
                       simulation_settings=sim_settings,
                       model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
                       pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']], 
                       population_samples=opt_object[['populationSample']], n_sub=sim_settings['n_subjects'],
                       gate_penalties = gate_penalty,
                       gamma = kurtosis,
                       verbose_level = verbose,
                       write_file = write_file_name,
                       dose_scaling = opt_object[['fun_doseScaling']],
                       fun_postprocessing=opt_object[['fun_postProcessing']],
                       lower=lower_bound, upper=upper_bound, 
                       method = use_method_name)
    }
    
    # save output
    save(optRes, file=paste0(output_file, method_name, "_mainfit.Rdata"))
    
    # save progress
    if(save_object){
      saveProject(opt_object)
    }
    
    if(return_optimization_result){
      return(optRes)
    }else{
      return(opt_object)
    }
  }
}
