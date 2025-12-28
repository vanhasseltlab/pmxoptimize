#GRID SEARCH---------------
#' Draw Grid Search Simulation Agenda
#'
#' Generates a grid of treatment parameter combinations for evaluation using with grid search.
#' Parameter values are created based on the lower and upper bounds defined in the optimization object's
#' \code{searchSpace}. The generated grid is saved to file and stored in the object.
#'
#' @import dplyr
#' @importFrom rlist list.rbind
#' @param opt_object An optimization task object as created by \code{initiateOptimizationTask}. Must include a defined \code{searchSpace}.
#' @param round_inputs A logical or integer vector indicating which treatment parameters should be rounded to whole numbers. If logical, \code{TRUE} values indicate rounding.
#' @param intersect_points Optional named vector. If provided, generates a sub-grid where one or more parameters are fixed at specified values (e.g. the identified optima). Parameters not mentioned in this vector are allowed to vary.
#' @param grid_resolution Integer. Number of points to generate per parameter range. Default is 50.
#'
#' @return The updated \code{opt_object} with a new \code{grid} element and status. The grid is also saved to \code{GridSearch.csv} in the project directory.
#' @export
drawSimulationGrid <- function(opt_object, round_inputs=c(), 
                               intersect_points=NULL, grid_resolution=50){
 
  # setup from opt object
  dose_inputs <- opt_object[['searchSpace']]
  
  # Create a list where each element is a sequence based on lb and ub
  seq_list <- lapply(seq_len(nrow(dose_inputs)), function(i) {
    seq(dose_inputs$lb[i], dose_inputs$ub[i], length.out = grid_resolution)
  })
  names(seq_list) <- rownames(dose_inputs)
  
  # expand grid to create simulation agenda
  par_values <- expand.grid(seq_list)
  
  # rounding based on input
  if(is.logical(round_inputs)) round_inputs <- which(round_inputs)
  for(i in round_inputs){
    par_values[,i] <- round(par_values[,i], 0)
  }
  par_values <- distinct(par_values)
  
  # select only intersects of the optima
  if(!is.null(intersect_points)){
    sim_agenda <- lapply(c(1:length(intersect_points)), function(xi){
      applied_par_values <- par_values
      applied_par_values[,names(intersect_points)[xi]] <- intersect_points[xi]
      applied_par_values <- distinct(applied_par_values)
      return(applied_par_values)
    }) %>% list.rbind() %>% distinct()
  }else{
    sim_agenda <- par_values
  }
  rownames(sim_agenda) <- NULL
  
  # append to object
  opt_object[['grid']] <- sim_agenda
  opt_object[['grid_run_status']] <- "Not run"
  
  # save grid
  write.csv(sim_agenda, paste0(opt_object[['WorkDirectory']], "GridSearch.csv"), row.names=F)
  
  # save object
  saveProject(opt_object)
  
  return(opt_object)
}

#' Manually Add Grid Search Simulation Agenda
#'
#' Receives a simulation agenda (table) for running Grid Search. Simulation agenda format needs to follow the format drawn by \code{drawSimulationGrid}.  
#'
#' @param opt_object An optimization task object as created by \code{initiateOptimizationTask}. 
#' @param sim_agenda Grid simulation agenda (table). See \code{drawSimulationGrid}.
#'
#' @return The updated \code{opt_object} with a new \code{grid} element and status. The grid is also saved to \code{GridSearch.csv} in the project directory.
#' @export
addManualSimulationGrid <- function(opt_object, sim_agenda){
  # append to object
  opt_object[['grid']] <- sim_agenda
  opt_object[['grid_run_status']] <- "Not run"
  
  # save grid
  write.csv(sim_agenda, paste0(opt_object[['WorkDirectory']], "GridSearch.csv"), row.names=F)
  
  # save object
  saveProject(opt_object)
  
  return(opt_object)
}

#' Run Grid Search
#'
#' Runs PKPD simulations over a grid of treatment parameter combinations stored in the optimization object.
#' Function \code{drawSimulationGrid} need to first be run before running this function to generate the grid table.
#' The results are evaluated using the objective function, and penalties are written to a result file.
#' Previously evaluated points are skipped unless \code{redo_run = TRUE}.
#' 
#' @import dplyr
#' @import rlist
#' @import parallel
#' @param opt_object An optimization task object that includes a valid simulation grid. The grid should be created using \code{drawSimulationGrid()}.
#' @param output_filename Optional character string specifying the name of the CSV file to save results. Default is \code{"GridSearchResult.csv"}.
#' @param n_cores Integer. Number of CPU cores to use for parallel simulation. If \code{n_cores > 1}, uses \code{mclapply()} for parallel execution. Default is 1 (sequential).
#' @param redo_run Logical. If TRUE, re-evaluates all grid points even if results already exist. Default is FALSE.
#' @param penalty_type string. Options are 'binary' (target failure rate) or 'continuous' (penalty value). Default is 'binary'. 
#'
#' @return The updated \code{opt_object} with updated \code{grid_run_status}. The results are appended to a CSV file in the working directory.
#' @export
gridSearchEvaluation <- function(opt_object, output_filename=NULL,
                                 n_cores=1, redo_run=F, penalty_type='binary'){
  if(is.null(opt_object[['grid']])){
    print("Grid missing! Please first draw the simulation grid using 'drawSimulationGrid'")
  }else{
    # setup
    sim_inputs <- opt_object[['simulationSettings']]
    n_subjects <- sim_inputs[['n_subjects']]
    grid_table <- opt_object[['grid']]
    
    # check previously generated values
    if(is.null(output_filename)){output_filename <- "GridSearchResult.csv"}
    output_file <- paste0(opt_object[['WorkDirectory']], output_filename)
    if(file.exists(output_file) & !redo_run){
      # read previously ran grid
      previous_grid <- read.csv(output_file, header=T) 
      
      if(penalty %in% colnames(previous_grid)){
        previous_grid <- dplyr::select(previous_grid, -penalty)
      }
      
      # summarize condition in grid
      previous_grid$condition <- apply(previous_grid, 1, function(x) paste(x, collapse='_'))
      
      # remove executed runs from the simulation agenda
      grid_table$condition <- apply(grid_table, 1, function(x) paste(x, collapse='_'))
      grid_table <- filter(grid_table, !condition %in% previous_grid$condition) %>% 
        dplyr::select(-condition)
    }else{
      # initiate output file
      empty_df <- data.frame(matrix(ncol = (ncol(grid_table)+
                                              nrow(opt_object[['treatmentObjectives']])), nrow = 0))
      obj_names <- apply(opt_object[['treatmentObjectives']], 1, function(x){
        paste("penalty", x['evaluation_column'], x['type'], sep="_")
      })
      colnames(empty_df) <- c(colnames(grid_table), obj_names)
      write.csv(empty_df, output_file, row.names = FALSE)
    }
    
    # status update
    opt_object[['grid_run_status']] <- "Run started . . . "
    opt_object[['grid_run_output_file']] <- output_filename
    saveProject(opt_object)
    print("Starting grid search . . . ")
    
    # run through the simulation agenda
    if(n_cores<=1){
      nullCallBack <- lapply(c(1:nrow(grid_table)), function(xi){
        penalty_value <- objFn(unlist(grid_table[xi,]), colnames(grid_table),
                               sim_inputs,
                               model=opt_object[['pkpdModel']], 
                               evaluation_matrix=opt_object[['treatmentObjectives']],
                               pkpd_pars=opt_object[['modelParameters']], 
                               fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']],
                               population_samples=opt_object[['populationSample']], n_sub=n_subjects,
                               fun_postprocessing=opt_object[['fun_postProcessing']],
                               verbose_level=0, write_file=NULL,
                               dose_scaling = opt_object[['fun_doseScaling']],
                               eval_format = penalty_type)
        
        penalty_value <- c(unlist(grid_table[xi,]), penalty_value)
        
        write.table(data.frame(t(penalty_value)), output_file,
                    sep=",", row.names=F, col.names=F, append=T)
        return(NULL)
      })
    }else{
      nullCallBack <- mclapply(c(1:nrow(grid_table)), function(xi){
        penalty_value <- objFn(unlist(grid_table[xi,]), colnames(grid_table),
                               sim_inputs,
                               model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
                               pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']],
                               population_samples=opt_object[['populationSample']], n_sub=n_subjects,
                               fun_postprocessing=opt_object[['fun_postProcessing']],
                               verbose_level=0,
                               write_file=NULL,
                               dose_scaling = opt_object[['fun_doseScaling']],
                               eval_format = penalty_type)
        penalty_value <- c(unlist(grid_table[xi,]), penalty_value)
        
        write.table(data.frame(t(penalty_value)), output_file,
                    sep=",", row.names=F, col.names=F, append=T)
        return(NULL)
      }, mc.cores=n_cores)
    }
  }
  
  # reporting
  opt_object[['grid_run_status']] <- "Completed"
  saveProject(opt_object)
  print("Grid search completed.")
  
  return(opt_object)
}

#PARETO FRONT APPROXIMATE FROM GRID SEARCH----------
#' Approximate Pareto Front from Grid Search Result
#'
#' Draws a list of Pareto-optimal points from grid search result. 
#' If a previous Pareto front approximate is available in the project directory, this result will be used unless redraw_front is TRUE.
#' Function returns error message if neither GridSearchResult.csv nor ParetoFrontApproximate.csv was found.
#'
#' Note: \code{\link{gridSearchEvaluation}} must be run beforehand to generate the grid from which the Pareto front will be approximated from.
#'
#' @import dplyr
#' @import rPref
#' @import purrr
#' @param opt_object An optimization task object as created by \code{initiateOptimizationTask}. Must include previously run optimization history.
#' @param write_output Boolean prompt to write the Pareto front approximate as output. Filename is 'ParetoFrontApproximate.csv' by default.
#' @param redraw_front Boolean prompt to manually require the function to re-draw Pareto front approximate from the latest grid search result.
#'
#' @return The updated \code{opt_object} with a new \code{paretoGrid} element and status. The grid is also saved to \code{ParetoGrid.csv} in the project directory.
#' @export
approximateParetoFront <- function(opt_object, write_output=T, redraw_front = F,
                                   grid_search_result_filename = NULL, output_file_name = NULL){
  # get grid run name
  if(is.null(grid_search_result_filename)){
    grid_search_result_filename <- ifelse(!is.null(opt_object[['grid_run_output_file']]),
                                          "GridSearchResult.csv", opt_object[['grid_run_output_file']])
  }
  
  # attempt to look for grid search result
  grid_result_file <- paste0(opt_object$WorkDirectory, "/", grid_search_result_filename)
  pareto_front_file <- paste0(opt_object$WorkDirectory, ifelse(is.null(output_file_name), 
                                                                       "/ParetoFrontApproximate.csv", paste0("/", output_file_name)))
  if(file.exists(grid_result_file) & (!file.exists(pareto_front_file) | redraw_front)){
    # read grid search result
    grid_search <- read.csv(grid_result_file, header=TRUE)
    
    # Build the expression: product of the low of all columns starting with "penalty_"
    penalty_cols <- names(grid_search)[grepl("^penalty_", names(grid_search))]
    pref_expr <- paste0("low(", penalty_cols, ")", collapse = " * ")
    
    # get list of pareto fronts
    pareto_fronts <- psel(grid_search, eval(parse(text = pref_expr)))
    
    # output pareto optimal points
    if(write_output){
      write.csv(pareto_fronts, pareto_front_file, row.names=F)
    }
  }else if(file.exists(pareto_front_file) & !redraw_front){
    # read previous pareto front approximate if output file exist and re-read prompt not given
    pareto_fronts <- read.csv(pareto_front_file, header=TRUE)
  }else{
    pareto_fronts <- "Error: Grid search result not found. Please run Grid Search Evaluation!"
  }
  
  # save result
  opt_object[['approx_pareto_front']] <- pareto_fronts
  saveProject(opt_object)
  
  return(opt_object)
}