
#' Check if a shiftestim came from model selection
#'
#' @param shiftestim a shiftestim object.
#'
#' @return logical.
check_selection <- function(shiftestim){
  str <- shiftestim$method
  grepl("with model selection", str)
}

#' Create a list from
#'
#' @param objective_value the objective value.
#' @param shifts_est the estimated shifts.
#' @param method the used method.
#' @param best_alphaOU the best alpha during the model selection phase.
#' @param best_lambda the best lambda during the model selection phase.
#'
#' @return a list.
create_listopt_pms <- function(objective_value, shifts_est, method,
                               best_alphaOU, best_lambda){
  list(par = shifts_est, value = objective_value,
       method = paste0(method, ", part of model selection"),
       better_parameters = c(better_alphaOU  = best_alphaOU,
                             better_lambda = best_lambda))
}

#' Extract models computed during selection
#'
#' @param shiftestim a shiftestim object.
#'
#' @return a list of shiftestim objects.
#' @export
extract_models <- function(shiftestim){

  if(!check_selection(shiftestim)){
    warning("No selection has been done during the computation of this model.")
    return(shiftestim)
  }

  df_selection <- shiftestim$optim_info$bic_selection

  list_models <- vector("list", nrow(df_selection))

  for(i in seq_len(nrow(df_selection))){
    listopt <-
      create_listopt_pms(objective_value = df_selection$objective_value[i],
                         shifts_est = df_selection$shifts_est[[i]],
                         method = sub(pattern = " with model selection",
                                      replacement = "", x = shiftestim$method),
                         best_alphaOU = shiftestim$alphaOU,
                         best_lambda = shiftestim$lambda)

    list_models[[i]] <- as_shiftestim(listopt = listopt, tree = shiftestim$tree,
                                      shiftestim$zscores_obs,
                                      alphaOU = df_selection$alphaOU[i],
                                      lambda = df_selection$lambda[i])

  }

  list_models

}
