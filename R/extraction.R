
#' Check if a shiftpunct came from model selection
#'
#' @param shiftpunct a shiftpunct object.
#'
#' @return logical.
check_selection <- function(shiftpunct){
  str <- shiftpunct$method
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
#' @param shiftpunct a shiftpunct object.
#'
#' @return a list of shiftpunct objects.
#' @export
extract_models <- function(shiftpunct){

  if(!check_selection(shiftpunct)){
    warning("No selection has been done during the computation of this model.")
    return(shiftpunct)
  }

  df_selection <- shiftpunct$optim_info$bic_selection

  list_models <- vector("list", nrow(df_selection))

  for(i in seq_len(nrow(df_selection))){
    listopt <-
      create_listopt_pms(objective_value = df_selection$objective_value[i],
                         shifts_est = df_selection$shifts_est[[i]],
                         method = sub(pattern = " with model selection",
                                      replacement = "", x = shiftpunct$method),
                         best_alphaOU = shiftpunct$alphaOU,
                         best_lambda = shiftpunct$lambda)

    list_models[[i]] <- as_shiftpunct(listopt = listopt, tree = shiftpunct$tree,
                                      shiftpunct$zscores_obs,
                                      alphaOU = df_selection$alphaOU[i],
                                      lambda = df_selection$lambda[i])

  }

  list_models

}
