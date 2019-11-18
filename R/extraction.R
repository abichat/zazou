
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
#' @param shift_est the estimated shifts.
#' @param method the used method.
#' @param best_alpha the best alpha during the model selection phase.
#' @param best_lambda the best lambda during the model selection phase.
#'
#' @return a list.
create_listopt_pms <- function(objective_value, shift_est, method,
                               best_alpha, best_lambda){
  list(par = shift_est, value = objective_value,
       method = paste0(method, ", part of model selection"),
       better_parameters = c(better_alpha  = best_alpha,
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
                         shift_est = df_selection$shift_est[[i]],
                         method = sub(pattern = " with model selection",
                                      replacement = "", x = shiftestim$method),
                         best_alpha = shiftestim$alpha,
                         best_lambda = shiftestim$lambda)

    list_models[[i]] <- as_shiftestim(listopt = listopt, tree = shiftestim$tree,
                                      shiftestim$zscores_obs,
                                      alpha = df_selection$alpha[i],
                                      lambda = df_selection$lambda[i])

  }

  list_models

}
