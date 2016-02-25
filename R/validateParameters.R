validateParametersCis <- function(glm, pars ,hmm_input){
  nint_glm <- glm
  lint_glm <- glm
  int_glm <- glm
  int_glm$coefficients <- pars[19:20]
  lint_glm$coefficients <- pars[16:17]
  nint_glm$coefficients <- pars[13:14]
  predict_min_nint <- predict(nint_glm, hmm_input[which.min(hmm_input[,2]),])
  predict_min_lint <- predict(lint_glm, hmm_input[which.min(hmm_input[,2]),])
  predict_min_int <- predict(int_glm, hmm_input[which.min(hmm_input[,2]),])
  predict_max_nint <- predict(nint_glm, hmm_input[which.max(hmm_input[,2]),])
  predict_max_lint <- predict(lint_glm, hmm_input[which.max(hmm_input[,2]),])
  predict_max_int <- predict(int_glm, hmm_input[which.max(hmm_input[,2]),])
  min_dist <- c(predict_min_nint, predict_min_lint, predict_min_int)
  max_dist <- c(predict_max_nint, predict_max_lint, predict_max_int)
  pars_rows <- c(13,16,19)
  order_min <- order(min_dist)
  order_max <- order(max_dist)
  if(all(order(min_dist) == 1:3) & all(order(max_dist) == 1:3)){
    return(list(pars=pars, valid= TRUE))
  } else if(all(order_min == order_max)){
    nint_rows <- pars_rows[order_min[1]]
    lint_rows <- pars_rows[order_min[2]]
    hint_rows <- pars_rows[order_min[3]]
    nint_pars <- pars[nint_rows:(nint_rows+2)]
    lint_pars <- pars[lint_rows:(lint_rows+2)]
    hint_pars <- pars[hint_rows:(hint_rows+2)]
    pars[19:21] <- hint_pars
    pars[16:18] <- lint_pars
    pars[13:15] <- nint_pars
    return(list(pars=pars, valid = TRUE))
  } else {
    return(list(pars=pars, valid = FALSE))
  }
}
