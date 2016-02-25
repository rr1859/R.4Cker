parameterEstimationCis = function(hmm_input,reps,trstart,respstart, instart,ineqfun){
  assign("ineqfun", ineqfun, envir = .GlobalEnv)
  print("Parameter estimation.....")
  ntimes <- rep(nrow(hmm_input)/reps,reps)
  glmformula <- formula(counts ~ distance)
  mod <- depmix(response = glmformula,
               data = hmm_input,
               nstates = 3,
               trstart = trstart,
               family = gaussian(link = "identity"),
               respstart = respstart,
               instart = instart,
               ntimes = ntimes
  )
  conr <- matrix(0,5,21)
  conr[1,16] <- -1
  conr[1,19] <- 1
  conr[2,13] <- -1
  conr[2,16] <- 1
  conr[3,14] <- 1
  conr[4,17] <- 1
  conr[5,20] <- 1
  f=file()
  sink(file=f)
  mod_fit <- fit(mod,verbose = FALSE,
                conrows = conr,
                conrows.lower = c(rep(0.1,2), rep(-Inf,3)),
                conrows.upper = c(rep(Inf,2), rep(0,3)),
                solnpcntrl = list(tol = 1e-4))
  sink()
  close(f)
  pars <- c(unlist(getpars(mod_fit)))
  if(unique(pars[13:21] == respstart) == TRUE)
    print("Warning: only 1 iteration. Using starting parameters")
  in_df <- data.frame(matrix(round(pars[1:3],2), ncol = 3))
  colnames(in_df) <- c("NI", "LI", "HI")
  tr_df <- data.frame(matrix(round(pars[4:12],2), ncol = 3, byrow = TRUE))
  colnames(tr_df) <- c("NI", "LI", "HI")
  rownames(tr_df) <- c("NI", "LI", "HI")
  em_df <- data.frame(matrix(round(pars[13:21],2), ncol = 3, byrow = TRUE))
  rownames(em_df) <- c("NI", "LI", "HI")
  colnames(em_df) <- c("B0-intercept", "Distance-intercept", "Standard deviation")
  return(list(pars=pars, mod_fit = mod_fit))
}


parameterEstimationTrans = function(hmm_input,reps,trstart,respstart, instart, ineqfun){
  assign("ineqfun", ineqfun, envir = .GlobalEnv)
  print("Parameter estimation.....")
  ntimes = rep(nrow(hmm_input)/reps,reps)
  glmformula = formula(counts ~ 1)
  max_value = max(hmm_input[,1])
  min_value = max_value*0.8
  mod = depmix(response = glmformula,
               data = hmm_input,
               nstates = 3,
               trstart = trstart,
               family = gaussian(),
               respstart = respstart,
               instart = instart,
               ntimes = ntimes
  )
  conr <- matrix(0,3,18)
  conr[1,15] = -1
  conr[1,17] = 1
  conr[2,13] = -1
  conr[2,15] = 1
  conr[3,17] = 1
  mod_fit = suppressMessages(fit(mod,verbose = FALSE,conrows = conr, conrows.lower = c(rep(0.05,2),min_value), conrows.upper = c(rep(Inf,2), max_value),
                solnpcntrl = list(tol = 1e-4)))
  pars <- c(unlist(getpars(mod_fit)))
  if(unique(pars[13:18] == respstart) == TRUE)
    print("Warning: only 1 iteration. Using starting parameters")
  in_df = data.frame(matrix(round(pars[1:3],2), ncol = 3))
  colnames(in_df) = c("NI", "LI", "HI")
  tr_df = data.frame(matrix(round(pars[4:12],2), ncol = 3, byrow = TRUE))
  colnames(tr_df) = rownames(tr_df) = c("NI", "LI", "HI")
  em_df = data.frame(matrix(round(pars[13:18],2), ncol = 2, byrow = TRUE))
  rownames(em_df) = c("NI", "LI", "HI")
  colnames(em_df) = c("B0-intercept",  "Standard deviation")
  return(list(pars=pars, mod_fit = mod_fit))
}


