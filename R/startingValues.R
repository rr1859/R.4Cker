startingValuesCis <- function(window_counts,num_windows,bait_coord, num_samples, norm_counts_log, dist_log, lower, upper){
  low <- 0.5
  high <- 0.8
  non_int <- data.frame()
  int_low <- data.frame()
  int <- data.frame()
  for(i in seq(1, (num_windows-32), by = 30)){
    data_sub <- cbind(c(norm_counts_log[i:(i+30),]), rep(dist_log[i:(i+30)],num_samples))
    qhigh <- quantile(data_sub[,1], high)
    qlow <- quantile(data_sub[,1], low)
    non_int <- rbind(non_int, data_sub[data_sub[,1] <= qlow,])
    int_low <- rbind(int_low, data_sub[data_sub[,1] > qlow & data_sub[,1] < qhigh,])
    int <- rbind(int, data_sub[data_sub[,1] >= qhigh,])
  }
  colnames(non_int) <- c("counts", "distance")
  colnames(int_low) <- c("counts", "distance")
  colnames(int) <- c("counts", "distance")
  non_int <- non_int[order(non_int[,2]),]
  int_low <- int_low[order(int_low[,2]),]
  int <- int[order(int[,2]),]
  glmformula <- formula(counts ~ distance)
  nint_par <- fitdistr(non_int[,1], "normal",lower=0.001)$estimate
  int_low_par <- fitdistr(int_low[,1], "normal",lower=0.001)$estimate
  int_par <- fitdistr(int[,1], "normal",lower=0.001)$estimate
  nint_glm <- glm(glmformula, data = non_int, family = gaussian(link = "identity"))
  int_low_glm <- glm(glmformula, data = int_low, family = gaussian(link = "identity"))
  int_glm <- glm(glmformula, data = int, family = gaussian(link = "identity"))
  int_low_glm[[1]][1] <- upper*int_glm[[1]][1]
  nint_glm[[1]][1] <- lower*int_glm[[1]][1]
  trstart <- c(0.5,0.25,0.25, 0.25,0.5,0.25,0.25,0.25,0.5)
  respstart <- c(nint_glm[[1]],nint_par[2]^2, int_low_glm[[1]], int_low_par[2]^2, int_glm[[1]], int_par[2]^2)
  instart <- c(0.5,0.5,0.5)
  instart_df <- data.frame(matrix(instart,ncol =3))
  colnames(instart_df) <- c("Non-interacting", "Low-interacting","Highly-interacting" )
  trstart_df <- data.frame(matrix(trstart, ncol = 3, byrow = TRUE))
  colnames(trstart_df) <-  c("NI", "LI", "HI")
  rownames(trstart_df) <- c("NI", "LI", "HI")
  em_df <- data.frame(matrix(respstart,ncol = 3, byrow=TRUE))
  colnames(em_df) <- c("B0-intercept", "Distance-intercept", "Standard deviation")
  rownames(em_df) <- c("NI", "LI", "HI")
  return(list(instart = instart, trstart = trstart, respstart = respstart, int_glm =int_glm))
}

startingValuesTrans = function(synth_hmm_input, lower, upper){
  non_int <- data.frame()
  int_low <- data.frame()
  int <- data.frame()
  qhigh <- quantile(synth_hmm_input[,1], upper)
  qlow <- quantile(synth_hmm_input[,1], lower)
  non_int <- data.frame(counts = synth_hmm_input[synth_hmm_input[,1] <= qlow,])
  int_low <- data.frame(counts = synth_hmm_input[synth_hmm_input[,1] > qlow & synth_hmm_input[,1] < qhigh,])
  int <- data.frame(counts = synth_hmm_input[synth_hmm_input[,1] >= qhigh ,])
  glmformula = formula(counts ~ 1)
  nint_par = fitdistr(non_int[,1], "normal",lower=0.001)$estimate
  int_low_par = fitdistr(int_low[,1], "normal",lower=0.001)$estimate
  int_par = fitdistr(int[,1], "normal",lower=0.001)$estimate
  nint_glm = glm(glmformula, data = non_int, family = gaussian())
  int_low_glm = glm(glmformula, data = int_low, family = gaussian())
  int_glm = glm(glmformula, data = int, family = gaussian())
  int_low_glm[[1]][1] = 0.8*int_glm[[1]][1]
  nint_glm[[1]][1] = 0.5*int_glm[[1]][1]
  trstart = c(0.5,0.25,0.25, 0.25,0.5,0.25,0.25,0.25,0.5)
  respstart =c(nint_glm[[1]],nint_par[2]^2, int_low_glm[[1]], int_low_par[2]^2, int_glm[[1]], int_par[2]^2)
  instart = c(0.5,0.5,0.5)
  instart_df = data.frame(matrix(instart,ncol=3))
  colnames(instart_df) = c("Non-interacting", "Low-interacting","Highly-interacting" )
  trstart_df = data.frame(matrix(trstart, ncol = 3, byrow = TRUE))
  colnames(trstart_df) = rownames(trstart_df) = c("NI", "LI", "HI")
  em_df = data.frame(matrix(respstart,ncol = 2, byrow=TRUE))
  colnames(em_df) = c("B0-intercept", "Standard deviation")
  rownames(em_df) = c("NI", "LI", "HI")
  return(list(instart = instart, trstart = trstart, respstart = respstart, int_glm =int_glm))
}
