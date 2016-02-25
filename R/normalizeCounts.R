div_gm <-function(x){
  return(x/geometric.mean(x))
}
normalizeCounts <- function(x){
  size_f <- t(apply(x, 1, div_gm))
  size_f_median <- colMedians(size_f)
  n_matrix <- t(t(x)/as.vector(size_f_median))
  n_matrix[is.infinite(n_matrix) | is.na(n_matrix)]=0
  return(n_matrix)
}
