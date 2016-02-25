removePCR <- function(data = data){
  if(length(data) > 1){
    q_75 <- quantile(data, 0.75)
    data[data > q_75] = round(q_75)
    return(sum(data))
  }
  else{
    return(1)
  }
}
