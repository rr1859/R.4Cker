#merge windows
merge_windows = function(windows){
  output = NULL
  chr=as.character(windows[1,1])
  start=windows[1,2]
  end=windows[1,3]
  for(i in 2:nrow(windows)){
    if(as.character(windows[i,1]) == chr & windows[i,2] <= end & windows[i,2] >= start){
      end=windows[i,3]
    }
    else{
      output=rbind(output, c(chr,start, end))
      chr=as.character(windows[i,1])
      start=windows[i,2]
      end=windows[i,3]
    }
  }
  return(data.frame(output))
}