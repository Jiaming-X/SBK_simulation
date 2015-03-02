####################################################################################################
####################################################################################################
## Code for processing the Solar data in Seattle

readFile<-function( afile = afile){
  con = file(afile, 'rt')
  aFile = readLines(con, warn = FALSE)
  theName = strsplit(aFile[1], split = "\t", perl = TRUE)[[1]][1:3]
  #theName = paste(theName[[1]][1:3], collapse = "#")
  #print(theName)
  aFile = aFile[-1]
  index = which( "" ==  aFile)
  #print(index)
  #print(length(index))
  if( length(index) > 0){
    aFile = aFile[-c(index)]
  }
 
  temp2 = lapply(aFile, function(i){
    temp = strsplit(i, split = "\t", perl = TRUE)
    return(as.numeric(temp[[1]][1:3]))
  })
 
  num = length(temp2)
  temp2 = matrix(unlist(temp2), byrow = TRUE, nrow = num)
  colnames( temp2 ) = theName
  
  close(con)
  return( temp2 )
}

main<-function(){
  setwd("~/Documents/RTG STATISTICS/solar_Data/Seattle/")
  somefiles = list.files("~/Documents/RTG STATISTICS/solar_Data/Seattle/", recursive = TRUE)
  
  tt = lapply(somefiles, readFile)
  names(tt) = somefiles

  return(tt)
}

FIND <- function(data = data, duration = duration, year = year, month = month){
  year = as.character( year %% 100)
  if( month >= 10){
    month = as.character( month)
  }else{
    month = paste0("0", as.character( month))
  }
  theName = paste0( duration, year, month)
  pattern = paste0("(SEP|SER)", theName)
  Find = grep(pattern = pattern, names(temp))
  
  if( length(Find) > 0 ){
    print( data[Find] )
  }else{
    print("The file does not exit")
  }
}

library(data.table)
setwd("~/Documents/RTG STATISTICS/solar_Data/Seattle/")
temp = main()
#FIND(temp, "Q", 2010, 12)
#save(temp, file = "Seattle.rda")


