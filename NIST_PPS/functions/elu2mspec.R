elu2mspec <- function(x){
  temp <- readLines(x)
  temp <- paste(temp, collapse = "") # collapse into one vector
  temp <- strsplit(gsub("(NAME)","~\\1",temp),"~")[[1]] # gsub("(NAME)","~\\1",test1) adds ~ to split on
  temp <- temp[-1] #remove heading blank line
  temp <- gsub("\\([0-9]*,[0-9]* [A-Z][0-9].[0-9]\\)", "", temp) #remove uncertain peaks
  temp <- gsub("NAME: .*SN", "NAME: SN", temp) # delete between NAME: and SN
  temp <- gsub("\\|.*RT", " RT", temp) # delete between | and RT
  temp <- gsub("\\|.*RI", " RI", temp) # delete between | and RI
  temp <- gsub("\\|.*NUM PEAKS:", "NUM PEAKS:", temp) #delete between | and NUM PEAKS:
  temp <- gsub("NAME: ", "\\\n\nNAME:", temp) #new line for NAME:
  temp <- gsub("SN", "\\\nSN: ", temp) #format for SN
  temp <- gsub("RT", "\\\nRT: ", temp) #format for RT
  temp <- gsub("RI", "\\\nRI: ", temp) #new line for Retention index
  temp[1] <- sub("\\\n\\\n", "", temp[1]) #remove leading 2 blank lines
  for (i in 1:length(temp)) {temp[i] <- sub("NAME:", paste("NAME: Peak",i), temp[i])} #add name to the compounds
  for (i in 1:length(temp)) {temp[i] <- sub("NUM PEAKS: [0-9]*", paste0("\\\nNUM PEAKS: ", length(gregexpr("\\(", temp[i])[[1]]), "\\\n"), temp[i])}
  temp <- gsub(" )", ")", temp)
  temp <- c(temp, "\n")
  y <- gsub(".*/", "", x)
  y <- sub("ELU", "MSPEC", y)
  y <- gsub(" ", "_", y)
  cat(temp, file = paste0("./data/", y))
  return(y)
  }
