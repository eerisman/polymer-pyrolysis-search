#library searching using package mssearchr
#library mspec must be "./libraries/pyrolysis.MSPEC"
#returns nested list of polymer.list as polymer result, peaks and pyrolyzate results, 
#query and lib as massearchr PreprocessedMassSpectra
mssearch <- function(x){
  query <- PreprocessMassSpectra(ReadMsp(x))  
  lib <- PreprocessMassSpectra(ReadMsp("./libraries/pyrolysis.MSPEC"))
  #signal to noise filter only take greater than 50
  sn <- NULL
  for (i in 1:length(query)){sn[i] <- as.numeric(query[[i]]$sn)}
  query <- query[which(sn >50)]
  
  #remove dup rts
  rt <- NULL
  for (i in 1:length(query)){rt[i] <- as.numeric(query[[i]]$rt)}
  query <- query[!duplicated(rt)]
  
  #library search
  results <- LibrarySearch(query, lib, n_hits = 20)
  
  #add query index to results
  for (i in 1:length(results)){results[[i]]$qindx <- i}
  
  #convert nested list to data frame
  results <- as.data.frame(do.call(rbind, results))
  results <- subset(results, results$mf>700)
  
  #modify mf delta RI
  for (i in 1:nrow(results)){results$dri[i] <- abs(as.numeric(query[[results$qindx[i]]]$ri)- as.numeric(lib[[results$idx[i]]]$ri))}
  #RI tolerance
  for (i in 1:nrow(results)){results$tol[i] <- as.numeric(query[[results$qindx[i]]]$ri)/100}
  # calc match factor penalty
  for (i in 1:nrow(results)){results$mfp[i] <- if (results$dri[i] < results$tol[i]){0}else{50*(results$dri[i]-results$tol[i])/results$tol[i]}}
  # calc mod mf
  for (i in 1:nrow(results)){results$modmf[i] <-results$mf[i] - results$mfp[i]}
  # ri for display
  for (i in 1:nrow(results)){results$ri[i] <- as.numeric(query[[results$qindx[i]]]$ri)}
  # rt for display
  for (i in 1:nrow(results)){results$rt[i] <- as.numeric(query[[results$qindx[i]]]$rt)}
  #subset results for modmf
  results <- subset(results, results$modmf >700)
  results <- results[order(results$qindx, -results$modmf),]
  write.csv(results, "./results/searchresults.csv", row.names = FALSE)
  #top hit
  tophit <- results[!duplicated(results$qindx),]
  #ordering by modmf
  tophit <- tophit[order(-tophit$modmf),] #added
  tophit <- tophit[!duplicated(tophit$idx),]
  #other hits
  otherhits <- results[duplicated(results$qindx),]
  otherhits <- subset(otherhits, !(otherhits$idx %in% tophit$idx))
  otherhits <- otherhits[order(-otherhits$modmf),] 
  otherhits <- otherhits[!duplicated(otherhits$idx),]
  #combine and sort tophit and other hits
  peaks <- rbind(tophit, otherhits)
  peaks <- peaks[order(peaks$qindx, -peaks$modmf),]
  peaks[,c(2,8,9,10,11)] <- round(peaks[,c(2,8,9,10,11)], 2)## rounding
  write.csv(peaks, "./results/peaksresults.csv", row.names = FALSE)
  for (i in 1:nrow(peaks)){peaks$polymer[i] <- lib[[peaks[i,3]]]$polymer}#polymer generate polymer list
  polymer.list <- unique(unlist(strsplit(peaks$polymer, "; ")))#unique polymers
  #polymer attribution
  count <- NULL
  for (i in 1:length(polymer.list)){count[i] <- length(grep(gsub("\\(.*", "",polymer.list[i]), peaks$polymer))}
  #libery pyrolyate count
  libcount <- NULL
  for (i in 1:length(polymer.list)){libcount[i] <- as.integer(gsub(".*\\((.+)\\)", "\\1",polymer.list[i]))}
  # average modmf
  meanmf <- NULL
  for (i in 1:length(polymer.list)){meanmf[i] <- mean(peaks$modmf[grep(gsub("\\(.*", "",polymer.list[i]),peaks$polymer)])}
  #combine polymer list
  polymer.list <- cbind.data.frame(polymer.list, count, libcount, meanmf)
  polymer.list$score <- polymer.list$meanmf*(polymer.list$count/polymer.list$libcount)
  polymer.list[,c(4,5)] <- round(polymer.list[,c(4,5)], 2)##rounding
  polymer.list <- polymer.list[order(-polymer.list$score),]
  write.csv(polymer.list, "./results/polymerlist.csv", row.names = FALSE)
  return(list(polymer.list,peaks, query, lib))
}