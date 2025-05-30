# RI calibration file generation by searching RI library, MF > 700, removing duplicates and detecting carbons
# x is path to Ri mspec file, RI library must be in ./libraries/RILIB.MSPEC output is ./AMIS/LIB/local.cal
rical <- function(x){
RI <- PreprocessMassSpectra(ReadMsp(x))
RIlib <- PreprocessMassSpectra(ReadMsp("./libraries/RILIB.MSPEC"))
sn <- NULL
for (i in 1:length(RI)){sn[i] <- as.numeric(RI[[i]]$sn)}
RI <- RI[which(sn > 50)]
results <- LibrarySearch(RI, RIlib, n_hits = 1)
mf <- NULL
for (i in 1:length(results)){mf[i] <- results[[i]]$mf}
results <- results[which(mf > 800)]
RI <- RI[which(mf > 800)]
mf <- mf[which(mf > 800)]
rt <- NULL
for (i in 1:length(RI)){rt[i] <- as.numeric(RI[[i]]$rt)}
results <- results[!duplicated(rt)]
rt <- rt[!duplicated(rt)]
mf <- mf[!duplicated(rt)]
best <- which.max(mf)
bestcarbon <- (as.numeric(results[[best]]$mw)-2)/14
cal <- cbind(rt, seq((bestcarbon-best+1)*100, length.out = length(rt),by = 100))
write.table(cal, "./AMDIS/LIB/Local.cal", row.names = FALSE, sep = ' ', col.names = FALSE)
return(cal)
}