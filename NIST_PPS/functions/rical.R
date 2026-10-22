# RI calibration file generation by searching RI library, MF > 700, removing duplicates and detecting carbons
# x is path to Ri mspec file, RI library must be in ./libraries/RILIB.MSPEC output is ./AMDIS/LIB/Local.cal
rical <- function(x){
  if (!file.exists(x)) return(NULL)
  RI <- PreprocessMassSpectra(ReadMsp(x))
  if (length(RI) == 0) return(NULL)
  RIlib <- PreprocessMassSpectra(ReadMsp("./libraries/RILIB.MSPEC"))
  sn <- sapply(RI, function(item) as.numeric(item$sn))
  RI <- RI[which(sn > 50)]
  if (length(RI) == 0) return(NULL)

  results <- LibrarySearch(RI, RIlib, n_hits = 1)
  if (length(results) == 0) return(NULL)

  mf <- sapply(results, function(item) item$mf)
  valid <- which(mf > 800)
  if (length(valid) == 0) return(NULL)

  results <- results[valid]
  RI <- RI[valid]
  mf <- mf[valid]

  rt <- sapply(RI, function(item) as.numeric(item$rt))
  dup <- duplicated(rt)
  results <- results[!dup]
  rt <- rt[!dup]
  mf <- mf[!dup]

  if (length(rt) == 0) return(NULL)

  best <- which.max(mf)
  bestcarbon <- (as.numeric(results[[best]]$mw)-2)/14
  cal <- cbind(rt, seq((bestcarbon-best+1)*100, length.out = length(rt), by = 100))
  dir.create("./AMDIS/LIB", showWarnings = FALSE, recursive = TRUE)
  write.table(cal, "./AMDIS/LIB/Local.cal", row.names = FALSE, sep = ' ', col.names = FALSE)
  return(cal)
}