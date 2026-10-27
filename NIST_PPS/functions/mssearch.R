# library searching using package mssearchr
# library mspec must be "./libraries/pyrolysis.MSPEC"
# returns nested list of polymer.list as polymer result, peaks and pyrolyzate results,
# query and lib as mssearchr PreprocessedMassSpectra
mssearch <- function(x, snc, mfc, rip, waive_out_of_range = TRUE){
  if (!file.exists(x)) return(NULL)
  query <- PreprocessMassSpectra(ReadMsp(x))  
  lib <- PreprocessMassSpectra(ReadMsp("./libraries/pyrolysis.MSPEC"))
  
  empty_polymer <- data.frame(
    Polymer = character(0),
    `Markers Found` = numeric(0),
    `Markers Expected` = numeric(0),
    `Mean Match Factor` = numeric(0),
    `Final Score` = numeric(0),
    check.names = FALSE, stringsAsFactors = FALSE
  )

  empty_peaks <- data.frame(
    Name = character(0), `Match Factor` = numeric(0), idx = integer(0),
    Formula = character(0), `Molecular Mass` = numeric(0), qindx = integer(0),
    `Delta RI` = numeric(0), `RI tolerance` = numeric(0), `Match Factor Penalty` = numeric(0),
    `Modified Match Factor` = numeric(0), `Retention Index (RI)` = numeric(0),
    `Retention Time` = numeric(0), polymer = character(0),
    check.names = FALSE, stringsAsFactors = FALSE
  )

  if (length(query) == 0) return(list(empty_polymer, empty_peaks, query, lib))

  # signal to noise filter only take greater than cutoff "snc"
  sn <- sapply(query, function(q) as.numeric(q$sn))
  query <- query[which(sn > snc)]
  if (length(query) == 0) return(list(empty_polymer, empty_peaks, query, lib))

  # remove duplicate RTs
  rt <- sapply(query, function(q) as.numeric(q$rt))
  query <- query[!duplicated(rt)]
  if (length(query) == 0) return(list(empty_polymer, empty_peaks, query, lib))
  
  # library search
  results <- LibrarySearch(query, lib, n_hits = 20, hitlist_columns = c("name", "formula", "mw"))
  if (length(results) == 0) return(list(empty_polymer, empty_peaks, query, lib))
  
  # add query index to results
  for (i in seq_along(results)){
    if (!is.null(results[[i]]) && is.data.frame(results[[i]]) && nrow(results[[i]]) > 0) {
      results[[i]]$qindx <- i
    }
  }
  
  results <- as.data.frame(do.call(rbind, results))
  if (is.null(results) || nrow(results) == 0) return(list(empty_polymer, empty_peaks, query, lib))

  results <- subset(results, results$mf > mfc)
  if (nrow(results) == 0) return(list(empty_polymer, empty_peaks, query, lib))
  
  # Vectorized calculation of dri, tol, mfp, modmf, ri, rt
  q_ri <- sapply(query, function(q) as.numeric(q$ri))
  q_rt <- sapply(query, function(q) as.numeric(q$rt))
  lib_ri <- sapply(lib, function(l) as.numeric(l$ri))

  query_ri_vec <- q_ri[results$qindx]
  query_rt_vec <- q_rt[results$qindx]
  lib_ri_vec <- lib_ri[results$idx]

  results$dri <- abs(query_ri_vec - lib_ri_vec)
  results$tol <- query_ri_vec / 100
  results$mfp <- ifelse(results$dri < results$tol, 0, rip * (results$dri - results$tol) / results$tol)

  # RI values outside calibration range are suspect and should not be penalized if selected
  if (isTRUE(waive_out_of_range) && file.exists("./AMDIS/LIB/Local.cal")) {
    cal <- tryCatch(read.csv("./AMDIS/LIB/Local.cal", sep = "", header = FALSE), error = function(e) NULL)
    if (!is.null(cal) && nrow(cal) > 0 && ncol(cal) >= 2) {
      min_ri <- min(cal[, 2], na.rm = TRUE)
      max_ri <- max(cal[, 2], na.rm = TRUE)
      outside_cal <- (query_ri_vec < min_ri | query_ri_vec > max_ri | lib_ri_vec < min_ri | lib_ri_vec > max_ri)
      outside_cal[is.na(outside_cal)] <- TRUE
      results$mfp[outside_cal] <- 0
    }
  }

  results$modmf <- results$mf - results$mfp
  results$ri <- query_ri_vec
  results$rt <- query_rt_vec

  results <- subset(results, results$modmf > mfc)
  if (nrow(results) == 0) return(list(empty_polymer, empty_peaks, query, lib))

  results <- results[order(results$qindx, -results$modmf),]

  # Setup unique parameter-and-timestamp-named subfolder in ./results/
  sample_name <- tools::file_path_sans_ext(basename(x))
  if (is.null(sample_name) || sample_name == "") sample_name <- "search"
  timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
  param_suffix <- sprintf("snc%s_mfc%s_rip%s_waive%s",
                          snc, mfc, rip,
                          ifelse(isTRUE(waive_out_of_range), "TRUE", "FALSE"))

  run_folder_name <- sprintf("%s_%s_%s", sample_name, param_suffix, timestamp_str)
  run_dir <- file.path("./results", run_folder_name)

  dir.create("./results", showWarnings = FALSE, recursive = TRUE)
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

  # Record parameter summary for repeating runs
  param_txt <- c(
    "NIST Polymer Pyrolysis Search - Run Parameters",
    "==============================================",
    paste("Date & Time:", Sys.time()),
    paste("Sample Name:", sample_name),
    paste("Input MSPEC File:", x),
    paste("S/N Cutoff (snc):", snc),
    paste("Match Factor Cutoff (mfc):", mfc),
    paste("RI Penalty Rate (rip):", rip),
    paste("Waive RI Penalty Outside Calibration Range (waive_out_of_range):", isTRUE(waive_out_of_range))
  )
  writeLines(param_txt, file.path(run_dir, sprintf("%s_search_parameters_%s_%s.txt", sample_name, param_suffix, timestamp_str)))

  # Write raw search results parameter- and time-stamped inside the unique subfolder
  write.csv(results, file.path(run_dir, sprintf("%s_rawsearchresults_%s_%s.csv", sample_name, param_suffix, timestamp_str)), row.names = FALSE)

  # top hit
  tophit <- results[!duplicated(results$qindx),]
  tophit <- tophit[order(-tophit$modmf),]
  tophit <- tophit[!duplicated(tophit$idx),]

  # other hits
  otherhits <- results[duplicated(results$qindx),]
  if (nrow(otherhits) > 0) {
    otherhits <- subset(otherhits, !(otherhits$idx %in% tophit$idx))
    if (nrow(otherhits) > 0) {
      otherhits <- otherhits[order(-otherhits$modmf),]
      otherhits <- otherhits[!duplicated(otherhits$idx),]
    }
  }

  peaks <- rbind(tophit, otherhits)
  if (nrow(peaks) == 0) return(list(empty_polymer, empty_peaks, query, lib))

  peaks <- peaks[order(peaks$qindx, -peaks$modmf),]
  peaks[, c(2, 7, 8, 9, 10)] <- round(peaks[, c(2, 7, 8, 9, 10)], 2)
  colnames(peaks) <- c("Name", "Match Factor", "idx", "Formula", "Molecular Mass", "qindx", "Delta RI", "RI tolerance", "Match Factor Penalty", "Modified Match Factor", "Retention Index (RI)", "Retention Time")

  # Write filtered search results parameter- and time-stamped inside the unique subfolder
  write.csv(peaks, file.path(run_dir, sprintf("%s_filteredsearchresults_%s_%s.csv", sample_name, param_suffix, timestamp_str)), row.names = FALSE)

  # polymer list generation
  peaks$polymer <- sapply(peaks$idx, function(idx_val) {
    p <- lib[[idx_val]]$polymer
    if (is.null(p)) "" else p
  })

  polymer_entries <- unlist(strsplit(peaks$polymer, "; "))
  polymer_entries <- polymer_entries[polymer_entries != ""]
  polymer.list <- unique(polymer_entries)

  if (length(polymer.list) == 0) {
    return(list(empty_polymer, peaks, query, lib))
  }

  base_names <- sub("\\s*\\([0-9]+\\)$", "", polymer.list)

  count <- numeric(length(polymer.list))
  libcount <- numeric(length(polymer.list))
  meanmf <- numeric(length(polymer.list))

  for (i in seq_along(polymer.list)) {
    base_name <- base_names[i]
    matched_rows <- which(sapply(peaks$polymer, function(p_str) grepl(base_name, p_str, fixed = TRUE)))
    count[i] <- length(matched_rows)

    m <- regmatches(polymer.list[i], regexec("\\(([0-9]+)\\)$", polymer.list[i]))[[1]]
    libcount[i] <- if (length(m) > 1) as.integer(m[2]) else 1

    if (count[i] > 0) {
      meanmf[i] <- mean(peaks$`Modified Match Factor`[matched_rows], na.rm = TRUE)
    } else {
      meanmf[i] <- 0
    }
  }

  polym_df <- data.frame(
    Polymer = base_names,
    `Markers Found` = count,
    `Markers Expected` = libcount,
    `Mean Match Factor` = meanmf,
    `Final Score` = meanmf * (count / libcount),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  polym_df[, c(4, 5)] <- round(polym_df[, c(4, 5)], 2)
  polym_df <- polym_df[order(-polym_df$`Final Score`), ]
  colnames(polym_df) <- c("Polymer", "Markers Found", "Markers Expected", "Mean Match Factor", "Final Score")

  # Write polymer list results parameter- and time-stamped inside the unique subfolder
  write.csv(polym_df, file.path(run_dir, sprintf("%s_polymerlist_%s_%s.csv", sample_name, param_suffix, timestamp_str)), row.names = FALSE)

  peaks$polymer <- sub("\\s*\\([0-9]+\\)$", "", peaks$polymer)

  return(list(polym_df, peaks, query, lib))
}
