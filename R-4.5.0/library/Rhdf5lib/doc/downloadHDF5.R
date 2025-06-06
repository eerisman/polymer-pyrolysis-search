## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(eval = FALSE)

## ----loadLibraries------------------------------------------------------------
# library(stringr)

## -----------------------------------------------------------------------------
# hdf5_source <- tempfile()
# download.file(url = "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.7/src/hdf5-1.10.7.tar.bz2", dest = hdf5_source)
# untar(tarfile = hdf5_source, exdir = tempdir())
# system2("mv", args = c(file.path(tempdir(), "hdf5-1.10.7"), file.path(tempdir(), "hdf5")))

## -----------------------------------------------------------------------------
# hdf5_dir <- file.path(tempdir(), "hdf5")

## -----------------------------------------------------------------------------
# unlink(x = file.path(hdf5_dir, c("examples", "fortran", "java",
#              "release_docs", "test", "testpar", "tools",
#              "c++/examples", "c++/test",
#              "hl/fortran", "hl/examples", "hl/tools", "hl/test",
#              "hl/c++/examples", "hl/c++/test")),
#        recursive = TRUE)

## -----------------------------------------------------------------------------
# configure_ac <- xfun::read_utf8(file.path(hdf5_dir, "configure.ac"))
# 
# ## modify list of build files
# start <- which(str_detect(configure_ac, pattern = "AC_CONFIG_FILES"))
# end <- which(str_detect(configure_ac[start:(length(configure_ac))], pattern = "\\)$"))[1] + start - 1
# config_files <- configure_ac[start:end]
# rm_idx <- which(str_detect(config_files, pattern = "test/|testpar/|tools/|examples/|fortran/|java/|h5c++/"))
# config_files <- config_files[-rm_idx]
# config_files[length(config_files)] <- paste0(tail(config_files, 1), "])")
# configure_ac[start] <- paste(config_files, collapse = "\n")
# configure_ac <- configure_ac[-((start+1):(end))]
# 
# ## remove reference to h5cc
# h5cc <- str_which(configure_ac, pattern = "chmod 755 [a-z/]*/h5cc")
# configure_ac <- configure_ac[-((h5cc):(h5cc+4))]
# 
# ## fortran headers
# fortran_inc <- str_which(configure_ac, pattern = "AC_CONFIG_HEADERS\\(\\[fortran/src/H5config_f\\.inc")
# configure_ac[fortran_inc:(fortran_inc+1)] <- paste("##", configure_ac[fortran_inc:(fortran_inc+1)])
# 
# ## write
# xfun::write_utf8(configure_ac, con = file.path(hdf5_dir, "configure.ac"))
# 
# ## C++ makefile
# make_cplusplus <- xfun::read_utf8(file.path(hdf5_dir, 'c++/Makefile.am'))
# idx <- str_which(make_cplusplus, "BUILD_CXX_CONDITIONAL")
# make_cplusplus[idx] <- "if BUILD_CXX_CONDITIONAL\n   SUBDIRS=src\nendif\nDIST_SUBDIRS = src"
# make_cplusplus <- make_cplusplus[-((idx+1):(length(make_cplusplus)-2))]
# xfun::write_utf8(make_cplusplus, con = file.path(hdf5_dir, "c++/Makefile.am"))
# 
# ## HL makefile
# make_hl <- xfun::read_utf8(file.path(hdf5_dir, 'hl/Makefile.am'))
# idx <- str_which(make_hl, "BUILD_HDF5_HL_CONDITIONAL")
# make_hl[idx] <- "if BUILD_HDF5_HL_CONDITIONAL\n   SUBDIRS=src $(CXX_DIR)\nendif\nDIST_SUBDIRS = src c++"
# make_hl <- make_hl[-((idx+1):(length(make_hl)-2))]
# xfun::write_utf8(make_hl, con = file.path(hdf5_dir, "hl/Makefile.am"))
# 
# ## HL C++ makefile
# make_hl_cpp <- xfun::read_utf8(file.path(hdf5_dir, 'hl/c++/Makefile.am'))
# idx <- str_which(make_hl_cpp, "^SUBDIRS=src")
# make_hl_cpp[idx] <- "SUBDIRS=src\nDIST_SUBDIRS=src"
# make_hl_cpp <- make_hl_cpp[-((idx+1):(length(make_hl_cpp)-2))]
# xfun::write_utf8(make_hl_cpp, con = file.path(hdf5_dir, "hl/c++/Makefile.am"))
# 
# ## Primary makefile
# make <- xfun::read_utf8(file.path(hdf5_dir, 'Makefile.am'))
# idx <- str_which(make, "SUBDIRS = src")[1]
# make[idx] <- "SUBDIRS = src . $(CXX_DIR) $(HDF5_HL_DIR)"
# make[idx+1] <- "DIST_SUBDIRS = src . c++ hl"
# make[idx+2]  <- ""
# idx <- str_which(make, "# Make all, tests, and \\(un\\)install")
# make[(idx+1):(idx+6)] <- paste0("##", make[(idx+1):(idx+6)])
# xfun::write_utf8(make, con = file.path(hdf5_dir, "Makefile.am"))

## -----------------------------------------------------------------------------
# code <- xfun::read_utf8(file.path(hdf5_dir, 'c++', 'src', 'H5Library.cpp'))
# code <- str_replace(code, '([ ]{1,})(exit\\()', replacement = '\\1std::\\2' )
# xfun::write_utf8(code, con = file.path(hdf5_dir, 'c++', 'src', 'H5Library.cpp'))

## -----------------------------------------------------------------------------
# system(command = paste0("cd ", hdf5_dir, " && autoconf"))
# system(command = paste0("cd ", hdf5_dir, " && aclocal"))
# system(command = paste0("cd ", hdf5_dir, " && automake"))
# unlink(file.path(hdf5_dir, "autom4te.cache"), recursive = TRUE)

## -----------------------------------------------------------------------------
# szip_source <- tempfile()
# download.file(url = "https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz",
#               dest = szip_source)
# untar(tarfile = szip_source, exdir = tempdir())
# system2("mv", args = c(file.path(tempdir(), "szip-2.1.1"), file.path(tempdir(), "szip")))

## -----------------------------------------------------------------------------
# szip_dir <- file.path(tempdir(), "szip")

## -----------------------------------------------------------------------------
# unlink(x = file.path(szip_dir, "test"),
#        recursive = TRUE)

## -----------------------------------------------------------------------------
# xfun::read_utf8(file.path(szip_dir, "configure.ac")) %>%
#   str_remove("test/Makefile") %>%
#   xfun::write_utf8(file.path(szip_dir, "configure.ac"))
# 
# xfun::read_utf8(file.path(szip_dir, "Makefile.am")) %>%
#   str_replace(pattern = "SUBDIRS=src test", replacement = "SUBDIRS=src") %>%
#   xfun::write_utf8(file.path(szip_dir, "Makefile.am"))

## -----------------------------------------------------------------------------
# system(command = paste0("cd ", szip_dir, " && autoconf"))
# system(command = paste0("cd ", szip_dir, " && aclocal"))
# system(command = paste0("cd ", szip_dir, " && automake"))
# unlink(file.path(szip_dir, "autom4te.cache"), recursive = TRUE)

## ----createTarball------------------------------------------------------------
# system2("mv", args = c(szip_dir, file.path(hdf5_dir, "szip")))
# system2("tar", args = c("-C", tempdir(), "-czf", file.path(tempdir(), "hdf5small_cxx_hl_1.10.7.tar.gz"), "hdf5"))

## -----------------------------------------------------------------------------
# if(file.exists("/tmp/hdf5small_cxx_hl_1.10.7.tar.gz")) { file.remove("/tmp/hdf5small_cxx_hl_1.10.7.tar.gz") }
# file.copy(file.path(tempdir(), "hdf5small_cxx_hl_1.10.7.tar.gz"), to = "/tmp/")

