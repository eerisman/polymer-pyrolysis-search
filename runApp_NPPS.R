# The purpose of this code is to clear working variables and run the shiny app

rm(list=ls())

shiny::runApp('./NIST_PPS',port=7777,launch.browser="TRUE")


