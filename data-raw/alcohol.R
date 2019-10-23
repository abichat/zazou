## code to prepare `alcohol` dataset goes here

library(StructFDR)

data("alcohol")

names(alcohol)[4] <- "taxonomy"

usethis::use_data(alcohol)
