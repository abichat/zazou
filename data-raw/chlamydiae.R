## code to prepare `chlamydia` dataset goes here

library(structSSI)

data("chlamydiae")

X <-
  chlamydiae@otu_table %>%
  as.data.frame() %>%
  as.matrix()

Y <- chlamydiae@sam_data$SampleType

tree <- chlamydiae@phy_tree

taxonomy <-
  chlamydiae@tax_table %>%
  as.data.frame() %>%
  as.matrix()

chlamydiae <- list(X = X, Y = Y, tree = tree, taxonomy = taxonomy)

usethis::use_data(chlamydiae)

