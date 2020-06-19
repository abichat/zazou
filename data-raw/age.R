## code to prepare `age` dataset goes here

library(curatedMetagenomicData)
library(tidyverse)

brito <-
  "BritoIL_2016.metaphlan_bugs_list.stool" %>%
  curatedMetagenomicData(dryrun = FALSE, counts = TRUE) %>%
  mergeData() %>%
  ExpressionSet2phyloseq(phylogenetictree = TRUE)

X <-
  brito@otu_table %>%
  as.data.frame() %>%
  rename_all(str_remove, "BritoIL_2016.metaphlan_bugs_list.stool:") %>%
  rename_all(str_remove, ".ST") %>%
  as.matrix()
rownames(X) <- str_remove(rownames(X), "s__")
X

Y <-
  brito@sam_data %>%
  unclass() %>%
  as.data.frame() %>%
  mutate(age_category = as.character(age_category)) %>%
  pull(age_category, name = subjectID)
Y

taxonomy <-
  brito@tax_table %>%
  unclass() %>%
  as_tibble(rownames = "taxonomy") %>%
  mutate(Strain = NULL,
         taxonomy = str_remove(taxonomy, "s__")) %>%
  column_to_rownames("taxonomy") %>%
  as.matrix()
taxonomy

tree <- brito@phy_tree
tree$node.label <- NULL
tree$tip.label <- str_remove(tree$tip.label, "s__")
tree

all(colnames(X) == names(Y))
all(rownames(X) == rownames(taxonomy))
all(rownames(X) == tree$tip.label)

age <- list(X = X, Y = Y, tree = tree, taxonomy = taxonomy)

usethis::use_data(age, overwrite = TRUE)
