library(devtools)
library(usethis)
library(testthat)


# use_build_ignore("dev_history.R")

# use_gpl3_license("Antoine Bichat")

# use_r("zscores")

# use_r("vectorized_tests")

# use_r("algebra")

# use_testthat()
# use_spell_check()

# use_test("algebra")

# use_github_links()

# use_r("optim_internal")
# use_test("unidirectionnal")

# use_r("optim_functions")
# use_test("compute_functions")

# use_readme_rmd()
# badgecreatr::badge_last_change()

# use_data_raw("alcohol")
# use_r("datasets")
# # Add Depends: R (>= 2.10) in decription

# use_test("vectorized_tests")
# use_test("zscores")

# chameleon::build_pkgdown()
# chameleon::open_pkgdown_function()

# badgecreatr::badge_packageversion()

# use_r("optim_global")

# use_r("incidence_matrix")
# use_test("incidence_matrix")

# use_vignette("simulations")

# use_r("covar")
# use_r("plot_shifts")

# use_r("shiftestim")

# use_r("tree_helpers")
# use_test("tree_helpers")

# use_r("simulations")

# use_vignette("comparisons")

# # Add biocViews: ggtree in Description

####

document()
load_all()
attachment::att_to_description()
use_tidy_description()


spell_check()
# spelling::update_wordlist()
test()

run_examples()

goodpractice::gp()
check()


# pkgdown::template_articles()
# pkgdown::template_reference()

####

install()
rmarkdown::render("README.Rmd", output_format = "md_document")
chameleon::build_pkgdown(yml = "_pkgdown.yml")
install()


