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
# use_test("model_selection")

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

# use_vignette("alcohol")

# use_data_raw("chlamydiae")

# use_vignette("chlamydiae")

# use_r("bic")

# use_test("shiftestim")

# # Change some test names

# use_r("extraction")

# use_r("clusters")
# use_test("clusters")

# use_r("optim_scaled")

# use_r("optim_desparsified")

# use_test("desparsified_internal")
# use_test("desparsified_midlevel")

# use_vignette("constraints")

# use_test("shiftestim")


####

document()
load_all()
attachment::att_to_description()
use_tidy_description()


spell_check()
# spelling::update_wordlist()
test()

run_examples()

covr::report()

check()
check(args = "--no-build-vignettes", build_args = "--no-build-vignettes")
goodpractice::gp()


# pkgdown::template_articles()
# pkgdown::template_reference()

####

install(upgrade = "never")
rmarkdown::render("README.Rmd")
file.remove("README.html")
chameleon::build_pkgdown(yml = "_pkgdown.yml")
install(upgrade = "never")


