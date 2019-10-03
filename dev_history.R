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



####

load_all()
document()
attachment::att_to_description()
use_tidy_description()


spell_check()
# spelling::update_wordlist()
test()

goodpractice::gp()
check()

