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


####

load_all()
document()
attachment::att_to_description()
use_tidy_description()

test()

goodpractice::gp()
check()

