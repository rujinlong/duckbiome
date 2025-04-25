golem::fill_desc(
  pkg_name = "duckbiome",
  pkg_title = "An R package for project xxx",
  pkg_description = "This is the package and shiny app for project duckbiome.",
  authors = person(given = "Jinlong", family = "Ru", email = "jinlong.ru@gmail.com", role = c("aut", "cre")),
  repo_url = "https://github.com/rujinlong/duckbiome",
  pkg_version = "0.0.1" # The Version of the package containing the App
)
golem::set_golem_options()
usethis::use_gpl3_license()
usethis::use_readme_rmd(open = FALSE)
usethis::use_code_of_conduct(contact = "jinlong.ru@gmail.com")
usethis::use_git()

# finishing setup git
# golem::use_recommended_tests()
# golem::use_favicon("inst/app/www/favicon.png")
golem::use_utils_ui(with_test = TRUE)
golem::use_utils_server(with_test = TRUE)
rstudioapi::navigateToFile("dev/02_dev.R")
