if(F){
  
  ## Change the oken
  #Sys.getenv("BEARER_NASA_TOKEN")
  #Sys.setenv(BEARER_NASA_TOKEN = "new_token_value")
  #Sys.setenv(NASA_USERNAME = "value")
  #Sys.setenv(NASA_PASSWORD = "value")

  setwd("~/Documents/Github/blackmarbler")
  
  roxygen2::roxygenise("~/Documents/Github/blackmarbler")
  
  pkgdown::build_favicons(pkg = "~/Documents/Github/blackmarbler", 
                          overwrite = FALSE)
  
  pkgdown::clean_site()
  
  usethis::use_pkgdown()
  pkgdown::build_site()
  pkgdown::build_site_github_pages()
  #usethis::use_pkgdown_github_pages()
  
  usethis::use_testthat(3)
  
  usethis::use_pkgdown()
  usethis::use_github_pages()
  usethis::use_pkgdown_github_pages() #####
  usethis::use_github_action_check_standard()
  
  ## Comand line code for building and checking package
  #R CMD build --as-cran "~/Documents/Github/googletraffic"
  #R CMD check --as-cran "~/Documents/Github/googletraffic/googletraffic_0.0.0.9000.tar.gz"
  
  devtools::check("~/Documents/Github/blackmarbler")
  
  devtools::check_win_devel("~/Documents/Github/blackmarbler")
  devtools::check_win_release("~/Documents/Github/blackmarbler")
  devtools::check_win_oldrelease("~/Documents/Github/blackmarbler")
  
  devtools::build("~/Documents/Github/blackmarbler")
}