if(F){
  setwd("~/Documents/Github/blackmarbler")
  
  roxygen2::roxygenise("~/Documents/Github/blackmarbler")
  
  pkgdown::clean_site()
  
  usethis::use_pkgdown()
  pkgdown::build_site()
  pkgdown::deploy_to_branch()
  
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