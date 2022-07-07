

# TODO --------------------------------------------------------------------


#TODO reverse copula direction
#TODO  - spearman's rho




# CRAN release ------------------------------------------------------------

#preparation
devtools::spell_check() 
# spelling::update_wordlist()

#check on different OS
devtools::check_rhub(interactive=FALSE, env_vars=c(R_COMPILE_AND_INSTALL_PACKAGES = "always"))
devtools::check_win_devel(quiet=TRUE)
# devtools::check_win_release()
# devtools::check_win_oldrelease()
beepr::beep()

codemetar::write_codemeta()


#release
devtools::release() #same with no questions = devtools::submit_cran()
#TODO accepter mail CRAN après release
#TODO update readme (commit, CRAN version)

beepr::beep()
