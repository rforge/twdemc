if(require("testthat", quietly=TRUE)) {
    # https://github.com/hadley/testthat/issues/86
    Sys.setenv("R_TESTS" = "")  
    test_package("blockDemac")
} else {
    warning("cannot run unit tests -- package testthat is not available")
}

