# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
library(roxygen2)
library(devtools)

#' Helper method for the main method.
#'
#' Tests to see if the input to the main function is valid.
#' Valid input is defined as 2-5 data frames, each with:
#' \itemize{
#'     \item The same number of columns.
#'     \item At least 2 unique values in the first column.
#' }
#'
#' @param input A list() of the data frames passed as input to the main function.
#' @return True if the input is valid, as defined above; false otherwise.
isInputValid = function(input){
  ## TO-DO: Method Stub
}
