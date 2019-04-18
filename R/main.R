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

#' Helper method which determines if a biomarker is normally distributed by group at alpha level 0.05.
#'
#' This function will be called by the implementations of bullet 1 and 2 on page 7 of the project slides.
#' Note that the shapiro.test() function limits sample size to 5000.
#'
#' @param dataFrame The data frame to test for normality by group.
#' @param marker The biomarker number to be tested; e.g., for the nth biomarker, input integer n.
#'                   Note that marker n corresponds to column n+1 in the data frame.
#' @return True if normality cannot be rejected for all groups; false otherwise.
isNormByGroup = function(dataFrame, marker){
  tbr = TRUE
  groupNames = levels(factor(dataFrame$group))
  # For each group, test for normality
  for(i in 1:length(groupNames)){
    group_i_indices = which(dataFrame$group == groupNames[i])
    group_i = dataFrame[group_i_indices, marker+1]
    # Test the values for biomarker n in group i for normality
    if(shapiro.test(group_i)$p.value < 0.05){
      # If p-value is less than alpha = 0.05, reject normality
      tbr = FALSE
      # No need to continue the loop, since one group isnt normal--break
      break
    }
  }
  return(tbr)
}

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
