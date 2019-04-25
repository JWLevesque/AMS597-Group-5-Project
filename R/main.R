#library(roxygen2)
#library(devtools)

#' Pools the p-values between data frames using the specified method.
#'
#' @param method The desired pooling method: "Fisher" for Fisher, "Stouffer" for Stouffer, "minP" for minimum p-value, or "maxP" for maximum p-value.
#' @param ... The data frames to pool the p-values from.  Must be between 2 and 5 data frames, all with the same number of columns, and at least 2 groups per data frame.
#' @return A vector containing pooled p-values for each biomarker.  For p biomarkers, the vector will be of length p.
#' @export
pooledPValues = function(method = "Fisher", ...){
  # Place the optional data frame arguments into a list
  frames = list(...)
  # Check the data frames for validity
  if(!isInputValid(frames)){
    # If invalid, thow an exception and stop execution
    stop("Input is invalid: invalid data frames.  Please see documentation.")
  }
  # Generate a list (by data frame) of vectors (by biomarker in data frame) of p-values for group difference in each biomarker.
  #     e.g., p-value for group difference in biomarker j from data frame i is given by unpooledPvalsList[[i]][j].
  unpooledPvalsList = list()
  for(i in 1:length(frames)){
    unpooledPvalsList[[i]] = getFramePvals(frames[[i]])
  }
  # Pool p-values according to the
  tbr = rep(NA, length(unpooledPvalsList[[1]]))
  if(method == "Stouffer"){
    # use Stouffer
    tbr = poolStouffer(unpooledPvalsList, frames)
  } else if(method == "minP"){
    # use minP
    tbr = poolMinP(unpooledPvalsList, frames)
  } else if(method == "maxP"){
    # use maxP
    tbr = poolMaxP(unpooledPvalsList, frames)
  } else if(method == "Fisher"){
    # use Fisher
    tbr = poolFisher(unpooledPvalsList, frames)
  } else {
    # Invalid input for method argument; throw exception
    stop("Invalid input for argument 'method'. Please see documentation.")
  }
  return(tbr)
}

#' Pools the p-values using the Stouffer method.
#'
#' @param unpooledPvalsList A list (by data frame) of vectors (by biomarker in data frame) of p-values for group difference in each biomarker.
#'                          E.g., p-value for group difference in biomarker j from data frame i is given by \code{unpooledPvalsList[[i]][j]}.
#' @param frameList A list of the original 2-5 data frames.  Note that data frame i is given by \code{frameList[[i]]}.
#' @return A vector containing pooled p-values for each biomarker.  For p biomarkers, the vector will be of length p.
#' @importFrom stats qnorm pnorm
poolStouffer = function(unpooledPvalsList, frameList){
  k<-length(frameList)
  Tstouffer<-rep(NA,length(unpooledPvalsList[[1]]))
  pooledPValues<-rep(NA,length(unpooledPvalsList[[1]]))
  for (i in 1: length(unpooledPvalsList[[1]])){
    Tstouffer[i] <-0
    for(j in 1:length(frameList)){
      Tstouffer[i]<-Tstouffer[i]+qnorm(unpooledPvalsList[[j]][i])/sqrt(k)
    }
    pooledPValues[i]<-pnorm(Tstouffer[i])
  }
  return(pooledPValues)
}

#' Pools the p-values using the Minimum p-value method.
#'
#' @param unpooledPvalsList A list (by data frame) of vectors (by biomarker in data frame) of p-values for group difference in each biomarker.
#'                          E.g., p-value for group difference in biomarker j from data frame i is given by \code{unpooledPvalsList[[i]][j]}.
#' @param frameList A list of the original 2-5 data frames.  Note that data frame i is given by \code{frameList[[i]]}.
#' @return A vector containing pooled p-values for each biomarker.  For p biomarkers, the vector will be of length p.
poolMinP = function(unpooledPvalsList, frameList){
  k<-length(frameList)
  pooledPValues<-rep(NA,length(unpooledPvalsList[[1]]))
  for (i in 1: length(unpooledPvalsList[[1]])){
    mincandidate<-rep(NA,length(frameList))
    for(j in 1:length(frameList)){
      mincandidate[j]<-unpooledPvalsList[[j]][i]
    }
    pooledPValues[i]<-min(mincandidate)
  }
  return(pooledPValues)
}

#' Pools the p-values using the Maximum p-value method.
#'
#' @param unpooledPvalsList A list (by data frame) of vectors (by biomarker in data frame) of p-values for group difference in each biomarker.
#'                          E.g., p-value for group difference in biomarker j from data frame i is given by \code{unpooledPvalsList[[i]][j]}.
#' @param frameList A list of the original 2-5 data frames.  Note that data frame i is given by \code{frameList[[i]]}.
#' @return A vector containing pooled p-values for each biomarker.  For p biomarkers, the vector will be of length p.
poolMaxP = function(unpooledPvalsList, frameList){
  k<-length(frameList)
  pooledPValues<-rep(NA,length(unpooledPvalsList[[1]]))
  for (i in 1: length(unpooledPvalsList[[1]])){
    maxcandidate<-rep(NA,length(frameList))
    for(j in 1:length(frameList)){
      maxcandidate[j]<-unpooledPvalsList[[j]][i]
    }
    pooledPValues[i]<-max(maxcandidate)
  }
  return(pooledPValues)
}

#' Pools the p-values using the Fisher method.
#'
#' @param unpooledPvalsList A list (by data frame) of vectors (by biomarker in data frame) of p-values for group difference in each biomarker.
#'                          E.g., p-value for group difference in biomarker j from data frame i is given by \code{unpooledPvalsList[[i]][j]}.
#' @param frameList A list of the original 2-5 data frames.  Note that data frame i is given by \code{frameList[[i]]}.
#' @return A vector containing pooled p-values for each biomarker.  For p biomarkers, the vector will be of length p.
#' @importFrom stats pchisq
poolFisher = function(unpooledPvalsList, frameList){
  k<-length(frameList)
  chisquare<-rep(NA,length(unpooledPvalsList[[1]]))
  pooledPValues<-rep(NA,length(unpooledPvalsList[[1]]))
  for (i in 1: length(unpooledPvalsList[[1]])){
    chisquare[i] <-0
    for(j in 1:length(frameList)){
      chisquare[i]<-chisquare[i]-2*log(unpooledPvalsList[[j]][i])
    }
    pooledPValues[i]<-pchisq(chisquare[i],lower.tail = FALSE,df=2*k)
  }
  return(pooledPValues)
}

#' Helper method which populates a vector with p-values for group difference in each biomarker from a single frame.
#'
#' This function tests for normality within groups, then applies the appropriate test for group difference.
#'
#' @param dataFrame The data frame to generate p-values from.
#' @return A vector containing the p-values for group difference for each biomarker in the specified data frame.  Note that if there are p biomarkers, this vector will be of length p.
#' @importFrom stats var.test t.test bartlett.test aov summary.aov oneway.test wilcox.test kruskal.test
getFramePvals = function(dataFrame){
  p = (length(dataFrame[1,]) - 1) # p = num_biomarkers
  tbr = rep(NA, p) # vector to populate with p-values
  numGroups = length(factor(dataFrame$group))
  groupNames = levels(factor(dataFrame$group))
  groupIndicesList = list() # To be populated with vectors of indices for each group.  groupIndicesList[[i]]= vector of indices of group i
  for(i in 1:length(groupNames)){
    groupIndicesList[[i]] = which(dataFrame$group == groupNames[i])
    #group_i = dataFrame[group_i_indices, marker+1]
  }
  # For each biomarker i in the data frame:  (Note that biomarker i corresponds to column i+1 in the data frame)
  for(i in 1:p){
    # Check for normality within groups
    isNorm = isNormByGroup(dataFrame,i)
    # If normal, then do either one-way ANOVA or two sample t-test according to quantity of groups
    if(isNorm){ # Normal between groups
      if(numGroups == 2){
        ### Do two-sample t-test
        groupOne = dataFrame[groupIndicesList[[1]], i+1]
        groupTwo = dataFrame[groupIndicesList[[2]], i+1]
        isHomoskedastic = ((var.test(groupOne,groupTwo)$p.value) >= 0.05) # Cannot reject equal variance if p-value >= alpha=0.05
        tbr[i] = (t.test(groupOne, groupTwo, var.equal = isHomoskedastic)$p.value) # add the p-value to the vector tbr
      } else {
        ### Do one-way ANOVA
        # Check equal variance assumption of ANOVA.  Since the data is normal by if(isNorm), Bartlett's test can be used instead of Levene:
        isHomoskedastic = ((bartlett.test(x = dataFrame[,i+1], g = dataFrame$group)$p.value) >= 0.05) # Cannot reject equal variance if p-value >= alpha=0.05
        # If equal variance assumption is plausible at alpha level 0.05, do ANOVA
        if(isHomoskedastic){
          aovRes = aov(dataFrame[,i+1]~factor(dataFrame$group))
          aovSummary = summary.aov(aovRes)
          tbr[i] = aovSummary[[1]][["Pr(>F)"]][[1]]
        }
        # If we reject equal variance between groups, use Welch ANOVA
        else {
          tbr[i] = (oneway.test(dataFrame[,i+1]~factor(dataFrame$group))$p.value)
        }
      }
    }
    # If not normal, then do either Kruskal Wallis test or the Wilcoxon rank sum test according to quantity of groups
    else { # Not normal between groups
      if(numGroups == 2){
        ### Do Wilcoxon rank sum test
        groupOne = dataFrame[groupIndicesList[[1]], i+1]
        groupTwo = dataFrame[groupIndicesList[[2]], i+1]
        tbr[i] = (wilcox.test(groupOne, groupTwo)$p.value) # add the p-value to the vector tbr
      } else {
        ### Do Kruskal Wallis test
        tbr[i] = (kruskal.test(x = dataFrame[,i+1], g = dataFrame$group)$p.value)
      }
    }
  }
  return(tbr)
}

#' Helper method which determines if a biomarker is normally distributed by group at alpha level 0.05.
#'
#' This function will be called by the implementations of bullet 1 and 2 on page 7 of the project slides.
#' Note that the \code{shapiro.test()} function limits sample size to \code{5000}.
#'
#' @param dataFrame The data frame to test for normality by group.
#' @param marker The biomarker number to be tested; e.g., for the nth biomarker, input integer n.
#'                   Note that marker n corresponds to column n+1 in the data frame.
#' @return True if normality cannot be rejected for all groups; false otherwise.
#' @importFrom stats shapiro.test
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
  tbr = TRUE
  # Check to see if there are between 2 and 5 data frames
  if((length(input)<2) || (length(input)>5)){
    tbr = FALSE #return(FALSE)
  }
  # Iterate through the 2-5 data frames
  numCols = length(input[[1]][1,])
  for(i in 1:length(input)){
    # Check that all data frames have the same number of columns
    if(length(input[[i]][1,]) != numCols){
      tbr = FALSE #return(FALSE)
      break
    }
    # Check that there are at least 2 unique values in the first column of the current data frame
    uniqueVals = length(factor(input[[i]][,1]))
    if(uniqueVals < 2){
      tbr = FALSE #return(FALSE)
      break
    }
  }
  # If this code is reached, there are between 2 and 5 data frames, all with the same number of columns,
  #     and at least 2 unique values in the first column.
  return(tbr) #return(TRUE)
}
