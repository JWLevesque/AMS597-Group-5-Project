## ---- include = FALSE----------------------------------------------------
#library(devtools)
devtools::load_all()
library(Group5Project)

## ------------------------------------------------------------------------
set.seed(123)
p<-100
data1 <- data.frame(group=sample(1:3,200,replace=TRUE),matrix(rnorm(p*200),ncol=p))
data2 <- data.frame(group=sample(1:2,150,replace=TRUE),matrix(rnorm(p*150),ncol=p))
pooledPValues(method = "Fisher",data1,data2)

## ------------------------------------------------------------------------
pooledPValues(method = "Stouffer",data1,data2)

## ------------------------------------------------------------------------
pooledPValues(method = "minP",data1,data2)

## ------------------------------------------------------------------------
pooledPValues(method = "maxP",data1,data2)

## ------------------------------------------------------------------------
set.seed(123)
p<-100
data1 <- data.frame(group=sample(1:3,200,replace=TRUE),matrix(rnorm(p*200),ncol=p))
getFramePvals(data1)

## ------------------------------------------------------------------------
set.seed(123)
p<-100
data1 <- data.frame(group=sample(1:3,200,replace=TRUE),matrix(rnorm(p*200),ncol=p))
data2 <- data.frame(group=sample(1:2,150,replace=TRUE),matrix(rnorm(p*150),ncol=p))
isInputValid(list(data1,data2))

## ------------------------------------------------------------------------
set.seed(123)
p<-100
data1 <- data.frame(group=sample(1:3,200,replace=TRUE),matrix(rnorm(p*200),ncol=p))
a<-list(data1,data1,data1,data1,data1,data1)
isInputValid(a)

## ------------------------------------------------------------------------
isNormByGroup(data1,3)

## ------------------------------------------------------------------------
set.seed(123)
p<-100
data3 <- data.frame(group=sample(1:3,200,replace=TRUE),matrix(rexp(p*200),ncol=p))
isNormByGroup(data3,4)

