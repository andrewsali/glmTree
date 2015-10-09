#' Generate a sequence of bagged trees
#' @param model.formula
#' @param input.data
bagged.trees <- function(model.formula,input.data) {
  N <- nrow(input.data)
  fitted.trees <- list()
  nUnique <- 0
  for (nn in 1:25) {
    print(nn)
    bagged.data <- input.data[sample(1:N,replace = TRUE),]
    fitted.trees[[nn]] <- rpart(model.formula,bagged.data,control = rpart.control(minsplit=round(N/nn),xval=0,cp=0),weights=weights)
    nUnique <- nUnique + length(unique(fitted.trees[[nn]]$where))
  }
  return(list(fitted.trees=fitted.trees,nUnique=nUnique))
}

predMatrix <- function(fitted.trees,input.data) {
  X <- matrix(0,nrow(input.data),fitted.trees$nUnique)
  currInd <- 1
  for (nn in 1:length(fitted.trees$fitted.trees)) {
    print(nn)
    predicted.class <- factor(rpart.predict.leaves(fitted.trees$fitted.trees[[nn]],input.data))
    X[,seq(currInd,currInd+length(levels(predicted.class))-1)] <- model.matrix(~-1+.,data=data.frame(class=predicted.class))
    currInd <- currInd + length(levels(predicted.class))
  }
  print(c(currInd,fitted.trees$nUnique))
  return(X)
}

predict.glmTree <- function(fitStruct,newdata){
  return(predict(fitStruct$model,newx=predMatrix(fitStruct$fitted.trees,newdata)))
}

glmTree <- function(model.formula,input.data,weights,sparse=FALSE) {
  fitted.trees <- bagged.trees(model.formula,input.data)

  y <- input.data[,as.character(model.formula)[2]]
  X <- predMatrix(fitted.trees,input.data)

  fitStruct <- list(model=cv.glmnet(X,y,intercept=TRUE,standardize=FALSE,alpha=.5,weights = weights),fitted.trees=fitted.trees)
  class(fitStruct) <- "glmTree"
  return(fitStruct)
}
