#' Generate a sequence of bagged trees
#' @param model.formula
#' @param input.data
bagged.trees <- function(model.formula,input.data,weights,nTrees=15) {
  N <- nrow(input.data)
  fitted.trees <- list()
  nUnique <- 0
  for (nn in 1:nTrees) {
    print(nn)
    sample.ind <- sample(1:N,replace = TRUE)
    bagged.data <- input.data[sample.ind,]
    weights.loc <<- weights[sample.ind]
    fitted.trees[[nn]] <- rpart(model.formula,bagged.data,control = rpart.control(minsplit=round(N/nn),xval=0,cp=0),weights=weights.loc)
    nUnique <- nUnique + length(unique(fitted.trees[[nn]]$where))
  }
  return(list(fitted.trees=fitted.trees,nUnique=nUnique))
}

predMatrix <- function(fitted.trees,input.data,sparse=FALSE) {
  print("Using sparse:")
  print(sparse)
  if (!sparse) {
    X <- matrix(0,nrow(input.data),fitted.trees$nUnique)
  } else {
    #X <- Matrix(0,nrow(input.data),fitted.trees$nUnique,sparse=TRUE)
    x <- rep(0,nrow(input.data) * length(fitted.trees$fitted.trees))
  }
  currInd <- 1
  for (nn in 1:length(fitted.trees$fitted.trees)) {
    print(proc.time())
    print(nn)
    predicted.class <- factor(rpart.predict.leaves(fitted.trees$fitted.trees[[nn]],input.data),levels=levels(factor(fitted.trees$fitted.trees[[nn]]$where)))
    if (!sparse) {
      X[,seq(currInd,currInd+length(levels(predicted.class))-1)] <- model.matrix(~-1+.,data=data.frame(class=predicted.class))
    } else {
        x[seq((nn-1)*nrow(input.data)+1,nn*nrow(input.data))] <- currInd+as.integer(predicted.class)-1
    #  X[,seq(currInd,currInd+length(levels(predicted.class))-1)] <- model.Matrix(~-1+.,data=data.frame(class=predicted.class),sparse=TRUE)
    }
    currInd <- currInd + length(levels(predicted.class))
  }
  print(c(currInd,fitted.trees$nUnique))
  if (sparse) {
    X <- sparseMatrix(rep(1:nrow(input.data),times=length(fitted.trees$fitted.trees)),x)
    print(dim(X))
  }
  return(X)
}

predict.glmTree <- function(fitStruct,newdata){
  return(predict(fitStruct$model,newx=predMatrix(fitStruct$fitted.trees,newdata)))
}

glmTree <- function(model.formula,input.data,weights,sparse=FALSE,nTrees=15) {
  fitted.trees <- bagged.trees(model.formula,input.data,weights = weights,nTrees = nTrees)

  y <- input.data[,as.character(model.formula)[2]]
  X <- predMatrix(fitted.trees,input.data,sparse=sparse)
  model.fit <- cv.glmnet(X,y,intercept=TRUE,standardize=FALSE,alpha=.5,weights = weights)
  plot(model.fit)
  fitStruct <- list(model=model.fit,fitted.trees=fitted.trees)
  class(fitStruct) <- "glmTree"
  return(fitStruct)
}
