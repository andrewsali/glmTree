#' Generate a sequence of bagged trees
#' @param model.formula
#' @param input.data
bagged.trees <- function(model.formula,input.data,weights,nTrees,log.base=2,bagged.trees = list(fitted.trees = list(),nUnique=0,nRuns=0)) {
  print(c(nTrees,log.base))
  N <- nrow(input.data)
  print("Fitting trees:")
  pb <- txtProgressBar(max=nTrees-1,style=3)
  for (nn in 1:nTrees) {
    setTxtProgressBar(pb,nn-1)
    sample.ind <- sample(1:N,replace = TRUE)
    bagged.data <- input.data[sample.ind,]
    weights.loc <<- weights[sample.ind]
    #rpart.contr <- rpart.control(minsplit=round(N/log.base^nn),xval=0,cp=0)
    rpart.contr <- rpart.control(maxdepth=nn,xval=0,cp=0)
    curr.tree <- rpart(model.formula,bagged.data,control = rpart.contr,weights=weights.loc,cost=runif(ncol(bagged.data)-1)^0)
    bagged.trees$fitted.trees[[length(bagged.trees$fitted.trees)+1]] <- curr.tree
    bagged.trees$nUnique <- bagged.trees$nUnique + length(unique(curr.tree$where))
  }
  close(pb)
  bagged.trees$nRuns <- bagged.trees$nRuns+1
  return(bagged.trees)
}

predMatrix <- function(bagged.trees,input.data,sparse=TRUE) {
  print("Using sparse:")
  print(sparse)
  if (!sparse) {
    X <- matrix(0,nrow(input.data),bagged.trees$nUnique)
  } else {
    #X <- Matrix(0,nrow(input.data),fitted.trees$nUnique,sparse=TRUE)
    x <- rep(0,nrow(input.data) * length(bagged.trees$fitted.trees))
  }
  currInd <- 1
  pb <- txtProgressBar(max=length(bagged.trees$fitted.trees)-1,style=3)
  for (nn in 1:length(bagged.trees$fitted.trees)) {
    setTxtProgressBar(pb,nn-1)
    predicted.class <- factor(rpart.predict.leaves(bagged.trees$fitted.trees[[nn]],input.data),levels=levels(factor(bagged.trees$fitted.trees[[nn]]$where)))
    if (!sparse) {
      X[,seq(currInd,currInd+length(levels(predicted.class))-1)] <- model.matrix(~-1+.,data=data.frame(class=predicted.class))
    } else {
        x[seq((nn-1)*nrow(input.data)+1,nn*nrow(input.data))] <- currInd+as.integer(predicted.class)-1
    #  X[,seq(currInd,currInd+length(levels(predicted.class))-1)] <- model.Matrix(~-1+.,data=data.frame(class=predicted.class),sparse=TRUE)
    }
    currInd <- currInd + length(levels(predicted.class))
  }
  close(pb)
  print(c(currInd,bagged.trees$nUnique))
  if (sparse) {
    X <- sparseMatrix(rep(1:nrow(input.data),times=length(bagged.trees$fitted.trees)),x,dims=c(nrow(input.data),bagged.trees$nUnique))
    print(dim(X))
  }
  return(X)
}

predict.glmTree <- function(fitStruct,newdata,s="lambda.1se"){
  return(predict(fitStruct$model,newx=predMatrix(fitStruct$bagged.trees,newdata,sparse=TRUE),s=s))
}

glmTree <- function(model.formula,input.data,weights,sparse=TRUE,log.base=1.5,nTrees=floor(logb(nrow(input.data),log.base)-logb(20,log.base)),seed=1,alpha=1,fitStruct=NULL,is.struct) {
  set.seed(seed)
  glmnet.control(eps=1e-9)

  if (!is.null(fitStruct)) {
    input.newdata <- input.data
    input.newdata$Sales <- input.newdata$Sales - predict(fitStruct,input.newdata,"lambda.min")[,1]*fitStruct$bagged.trees$nRuns/(fitStruct$bagged.trees$nRuns+1)
    bagged.trees <- bagged.trees(model.formula,input.newdata,weights = weights,nTrees = nTrees,log.base=log.base,fitStruct$bagged.trees)
  } else {
    bagged.trees <- bagged.trees(model.formula,input.data,weights = weights,nTrees = nTrees,log.base=log.base)
  }

  y <- input.data[,as.character(model.formula)[2]]
  X <- predMatrix(bagged.trees,input.data,sparse=sparse)

  model.fit <- cv.glmnet(X,y,intercept=TRUE,standardize=FALSE,alpha=alpha,weights = weights,lambda.min.ratio=1e-9)
  plot(model.fit)
  fitStruct <- list(model=model.fit,bagged.trees=bagged.trees)
  class(fitStruct) <- "glmTree"
  return(fitStruct)
}
