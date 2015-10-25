#' Generate a sequence of bagged trees
#' @param model.formula
#' @param input.data
bagged.trees <- function(model.formula,input.data,nTrees,log.base=2,bagged.trees = list(fitted.trees = list(),nUnique=0,nRuns=0)) {
  print(c(nTrees,log.base))
  N <- nrow(input.data)
  print("Fitting trees:")
  pb <- txtProgressBar(max=nTrees-1,style=3)
  for (nn in 1:nTrees) {
    setTxtProgressBar(pb,nn-1)
    #sample.ind <- sample(1:N,replace = FALSE)
    #bagged.data <- input.data[sample.ind,]
    rpart.contr <- rpart.control(minsplit=round(N/log.base^nn),xval=0,cp=0,maxsurrogate=6)
    #rpart.contr <- rpart.control(maxdepth=nn,xval=0,cp=0)
    curr.tree <- rpart(model.formula,input.data,control = rpart.contr,weights=w,cost=runif(ncol(input.data)-2))
    bagged.trees$fitted.trees[[length(bagged.trees$fitted.trees)+1]] <- curr.tree
    #bagged.trees$nUnique <- bagged.trees$nUnique + length(unique(curr.tree$where))
    bagged.trees$nUnique <- bagged.trees$nUnique + nrow(curr.tree$frame)
  }
  close(pb)
  bagged.trees$nRuns <- bagged.trees$nRuns+1
  return(bagged.trees)
}

predMatrix <- function(bagged.trees,input.data,sparse=TRUE,isNewData=FALSE) {
  print("Using sparse:")
  print(sparse)

  x <- rep(0,nrow(input.data) * length(bagged.trees$fitted.trees))
  currInd <- 1
  pb <- txtProgressBar(max=length(bagged.trees$fitted.trees)-1,style=3)
  for (nn in 1:length(bagged.trees$fitted.trees)) {
    setTxtProgressBar(pb,nn-1)
    if (isNewData) {
      #predicted.class <- factor(rpart.predict.leaves(bagged.trees$fitted.trees[[nn]],input.data),levels=levels(factor(bagged.trees$fitted.trees[[nn]]$where)))
      predicted.class <- factor(rpart.predict.leaves(bagged.trees$fitted.trees[[nn]],input.data),levels=1:nrow(bagged.trees$fitted.trees[[nn]]$frame))
    }  else {
      predicted.class <- factor(bagged.trees$fitted.trees[[nn]]$where)
    }
    if (sum(is.na(predicted.class))>0) { browser()}
    x[seq((nn-1)*nrow(input.data)+1,nn*nrow(input.data))] <- currInd+as.integer(predicted.class)-1
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

glmTree <- function(model.formula,input.data,weights,sparse=TRUE,log.base=1.5,nTrees=floor(logb(nrow(input.data),log.base)-logb(20,log.base)),seed=1,alpha=0,fitStruct=list(nRuns=0),pred.data,s="lambda.min") {
  set.seed(seed)
  cat("Doing iteration:",fitStruct$nRuns+1)
  glmnet.control(eps=1e-9)

  struct.int <- sample(1:nrow(input.data),size=round(nrow(input.data)/4),replace = FALSE)
  bagged.trees <- bagged.trees(model.formula,input.data[struct.int,],nTrees = nTrees,log.base=log.base)

  y <- input.data[-struct.int,as.character(model.formula)[2]]
  X <- predMatrix(bagged.trees,input.data[-struct.int,],sparse=TRUE,isNewData = TRUE)
  X.new <- predMatrix(bagged.trees,pred.data,sparse=TRUE,isNewData = TRUE)

  if (!is.null(fitStruct$X)) {
    X <- cBind(fitStruct$X,X)
    X.new <- cBind(fitStruct$X.new,X.new)
  }
  model.fit <- cv.glmnet(X,y,intercept=TRUE,standardize=FALSE,alpha=alpha,weights = input.data$w[-struct.int],lambda.min.ratio=1e-9,thresh=1e-6,nlambda=20)
  plot(model.fit)

  # creating predictions
  fitStruct <- list(input.predict = predict(model.fit,newx=X,s=s)[,1], new.predict = predict(model.fit,newx=X.new,s=s)[,1],nRuns = fitStruct$nRuns+1)
  class(fitStruct) <- "glmTree"

  return(fitStruct)
}
