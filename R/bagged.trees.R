#' Generate a sequence of bagged trees
#' @param model.formula
#' @param input.data
bagged.trees <- function(model.formula,input.data) {
  N <- nrow(input.data)
  fitted.trees <- list()

  for (nn in 1:round(sqrt(nrow(input.data)))) {
    bagged.data <- input.data[sample(1:N,replace = TRUE),]
    fitted.trees[[nn]] <- rpart(model.formula,bagged.data,control = rpart.control(minsplit=round(N/sqrt(nn)),xval=0,cp=0))
  }
  return(fitted.trees)
}

predMatrix <- function(fitted.trees,input.data) {
  X <- NULL
  for (nn in 1:length(fitted.trees)) {
    X <- cbind(X,model.matrix(~-1+.,data=data.frame(class=factor(rpart.predict.leaves(fitted.trees[[nn]],input.data)))))
  }
  return(X)
}

glmTree <- function(model.formula,input.data) {
  fitted.trees <- bagged.trees(model.formula,input.data)

  y <- input.data[,as.character(model.formula)[2]]
  X <- predMatrix(fitted.trees,input.data)

  return(list(model=cv.glmnet(X,y,intercept=TRUE,standardize=FALSE,alpha=.5),fitted.trees=fitted.trees))
}
