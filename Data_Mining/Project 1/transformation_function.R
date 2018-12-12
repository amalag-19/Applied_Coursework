# Function: trsfrm.gen 
#
# Purpose: Generate transformation functions. Normalization during 
# cross validation needs to be based on training data only. So, since
# sample means and sds from the training data will change depending
# on which partition of the data is used for training, different
# transformations will be used for each iteration.
#
# Arguments: 
# train.data - training data matrix with column names (should all
#             be numeric). 
# Returns:
# trsfrm - function for log transformation and normalization based
# on sample means and sds from trainin data. trsfrm() takes a matrix
# with the same column names passed to trsfrm.gen (can be both 
# training and test data)
#
# Example:
# tr.fun <- trsfrm.gen(train.data)
# tr.fun(test.data)


trsfrm.gen <- function(train.data){
  ncols <- dim(train.data)[2]
  for(i in 1:ncols){
    if(names(train.data)[i] %in% c("BrdIndx", "Area", 
                                   "Compact", "ShpIndx","SD_NIR", "LW")){
      train.data[,i] <- log(train.data[,i])
    }
  }
  train.data <- as.matrix(train.data)
  mns <- apply(X = train.data, MARGIN = 2, FUN = mean)
  sds <- apply(X = train.data, MARGIN = 2, FUN = mean)
  trsfrm <- function(x, mn = mns, sd = sds){
    ncols <- dim(x)[2]
    for(i in 1:ncols){
      if(names(x)[i] %in% c("BrdIndx", "Area", 
                            "Compact", "ShpIndx","SD_NIR", "LW")){
        x[,i] <- log(x[,i])
      }
    }
    x <- as.matrix(x)
    nrml.data <- matrix(NA, ncol = ncol(x), nrow = nrow(x))
    xcols <- dim(x)[2]
    for(i in 1:xcols){
      nrml.data[,i] <- (x[,i] - mn[i])/sd[i]
    }
    return(nrml.data)
  }
  return(trsfrm)
}

