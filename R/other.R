

isWholeNumber <- function(x, tol = .Machine$double.eps^0.5) {
  #Checks whether a value is a whole number (used in nnhistMC reporting)
  abs(x - round(x)) < tol
}

# Check for odd number
isOdd <- function(x){ x %% 2 != 0 }

## 

# Find-Replace function
mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}
##
