
#' @export
ngrResolve <- function(ngr, centre=FALSE, grid="British"){
  if (!is.na(ngr)){
      ##Convert UK grid references to full coordinates
      ngr <- gsub(pattern="[[:space:]]+", replacement="", ngr)
      lets <- substr(ngr,1,2)
      if (grid=="British"){
          nums <- substr(ngr,3,nchar(ngr))
          prefixes <- data.frame(L=c("HP", "HT", "HU", "HY", "HZ", "NA", "NB", "NC", "ND", "NF", "NG", "NH", "NJ", "NK", "NL", "NM", "NN", "NO", "NR", "NS", "NT", "NU", "NW", "NX", "NY", "NZ", "SC", "SD", "SE", "SH", "SJ", "SK", "SM", "SN", "SO", "SP", "SR", "SS", "ST", "SU", "SV", "SW", "SX", "SY", "SZ", "TA", "TF", "TG", "TL", "TM", "TQ", "TR", "TV"), X=c("4", "3", "4", "3", "4", "0", "1", "2", "3", "0", "1", "2", "3", "4", "0", "1", "2", "3", "1", "2", "3", "4", "1", "2", "3", "4", "2", "3", "4", "2", "3", "4", "1", "2", "3", "4", "1", "2", "3", "4", "0", "1", "2", "3", "4", "5", "5", "6", "5", "6", "5", "6", "5"), Y=c("12", "11", "11", "10", "10", "9", "9", "9", "9", "8", "8", "8", "8", "8", "7", "7", "7", "7", "6", "6", "6", "6", "5", "5", "5", "5", "4", "4", "4", "3", "3", "3", "2", "2", "2", "2", "1", "1", "1", "1", "0", "0", "0", "0", "0", "4", "3", "3", "2", "2", "1", "1", "0"))
      } else if (grid=="Irish"){
          nums <- substr(ngr,2,nchar(ngr))
          nums[isOdd(nchar(nums))] <- paste(0,nums[isOdd(nchar(nums))],sep="")
          northings <- data.frame(L=c("A","B","C","D","E","F","G","H","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"),Y=rep(c(4:0),each=5))
          eastings <- data.frame(L=c("A","B","C","D","E","F","G","H","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"),X=rep(c(0:4),times=5))
          prefixes <- merge(eastings, northings, by="L")
      }
      epre <- as.character(prefixes$X[prefixes$L==lets])
      npre <- as.character(prefixes$Y[prefixes$L==lets])
      w <- nchar(nums) / 2
      e <- substr(nums,1,w)
      n <- substr(nums, w+1, nchar(nums))
      e <- paste(epre,e,sep="")
      n <- paste(npre,n,sep="")
      if (centre) mv <- "5" else mv <- "0"
      if (w==1){
          e <- as.numeric(paste(e,mv,"000",sep=""))
          n <- as.numeric(paste(n,mv,"000",sep=""))
      } else if (w==2){
          e <- as.numeric(paste(e,mv,"00",sep=""))
          n <- as.numeric(paste(n,mv,"00",sep=""))
      } else if (w==3){
          # six fig ref
          e <- as.numeric(paste(e,mv,"0",sep=""))
          n <- as.numeric(paste(n,mv,"0",sep=""))
      } else if (w==4){
          e <- as.numeric(paste(e,mv,sep=""))
          n <- as.numeric(paste(n,mv,sep=""))
      } else if (w==5){
          e <- as.numeric(paste(e,mv,sep=""))/10
          n <- as.numeric(paste(n,mv,sep=""))/10
      } else if (w==6){
          e <- as.numeric(paste(e,mv,sep=""))/100
          n <- as.numeric(paste(n,mv,sep=""))/100
      } else if (w==7){
          e <- as.numeric(paste(e,mv,sep=""))/1000
          n <- as.numeric(paste(n,mv,sep=""))/1000
      } else if (w==8){
          e <- as.numeric(paste(e,mv,sep=""))/10000
          n <- as.numeric(paste(n,mv,sep=""))/10000
      } else {
          e <- NA
          n <- NA
      }
      return(c(e,n))
  } else {
      return(NA)
  }
}

##
 
#' @title Convert BP dates to BC/AD format 
#' @description Converts calibrated BP dates to BC/AD dates, omitting 'year 0' 
#' @param x A numerical vector (currently no checks that these numbers are in a sensible range). 
#' @return A vector with BC/BCE dates expressed as negative numbers and AD/CE dates as positive ones.
#' @examples
#' BPtoBCAD(4200)
#' @export
BPtoBCAD <- function(x){
    res <- matrix(c(x, rep(NA,length(x))), ncol=2)
    res[x < 1950,2] <- 1950-res[x < 1950,1]
    res[x >= 1950,2] <- 1949-res[x >= 1950,1]
    return(res[,2])
}

#' @title Convert BC/AD dates to BP format
#' @description Converts BC/AD dates to BP fromat while handling the absence of 'year 0' 
#' @param x A numerical vector (currently no checks that these numbers are in a sensible range).
#' @return A vector with BC/BCE dates expressed as negative numbers and AD/CE dates as positive ones.
#' @examples
#' BCADtoBP(-1268)
#' @export
BCADtoBP <- function(x){
    res <- matrix(c(x, rep(NA,length(x))), ncol=2)
    res[x > 0,2] <- abs(res[x > 0,1] - 1950)
    res[x < 0,2] <- abs(res[x < 0,1] - 1949)
    return(res[,2])
}
