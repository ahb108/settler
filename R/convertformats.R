
#' @title Resolve a Britsh or Irish National Grid Reference. 
#' @description Convert a Britsh or Irish National Grid Reference to a coordinate pair.
#' @param ngr A character string with a British National Grid Reference
#' @param centre Return the coordinate of the origin of the local square (default) or the centre.
#' @param grid Input grid system (currently British or Irish)
#' @return A numeric coordinate pair
#' @examples
#' ngrResolve("SJ456782")
#' @export
ngrResolve <- function(ngr, centre=FALSE, grid="British"){
  if (!is.na(ngr)){
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

#' @title Change the resolution of a British National Grid Reference. 
#' @description Increase of decrease the resolution of a British National Grid Reference. 
#' @param x A character string with a British National Grid Reference
#' @param figures The desired precision (an even whole number)
#' @return A character string with a new British National Grid Reference
#' @examples
#' x <- "SJ456782"
#' ngrResolution(x,4)
#' ngrResolution(x,8)
#' ngrResolution(x,7)
#' @export
ngrResolution <- function(x, figures){

    if (figures%%1!=0 | figures < 0 | figures %% 2 != 0){
        stop("The figures argument must be an even, whole number greater than or equal to 0.")
    }
    prefix <- substr(x,1,2)
    suffix <- substr(x,3,nchar(x))
    if (nchar(suffix) %% 2 == 0){
        xcomp <- substr(suffix,1,nchar(suffix)/2)
        ycomp <- substr(suffix,(nchar(suffix)/2)+1,nchar(suffix))
        if (figures <= nchar(suffix)){
            newxcomp <- substr(xcomp,1, figures/2)
            newycomp <- substr(ycomp,1, figures/2)
            newngr <- paste(prefix, newxcomp, newycomp, sep="")
        } else {
            add0s <- (figures - nchar(suffix))/2
            add0s <- paste0(rep("0",add0s), collapse="")
            newngr <- paste(prefix, xcomp, add0s, ycomp, add0s, sep="")
        }
    } else {
        stop("Grid reference has an incorrect (odd) number of figures.")
    }
    return(newngr)
}

#' @title Convert a coordinate pair to a different coordinate system 
#' @description Converts a single coordinate pair to a different coordinate system  
#' @param coords A vector of length 2 with the coordinates
#' @param proj1 The CRS of the input coordinates 
#' @param proj2 The target CRS of the output coordinates
#' @return A vector of length 2 with the converted coordinate pair.
#' @examples
#' bng <- CRS("+init=epsg:27700") # British National Grid
#' lonlat <- CRS("+init=epsg:4326") # LatLon WGS84
#' xyTransform(c(345600,378200), bng, lonlat)
#' @export
xyTransform <- function(coords, proj1, proj2){
    if (as.character(class(proj1))!="CRS" | as.character(class(proj2))!="CRS"){
        stop("One or both of proj1 and proj2 are not correct CRS class definitions.")
    }
    tmp <- SpatialPoints(t(as.matrix(coords)), proj4string=proj1)
    tmp <- spTransform(tmp,proj2)
    coordinates(tmp)
    res <- as.numeric(coordinates(tmp))
    return(res)
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
