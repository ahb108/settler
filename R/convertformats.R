

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


ah2sp <- function(x, increment=360, rnd=10, proj4string=CRS(as.character(NA))){
 require(alphahull)
 require(maptools)
 if (class(x) != "ahull"){
   stop("x needs to be an ahull class object")
 }
 # Extract the edges from the ahull object as a dataframe
 xdf <- as.data.frame(x$arcs)
 # Remove all cases where the coordinates are all the same      
 xdf <- subset(xdf,xdf$r > 0)
 res <- NULL
 if (nrow(xdf) > 0){
   # Convert each arc to a line segment
   linesj <- list()
   prevx<-NULL
   prevy<-NULL
   j<-1
   for(i in 1:nrow(xdf)){
     rowi <- xdf[i,]
     v <- c(rowi$v.x, rowi$v.y)
     theta <- rowi$theta
     r <- rowi$r
     cc <- c(rowi$c1, rowi$c2)
     # Arcs need to be redefined as strings of points. Work out the number of points to allocate in this arc segment.
     ipoints <- 2 + round(increment * (rowi$theta / 2),0)
     # Calculate coordinates from arc() description for ipoints along the arc.
     angles <- anglesArc(v, theta)
     seqang <- seq(angles[1], angles[2], length = ipoints)
     x <- round(cc[1] + r * cos(seqang),rnd)
     y <- round(cc[2] + r * sin(seqang),rnd)
     # Check for line segments that should be joined up and combine their coordinates
     if (is.null(prevx)){
       prevx<-x
       prevy<-y
     } else if (x[1] == round(prevx[length(prevx)],rnd) && y[1] == round(prevy[length(prevy)],rnd)){
         if (i == nrow(xdf)){
         #We have got to the end of the dataset
           prevx<-append(prevx,x[2:ipoints])
           prevy<-append(prevy,y[2:ipoints])
           prevx[length(prevx)]<-prevx[1]
           prevy[length(prevy)]<-prevy[1]
           coordsj<-cbind(prevx,prevy)
           colnames(coordsj)<-NULL
           # Build as Line and then Lines class
           linej <- Line(coordsj)
           linesj[[j]] <- Lines(linej, ID = as.character(j))
         } else {
           prevx<-append(prevx,x[2:ipoints])
           prevy<-append(prevy,y[2:ipoints])
         }
       } else {
     # We have got to the end of a set of lines, and there are several such sets, so convert the whole of this one to a line segment and reset.
         prevx[length(prevx)]<-prevx[1]
         prevy[length(prevy)]<-prevy[1]
         coordsj<-cbind(prevx,prevy)
         colnames(coordsj)<-NULL
     # Build as Line and then Lines class
         linej <- Line(coordsj)
         linesj[[j]] <- Lines(linej, ID = as.character(j))
         j<-j+1
         prevx<-NULL
         prevy<-NULL
       }
   }
   # Promote to SpatialLines
   lspl <- SpatialLines(linesj)
   # Convert lines to polygons
   # Pull out Lines slot and check which lines have start and end points that are the same
   lns <- slot(lspl, "lines")
   polys <- sapply(lns, function(x) { 
     crds <- slot(slot(x, "Lines")[[1]], "coords")
     identical(crds[1, ], crds[nrow(crds), ])
   }) 
   # Select those that do and convert to SpatialPolygons
   polyssl <- lspl[polys]
   list_of_Lines <- slot(polyssl, "lines")
   sppolys <- SpatialPolygons(list(Polygons(lapply(list_of_Lines, function(x) { Polygon(slot(slot(x, "Lines")[[1]], "coords")) }), ID = "1")), proj4string=proj4string)
   # Create a set of ids in a dataframe, then promote to SpatialPolygonsDataFrame
   hid <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "ID"))
   areas <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "area"))
   df <- data.frame(hid,areas)
   names(df) <- c("HID","Area")
   rownames(df) <- df$HID
   res <- SpatialPolygonsDataFrame(sppolys, data=df)
   res <- res[which(res@data$Area > 0),]
 }  
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
