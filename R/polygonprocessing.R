

removeIslands <- function(x){
  # Remove holes form polygon (based on function in wild1 library)
  require(sp)
  isles <-  unlist(lapply(x@Polygons, function(p) p@hole))
  p <- Polygons(x@Polygons[!isles], ID=x@ID)
  return(p)
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
  
spdfClip <- function(toclip, clipby){

  #Class checks
  if(!class(toclip)[1] %in% c("SpatialLines","SpatialPolygons","SpatialLinesDataFrame","SpatialPolygonsDataFrame")){stop("'toclip' is not a SpatialLines* or SpatialPolygons* objects.")}
  if(!inherits(clipby, "SpatialPolygons")){stop("'clipby' is not a SpatialPolygons* object.")}
  #Projection compatibility
  if (!identical(proj4string(toclip),proj4string(clipby))){stop("'toclip' and 'clipby' are not the same CRS.")}
  #Set up up dataframe to fill if necessary
  if (class(toclip)[1]=="SpatialLinesDataFrame"){
    tc <- as(toclip,"SpatialLines")
    df <- toclip@data[1, ,drop=FALSE]
    df <- df[-1, ,drop=FALSE]
    df <- cbind(df,data.frame(RN=numeric(0)))
  } else if (class(toclip)[1]=="SpatialPolygonsDataFrame"){
    tc <- as(toclip,"SpatialPolygons")
    #df <- toclip@data[1,]
    df <- toclip@data[1, ,drop=FALSE]
    df <- df[-1, ,drop=FALSE]
    df <- cbind(df,data.frame(RN=numeric(0)))
  } else {
    tc <- toclip
  }
  #Loop through and clip (preserving dataframe links if necessary)
  clippedlist <- vector("list", length(tc))
  dflist <- vector("list", length(tc))
  for (i in 1:length(tc)){
    clippedtmp <- gIntersection(tc[i,],clipby)
    if (!is.null(clippedtmp)){
      if (class(tc)[1]=="SpatialLines"){
        RN <- lapply(slot(tc[i,],"lines"), function(x) {slot(x, "ID")})[[1]][1]
      } else {
        RN <- lapply(slot(tc[i,],"polygons"), function(x) {slot(x, "ID")})[[1]][1]
      }
      clippedtmp <- spChFIDs(clippedtmp, RN)
      clippedlist[[i]] <- clippedtmp      
      if (class(toclip)[1]=="SpatialLinesDataFrame"||class(toclip)[1]=="SpatialPolygonsDataFrame"){
        dfi <- toclip@data[i, ,drop=FALSE]
        dfi$RN <- RN
        dflist[[i]] <- dfi
      }
    }
  }
  clippedlist <- clippedlist[!sapply(clippedlist, is.null)]
  clippedall <- do.call("rbind", clippedlist)
  dflist <- dflist[!sapply(dflist, is.null)]
  df <- do.call("rbind", dflist) 
  #Recombine the results with the dataframe if necessary
  if (!is.null(clippedall)){
    if (class(toclip)[1]=="SpatialLinesDataFrame"){
      rownames(df) <- df$RN
      clippedall <- SpatialLinesDataFrame(clippedall, data=df)
    } else if (class(toclip)[1]=="SpatialPolygonsDataFrame"){
      rownames(df) <- df$RN
      clippedall <- SpatialPolygonsDataFrame(clippedall, data=df)
    }
  }
  return(clippedall)
}

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

gExplode <- function(x) {
  # a modified version of something by Josh O'Brien
    a <- deparse(substitute(x))
    if(class(x)[[1]] == "SpatialLines" | class(x)[[1]] == "SpatialLinesDataFrame"){
        b <- unlist(lapply(x@lines, function(c) c@Lines))
        SpatialLines(lapply(seq_along(b), function(c) Lines(b[c], ID=paste0(a,c))), proj4string=CRS(proj4string(x)))
    } else if (class(x)[[1]] == "SpatialPolygons" | class(x)[[1]] == "SpatialPolygonsDataFrame"){
        b <- unlist(lapply(x@polygons, function(c) c@Polygons))
        SpatialPolygons(lapply(seq_along(b), function(c) Polygons(b[c], ID=paste0(a,c))), proj4string=CRS(proj4string(x)))
    } else {
        stop("Input must be of class SpatialLines* or SpatialPolygons*")
    }
}

##

buffers <- function(x, bands, xIds=NULL, rings=TRUE, bMerge=FALSE, ...){
    for(i in 1:length(bands)){
        if (bMerge){
            pids <- paste("",bands[i], sep="_")
        } else if (!is.null(xIds)){
            pids <- paste(xIds,bands[i],sep="_")
        } else {
            pids <- paste(1:length(x),bands[i],sep="_")
        }
        b <- gBuffer(x, byid=TRUE,id=NULL, width=bands[i], ...)
        if (bMerge){ b <- gUnaryUnion(b) }
        if (i == 1){
            b <- spChFIDs(b, pids)
            bufs <- b
            a <- b
        } else {
            d <- b
            # clip out previous buffer
            if (rings){
                d <- gDifference(d, a, byid=TRUE, id=NULL)
                if (!bMerge){
                    d <- d[seq(1,length(d),sqrt(length(d))+1)]
                }
            }
            d <- spChFIDs(d, pids)
            bufs <- spRbind(bufs,d)
            a <- b
        }    
    }
    if (!is.null(xIds)){
        # Convert to SpatialPolygonsDataFrame
        pids <- sapply(slot(bufs, "polygons"), function(x) slot(x, "ID"))
        df <- data.frame(PID=pids, row.names=pids)
        bufs <- SpatialPolygonsDataFrame(bufs, data=df)
        if (!bMerge){
            bufs$xID <- sapply(strsplit(as.character(bufs$PID),"_"), "[", 1)
        }
        bufs$Buffer <- as.numeric(sapply(strsplit(as.character(bufs$PID),"_"), "[", 2))
    }
    return(bufs)
}

##

voronoi <- function(x,xID,win=NULL, ...){
    # wraps deldir for SpatialPolygons
    if (is.null(win)){ rw1 <- NULL } else { rw1 <- as.vector(t(bbox(win))) }
    a <- deldir(coordinates(x)[,1], coordinates(x)[,2], rw=rw1, ...)
    w <- tile.list(a)
    polys <- vector(mode='list', length=length(w))
    for (i in seq(along=polys)) {
        pcrds <- cbind(w[[i]]$x, w[[i]]$y)
        pcrds <- rbind(pcrds, pcrds[1,])
        polys[[i]] <- Polygons(list(Polygon(pcrds)), ID=as.character(i))
    }
    sp1 <- SpatialPolygons(polys)
    spdf1 <- SpatialPolygonsDataFrame(sp1, data=data.frame(PtID=xID, row.names=sapply(slot(sp1, 'polygons'), function(x) slot(x, 'ID'))))
    if (!is.null(win)){
        proj4string(spdf1) <- CRS(proj4string(win))
        spdf1 <- spdfClip(spdf1, win)
    }
    return(spdf1)
}

##

mwvoronoi <- function(pts, ptID, ptWeights=rep(1,nrow(pts)), win=NULL, maxdist=NULL, nsteps=36, ...){

    # Multiplicatively weighted Voronoi tesselation
    for (a in 1:nrow(pts)){
        da <- NULL
        for (b in 1:nrow(pts)){
            if (a != b){
                if (ptWeights[a] <= ptWeights[b]){
                    p1 <- coordinates(pts)[a,]
                    p2 <- coordinates(pts)[b,]
                    w1 <- ptWeights[a]
                    w2 <- ptWeights[b]
                } else {
                    p1 <- coordinates(pts)[b,]
                    p2 <- coordinates(pts)[a,]
                    w1 <- ptWeights[b]
                    w2 <- ptWeights[a]
                }
                d12 <- sqrt(sum((p1 - p2) ^ 2)) #euclidean distance
                c1 <- (w2^2*p1-w1^2*p2) / (w2^2-w1^2)
                r1 <- (w1*w2*d12) / (w1^2-w2^2)
                if (ptWeights[a] == ptWeights[b]){
                    tmpts <- rbind(pts[a,],pts[b,])
                    ac <- voronoi(tmpts,c(w1,w2),win=win,nsteps=nsteps)
                    ac <- as(ac[1,],"SpatialPolygons")
                } else {
                    ac <- spEllipse(c1[1],c1[2],r1, nsteps=nsteps)
                }
                if (ptWeights[a] > ptWeights[b]){
                    ac <- gDifference(win,ac)
                }
                if (is.null(da)){
                    da <- gIntersection(ac,win)
                } else {
                    da <- gIntersection(ac,da)
                }
            }
        }
        if (!is.null(maxdist)){
            if (!is.numeric(maxdist) | maxdist <= 0){
                stop("Maximum distance must be a postive number.")
            } else {
                b <- gBuffer(pts[a,], width=maxdist, quadsegs=round(nsteps/4,0))
                da <- gIntersection(da,b)
            }
        }
        da <- spChFIDs(da, paste(pts$ptID[a],sep=""))
        da <- SpatialPolygonsDataFrame(da, data=data.frame(ptID=pts$ptID[a], weight=ptWeights[a], row.names=paste(pts$ptID[a],sep=""))) # convert to SpatialPolygonsDataFrame
        if (a==1){
            res <- da
        } else {
            res <- spRbind(res,da)
        }
    }
    proj4string(res) <- CRS(proj4string(pts))
    return(res)
}

##
