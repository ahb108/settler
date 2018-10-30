#' Produce a very simple map scale-bar
#'
#' This function draw a simple scale-bar in the form of a box with a certain width and height (in map units) and hung on its top right corner.
#'
#' @param x Top left x coordinate
#' @param y Top left y coordinate
#' @return a plotted polygon
#' @export
scaleBar <- function(x, y, width, height,...){
  polygon(x=c(x, (x+width), (x+width), x), y=c(y, y, (y-height), (y-height)), ...)
}

#' Produce an circle or ellipse as a SpatialPolygons object
#'
#' This function draw a simple circle or ellipse, centred on a certain coordinate pair and of a certain size (the circle's radius or height, width and azimuth for an ellipse). The smooothness of the circle can be controlled by the chosen number of vertices placed along the perimeter, and the circle can be assigned a projection system if desired. 
#'
#' @param xc x coordinate of circle or ellipse centre
#' @param yc y coordinate of circle or ellipse centre
#' @param a the size of the circle radius in map units or the length of the longer axis of the ellipse
#' @param b the size of the shorter axis of the ellipse, if necessary (default is circular with b=a).
#' @param angle the map orientation of the shorter axis of the ellipse, if necessary (default is longer axis is oriented north).
#' @param nsteps the number of vertices used to describe the circle of the ellipse, if necessary (i.e. the smoothness of the plotted circle, default is 360)
#' 
#' @return a plotted polygon
#' @export
spEllipse <- function(xc, yc, a, b=a, angle=0, nsteps=360, proj4string=CRS(as.character(NA)), noSp=FALSE, ...){

  theta <- seq(0,(2*pi),len=nsteps)
  angle <- (90 - angle) * (pi /180)

  x <- xc + ((a * cos(theta) * cos(angle)) - (b * sin(theta) * sin(angle)))
  y <- yc + ((a * cos(theta) * sin(angle)) + (b * sin(theta) * cos(angle)))

  cpol<-cbind(x,y)
  cpol <- rbind(cpol, cpol[1,])
  if (noSp){
    return(cpol)
  } else {
    cpolp <- Polygons(list(Polygon(cpol)), ID="1")
    cpolsp <- SpatialPolygons(list(cpolp), proj4string=proj4string)
    return(cpolsp)
  }
}


#' Produce a simple north arrow
#'
#' This function draw a simple north arrow of a certain size and at a certain map location.
#'
#' @param x x coordinate for the centre of the bounding circle of the north arrow
#' @param y y coordinate for the centre of the bounding circle of the north arrow
#' @param r Radius of the bounding circle of the north arrow
#' 
#' @return a plotted polygon
#' @export

northArrow <- function(x, y, r, type="simple", fill="black", bkgrd=NA, ...){
  
  if (type=="simple"){
    polygon(spEllipse(x,y,r, noSp=TRUE), col=bkgrd, ...)
    bl <- c((x + r * sin(210*pi/180)),(y + r * cos(210*pi/180)))
    br <- c((x + r * sin(150*pi/180)),(y + r * cos(150*pi/180)))
    polygon(x=c(x, br[1], x,  bl[1],x),y=c(y+r, br[2], (y-(0.5 * r)), bl[2],y+r), col=fill, ...)
  } else {stop("This north arrow type does not exist. Try leaving `type' blank.")}
}

##

# Create a graduated symbol function
gradSymbols <- function(x, mincex, maxcex, valrange=NULL, sqrttrans=FALSE,...){
  if (sqrttrans) { x <- sqrt(x) }
  if (is.null(valrange)){
    minx <- min(x)
    maxx <- max(x)
  } else if (is.vector(valrange) & length(valrange)==2){
    minx <- valrange[1]
    maxx <- valrange[2]
  } else { stop("valrange must be NULL or a vector of length 2.\n") }
  y <- mincex + (((x-minx)*maxcex - (x-minx)*mincex) / (maxx-minx))
  return(y)
}


#' Addition of both a north arrow and scalebar to map
#'
#' This is a helper function to add both north arrow and scalebar to a map.
#'
#' @param x x coordinate for the centre of the added elements (the centre of the text label)
#' @param y y coordinate for the centre of the added elements (the centre of the text label)
#' @param scalesize width of the scalebar in map units (this function assumes projcted units but will not throw a warning if unprojected units are used.
#' @param scalelabel text label for the state width of the scalebar (it requires the user to enter this and is not automatically generated from scalesize)
#' 
#' @return a plotted north arrow and scale-bar
#' @export
#' 
mapElements <- function(x, y, scalesize, scalelabel, cex.label=0.75, type.north="simple", lwd.north=0.5, offset.north=height.scale*1.5, r.north=scalesize/5, bkgrd.north="white", fill.north="black", offset.scale=height.scale, height.scale=scalesize/8, border.scale="black", col.scale="white", lwd.scale=0.5){
    
    northArrow(x, y+(r.north+offset.north), r=r.north, type=type.north, lwd=lwd.north, fill=fill.north, bkgrd=bkgrd.north)
    scaleBar(x-(scalesize/2), y-(offset.scale), width=scalesize, height=height.scale, lwd=lwd.scale, col=col.scale, border=border.scale)
    text(x, y, labels=scalelabel, cex=cex.label)

}

#' @export
#' 
pieSymbols <- function (x, y, values, sizes=NULL, labels=names(x), edges=360, clockwise=FALSE, init.angle=if(clockwise) 90 else 0, density=NULL, angle=45, col=NULL, border=NULL, lty=NULL, main=NULL, ...){
    ## This is a modification of the basic pie() function to allow overplotting on map data
    t2xy <- function(t){
        t2p <- twopi * t + init.angle * pi/180
        list(x=r * cos(t2p), y=r * sin(t2p))
    }
    for (a in 1:nrow(values)){
        v <- values[a,]
        if (is.null(sizes)){
            sizes <- rep((max(y)-min(y))/3,length(sizes))
        } else if (length(sizes)==1){
            r <- rep(sizes,length(sizes))
        } else {
            r <- sizes[a]
        }
        if (!is.numeric(v) || any(is.na(v) | v < 0)){ 
            stop("values must be positive and numeric.")
        }
        v <- c(0, cumsum(v)/sum(v))
        dv <- diff(v)
        nv <- length(dv)
        if (!is.null(border)){ border <- rep_len(border, nv) }
        if (!is.null(lty)){ 
            lty <- rep_len(lty, nv)
            angle <- rep(angle, nv)
        }
        if (!is.null(density)){ density <- rep_len(density, nv) }
        if (clockwise){ twopi <- -2 * pi } else { twopi <- 2 * pi }
        for (i in 1L:nv) {
            n <- max(2, floor(edges * dv[i]))
            P <- t2xy(seq.int(v[i], v[i + 1], length.out=n))
            polygon(c(x[a]+P$x,x[a]), c(y[a]+P$y,y[a]), density=density[i], angle=angle[i], 
                    border=border[i], col=col[i], lty=lty[i])
        }
    }
}
