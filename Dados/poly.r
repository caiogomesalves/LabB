##
##  poly.R -- R version of poly.S -- Brad conversion 10/23/05
##
###################
# S-plus polygon extractor v. 2.1  Lance Waller 8/22/01 
# (including some of Brad Carlin's output formatting)

# This converts Splus map() format boundary files into
# a list suitable for input to GeoBUGS v. 1.4

# Uses mapgetl, mapgetg, and makepoly functions.

#  first:  source("poly.S")
#  then (e.g.): mkpoly("idaho")
#  This will generate idaho.txt, suitable for reading into GeoBUGS.
###################

mkpoly <- function(state)
{
  toradians <- atan(1)/45
  radiusearth <- 0.5*(6378.2+6356.7)
  sine51 <- sin( 51.5*toradians )

  ## in lieu of a call to X11():
  {
  outfile2 <- paste(state,".ps",sep="")
  postscript(outfile2,horizontal=F)
  }
  ##

  # First, list state for which county boundaries are desired:
#  state <- "idaho"
  #state <- "minnesota"
  #state <- "connecticut"

  library(maps)

  outfile <- paste(state,".txt",sep="")
  # mnpolys <- map("county",state,fill=T)
  namesvec <- map("county",state,namesonly=T)

  write(paste("map:",length(namesvec),"\n"),outfile)

  for (i in 1:length(namesvec)) {
    write(paste( i, paste("grid",i,sep="")),outfile,append=T)
  }

  for (i in 1:length(namesvec)) {
    exact <- T
    gon <- mapname("county",namesvec[i],exact)
    line <- mapgetg("county", gon, fill=T)
    coord <- mapgetl("county", line$number)
  
    gonsize <- line$size
#    Brad tweaks, next 3 lines:
    color <- 1
    color <- rep(color, length = length(gonsize))
    keep <- !is.na(color)
    coord[c("x", "y")] <- makepoly(coord, gonsize, keep)

  # Add repeat of first point to close polygon
  #    Brad tweak:  Comment these out; not needed by GeoBUGS!
  #  coord$x <- c(coord$x,coord$x[1])
  #  coord$y <- c(coord$y,coord$y[1])
  
  # Lance rough conversion of US lat/long to km (used by GeoBUGS):
  #   (see also forum.swarthmore.edu/dr.math/problems/longandlat.html)
# radius of earth: 
# r = 3963.34 (equatorial) or 3949.99 (polar) mi
#   = 6378.2 or 6356.7 km
# which implies: km per mile  = 1.609299 or 1.609295
# a change of 1 degree of latitude corresponds to the same number
# of km, regardless of longitude.  arclength=r*theta, so the
# multiplier for coord$y should probably be just the radius of
# earth.  
# On the other hand, a change of 1 degree in longitude corresponds
# to a different distance, depending on latitude.  (at N pole,
# the change is essentially 0.  at the equator, use equatorial
# radius.
# Perhaps for U.S., might use an "average" latitude, 30 deg is
# roughly Houston, 49deg is most of N bdry of continental 48 states.
# 0.5(30+49)=39.5 deg.  so use r approx 6378.2*sin(51.5)

##@@  Let's try these new multipliers:
#    coord$x <- (coord$x*atan(1)/45)*6416.6
#    coord$y <- (coord$y*atan(1)/45)*6416.6

    coord$x <- (coord$x*toradians)*radiusearth*sine51
    coord$y <- (coord$y*toradians)*radiusearth
##@@

    coordmat <- cbind( rep(paste("grid",i,sep=""),length(coord$x)) ,
                     round(coord$x,5), 
                     round(coord$y,5) )
                     
    if (i == 1) {
        write("",outfile,append=T)}
    else {
      write(c(NA,NA,NA),outfile,append=T)
    }
    write(t(coordmat),outfile,append=T,ncol=3)                            

    # Uncomment if you want the "movie" of county maps...               
      plot(coord,type="l")
      polygon(coord)
    
    # A counter for the loop, comment out if unnecessary.
      print( paste("i = ",i) )
  } ## for i

  write("END",outfile,append=T)
}


