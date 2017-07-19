require(rhdf5)
require(bigh5)


#' Collective polarization order
#'
#' @param frameData A data frame with the columns R, C, u, and V specifying optical flow properties
#'
#' @return A scalar quantifying the polarization order (as in Couzin et al. 2002)
#' @export
#'
#' @examples
measurePolarization = function(frameData) {
  colnames(frameData) = c("R", "C", "u", "v")
  frameData$speed = sqrt(frameData$u^2+frameData$v^2)
  frameData = frameData[frameData$speed > 0,]
  frameData$unitu = frameData$u/frameData$speed
  frameData$unitv = frameData$v/frameData$speed
  sumu = sum(frameData$unitu)
  sumv = sum(frameData$unitv)
  polarization = sqrt(sumu^2 + sumv^2)/(nrow(frameData))
  polarization = frameData %>% summarise(sumu = sum(unitu), sumv = sum(unitv),count=n()) %>% mutate(polarization=sqrt(sumu^2+sumv^2)/count)
  return(polarization$polarization)
}




#' Collective rotation order
#'
#' @param frameData A data frame with the columns R, C, u, and V specifying optical flow properties
#'
#' @return A scalar quantifying the rotation order (as in Couzin et al. 2002)
#' @export
#'
#' @examples
measureRotation = function(frameData) {
  currFrame = frameData
  colnames(currFrame) = c("R", "C", "u", "v")
  currFrame$centerR = mean(currFrame$R) # What is the center of mass?
  currFrame$centerC = mean(currFrame$C)

  # Get unit vector velocities
  currFrame$speed = sqrt(currFrame$u^2+currFrame$v^2)
  currFrame = currFrame[currFrame$speed > 0,] # vectors of magnitude 0 will have undefined unit vectors; remove.
  currFrame$unitu = currFrame$u/currFrame$speed # Define the unit vectors, and the relative position of particles from center.
  currFrame$unitv = currFrame$v/currFrame$speed

  # Get unit position vectors from center.
  currFrame$relR = currFrame$R - currFrame$centerR # Define relative position and unit vector of relative position.
  currFrame$relC = currFrame$C - currFrame$centerC
  currFrame$distFromCenter = sqrt(currFrame$relR^2 + currFrame$relC^2)
  currFrame = currFrame[currFrame$distFromCenter > 0,]
  currFrame$unitR = currFrame$relR/currFrame$distFromCenter
  currFrame$unitC = currFrame$relC/currFrame$distFromCenter

  currFrame$rotationContribution = currFrame$unitR * currFrame$unitv - currFrame$unitC * currFrame$unitu # Cross product of two vectors.
  rotation = mean(currFrame$rotationContribution)
  return(rotation)
}




#' Collective dilatation order
#'
#' This is based off of Attanasi et al. 2014 paper on swarmes of midges where they define an order parameter that measures the extent to which
#' particles are moving outward (positive) or inward (negative) toward the center of mass of the collective.
#'
#' @param frameData A data frame with the columns R, C, u, and V specifying optical flow properties
#'
#' @return A scalar quantifying the dilatation order (as in Couzin et al. 2002)
#' @export
#'
#' @examples
measureDilatation = function(frameData) {
  currFrame = frameData
  colnames(currFrame) = c("R", "C", "u", "v")
  currFrame$centerR = mean(currFrame$R) # What is the center of mass?
  currFrame$centerC = mean(currFrame$C)

  currFrame$speed = sqrt(currFrame$u^2+currFrame$v^2)
  currFrame = currFrame[currFrame$speed > 0,] # vectors of magnitude 0 will have undefined unit vectors; remove.
  currFrame$unitu = currFrame$u/currFrame$speed # Define the unit vectors, and the relative position of particles from center.
  currFrame$unitv = currFrame$v/currFrame$speed

  currFrame$relR = currFrame$R - currFrame$centerR # Define relative position and unit vector of relative position.
  currFrame$relC = currFrame$C - currFrame$centerC
  currFrame$distFromCenter = sqrt(currFrame$relR^2 + currFrame$relC^2)
  currFrame = currFrame[currFrame$distFromCenter > 0,]
  currFrame$unitR = currFrame$relR/currFrame$distFromCenter
  currFrame$unitC = currFrame$relC/currFrame$distFromCenter

  currFrame$dilatationContribution = currFrame$unitC * currFrame$unitv + currFrame$unitR * currFrame$unitu # Dot product of two vectors.
  dilatation = mean(currFrame$dilatationContribution)
  return(dilatation)
}



#' Plot of velocity field
#'
#'This function takes in a 4-column data frame that contains the requisite data on a velocity field, where each row is
#'a vector in space which has an x and y position and an x and y velocity component. It returns a ggplot object.
#'
#' @param frameData A 4-column data frame consisting of (from left to right) the x position, y position, x-component, and y-component of the velocity.
#' @param arrowl A scalar multiple to apply to the velocity components, for visualization.
#' @param colormapped If F, plot as arrow line segments. Else, plot as color heatmap.
#'
#' @return A ggplot2 object that is the plot of the velocity field.
#' @export
#'
#' @examples
quiverPlot <- function(frameData,arrowl, colormapped=F) {
  require(ggplot2)
  colnames(frameData) = c("x","y","u","v")
  frameData = mutate(frameData,v=v*arrowl, u=u*arrowl)

  if (colormapped) {
    frameData = frameData %>% mutate(speed = sqrt(v^2 + u^2), angle = atan2(x = u, y = v))
    cbPalette = c("magenta","red","yellow","green","cyan","blue","magenta")
    p = ggplot(frameData, aes(x, y, fill=angle, alpha=speed)) + geom_tile() +
      scale_color_gradientn(name = expression(theta), colours = cbPalette, limits = c(-pi,pi), breaks=c(-pi, 0, pi),
                            labels=c(expression(paste("-",pi,sep="")), 0, expression(paste(pi)))) +
      scale_fill_gradientn(name = expression(theta), colours = cbPalette, limits = c(-pi,pi), breaks=c(-pi, 0, pi),
                           labels=c(expression(paste("-",pi,sep="")), 0, expression(paste(pi))))

  } else {
    p = ggplot(frameData, aes(x, y, xend=x+u, yend=y+v)) +
      geom_segment(arrow=arrow(angle=20,length=unit(0.2,"cm")))
  }
  p = p + coord_fixed(ratio = 1) + theme(axis.title = element_blank())
  return(p)
}
