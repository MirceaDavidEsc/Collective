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
  
  currFrame$rotationContribution = currFrame$unitC * currFrame$unitv - currFrame$unitR * currFrame$unitu # Cross product of two vectors.
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
  
  currFrame$dilatationContribution = currFrame$unitR * currFrame$unitv + currFrame$unitC * currFrame$unitu # Dot product of two vectors.
  dilatation = mean(currFrame$rotationContribution)
  return(dilatation)
}



#' Measure the descriptive statistical properties of a collective.
#'
#' @param dataFolder The folder containing the collective order dataset to use
#' @param frameChunkSize The number of frames to load into memory at any given time. Useful for batch read-write (memory-dependent)
#' @param dplyr Is the dplyr package available?Can use it to speed up processing.
#'
#' @return A data frame with the frame number, polarization, rotation, mean speed, mean v, mean u, and the max speed for each frame.
#' @export
#'
#' @examples
measureInternalOrder = function(dataFolder, frameChunkSize=3600, dplyr=FALSE) {
  require(doParallel)
  #Read the frames sequentially.
  #Produce the order data s doen in CharacterizeVectorFields.m
  #1. Get the number of frames to search

  frameDataIndices = h5read(paste(dataFolder,"longformOpticalFlow.hdf5",sep="/"),"/longformFrameIndices")
  frameDataIndices = as.data.frame(frameDataIndices) %>% filter(V1 > 0)
  colnames(frameDataIndices) = c("startFrame","dataIndexStart","dataIndexEnd")
  numFrames = dim(frameDataIndices)[1]
  H5close()

  #2. Produce a data startFrame to store the aggregate data.
  frameStartPoints = seq(from = 1, by = frameChunkSize, to = numFrames)
  if (file.exists("OrderMeasures.csv")) {
    orderMeasuresDataFrame = read.csv("OrderMeasures.csv")
    orderMeasuresDataFrame = filter(orderMeasuresDataFrame,Frame > 0) %>% select(Frame,Polarization,Rotation,AverageSpeed,averageVelocityX,averageVelocityY)
    frameStartPoints = frameStartPoints[frameStartPoints > max(orderMeasuresDataFrame$Frame)]
  } else {
    orderMeasuresDataFrame = as.data.frame(matrix(data = 0, nrow = length(frameDataIndices$startFrame), ncol = 7))
    colnames(orderMeasuresDataFrame) = c("Frame","Rotation","Polarization","AverageSpeed","averageVelocityX","averageVelocityY","maxSpeed")
  }

  #3. For all frames, load the frames and calculate the order measures
  for (startFrame in frameStartPoints) {
    print(paste("Measuring chunk",startFrame))
    chunkEnd = min(startFrame+frameChunkSize,numFrames)-1
    theseFrames = startFrame:chunkEnd
    #Get the startFrame data
    print("Reading in data")
    frameData = readVelocityFieldAtFrame("longformOpticalFlow.hdf5",theseFrames)
    H5close()

    print("Processing measurements...")
    cl = makeCluster(7)
    print("Rotation...")
    rotation = foreach(thisFrame = frameData, .combine=c, .packages="dplyr") %dopar% {
      measureRotation(thisFrame)
    }
    print("Polarization...")
    polarization = foreach(thisFrame = frameData, .combine=c, .packages="dplyr") %dopar% {
      measurePolarization(thisFrame)
    }
    print("Mean speed...")
    meanSpeed = foreach(thisFrame = frameData, .combine=c, .packages="dplyr") %dopar% {
      mean(sqrt(thisFrame$u^2+thisFrame$v^2))
    }

    meanU = foreach(thisFrame = frameData, .combine=c, .packages="dplyr") %dopar% {
      mean(thisFrame$u)
    }

    meanV = foreach(thisFrame = frameData, .combine=c, .packages="dplyr") %dopar% {
      mean(thisFrame$v)
    }

    maxSpeed = foreach(thisFrame = frameData, .combine=c, .packages="dplyr") %dopar% {
      max(sqrt(thisFrame$v^2+thisFrame$u^2))
    }

    thisChunkData = cbind(theseFrames,polarization,rotation,meanSpeed,meanU,meanV,maxSpeed)
    colnames(thisChunkData) = colnames(orderMeasuresDataFrame)
    orderMeasuresDataFrame[theseFrames,] = thisChunkData

    rm(list="frameData")
    gc()
  }
  write.csv(orderMeasuresDataFrame, file="OrderMeasures.csv",row.names=F)
}
