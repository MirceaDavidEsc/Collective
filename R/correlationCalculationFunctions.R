require(bigsplines)
require(dplyr)




#' Get correlation lengths from profiles
#'
#' This function takes a data frame (with grouping) and determines the zero-intercept (correlation length) for each profile.
#' It is dependent on nesting ahd purrr map functions to iterate over all the profiles.
#'
#' @param correlationProfiles
#'
#' @return
#' @export
#'
#' @examples
calcCorrelationLengths = function(groupedProfiles) {

  # Nest the profiles appropriately so that it can work with my functions.
  nestedDist = groupedProfiles %>% select(domain) %>% nest(.key = "domain")
  nestedV = groupedProfiles %>% select(vCorr) %>% nest(.key = "mVCorr")
  nestedD = groupedProfiles %>% select(dCorr) %>% nest(.key = "mDCorr")
  nestedS = groupedProfiles %>% select(sCorr) %>% nest(.key = "mSCorr")
  nestedProfiles = inner_join(nestedV, nestedD) %>% inner_join(nestedS) %>% inner_join(nestedDist)

  # Get the zero-crossing of all the correlation profiles
  nestedProfiles = nestedProfiles %>% mutate(vZero = map2(mVCorr, .y = domain, zerosFromDF),
                                             dZero = map2(mDCorr, .y = domain, zerosFromDF),
                                             sZero = map2(mSCorr, .y = domain, zerosFromDF))

  corrLengths = nestedProfiles %>% select(-domain, -mVCorr, -mDCorr, -mSCorr) %>% unnest()
  return(corrLengths)
}


# The function below is used so I don't need to store all pairwise correlations, which is memory-intensive.
#' Calculate fluctuation correlations in a full velocity field.
#'
#' @param velocityField A data frame containing all the vectors in a velocity field (full velocities)
#' @param sampleSize How many vectors to use from the original field for pair-wise correlation comparisons.
#'
#' @return The correlation functions/profiles for the vector field.
#' @export
#'
#' @examples
calcCorrelationFunction = function(velocityField, sampleSize = 2000) {
  if (dim(velocityField)[1] > sampleSize) {
    sampleRows = sample(1:dim(velocityField)[1], sampleSize)
  } else {
    sampleRows = 1:dim(velocityField)[1]
  }


  # Get pairwise indices for a manageable subset of vectors, and get their distances and speed correlations
  fluctuationField = calculateFluctuationField(velocityField)
  forCalculation = velocityField %>% slice(sampleRows)
  vectorFluctuations = fluctuationField %>% slice(sampleRows)


  distances = getPairwiseDistances(forCalculation[, 1:2])
  speedCorrelations = getSpeedCorrelation(forCalculation[, 3:4])
  velocityCorrelations = getVelocityCorrelation(vectorFluctuations[, 3:4])
  directionCorrelations = getDirectionCorrelation(vectorFluctuations[, 3:4])

  pairwisecorrelations = data.frame(distances, velocityCorrelations, directionCorrelations, speedCorrelations)

  # From the pair-wise correlations, calculate the correlation profiles
  vCorrOut = normCorrFunction(pairwisecorrelations$distances, pairwisecorrelations$velocityCorrelations)
  dCorrOut = normCorrFunction(pairwisecorrelations$distances, pairwisecorrelations$directionCorrelations)
  sCorrOut = normCorrFunction(pairwisecorrelations$distances, pairwisecorrelations$speedCorrelations)

  corrOut = vCorrOut %>% mutate(directionCorr = dCorrOut[, -1], speedCorr = sCorrOut[, -1])
  colnames(corrOut) = c("domain", "vCorr", "dCorr", "sCorr")
  return(corrOut)
}



#' Isolate individual velocity fluctuations by removing collective
#'
#' This function takes in
#'
#' @param frameVectorField Takes in a N-by-4 velocity field where each vector consists of a point and its velocity
#' defined as (x, y, vx, vy).
#'
#' @return A N-by-4 velocity field that provides the relative position and fluctuation velocity of each point..
#'
#' @examples
#' velocityField = data.frame(x = runif(0,10,100), y = runif(0,10,100), vx = rnorm(100, 0, 2), vy = rnorm(100, 0, 2))
#' fluctField = calculateFluctuationField(velocityField)
calculateFluctuationField = function(frameVectorField) {
  # Function: takes a vector field. Determine each point's next position based on its respective vector.
  # Then, center both current and future positions at their respective centroid.
  colnames(frameVectorField) = c("X", "Y", "vX", "vY")
  positionsInterp = frameVectorField %>% mutate(X2 = X + vX, Y2 = Y + vY) %>%
    mutate(relX = X - mean(X), relY = Y - mean(Y), relX2 = X2 - mean(X2), relY2 = Y2 - mean(Y2)) %>%
    select(relX, relY, relX2, relY2)

  # For the centered points from now and future position, what is optimal transformation to reduce distance?
  # The start/null hypothesis is that there is no rotation (0) and no dilatation (1).
  transformParams = optim(c(0,1), measureDeviation, data = positionsInterp)
  rotationMatrix = transformParams$par[[1]]
  dilation = transformParams$par[[2]]

  # Apply optimal rotation and dilation on future points to get best alighnment, what is left is fluctuations.
  realignedPoints = affineTransform(positionsInterp[,3:4], rotate = rotationMatrix, dilate = dilation)

  # Subtract past positions from future positions post-transformation to get fluctuation vectors.
  alignedPositions = cbind(positionsInterp, realignedPoints) %>%
    mutate(fluctX = newx - relX, fluctY = newy - relY) %>%
    select(relX, relY, fluctX, fluctY)
  return(alignedPositions)
}


calculateDilatationalOrder = function(frameVectorField) {
  require(dplyr)
  center = colMeans(frameVectorField[,1:2])
  transformedField = frameVectorField %>% mutate(relX = X - center[1], relY = Y - center[2], speed = sqrt(vX^2 + vY^2)) %>%
    mutate(hypotCenter = sqrt(relX^2 + relY^2)) %>%
    mutate(unitX = relX/hypotCenter, unitY = relY/hypotCenter, uVX = vX/speed, uVY = vY/speed)

  dilatation = sum((transformedField$unitX * transformedField$uVX + transformedField$unitY * transformedField$uVY))/dim(transformedField)[1]
  return(dilatation)
}


#' Title
#'
#' @param transformParams
#' @param data
#'
#' @return
#' @export
#'
#' @examples
measureDeviation = function(transformParams, data) {
  transPoints = affineTransform(data[,3:4], transformParams[1], transformParams[2])

  sum(sqrt((data[[1]] - transPoints[[1]])^2 + (data[[2]] - transPoints[[2]])^2))
}



#' Title
#'
#' @param points
#' @param rotate
#' @param dilate
#'
#' @return
#' @export
#'
#' @examples
affineTransform = function(points,rotate, dilate) {
  points = data.frame(newx = points[[1]]*cos(rotate) + points[[2]]*sin(rotate), newy = points[[1]]*-sin(rotate)+points[[2]]*cos(rotate))
  points = points*dilate
  return(points)
}

#' Title
#'
#' @param thisFrame
#'
#' @return
#' @export
#'
#' @examples
removeTranslationalComponent = function(thisFrame) {
  ###FUNCTION: Taking a vector field as input, this function removes any mean translational motion captured by the vector field (u,v) and produces (relu,relv) and places the coordinates of the vector field (R,C) in a center-of-mass reference frame (relR,relC).

  require(dplyr)
  returnFrame = mutate(thisFrame,relR=R-mean(R),relC=C-mean(C),relu=u-mean(u),relv=v-mean(v))
  return(returnFrame)
}

#' Title
#'
#' @param thisFrame
#'
#' @return
#' @export
#'
#' @examples
removeRotationalComponent = function(thisFrame) {
  ###FUNCTION: This function removes the rotational movement (from the center-of-mass reference frame) of a vector field. This function assumes the object(s) represented by the vector field rotates as a solid object, such that the angular speed of rotation is linearly proportional to the radial distance from the center of mass.
  require(dplyr)
  #For each position, find the projected velocity that is perpendicular to radial vector
  #Calculate the speed/magnitude of this perpendicular vector.
  returnedFrame = mutate(thisFrame,mag = sqrt(relv^2+relu^2)) %>%
    mutate(magFromCenter = sqrt(relR^2+relC^2)) %>%
    mutate(unitrelR = relR/magFromCenter,unitrelC = relC/magFromCenter) %>%
    mutate(compS = (relR*relu+(-relC)*relv)/magFromCenter) %>%
    #Using the speed and radial distance, calculate the angular velocity for that point.
    mutate(angularS = compS/magFromCenter)

  #Determine the average angular velocity over all vectors.
  averageAngularS = summarise(returnedFrame,meanS=mean(angularS,na.rm = T))
  #Using the average angular velocity and radial distance, calculate the orthogonal vector to be subtracted from each vector.
  finalFrame = mutate(returnedFrame,tosubtu = averageAngularS$meanS*relR,tosubtv = -relC*averageAngularS$meanS) %>%

    mutate(fluctu = relu - tosubtu, fluctv = relv - tosubtv)
  return(finalFrame)
}



correlationSlope = function(distances, correlations) {
  # Calculates the slope of the correlation profile at the zero intercept, for the rescaled profile.
  # It does this using a very robust method: take closest five points, fit a quadratic, then derivative.
  # Based off of advice from http://www.theanalysisfactor.com/r-tutorial-4/
  quadraticFit = lm(correlations ~ distances + I(distances^2))
  definePoly = as.polynomial(quadraticFit$coefficients)
  derivativePoly = deriv(definePoly)
  slopeAtZero = predict(derivativePoly, 1)
  return(slopeAtZero)
}




#' Calculate susceptibility (maximum cumulative correlation)
#' Uses the method developed in Attanasi et al. 2014 to calculate the finite-size susceptibility or maximal cumulative correlation
#' of a collectively moving system. The cumulative correlation is calculated using the trapezoidal (Simpson's) rule for numerical integration.
#' This function filters the correlation profile to only consider values before the first zero-crossing if it detects any negative correlation values.
#'
#' @param distance The domain over which to integrate. If not provided, daata points are assumed to be evenly spaced with distance 1.
#' @param correlation The correlation at the corresponding distance.
#'
#' @return
#' @export
#'
#' @examples
calculateSusceptibility = function(distance = 1:length(correlation), correlation) {
  require(pracma)
  # Calculate susceptibility using trapezoidal rule of the curve.
  # Filter at the first zero crossing if not done already
  if (any(correlation < 0)) {
    zeroCrossing = getFirstZeroCrossing(correlation)
    distance = distance[1:zeroCrossing]
    correlation = correlation[1:zeroCrossing]
  }

  susceptibility = trapz(distance, correlation)
  return(susceptibility)
}



#' Title
#'
#' Calculates the pair-wise correlations between all vectors in a velocity field.
#' A vector can be provided to allow within- and between-group comparisons (not implemented).
#' Vectors can be subsampled to deal with quadratic memory complexity of calculating pair-wise correlations.
#'
#' @param frameData a N-by-5 data frame storing N velocity vectors with column names R, C, u, v, and Frame.
#' @param memoryMaxVectors the number of vectors to subset the data frame by in order to avoid memory issues.
#'
#' @return The pairwise velocity, direction, and speed correlation for all pairs of vectors found in frameData
#' @export
#'
#' @examples
getPairwiseCorrelations = function(frameData, memoryMaxVectors = 3000) {
  # From positions and full velocities, get the positions of each vector relative to the center and
  # velocities relative to the mean.
  colnames(frameData) = c("X", "Y", "vX", "vY")
  velocityFrame = frameData %>% mutate(speed = sqrt(vX^2+vY^2))
  velocityFrame = cbind(velocityFrame, calculateFluctuationField(velocityFrame[, 1:4]))

  velocityFrame = velocityFrame %>% select(relR,relC,speed,fluctu,fluctv) %>%
    mutate(unitfluctu = fluctu/sqrt(fluctu^2+fluctv^2), unitfluctv = fluctv/sqrt(fluctu^2+fluctv^2)) %>%
    mutate(speedFluctuation = speed - mean(speed))

  # Produce all index pairs, excluding repeat entries by only getting upper-triangle.
  numVectors = dim(velocityFrame)[1]
  pairwiseInds = expand.grid(1:numVectors, 1:numVectors)
  inds1 = matrix(pairwiseInds[,1], nrow = numVectors)
  inds2 = matrix(pairwiseInds[,2], nrow = numVectors)
  indsi = inds1[upper.tri(inds1,T)]
  indsj = inds2[upper.tri(inds2,T)]
  pairwiseInds = cbind(indsi,indsj)

  # Calculate Euclidean distances
  coordinates = cbind(velocityFrame[pairwiseInds[,1],1:2], velocityFrame[pairwiseInds[,2],1:2])
  distances = sqrt((coordinates[,1] - coordinates[,3])^2 + (coordinates[,2] - coordinates[,4])^2)

  velocityCorrelation = velocityFrame[pairwiseInds[,1], 4]*velocityFrame[pairwiseInds[,2], 4] + velocityFrame[pairwiseInds[,1], 5]*velocityFrame[pairwiseInds[,2], 5]

  directionCorrelation = velocityFrame[pairwiseInds[,1], 6]*velocityFrame[pairwiseInds[,2], 6] + velocityFrame[pairwiseInds[,1], 7]*velocityFrame[pairwiseInds[,2], 7]

  speedCorrelation = velocityFrame[pairwiseInds[,1], 8]*velocityFrame[pairwiseInds[,2], 8]

  pairwiseCorrelations = data.frame(distance = distances, velocityCorrelation = velocityCorrelation, directionalCorrelation = directionCorrelation, speedCorrelation = speedCorrelation)
  return(pairwiseCorrelations)
}

getPairwiseDistances = function(coordinates) {
  coordinates = as.matrix(coordinates)
  distances = as.matrix(dist(coordinates, diag=T))
  distances = distances[upper.tri(distances, diag=T)]
  return(distances)
}

getVelocityCorrelation = function(fluctuationVectors) {
  fluctuationVectors = as.matrix(fluctuationVectors)
  tFluctuationVectors = t(fluctuationVectors)
  velocityCorrelations = fluctuationVectors %*% tFluctuationVectors
  velocityCorrelations = velocityCorrelations[upper.tri(velocityCorrelations, diag=T)]
  return(velocityCorrelations)
}

getDirectionCorrelation = function(fluctuationVectors) {
  colnames(fluctuationVectors) = c("vX", "vY")
  unitVectors = fluctuationVectors %>% mutate(modulus = sqrt(vX^2 + vY^2)) %>%
    mutate(uvX = vX/modulus, uvY = vY/modulus) %>%
    select(uvX, uvY) %>% as.matrix()
  tUnitVectors = t(unitVectors)
  directionCorrelations = unitVectors %*% tUnitVectors
  directionCorrelations = directionCorrelations[upper.tri(directionCorrelations, diag=T)]
  return(directionCorrelations)
}

getSpeedCorrelation = function(velocityVectors) {
  colnames(velocityVectors) = c("vX", "vY")
  speedFlucts = velocityVectors %>% mutate(speed = sqrt(vX^2 + vY^2)) %>%
    mutate(speedFluct = speed - mean(speed)) %>% select(speedFluct) %>% as.matrix
  tSpeedFlucts = t(speedFlucts)
  speedCorrelations = speedFlucts %*% tSpeedFlucts
  speedCorrelations = speedCorrelations[upper.tri(speedCorrelations, diag=T)]
  return(speedCorrelations)
}


getPairwiseIndices = function(numVectors) {
  # Produce all index pairs, excluding repeat entries by only getting upper-triangle.
  pairwiseInds = expand.grid(1:numVectors, 1:numVectors)
  inds1 = matrix(pairwiseInds[,1], nrow = numVectors)
  inds2 = matrix(pairwiseInds[,2], nrow = numVectors)
  indsi = inds1[upper.tri(inds1,T)]
  indsj = inds2[upper.tri(inds2,T)]
  pairwiseInds = cbind(indsi,indsj)
}


#' Title
#'
#' @param frameData
#'
#' @return
#' @export
#'
#' @examples
normCorrFunction = function(domain,correlations, resolution=1, nknots = 50) {
  splineVal = bigspline(domain, correlations, type="cub", nknots)

  # Define prediction range and intervals
  maxDist = max(domain)
  predictDomain = seq(from = 0, to = maxDist, by = resolution)
  predictOut = predict(splineVal, predictDomain)
  normalizedOut = predictOut/predictOut[[1]]
  # plot(predictDomain, normalizedOut)
  # head(predictDomain[normalizedOut < 0], 1)
  # head(roughPredict$domain[roughPredict$meanCorr < 0], 1)
  corrDF = data.frame(domain = predictDomain, out = normalizedOut)
  return(corrDF)
}



#' Title
#'
#' @param thisVector
#' @param distMapping
#'
#' @return
#' @export
#'
#' @examples
getFirstZeroCrossing = function(thisVector,distMapping=1:length(thisVector)) {
  signDiff = diff(sign(thisVector))
  crossing = first(which(signDiff != 0))
  if (is.na(crossing))
    crossing=NA
  else {
    distances = distMapping[(crossing):(crossing+1)]
    values = thisVector[crossing:(crossing+1)]
    crossing = predict(lm(distances ~ values),data.frame(values=c(0)))
  }
  #Fit a linear model between the crossing point data line and interpolate the crossing point.
  return(crossing)
}



dilatationalOrderParameter = function(velocityFieldDF) {

}



identifyBoundary = function(X, Y, marginPercent, tryAlpha) {
  centered = cbind(X - mean(X), Y - mean(Y))
  jitterPoints = as.data.frame(apply(centered,2, jitter))

  hull = ahull(jitterPoints$V1, jitterPoints$V2, alpha=tryAlpha)

  inahull(hull, cbind(jitterPoints$V1/(1-marginPercent), jitterPoints$V2/(1-marginPercent)))
}

#' Title
#'
#' @param frameProfiles
#' @param offset
#'
#' @return
#' @export
#'
#' @examples
getCorrelationLengths = function(frameProfiles,offset=0) {
  # Takes in a three-column data frame consisting of a frame number, a correlation value, and a distance value.
  colnames(frameProfiles) = c("Frame","correlation","distance")
  allFrames = unique(frameProfiles$Frame)
  storage = rep(0,length(allFrames))
  zeroCrossings = data.frame(Frame=storage,CorrelationLength=storage)

  for (thisFrameIndex in 1:length(allFrames)) {
    thisFrame = allFrames[thisFrameIndex]
    zeroCrossings$Frame[thisFrameIndex] = thisFrame
    # Get the 0 crossing of that frame and save in a data frame.
    thisFrameCorrelations = filter(frameProfiles,Frame==thisFrame)
    zeroCrossings$CorrelationLength[thisFrameIndex] = getFirstZeroCrossing(thisFrameCorrelations$correlation,distMapping = thisFrameCorrelations$distance)
  }
  # Return the data frame.
  return(zeroCrossings)
}


loadSourceVideoTypes = function(masterDir) {
  coolsnapVideos = dirname(list.files(path = masterDir, pattern="^movie\\.avi$",recursive = T,full.names = T))
  hammamatsuVideos = dirname(list.files(path=masterDir, pattern="^hamm-movie\\.avi$", recursive=T, full.names=T))
  numVideos = length(coolsnapVideos) + length(hammamatsuVideos)
  videoSources = data.frame(folder=c(rep("",numVideos)),hammVideo=F, stringsAsFactors = F)
  videoSources$folder[1:length(coolsnapVideos)] = coolsnapVideos
  videoSources$folder[(length(coolsnapVideos)+1):numVideos] = hammamatsuVideos
  videoSources$hammVideo[(length(coolsnapVideos)+1):numVideos] = T
  return(videoSources)
}


loadDataByFolder = function(masterDir,filepattern) {
  filesList = list.files(path=masterDir, pattern=filepattern, recursive = T, full.names = T)
  for (thisFile in filesList) {
    tempFile = read.csv(thisFile)
    tempFile$folder = thisFile
    if (thisFile == filesList[1]) {
      returnData = tempFile
    } else {
      returnData = rbind(returnData,tempFile)
    }
    colnames(returnData) = colnames(tempFile)
  }
  return(returnData)
}


getCorrelationLengthTimeSeries = function(hdf5Data,forFrames) {
  require(doParallel)
  numFrames = length(forFrames)
  correlationLengths = data.frame(Frame = forFrames, velocityCorrelationLength = rep(NA,numFrames), directionalCorrelationLength = rep(NA,numFrames), speedCorrelationLength = rep(NA,numFrames))

  allVelocityFields = readVelocityFieldAtFrame(hdf5Identifier = hdf5Data,frameNumbers = forFrames)

  # Set up parallel cluster environment.
  cl = makeCluster(7)
  registerDoParallel(cl)
  clusterEvalQ(cl,expr = "library(rhdf5);library(dplyr)")
  listOfExportFunctions = c("calculateFluctuationField", "measureDeviation", "affineTransform", "getPairwiseCorrelations")
  clusterExport(cl,varlist = listOfExportFunctions)
  correlationProfiles = clusterApply(cl, x = allVelocityFields, fun = calculateCorrelationProfile)
  stopCluster(cl)

  allFrameCorrelation = bind_rows(correlationProfiles)
  vCorrLength = select(allFrameCorrelation,Frame,meanVelocityCorrelation,meanDist) %>% getCorrelationLengths(offset = 0)
  colnames(vCorrLength) = c("Frame","velocityCorrelationLength")
  dCorrLength = select(allFrameCorrelation,Frame,meanDirectionalCorrelation,meanDist) %>% getCorrelationLengths(offset = 0)
  colnames(dCorrLength) = c("Frame","directionCorrelationLength")
  sCorrLength = select(allFrameCorrelation,Frame,meanSpeedCorrelation,meanDist) %>% getCorrelationLengths(offset = 0)
  colnames(sCorrLength) = c("Frame","speedCorrelationLength")
  theseCorrelationLengths = inner_join(vCorrLength,dCorrLength,by="Frame") %>% inner_join(sCorrLength,by="Frame")
  return(list(theseCorrelationLengths,allFrameCorrelation))
}

getInactiveCells = function(currFrame, speedThreshold=0.01) {
  currFrame = currFrame %>% mutate(speed = sqrt(u^2+v^2)) %>%
    mutate(inactiveCell = speed < speedThreshold)
}

