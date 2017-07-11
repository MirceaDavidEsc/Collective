# This is a series of tests to use to make sure the measurements are working correctly.
library(readr)
library(dplyr)
testRotation = function() {
  rotationField = read_rds("rotate.rds")
  rotation = measureRotation(rotationField)
  polarization = measurePolarization(rotationField)
  dilatation = measureDilatation(rotationField)
  quiverPlot(rotationField,1)
  return(abs(rotation) - polarization - abs(dilatation) > 0.8)
}


testPolarization = function() {
  polarizationField = read_rds("translate.rds")
  rotation = measureRotation(polarizationField)
  polarization = measurePolarization(polarizationField)
  dilatation = measureDilatation(polarizationField)
  quiverPlot(polarizationField,1)
  return(polarization - abs(rotation) - abs(dilatation) > 0.95)
}


testDilatation = function() {
  dilationField = read_rds("dilate.rds")
  rotation = measureRotation(dilationField)
  polarization = measurePolarization(dilationField)
  dilatation = measureDilatation(dilationField)
  quiverPlot(dilationField,1)
  return(abs(dilatation) - abs(rotation) - polarization > 0.95)
}


testReal = function() {
  realField = read_rds("real.rds") %>% select(-Frame)
  rotation = measureRotation(realField)
  polarization = measurePolarization(realField)
  dilatation = measureDilatation(realField)
  quiverPlot(realField,1)
  return(c(rotation,polarization,dilatation))
}

testRotation()
testPolarization()
testDilatation()
testReal()
