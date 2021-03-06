% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InternalOrderMeasureFunctions.R
\name{measureInternalOrder}
\alias{measureInternalOrder}
\title{Measure the descriptive statistical properties of a collective.}
\usage{
measureInternalOrder(dataFolder, frameChunkSize = 3600, dplyr = FALSE)
}
\arguments{
\item{dataFolder}{The folder containing the collective order dataset to use}

\item{frameChunkSize}{The number of frames to load into memory at any given time. Useful for batch read-write (memory-dependent)}

\item{dplyr}{Is the dplyr package available?Can use it to speed up processing.}
}
\value{
A data frame with the frame number, polarization, rotation, mean speed, mean v, mean u, and the max speed for each frame.
}
\description{
Measure the descriptive statistical properties of a collective.
}
