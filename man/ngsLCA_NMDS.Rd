% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ngsLCA_NMDS.R
\name{ngsLCA_NMDS}
\alias{ngsLCA_NMDS}
\title{Taxa profiles NMDS}
\usage{
ngsLCA_NMDS(path, run = "run01", dimension = 3, trymax = 100)
}
\arguments{
\item{path}{working directory, same to 'ngsLCA_profile'.}

\item{run}{name of the run, default is "run01".}

\item{dimension}{dimension of the NMDS; default is 3.}

\item{trymax}{maximum number of random starts in search of convergent solutions for NMDS; default is 100.}
}
\value{

}
\description{
Perform Non-metric Multidimensional Scaling (NMDS) on taxa profiles. Results will be in 'path/run/NMDS/'. See vegan::metaMDS for details.
}
\examples{
ngsLCA_NMDS(path="working_directory/",
            run="run01",
            dimension=3,
            trymax=1000)


## This will perform NMDS on the complete taxa
## profile as well as taxa groups and ranks
## (if available) of run01'. Results will be in
## 'working_directory/run01/NMDS/'.


}