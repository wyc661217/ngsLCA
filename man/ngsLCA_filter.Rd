% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ngsLCA_filter.R
\name{ngsLCA_filter}
\alias{ngsLCA_filter}
\title{Filter the Combined Taxa Profile}
\usage{
ngsLCA_filter(
  path,
  run = "run01",
  remove.taxa = NULL,
  threshold.1 = 2,
  threshold.2 = 0,
  threshold.3 = 5
)
}
\arguments{
\item{path}{working directory, same to \code{\link{ngsLCA_profile}}.}

\item{run}{name of the run, default is "run01".}

\item{remove.taxa}{a list of NCBI taxaID indicating taxa that will be removed. This can be a comma separated vector or the path to a text file listing one taxaID in each line.}

\item{threshold.1}{minimum reads number required for confirming a taxon in each sample; default is 2.}

\item{threshold.2}{minimum read percentage (to the total reads number of a given sample) required for confirming a taxon in each sample, ranging from 0 to 1; default is 0.}

\item{threshold.3}{minimum sum of reads across all samples required for confirming a taxon in the combined taxa profile; default is 5.}
}
\value{
An updated complete_profile.txt (after taxa filters).
}
\description{
Filter the combined taxa profile. Running this function will update the existing "complete_profile.txt".
}
\examples{
ngsLCA_filter(path=system.file("extdata","lca_files",package="ngsLCA"),
              run="run01",
              remove.taxa="1115501,10114",
              threshold.1=2,
              threshold.2=0,
              threshold.3=2)


## This will filter and update the combined taxa profile:
## "path/run01/taxonomic_profiles/complete_profile.txt"
## by the 3 thresholds supplied, and by removing taxa
## with taxaID 1115501 and 10114.


}
