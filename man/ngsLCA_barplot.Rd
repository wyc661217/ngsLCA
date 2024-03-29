% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ngsLCA_barplot.R
\name{ngsLCA_barplot}
\alias{ngsLCA_barplot}
\title{Generate Taxa Abundance Barplots}
\usage{
ngsLCA_barplot(path, run = "run01", taxa.number = 10)
}
\arguments{
\item{path}{working directory, same to \code{\link{ngsLCA_profile}}.}

\item{run}{name of the run, default is "run01".}

\item{taxa.number}{maximum number of taxa will be shown in barplots; default is 10.}
}
\value{
Taxa abundance barplots.
}
\description{
Generate barplots for taxa profiles into "path/run/barplot/" to show the reads abundance.
}
\examples{
ngsLCA_barplot(path=system.file("extdata","lca_files",package="ngsLCA"),
               run="run01",
               taxa.number=10)


## This will generate barplots for the complete taxa profile
## as well as taxa groups and ranks (if available) of "run01"
## in "path/run01/barplot/", with showing the 10 most abundant
## taxa in each profile.


}
