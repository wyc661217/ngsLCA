% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ngsLCA_rank.R
\name{ngsLCA_rank}
\alias{ngsLCA_rank}
\title{Classify Taxa to Taxonomic Ranks}
\usage{
ngsLCA_rank(path, run = "run01", rank.name = "species,genus,family")
}
\arguments{
\item{path}{working directory, same to \code{\link{ngsLCA_profile}}.}

\item{run}{name of the run, default is "run01".}

\item{rank.name}{a comma separated vector listing the taxonomic ranks that will be used for classifying taxa profiles; default is "species,genus,family"}
}
\value{
Taxa profiles clustered into taxa ranks.
}
\description{
Classify the combined taxa profile (and grouped taxa profiles generated by \code{\link{ngsLCA_group}} if available) into user-defined taxonomic ranks. Results will be in "path/run/taxonomic_profiles/taxa_ranks/".
}
\examples{
ngsLCA_rank(path=system.file("extdata","lca_files",package="ngsLCA"),
            run="run01",
            rank.name="species,genus")


## This will classify the combined taxa profile (and
## grouped taxa profiles if available) of "run01" into
## species and genus, by merging all taxa below a species
## into that species, and all taxa below a genus into
## that genus. Generated files will be in
## "path/run01/taxonomic_profiles/taxa_ranks/".


}
