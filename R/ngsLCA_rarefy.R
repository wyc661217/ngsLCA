#'
#' @title Rarefaction taxa reads abundance
#'
#' @description Generate rarefaction curves for taxa profiles in 'path/run/rarefaction/'. See vegan::rarefy for details.
#'
#' @param path working directory, same to \code{\link{ngsLCA_profile}}.
#' @param run name of the run, default is "run01".
#'
#' @importFrom vegan rarefy specnumber rarecurve
#' @importFrom utils read.csv
#' @importFrom graphics abline
#' @export
#'
#' @examples
#' ngsLCA_rarefy(path=system.file("extdata","lca_files",package="ngsLCA"),
#'               run="run01")
#'
#' ## This will generate rarefaction curves for the
#' ## complete taxa profile as well as taxa groups
#' ## and ranks (if available) of 'run01'. Results
#' ## will be in 'extdata/lca_files/run01/rarefaction/'.
#'
#'
ngsLCA_rarefy=function(path,
                       run="run01"){


  #modify input path
  if (!grepl('/$', path)) {
    path = paste(path,"/",sep="")
  }


  #whether ngsLCA_profile performed
  if (length(dir(paste(path, run, "/intermediate/", sep=""), pattern =  "taxa_profile_v1.txt")) == 0) {
    cat(paste("\n\n\t-> Error: required input file not found under '", run, "' folder, please run 'ngsLCA_profile' first.\n\n",
              sep = ""))
    stop()
  }else{
    cat("\n\n\t-> Perform rarefaction.\n\n")
  }


  #Function for rarefaction
  Raref = function(DF, OutName){

    if (pmin(dim(DF)[1],dim(DF)[2]) >= 2) {

      DF1 = t(DF)
      S = specnumber(DF1) # observed number of taxa
      raremax = min(rowSums(DF1))
      Srare = rarefy(DF1, raremax)
      #
      pdf(paste(path, run, "/rarefaction/", OutName, "_rarefy_scaling.pdf", sep=""), width=12, height=8)
      plot(S, Srare, xlab = "Observed No. of taxa", ylab = "Rarefied No. of taxa")
      abline(0, 1)
      dev.off()
      #
      pdf(paste(path, run, "/rarefaction/", OutName, "_rarefy_rarecurve.pdf", sep=""), width=12, height=8)
      rarecurve(DF1, step = 20, sample = raremax, col = "blue", cex = 1.0, ylab = "Observed taxa")
      dev.off()
    } else {
      cat(paste("\n\n\t-> For ", OutName, ", number of samples or taxa is less than 2, insufficient for rarefaction.\n\n",sep = ""))
    }
  }


  #function for preparing data
  dataPrep = function(DF){

    if (length(which(duplicated(DF[,2])))>0) {
      N=which(DF[,2] %in% DF[which(duplicated(DF[,2])),2])
      DF[N,2] = paste(DF[N,2],DF[N,3],sep="_")
      if (length(which(duplicated(DF[,2])))>0) {
        N=which(DF[,2] %in% DF[which(duplicated(DF[,2])),2])
        DF[N,2] = paste(DF[N,1],DF[N,2],sep="_")
      }
    }

    row.names(DF) = DF[,2]
    DF2 = as.data.frame(DF[,-c(1:3)])
    colnames(DF2) = colnames(DF)[-c(1:3)]
    rownames(DF2) = rownames(DF)

    return(DF2)
  }


  #create folder
  if (length(dir(paste(path, run, sep = ""), pattern = "rarefaction")) > 0) {
    cat(paste("\n\n\t-> '", path, run, "/rarefaction' already exists; files inside will be overwritten.\n\n", sep = ""))
  }else{
    dir.create(paste(path, run, "/rarefaction", sep=""))
  }


  #rarefy complete taxa profile
  X1 = read.csv(paste(path, run, "/taxonomic_profiles/complete_profile.txt", sep=""),
                sep="\t", quote="", check.names=F, stringsAsFactors=F)
  X2 = dataPrep(X1)
  Raref(DF=X2, OutName="all_taxa")


  #taxa ranks
  if (length(dir(paste(path, run, "/taxonomic_profiles/taxa_ranks",sep=""), pattern = ".txt"))>0){

    file.list = dir(paste(path, run, "/taxonomic_profiles/taxa_ranks",sep=""), pattern = ".txt")

    for (i in 1:length(file.list)) {
      X1 = read.csv(paste(path, run, "/taxonomic_profiles/taxa_ranks/", file.list[i], sep=""),
                    sep="\t", quote="",check.names=F, stringsAsFactors=F)
      X2 = dataPrep(X1)

      Name = sub(".txt", "", file.list[i])
      Raref(DF=X2, OutName=Name)
    }
  }
}
