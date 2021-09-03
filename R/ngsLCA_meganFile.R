#'
#' @title Generate Megan files
#'
#' @description Generate files that can be inputted into Megan for taxa profiling. Results will be in 'path/run/megan/'
#'
#' @param path working directory, same to \code{\link{ngsLCA_profile}}.
#' @param run name of the run, default is "run01".
#'
#' @importFrom utils read.csv write.table
#' @export
#'
#' @examples
#' ngsLCA_meganFile(path=system.file("extdata","lca_files",package="ngsLCA"),
#'                  run="run01")
#'
#'
#' ## This will generate Megan format taxa profiles for
#' ## the complete taxa profile as well as taxa groups
#' ## and ranks (if available) of 'run01'. Resutls will
#' ## be in extdata/lca_files/run01/megan/'. These files
#' ## can be directly used as input csv files in Megan
#' ## for taxa profiling.
#'
#'
ngsLCA_meganFile=function(path,
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
    cat("\n\n\t-> Generate Megan output.\n\n")
  }


  #create folder
  if (length(dir(paste(path, run, sep = ""), pattern = "megan")) > 0) {
    cat(paste("\n\n\t-> '", path, run, "/megan' already exists; files inside will be overwritten.\n\n", sep = ""))
  }else{
    dir.create(paste(path, run, "/megan", sep=""))
  }


  #complete profile
  X1 = read.csv(paste(path, run, "/taxonomic_profiles/complete_profile.txt", sep=""),
                sep="\t", quote="", check.names=F, stringsAsFactors=F)

  X1 = X1[,-c(2,3)]
  colnames(X1)[1] = "#dataset"
  write.table(X1, file = paste(path, run, "/megan/all_taxa.txt", sep=""), quote=F, row.names=F, col.names=T, sep=",")


  #taxa groups
  if (length(dir(paste(path, run, "/taxonomic_profiles/taxa_groups",sep=""), pattern = ".txt"))>0){

    file.list = dir(paste(path, run, "/taxonomic_profiles/taxa_groups",sep=""), pattern = ".txt")

    for (i in 1:length(file.list)) {
      X1 = read.csv(paste(path, run, "/taxonomic_profiles/taxa_groups/", file.list[i], sep=""),
                    sep="\t", quote="",check.names=F,stringsAsFactors=F)
      X1 = X1[,-c(2,3)]
      colnames(X1)[1] = "#dataset"
      write.table(X1, file = paste(path, run, "/megan/", file.list[i], sep=""), quote=F, row.names=F, col.names=T, sep=",")
    }
  }


  #taxa ranks
  if (length(dir(paste(path, run, "/taxonomic_profiles/taxa_ranks",sep=""), pattern = ".txt"))>0){

    file.list = dir(paste(path, run, "/taxonomic_profiles/taxa_ranks",sep=""), pattern = ".txt")
    for (i in 1:length(file.list)) {
      X1 = read.csv(paste(path, run, "/taxonomic_profiles/taxa_ranks/", file.list[i], sep=""),
                    sep="\t", quote="",check.names=F,stringsAsFactors=F)
      X1 = X1[,-c(2,3)]
      colnames(X1)[1] = "#dataset"
      write.table(X1, file = paste(path, run, "/megan/", file.list[i], sep=""), quote=F, row.names=F, col.names=T, sep=",")
    }
  }
}
