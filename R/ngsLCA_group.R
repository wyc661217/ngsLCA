#'
#' @title Group taxa in the combined taxa profile
#'
#' @description Group the combined taxa profile into user-defined taxa units. Results will be in 'pth/run/taxonomic_profiles/taxa_groups/'.
#'
#' @param path working directory, same to \code{\link{ngsLCA_profile}}.
#' @param run name of the run, default is "run01".
#' @param group.name a comma seprated vector listing the taxonomic units that will be used for grouping taxa, need to be as the scientific names of NCBI taxonomy (https://www.ncbi.nlm.nih.gov/taxonomy); default is "Viruses,Archaea,Bacteria,Fungi,Viridiplantae,Metazoa".
#' @param threshold.perGroup minimum reads percentage (to the total reads number of each group) required for confirming a taxon in each group of each sample, ranging from 0 to 1.
#'
#' @importFrom tidyr separate %>%
#' @importFrom utils read.csv write.table
#' @export
#'
#' @examples
#' ngsLCA_group(path=system.file("extdata","lca_files",package="ngsLCA"),
#'              run="run01",
#'              group.name="Viridiplantae,Metazoa",
#'              threshold.perGroup=0.01)
#'
#' ## This will extract all taxa under Viridiplantae and Metazoa
#' ## from the combined taxa profile of 'run01', and filter
#' ## these 2 taxa profiles by threshold.perGroup to generate 2
#' ## files named "Viridiplantae.txt" and "Metazoa.txt" in
#' ## "extdata/lca_files/run01/taxonomic_profiles/taxa_groups/".
#'
#'
ngsLCA_group = function(path,
                        run="run01",
                        group.name="Viruses,Archaea,Bacteria,Fungi,Viridiplantae,Metazoa",
                        threshold.perGroup = 0){


  #modify input path
  if (!grepl('/$', path)) {
    path = paste(path,"/",sep="")
  }


  #local variable
  taxa = NULL


  #whether ngsLCA_profile performed
  if (length(dir(paste(path, run, "/intermediate/", sep=""),pattern =  "taxa_profile_v1.txt")) == 0) {
    cat(paste("\n\n\t-> Error: required input file not found under '", run, "' folder, please run 'ngsLCA_profile' first.\n\n",
              sep = ""))
    stop()
  }else{
    cat("\n\n\t-> Group taxa.\n\n")
  }


  #generate taxa branch
  FileList = dir(path, pattern = ".lca$")

  if(length(FileList)==0){
    cat("\n\n\t-> Error: no lca files found in the appointed working directory.\n\n")
    stop()
  }

  PATH = read.csv(paste(path, FileList[1], sep=""), header=F, sep="\t", stringsAsFactors=F,
                  fill=T, col.names = paste0("V",seq_len(60)),comment.char = "#")
  PATH = PATH[,-1]

  if (length(FileList)>1) {
    for (i in 2:length(FileList)) {
      PATH.1 = read.csv(paste(path, FileList[i], sep=""), header=F, sep="\t", stringsAsFactors=F,
                        fill=T, col.names = paste0("V",seq_len(60)),comment.char = "#")
      PATH.1 = PATH.1[,-1]
      PATH = rbind(PATH,PATH.1)
      PATH = PATH[!duplicated(PATH[,1]),]
    }
  }else{
    PATH = PATH[!duplicated(PATH[,1]),]
  }

  write.table(PATH, file = paste(path, run, "/intermediate/", "taxa_branch.txt", sep=""), quote=F,
              row.names=F, col.names=T, sep="\t")


  #group and write files
  TaxaUnit = strsplit(group.name,",")[[1]]
  thrG=as.numeric(threshold.perGroup)
  tax.branch = PATH
  DF = read.csv(paste(path, run,"/intermediate/", "taxa_profile_final.txt", sep=""),
                 sep="\t", quote="", check.names=F, stringsAsFactors=F)

  if (length(dir(paste(path, run, "/taxonomic_profiles/", sep = ""), pattern = "taxa_groups")) > 0) {
    cat(paste("\n\n\t-> '", path, run, "/taxonomic_profiles/taxa_groups' already exists; files inside will be overwritten.\n\n",sep = ""))
  }else{
    dir.create(paste(path, run, "/taxonomic_profiles/taxa_groups", sep=""))
    dir.create(paste(path, run, "/intermediate/taxa_groups", sep=""))
  }

  L1 = NULL
  L1[TaxaUnit] = list(NULL)

  for (i in 1:dim(DF)[1]) {
    for (j in 1:length(TaxaUnit)) {
      if (length(grep(TaxaUnit[j], tax.branch[which(tax.branch[,1] %in%  DF[i,1]),])) != 0) {
        L1[[j]] = c(L1[[j]],i)
      }
    }
  }

  name = data.frame(taxa=names(L1))

  for (i in 1:length(L1)) {
    if (length(L1[[i]]) != 0) {

      DF1 = DF[L1[[i]],]

      if (dim(DF1)[2]>2) {
        for (j in 1:dim(DF1)[1]) {
          DF1[j, which(DF1[j,-1] <= colSums(DF1[,-1]) * thrG)+1] = 0
        }
      }else{
        for (j in 1:dim(DF1)[1]) {
          DF1[j, which(DF1[j,-1] <= sum(DF1[,-1]) * thrG)+1] = 0
        }
      }

      DF2 = DF1[rowSums(as.data.frame(DF1[,-1])) != 0,]
      write.table(DF2, file = paste(path, run, "/intermediate/taxa_groups/", name$taxa[i], ".txt", sep=""), quote=F,
                  row.names=F, col.names=T, sep="\t")

      DF2.1 = DF2 %>% separate(taxa, sep=":", c("taxa_ID","taxa_name","taxonomic_rank"))
      write.table(DF2.1, file = paste(path, run, "/taxonomic_profiles/taxa_groups/", name$taxa[i], ".txt", sep=""),
                  quote=F,row.names=F, col.names=T, sep="\t")
    }else{
      cat(paste("\n\n\t->  No taxa grouped into ", name$taxa[i], ".\n\n", sep = ""))
    }
  }
}
