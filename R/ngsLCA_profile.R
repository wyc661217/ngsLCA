#'
#' @title Generate a combined taxa profile
#'
#' @description Read lca files in the appointed directory to generate a combined taxa profile named 'complete_profile.txt' in folder 'path/run/taxonomic_profiles/'.
#'
#' @param path working directory, which should contain lca files that will be processed.
#' @param run run name, a folder that will be created for storing generated files and results; default is "run01".
#' @param remove.sample a comma separated vector listing file names indicating samples that will NOT be included for analysis, optional.
#' @param metadata path to sample metadata, an example at github.com/miwipe/ngsLCA/blob/master/R_script/metadata.txt, optional.
#'
#' @importFrom tidyr separate %>%
#' @importFrom stats aggregate
#' @importFrom utils read.csv write.table
#' @export
#'
#' @examples
#' ngsLCA_profile(path="working_directory/",
#'                run="run01",
#'                remove.sample="sample1.lca,sample2",
#'                metadata="path_to_metadata/metadata.txt")
#'
#' ## This will combine all lca files in "working_directory/"
#' ## to generate a file named 'complete_profile.txt' in
#' ## 'working_directory/run01/taxonomic_profiles/', with
#' ## replacing the file names by sample names supplied in
#' ## metadata.txt, and removing sample1 and sample2.
#'
#'
ngsLCA_profile = function(path,
                          run="run01",
                          remove.sample=NULL,
                          metadata=NULL){


  #modify input path
  if (!grepl('/$', path)) {
    path = paste(path,"/",sep="")
  }


  #local variable
  taxa = NULL


  #list lca files
  FileList = dir(path, pattern = ".lca$")

  if(length(FileList)==0){
    cat("\n\n\t-> Error: no lca files found in the appointed working directory.\n\n")
    stop()
  } else{
    cat("\n\n\t-> Data read-in and pre-process, may take some time depending on number and size of the input files.\n\n")
  }


  #read in lca files
  DF1 = data.frame(taxa=character(), stringsAsFactors=F)

  for (i in 1:length(FileList)) {

    DF2.1 =  read.csv(paste(path, FileList[i], sep=""), header=F, sep="\t", stringsAsFactors=F, fill=T,
                     col.names = paste0("V",seq_len(60)), skip=1)

    if(dim(DF2.1)[1]>0){

      DF2.2 = data.frame(taxa=DF2.1[,2],
                        count=rep(1, dim(DF2.1)[1]),
                        stringsAsFactors = F)

      DF2.3 = aggregate(DF2.2[,2]~DF2.2[,1], data=DF2.2, FUN=sum)
      colnames(DF2.3) = c("taxa",sub(".lca","",FileList[i]))
      DF1 = merge(DF1, DF2.3, by="taxa", all=T)

    }else{
      paste("\n\n\t-> ", FileList[i], " contains no taxa, skipped.\n\n", sep = "")
    }
  }

  DF1[is.na(DF1)] = 0


  #removes samples in the sample removing list
  if (!is.null(remove.sample)) {
    remove.sample = strsplit(remove.sample,",")[[1]]
    remove.sample = sub(".lca", "", remove.sample)
    DF1 = DF1[,!(colnames(DF1) %in% remove.sample)]
  }


  #replacing the file names with metadata
  if(is.null(metadata)){
    cat("\n\n\t-> Metadata file not supplied; file names will be illustrated.\n\n")
  } else {
    cat("\n\n\t-> Metadata supplied and will be used.\n\n")

    Mdata = read.csv(metadata, quote="", stringsAsFactors=F, header=F,sep = "\t",comment.char = "#")
    Mdata$V1 = sub(".lca","",Mdata$V1)
    Mdata = Mdata[order(Mdata$V3),]

    N = 0
    for (i in 1:dim(Mdata)[1]) {
      M = which(colnames(DF1)[-1] %in%  Mdata$V1[i])
      if (length(M)>0) {
        colnames(DF1)[M+1] = Mdata$V2[i]
        N = c(N,M+1)
      }
    }

    N = N[-1]
    M = which(!2:dim(DF1)[2] %in% N) + 1
    DF1 = DF1[,c(1,N,M)]
  }


  #write taxa profile
  if (dim(DF1)[1]>0) {

    if (length(dir(path,pattern = run)) > 0) {
      cat(paste("\n\n\t-> '", path, run, "' already exists; files inside will be overwritten.\n\n",sep = ""))
    }else{
      dir.create(paste(path, run, sep=""))
      dir.create(paste(path, run, "/intermediate", sep=""))
      dir.create(paste(path, run, "/taxonomic_profiles", sep=""))
    }

    write.table(DF1, file = paste(path, run, "/intermediate/", "taxa_profile_v1.txt", sep=""), quote=F,
                row.names=F, col.names=T, sep="\t")
    write.table(DF1, file = paste(path, run, "/intermediate/", "taxa_profile_final.txt", sep=""), quote=F,
                row.names=F, col.names=T, sep="\t")

    DF3 = DF1 %>% separate(taxa, sep=":", c("taxa_ID","taxa_name","taxonomic_rank"))
    write.table(DF3, file = paste(path, run, "/taxonomic_profiles/", "complete_profile.txt", sep=""), quote=F,
                row.names=F, col.names=T, sep="\t")

  }else{
    cat("\n\n\t-> ERROR: no alignments appear in the input lca files.\n\n")
    stop()
  }
}
