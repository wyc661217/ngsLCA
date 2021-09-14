#'
#' @title Generate Taxa Abundance Barplots
#'
#' @description Generate barplots for taxa profiles into "path/run/barplot/" to show the reads abundance.
#'
#' @param path working directory, same to \code{\link{ngsLCA_profile}}.
#' @param run name of the run, default is "run01".
#' @param taxa.number maximum number of taxa will be shown in barplots; default is 10.
#'
#' @import ggplot2
#' @importFrom reshape melt
#' @importFrom utils read.csv
#' @export
#'
#' @examples
#' ngsLCA_barplot(path=system.file("extdata","lca_files",package="ngsLCA"),
#'                run="run01",
#'                taxa.number=10)
#'
#'
#' ## This will generate barplots for the complete taxa profile
#' ## as well as taxa groups and ranks (if available) of "run01"
#' ## in "path/run01/barplot/", with showing the 10 most abundant
#' ## taxa in each profile.
#'
#'
ngsLCA_barplot=function(path,
                        run="run01",
                        taxa.number=10){


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
    cat("\n\n\t-> Generate barplot.\n\n")
  }


  #create folder
  if (length(dir(paste(path, run, sep = ""), pattern = "barplot")) > 0) {
    cat(paste("\n\n\t-> '", path, run, "/barplot' already exists; files inside will be overwritten.\n\n", sep = ""))
  }else{
    dir.create(paste(path, run, "/barplot", sep=""))
  }


  #local variables
  variable = value = taxa = NULL


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


  #complete taxa profile
  X1 = read.csv(paste(path, run, "/taxonomic_profiles/complete_profile.txt", sep=""),
                sep="\t", quote="", check.names=F, stringsAsFactors=F)
  X2 = dataPrep(X1)

  for (i in 1:dim(X2)[2]) {
    if (sum(X2[,i])>0) {
      X2[,i] = X2[,i]/sum(X2[,i])
    }
  }

  X2$sum = rowSums(X2)
  X2 = X2[order(-X2$sum),]

  if (dim(X2)[1] > taxa.number) {
    X2 = X2[1:taxa.number,-dim(X2)[2]]
  } else{
    X2 = X2[,-dim(X2)[2]]
  }

  sample_list = colnames(X2)
  X2$taxa = rownames(X2)
  X2 = melt(X2,"taxa")
  X2$variable = factor(X2$variable,sample_list,ordered = T)

  p0=ggplot(X2, aes(x=variable, y=value,fill =taxa)) +
    geom_bar(stat="identity")+
    labs(y="Proportion (%)")+
    theme(axis.title.x = element_blank(),
          plot.title = element_text(size=17, face="bold",hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          axis.title.y = element_text(size = 17),
          axis.text.x = element_text(size = 15, colour = "black" ,angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(size = 15, colour = "black"),
          panel.grid.major = element_blank(),
          axis.ticks.length=unit(.25, "cm"),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black",fill = "transparent",size=1.5),
          legend.position="bottom")

  ggsave(plot=p0, height=10, width=15, dpi=300,
         filename=paste(path, run, "/barplot/all_taxa_barplot.pdf", sep=""),
         useDingbats=FALSE)


  #taxa ranks
  if (length(dir(paste(path, run, "/taxonomic_profiles/taxa_ranks",sep=""), pattern = ".txt"))>0){

    file.list = dir(paste(path, run, "/taxonomic_profiles/taxa_ranks",sep=""), pattern = ".txt")
    for (i in 1:length(file.list)) {
      X1 = read.csv(paste(path, run, "/taxonomic_profiles/taxa_ranks/", file.list[i], sep=""),
                    sep="\t", quote="",check.names=F,stringsAsFactors=F)
      X2 = dataPrep(X1)

      for (j in 1:dim(X2)[2]) {
        if (sum(X2[,j])>0) {
          X2[,j] = X2[,j]/sum(X2[,j])
        }
      }

      X2$sum = rowSums(X2)
      X2 = X2[order(-X2$sum),]

      if (dim(X2)[1] > taxa.number) {
        X2 = X2[1:taxa.number,-dim(X2)[2]]
      } else{
        X2 = X2[,-dim(X2)[2]]
      }

      sample_list = colnames(X2)
      X2$taxa = rownames(X2)
      X2 = melt(X2,"taxa")
      X2$variable = factor(X2$variable,sample_list,ordered = T)

      Name = sub(".txt", "", file.list[i])

      p1 = ggplot(X2, aes(x=variable, y=value,fill =taxa)) +
        geom_bar(stat="identity")+
        labs(y="Proportion (%)")+
        theme(axis.title.x = element_blank(),
              plot.title = element_text(size=17, face="bold",hjust = 0.5),
              legend.title = element_blank(),
              legend.text = element_text(size = 15),
              axis.title.y = element_text(size = 17),
              axis.text.x = element_text(size = 15, colour = "black" ,angle = 45,hjust = 1,vjust = 1),
              axis.text.y = element_text(size = 15, colour = "black"),
              panel.grid.major = element_blank(),
              axis.ticks.length=unit(.25, "cm"),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black",fill = "transparent",size=1.5),
              legend.position="bottom")

      ggsave(plot=p1, height=10, width=15, dpi=300,
             filename=paste(path, run, "/barplot/", Name, "_barplot.pdf",sep=""),
             useDingbats=FALSE)
    }
  }
}
