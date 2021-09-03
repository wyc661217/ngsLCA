#'
#' @title Taxa profiles NMDS
#'
#' @description Perform Non-metric Multidimensional Scaling (NMDS) on taxa profiles. Results will be in 'path/run/NMDS/'. See vegan::metaMDS for details.
#'
#' @param path working directory, same to \code{\link{ngsLCA_profile}}.
#' @param run name of the run, default is "run01".
#' @param dimension dimension of the NMDS; default is 3.
#' @param trymax maximum number of random starts in search of convergent solutions for NMDS; default is 100.
#'
#' @importFrom vegan metaMDS stressplot
#' @importFrom  grDevices pdf dev.off
#' @importFrom utils read.csv
#' @import ggplot2
#' @export
#'
#' @examples
#' ngsLCA_NMDS(path=system.file("extdata","lca_files",package="ngsLCA"),
#'             run="run01",
#'             dimension=3,
#'             trymax=1000)
#'
#'
#' ## This will perform NMDS on the complete taxa
#' ## profile as well as taxa groups and ranks
#' ## (if available) of run01'. Results will be in
#' ## 'extdata/lca_files/run01/NMDS/'.
#'
#'
ngsLCA_NMDS=function(path,
                     run="run01",
                     dimension = 3,
                     trymax = 100){

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
    cat("\n\n\t-> Perform NMDS; will take some time depending on input file size and the trymax.\n\n")
  }


  if (dimension < 2){
    cat("\n\n\t-> Error: dimension must be equal or greater than 2.\n\n")
    stop()
  }


  #create folder
  if (length(dir(paste(path, run, sep = ""), pattern = "NMDS")) > 0) {
    cat(paste("\n\n\t-> '", path, run, "/NMDS' already exists; files inside will be overwritten.\n\n", sep = ""))
  }else{
    dir.create(paste(path, run, "/NMDS", sep=""))
  }


  #function for NMDS
  NMDSc = function(DF, OutName){

    if (dim(DF)[1] <= dimension | dim(DF)[1] <= dimension) {
      cat(paste("\n\n\t-> For ", OutName, ", number of samples or taxa is less than or equal to the NMDS dimensions; NMDS not performed.\n\n",sep = ""))
    }else{

      MDS1 = MDS2 = NULL
      NMDS = metaMDS(DF, k=dimension, trymax=trymax) #NMDS cluster
      saveRDS(NMDS, file = paste(path, run, "/NMDS/", OutName, "_NMDS_data.rda", sep=""))

      #output NMDS_stressplot figure
      pdf(paste(path, run, "/NMDS/", OutName, "_NMDS_stressplot.pdf", sep=""), width=12, height=12)
      stressplot(NMDS)
      dev.off()

      #output NMDS figure
      X3 = as.data.frame(NMDS$points)[,c(1,2)]
      X4 = as.data.frame(NMDS$species)[,c(1,2)]

      p1 = ggplot(X4, aes(MDS1, MDS2))+
        geom_point(size=4,colour="#7F7F7F")+
        geom_text(aes(label = rownames(X4)),colour="black",size=5)+
        labs(title=paste(OutName, "_NMDS", sep=""),
             y="NMDS2", x="NMDS1") +
        theme(axis.title.x = element_text(size = 17),
              title = element_text(size = 19, face="bold"),
              axis.title.y = element_text(size = 17),
              axis.text.x = element_text(size = 13, colour = "black"),
              axis.text.y = element_text(size = 13, colour = "black"),
              legend.title = element_blank(),
              legend.text = element_text(size = 17),
              panel.grid.major = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black",size=0.5),
              panel.border = element_rect(colour = "black",fill = "transparent",size=1.8))+
        geom_density_2d()+
        geom_point(data = X3,
                   mapping = aes(x = MDS1, y = MDS2),size=1.5,pch=17,colour="#FFE1FF")

      ggsave(plot=p1,height=8,width=10,dpi=500, filename=paste(path, run, "/NMDS/", OutName, "_NMDS.pdf", sep=""), useDingbats=FALSE)
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


  #NMDS on complete taxa profiles
  X1 = read.csv(paste(path, run, "/taxonomic_profiles/complete_profile.txt", sep=""),
                sep="\t", quote="", check.names=F, stringsAsFactors=F)
  X2 = dataPrep(X1)
  NMDSc(DF = X2, OutName = "all_taxa")


  #NMDS on ranks
  if (length(dir(paste(path, run, "/taxonomic_profiles/taxa_ranks",sep=""), pattern = ".txt"))>0){

    file.list = dir(paste(path, run, "/taxonomic_profiles/taxa_ranks",sep=""), pattern = ".txt")
    for (i in 1:length(file.list)) {
      X1 = read.csv(paste(path, run, "/taxonomic_profiles/taxa_ranks/", file.list[i], sep=""),
                    sep="\t", quote="",check.names=F, stringsAsFactors=F)
      X2 = dataPrep(X1)

      Name = sub(".txt", "", file.list[i])
      NMDSc(DF=X2, OutName=Name)
    }
  }
}
