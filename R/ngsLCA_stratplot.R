#'
#' @title Generate Taxa Abundance Stratplots
#'
#' @description Generate stratplots for taxa profiles into "path/run/stratplot/" to show the reads abundance.
#'
#' @param path working directory, same to \code{\link{ngsLCA_profile}}.
#' @param run name of the run, default is "run01".
#' @param taxa.number maximum number of taxa will be shown in stratplots; default is 10.
#'
#' @return Taxa abundance stratplots.
#' @import ggplot2
#' @importFrom ggpubr ggarrange
#' @importFrom utils read.csv
#' @importFrom ggthemes theme_tufte
#' @importFrom cowplot save_plot
#' @export
#'
#' @examples
#' ngsLCA_stratplot(path=system.file("extdata","lca_files",package="ngsLCA"),
#'                  run="run01",
#'                  taxa.number=10)
#'
#'
#' ## This will generate stratplots for the complete taxa profile
#' ## as well as taxa groups and ranks (if available) of "run01",
#' ## with showing the 10 most abundant taxa in each profile.
#' ## Results will be in "path/run01/stratplot/".
#'
#'
ngsLCA_stratplot=function(path,
                          run="run01",
                          taxa.number=10){


  #modify input path
  if (!grepl('/$', path)) {
    path = paste(path,"/",sep="")
  }


  #local variables
  taxa = name = plot1 = plot2 = plot3 = plot4 = plot5 = plot6 = plot7 = plot8 = plot9 = plot10 = NULL


  #whether ngsLCA_profile performed
  if (length(dir(paste(path, run, "/intermediate/", sep=""), pattern =  "taxa_profile_v1.txt")) == 0) {
    cat(paste("\n\n\t-> Error: required input file not found under '", run, "' folder, please run 'ngsLCA_profile' first.\n\n",
              sep = ""))
    stop()
  }else{
    cat("\n\n\t-> Generate stratplot.\n\n")
  }


  #create folder
  if (length(dir(paste(path, run, sep = ""), pattern = "stratplot")) > 0) {
    cat(paste("\n\n\t-> '", path, run, "/stratplot' already exists; files inside will be overwritten.\n\n", sep = ""))
  }else{
    dir.create(paste(path, run, "/stratplot", sep=""))
  }


  #stratplots function
  stratplots <- function(DF,plotname){

    plots =list()

    if (any(is.na(suppressWarnings(as.numeric(rownames(DF)))))){

      for (i in 1:dim(DF)[2]) {

        Name = colnames(DF)[i]
        d3 = data.frame(name=rownames(DF),
                        taxa = DF[,i],
                        stringsAsFactors = F)
        d3$order = 1:dim(d3)[1]

        if(i==1){
          d5 <- ggplot(d3, aes(order, taxa)) +
            geom_ribbon(aes(ymin = 0, ymax= taxa), fill='#4682b4') +
            theme_tufte() + coord_flip() + ggtitle(Name) +
            scale_x_continuous(breaks=d3$order,labels=d3$name,trans = "reverse") +
            theme(axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  axis.text.x = element_text(colour = "black"),
                  axis.text.y = element_text(colour = "black"),
                  axis.line.x = element_line(colour = "black",size=0.25),
                  axis.line.y = element_line(colour = "black",size=0.25),
                  axis.ticks.length=unit(.25, "cm"))
          plots[[i]] = d5
        }else{
          d5 <- ggplot(d3, aes(order, taxa)) +
            geom_ribbon(aes(ymin = 0, ymax= taxa), fill='#4682b4') +
            theme_tufte() + coord_flip() + ggtitle(Name) +
            scale_x_continuous(breaks=d3$order,labels=d3$name,trans = "reverse") +
            theme(axis.title.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.text.x = element_text(colour = "black"),
                  axis.line.x = element_line(colour = "black",size=0.25),
                  axis.ticks.length=unit(.25, "cm"))
          plots[[i]] = d5
        }
      }
    }


    if (all(!is.na(suppressWarnings(as.numeric(rownames(DF)))))){

      for (i in 1:dim(DF)[2]) {

        Name = colnames(DF)[i]
        d3 = data.frame(name=rownames(DF),
                        taxa = DF[,i],
                        stringsAsFactors = F)
        d3$name = as.numeric(d3$name)

        if(i==1){

          d5 <- ggplot(d3, aes(name, taxa)) +
            geom_ribbon(aes(ymin = 0, ymax= taxa), fill='#4682b4') +
            theme_tufte() + coord_flip() + ggtitle(Name) +
            theme(axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  axis.text.x = element_text(colour = "black"),
                  axis.text.y = element_text(colour = "black"),
                  axis.line.x = element_line(colour = "black",size=0.25),
                  axis.line.y = element_line(colour = "black",size=0.25),
                  axis.ticks.length=unit(.25, "cm"))
          plots[[i]] = d5
        }else{
          d5 <- ggplot(d3, aes(name, taxa)) +
            geom_ribbon(aes(ymin = 0, ymax= taxa), fill='#4682b4') +
            theme_tufte() + coord_flip() + ggtitle(Name) +
            theme(axis.title.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.text.x = element_text(colour = "black"),
                  axis.line.x = element_line(colour = "black",size=0.25),
                  axis.ticks.length=unit(.25, "cm"))
          plots[[i]] = d5
        }
      }
    }

    out_plot = ggarrange(plotlist=plots, ncol=i)
    save_plot(out_plot, filename = paste(path, run, "/stratplot/",plotname,".pdf", sep=""), base_width = (i*3), base_height = 8)
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


  #complete taxa profile
  X1 = read.csv(paste(path, run, "/taxonomic_profiles/complete_profile.txt", sep=""),
                sep="\t", quote="", check.names=F,stringsAsFactors=F)
  X2 = dataPrep(X1)

  for (i in 1:dim(X2)[2]){
    if (sum(X2[,i])>0) {
      X2[,i] = X2[,i]/sum(X2[,i])
    }
  }

  if (dim(X2)[1]>taxa.number) {
    X2 = X2[order(rowSums(X2),decreasing = T)[1:taxa.number],]
  }
  X2 = as.data.frame(t(X2),stringsAsFactors = F)

  if (dim(X2)[2] < 11) {
    x3 = X2
    stratplots(x3,"all_taxa")
  }


  if (dim(X2)[2] > 10) {

    M=as.integer(dim(X2)[2]/10)

    if (dim(X2)[2]/10 == M) {

      for (n in 1:M) {
        x3 = X2[,(10*n-9):(10*n)]
        stratplots(x3,paste("all_taxa_",n,sep = ""))
      }
    }

    if (dim(X2)[2]/10 != M) {

      for (n in 1:M) {
        x3 = X2[,(10*n-9):(10*n)]
        stratplots(x3,paste("all_taxa_",n,sep = ""))
      }

      x3 = as.data.frame(X2[,(M*10+1):dim(X2)[2]])
      rownames(x3) = rownames(X2)
      colnames(x3) = colnames(X2)[(M*10+1):dim(X2)[2]]
      stratplots(x3,paste("all_taxa_",n+1,sep = ""))

    }
  }


  #taxa ranks
  if (length(dir(paste(path, run, "/taxonomic_profiles/taxa_ranks",sep=""), pattern = ".txt"))>0){

    file.list = dir(paste(path, run, "/taxonomic_profiles/taxa_ranks",sep=""), pattern = ".txt")
    NAMES = sub(".txt","",file.list)

    for (i in 1:length(file.list)) {

      X1 = read.csv(paste(path, run, "/taxonomic_profiles/taxa_ranks/", file.list[i], sep=""),
                    sep="\t", quote="",check.names=F,stringsAsFactors=F)
      X2 = dataPrep(X1)

      for (m in 1:dim(X2)[2]){
        if (sum(X2[,m])>0) {
          X2[,m] = X2[,m]/sum(X2[,m])
        }
      }

      if (dim(X2)[1]>taxa.number) {
        X2 = X2[order(rowSums(X2),decreasing = T)[1:taxa.number],]
      }
      X2 = as.data.frame(t(X2),stringsAsFactors = F)

      if (dim(X2)[2] < 11) {
        x3 = X2
        stratplots(x3,file.list[i])
      }

      if (dim(X2)[2] > 10) {

        M=as.integer(dim(X2)[2]/10)

        if (dim(X2)[2]/10 == M) {

          for (n in 1:M) {
            x3 = X2[,(10*n-9):(10*n)]
            stratplots(x3,paste(NAMES[i],"_",n,sep = ""))
          }
        }

        if (dim(X2)[2]/10 != M) {

          for (n in 1:M) {
            x3 = X2[,(10*n-9):(10*n)]
            stratplots(x3,paste(NAMES[i],"_",n,sep = ""))
          }

          x3 = as.data.frame(X2[,(M*10+1):dim(X2)[2]])
          rownames(x3) = rownames(X2)
          colnames(x3) = colnames(X2)[(M*10+1):dim(X2)[2]]
          stratplots(x3,paste(NAMES[i],"_",n+1,sep = ""))

        }
      }
    }
  }
}
