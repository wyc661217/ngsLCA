#'
#' @title Count Reads and Taxa Numbers
#'
#' @description Count the reads number and taxa number after each filtering and de-contaminating, as well as in each taxa group. Results will be generated into "path/run/counts/"
#'
#' @param path working directory, same to \code{\link{ngsLCA_profile}}.
#' @param run name of the run, default is "run01".
#'
#' @return Text files and figures showing reads number and taxa number.
#' @import ggplot2
#' @importFrom utils read.csv read.delim write.table
#' @export
#'
#' @examples
#' ngsLCA_count(path=system.file("extdata","lca_files",package="ngsLCA"),
#'              run="run01")
#'
#'
#' ## This will generate tables and figures showing
#' ## the numbers of reads and taxa for "run01". Results
#' ## will be in "path/run01/counts/".
#'
#'
ngsLCA_count=function(path,
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
    cat("\n\n\t-> Count reads and taxa numbers.\n\n")
  }


  #local variables
  step = number = NULL


  #original taxa profile
  DF1 = read.delim(paste(path, run, "/intermediate/", "taxa_profile_v1.txt", sep=""),
                   stringsAsFactors=FALSE,header = T,check.names = F)
  ReadNO = data.frame(sample=c(colnames(DF1)[-1]),
                      original_Reads=colSums(DF1[-1]),
                      stringsAsFactors = F)

  for (i in 2:dim(DF1)[2]) {
    DF1[which(DF1[,i] > 0),i] = 1
  }
  TaxaNO =  data.frame(sample=c(colnames(DF1)[-1]),
                       taxa_in_original_lca=colSums(DF1[,-1]),
                       stringsAsFactors = F)


  #filtered taxa profile
  if (length(dir(paste(path, run, "/intermediate/", sep=""), pattern =  "taxa_profile_v2.3.txt")) > 0) {

    DF2.1 = read.delim(paste(path, run, "/intermediate/", "taxa_profile_v2.1.txt", sep=""),
                       stringsAsFactors=FALSE,header = T,check.names = F)
    DF2.2 = read.delim(paste(path, run, "/intermediate/", "taxa_profile_v2.2.txt", sep=""),
                       stringsAsFactors=FALSE,header = T,check.names = F)
    DF2.3 = read.delim(paste(path, run, "/intermediate/", "taxa_profile_v2.3.txt", sep=""),
                       stringsAsFactors=FALSE,header = T,check.names = F)

    ReadNO$passed_threshold.1 = colSums(DF2.1[,-1])
    ReadNO$passed_threshold.2 = colSums(DF2.2[,-1])
    ReadNO$passed_threshold.3 = colSums(DF2.3[,-1])

    for (i in 2:dim(DF1)[2]) {
      DF2.1[which(DF2.1[,i] > 0),i] = 1
      DF2.2[which(DF2.2[,i] > 0),i] = 1
      DF2.3[which(DF2.3[,i] > 0),i] = 1
    }

    TaxaNO$passed_threshold.1 = colSums(DF2.1[,-1])
    TaxaNO$passed_threshold.2 = colSums(DF2.2[,-1])
    TaxaNO$passed_threshold.3 = colSums(DF2.3[,-1])
  }


  #de-contaminated taxa profile
  if (length(dir(paste(path, run, "/intermediate/", sep=""), pattern =  "taxa_profile_v3.txt")) > 0) {

    DF3 = read.delim(paste(path, run, "/intermediate/", "taxa_profile_v3.txt", sep=""),
                     stringsAsFactors=FALSE,header = T,check.names = F)
    ReadNO$contaminate_removed = colSums(DF3[,-1])

    for (i in 2:dim(DF1)[2]) {
      DF3[which(DF3[,i] > 0),i] = 1
    }
    TaxaNO$contaminate_removed = colSums(DF3[,-1])
  }


  #order columns
  if (dim(ReadNO)[2] == 6) {
    ReadNO = ReadNO[,c(1,order(ReadNO[1,-1], decreasing = T) + 1)]
    TaxaNO = TaxaNO[,c(1,order(ReadNO[1,-1], decreasing = T) + 1)]
  }


  #groups
  if (length(dir(paste(path, run, "/intermediate/taxa_groups",sep=""), pattern = ".txt"))>0){

    file.list = dir(paste(path, run, "/intermediate/taxa_groups/", sep=""), pattern = ".txt")

    for (i in 1:length(file.list)) {

      DF = read.csv(paste(path, run, "/intermediate/taxa_groups/", file.list[i], sep=""),sep="\t",
                    quote="", check.names=F,stringsAsFactors=F)
      ReadNO$new = colSums(DF[,-1])
      colnames(ReadNO)[dim(ReadNO)[2]] = sub(".txt", "", file.list[i])

      for (j in 2:dim(DF)[2]) {
        DF[which(DF[,j] > 0),j] = 1
      }

      TaxaNO$new = colSums(DF[,-1])
      colnames(TaxaNO)[dim(TaxaNO)[2]] = sub(".txt", "", file.list[i])
    }
  }


  #create folder and write files
  if (length(dir(paste(path, run, sep = ""), pattern = "counts")) > 0) {
    cat(paste("\n\n\t-> '", path, run, "/counts' already exists; files inside will be overwritten.\n\n", sep = ""))
  }else{
    dir.create(paste(path, run, "/counts", sep=""))
    dir.create(paste(path, run, "/counts/readsNO", sep=""))
    dir.create(paste(path, run, "/counts/taxaNO", sep=""))
  }

  write.table(ReadNO,quote=F, row.names=F, col.names=T, sep="\t",
              file = paste(path, run, "/counts/","ReadNO.txt", sep=""))
  write.table(TaxaNO,quote=F, row.names=F, col.names=T, sep="\t",
              file = paste(path, run, "/counts/","TaxaNO.txt", sep=""))


  #generate figures
  if (dim(ReadNO)[2] > 2) {

    for (i in 1:dim(ReadNO)[1]) {

      D1 = data.frame(step=colnames(ReadNO)[-1],
                      number=as.numeric(ReadNO[i,-1]),
                      stringsAsFactors = F)
      D1$step = factor(D1$step, as.character(D1$step),ordered = T)

      D2 = data.frame(step=colnames(TaxaNO)[-1],
                      number=as.numeric(TaxaNO[i,-1]),
                      stringsAsFactors = F)
      D2$step = factor(D2$step, as.character(D2$step),ordered = T)

      Name = ReadNO$sample[i]


      p1=ggplot(data=D1, aes(x=step, y=number)) +
        geom_bar(stat="identity")+
        ggtitle(ReadNO$sample[i])+
        labs(y = "Reads number") +
        theme(axis.title.x = element_text(size = 13),
              axis.title.y = element_text(size = 13),
              axis.text.x = element_text(size = 11, colour = "black" ,angle = 45,hjust = 1,vjust = 1),
              axis.text.y = element_text(size = 11, colour = "black"),
              legend.text = element_text(size = 15),
              panel.grid.major = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              axis.line.x = element_line(colour = "black",size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black",size = 0.5, linetype = "solid"),
              axis.ticks.length = unit(0.2, "cm"))

      ggsave(plot=p1, height=8, width=6, dpi=300,
             filename=paste(path, run, "/counts/readsNO/",ReadNO$sample[i],"_ReadNO.pdf", sep=""),
             useDingbats=FALSE)

      p2=ggplot(data=D2, aes(x=step, y=number)) +
        geom_bar(stat="identity")+
        ggtitle(ReadNO$sample[i])+
        labs(y = "Reads number") +
        theme(axis.title.x = element_text(size = 13),
              axis.title.y = element_text(size = 13),
              axis.text.x = element_text(size = 11, colour = "black" ,angle = 45,hjust = 1,vjust = 1),
              axis.text.y = element_text(size = 11, colour = "black"),
              legend.text = element_text(size = 15),
              panel.grid.major = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              axis.line.x = element_line(colour = "black",size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black",size = 0.5, linetype = "solid"),
              axis.ticks.length = unit(0.2, "cm"))

      ggsave(plot=p2, height=8, width=6, dpi=300,
             filename=paste(path, run, "/counts/taxaNO/",ReadNO$sample[i],"_TaxaNO.pdf", sep=""),
             useDingbats=FALSE)
    }
  }
}
