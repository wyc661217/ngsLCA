#'
#' @title Generate taxa abundance heatmaps
#'
#' @description Generate heatmaps for taxa profiles into 'path/run/heatmap/' to show the reads abundance of each taxon.
#'
#' @param path working directory, same to 'ngsLCA_profile'.
#' @param run name of the run, default is "run01".
#' @param taxa.number maximum number of taxa will be shown in heatmaps; default is 30.
#'
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar unit
#' @importFrom utils read.csv
#' @export
#'
#' @examples
#' ngsLCA_heatmap(path="working_directory/",
#'                run="run01",
#'                taxa.number=20)
#'
#'
#' ## This will generate heatmaps for the complete taxa profile
#' ## as well as taxa groups and ranks (if available) of 'run01'
#' ## in 'working_directory/run01/heatmap/', with showing the 20
#' ## most abundant taxa in each profile.
#'
#'
ngsLCA_heatmap=function(path,
                        run="run01",
                        taxa.number=30){


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
    cat("\n\n\t-> Generate heatmaps.\n\n")
  }


  #function for heatmap
  HeatMap = function(DF){
    #transfer data into percentage and subset
    for (i in 1:dim(DF)[2]) {
      if(sum(DF[,i]) != 0){
        DF[,i] = DF[,i]/sum(DF[,i])
      }
    }
    DF$sum = rowSums(DF)
    DF = DF[order(-DF$sum),]

    if (dim(DF)[1] > taxa.number) {
      DF = DF[1:taxa.number,-dim(DF)[2]]
    } else{
      DF = DF[,-dim(DF)[2]]
    }

    DF = as.matrix(DF)
    if (length(which(colSums(DF) == 0))>0) {
      F1 = Heatmap(DF, column_dend_height = unit(1.5, "cm"), row_dend_width = unit(3, "cm"), show_row_names = T, show_column_names = T,
                   row_names_gp = gpar(cex=0.8), column_names_gp = gpar(cex=0.8), cluster_rows = T, cluster_columns= F,
                   clustering_distance_rows = "pearson", clustering_distance_columns = "euclidean",
                   col = colorRamp2(c(0, 0.001, 0.01, 0.03, 0.06, 0.2),
                                    c("white", "cornflowerblue", "yellow", "#FD8D3C","#E31A1C","#B10026")),
                   rect_gp = gpar(col = "gray12", lty = 1, lwd = 0.5),
                   heatmap_legend_param = list(title = "abundance", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8),
                                               legend_height = unit(10, "cm"), color_bar = "continous"))
    }else{
      F1 = Heatmap(DF, column_dend_height = unit(1.5, "cm"), row_dend_width = unit(3, "cm"), show_row_names = T, show_column_names = T,
                   row_names_gp = gpar(cex=0.8), column_names_gp = gpar(cex=0.8), cluster_rows = T, cluster_columns= T,
                   clustering_distance_rows = "pearson",clustering_distance_columns = "euclidean",
                   col = colorRamp2(c(0, 0.001, 0.01, 0.03, 0.06, 0.2),
                                    c("white", "cornflowerblue", "yellow", "#FD8D3C","#E31A1C","#B10026")),
                   rect_gp = gpar(col = "gray12", lty = 1, lwd = 0.5),
                   heatmap_legend_param = list(title = "abundance", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8),
                                               legend_height = unit(10, "cm"), color_bar = "continous"))
    }
    return(F1)
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
  if (length(dir(paste(path, run, sep = ""), pattern = "heatmap")) > 0) {
    cat(paste("\n\n\t-> '", path, run, "/heatmap' already exists; files inside will be overwritten.\n\n", sep = ""))
  }else{
    dir.create(paste(path, run, "/heatmap", sep=""))
  }


  #generate heapmaps
  X1 = read.csv(paste(path, run, "/taxonomic_profiles/complete_profile.txt", sep=""),
                sep="\t", quote="", check.names=F,stringsAsFactors=F)
  X2 = dataPrep(X1)

  if (dim(X2)[2]>1) {
    pdf(paste(path, run, "/heatmap/complete_profile_heatmap.pdf", sep=""), width=8, height=13)
    print({HeatMap(X2)})
    dev.off()
  }else{
    cat("\n\n\t-> Only one sample detected, heatmap cannot be generated.\n\n")
  }

  if (length(dir(paste(path, run, "/taxonomic_profiles/taxa_ranks",sep=""), pattern = ".txt"))>0){

    file.list = dir(paste(path, run, "/taxonomic_profiles/taxa_ranks",sep=""), pattern = ".txt")

    for (i in 1:length(file.list)) {

      X1 = read.csv(paste(path, run, "/taxonomic_profiles/taxa_ranks/", file.list[i], sep=""),
                    sep="\t", quote="",check.names=F,stringsAsFactors=F)
      X2 = dataPrep(X1)

      Name = sub(".txt", "", file.list[i])

      if (dim(X2)[2]>1) {
        pdf(paste(path, run, "/heatmap/",Name,"_heatmap.pdf", sep=""), width=8, height=13)
        print({HeatMap(X2)})
        dev.off()
      } else {
        cat(paste("\n\n\t-> Only one sample detected in ", file.list[i],", heatmap not generated for it.\n\n",sep = ""))
      }
    }
  }
}
