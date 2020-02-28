##################################################
## R script for miRNet
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# init resources for analysis
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataType PARAM_DESCRIPTION
#' @param analType PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Init.Data
#' @export 
Init.Data<-function(dataType, analType){
    dataSet <<- list(data=list());
    data.type <<- dataType; # mir or xeno.mir
    anal.type <<- analType;
    current.msg <<- "";
    msg.vec <<- vector(mode="character");
    module.count <<- 0;
    lib.path <<- "../../data/libs/";

    if(file.exists("/home/glassfish/sqlite/")){ #public server
        sqlite.path <<- "/home/glassfish/sqlite/mirnet/";
        sqlite.geneid.path <<- "/home/glassfish/sqlite/";
        library(BiocParallel);
        register(SerialParam());
    }else if(file.exists("/Users/xia/Dropbox/sqlite/")){# xia local
        sqlite.path <<- "/Users/xia/Dropbox/sqlite/mirnet/";
        sqlite.geneid.path <<- "/Users/xia/Dropbox/sqlite/";
    }else if(file.exists("/home/le/sqlite/gene-id-mapping/")){# le local
        sqlite.path <<- "/home/le/sqlite/mirnet/";
        sqlite.geneid.path <<- "/home/le/sqlite/gene-id-mapping/";
    }else if(file.exists("//Users/soufanom/Documents/Postdoc/Projects/gene-id-mapping")){# Othman local
        sqlite.path <<- "/Users/soufanom/Documents/Postdoc/Projects/mirNet/";
        sqlite.geneid.path <<- "/Users/soufanom/Documents/Postdoc/Projects/gene-id-mapping";
    }

    # preload some general package
    library("RSQLite");
    library('Cairo');
    CairoFonts("Arial:style=Regular","Arial:style=Bold","Arial:style=Italic","Helvetica","Symbol")
    print("miRNet init done!");
}

# read tab delimited file
# stored in dataSet list object
# can have many classes, stored in meta.info
# type: array, count, qpcr
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataName PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ReadTabExpressData
#' @export 
ReadTabExpressData <- function(dataName) {

    dataSet <- ReadTabData(dataName);

    # rename data to data.orig
    int.mat <- dataSet$data;
    dataSet$cls <- dataSet$meta.info[,1];
    dataSet$data <- NULL;
    dataSet$listData <- FALSE;

    msg <- paste("a total of ", ncol(int.mat), " samples and ", nrow(int.mat), " features were found. ");

    # remove NA, null
    row.nas <- apply(is.na(int.mat)|is.null(int.mat), 1, sum);
    good.inx<- row.nas/ncol(int.mat) < 0.5;
    if(sum(!good.inx) > 0){
        int.mat <- int.mat[good.inx,];
        msg <- c(msg, paste("removed ", sum(!good.inx), " features with over 50% missing values"));
    }
    # remove constant values
    filter.val <- apply(int.mat, 1, IQR, na.rm=T);
    good.inx2 <- filter.val > 0;
    if(sum(!good.inx2) > 0){
        int.mat <- int.mat[good.inx2,];
        msg <- c(msg, paste("removed ", sum(!good.inx2), " features with constant values"));
    }

    if(nrow(int.mat) > 2000){
        filter.val <- filter.val[good.inx2];
        rk <- rank(-filter.val, ties.method='random');

        var.num <- nrow(int.mat);
        kept.num <- 0.95*var.num;
        int.mat <- int.mat[rk < kept.num, ];
        # msg <- c(msg, paste("removed 5% features with near-constant values"));
    }

    minVal <- min(int.mat, na.rm=T);
    na.inx <- is.na(int.mat);
    if(sum(na.inx) > 0){
        int.mat[na.inx] <- minVal/2;
        # msg <- c(msg, "the remaining", sum(na.inx), "missing variables were replaced with data min");
    }
    current.msg <<- paste(msg, collapse="; ");
    #dataSet$data.proc <- int.mat;
    saveRDS(int.mat, file="data.proc");
    dataSet <<- dataSet;
    return (1);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetClassInfo
#' @export 
GetClassInfo <- function(){
    return(levels(dataSet$cls));
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetAnotNames
#' @export 
GetAnotNames<-function(){
    return(rownames(dataSet$data.anot));
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mirs PARAM_DESCRIPTION
#' @param orgType PARAM_DESCRIPTION
#' @param idType PARAM_DESCRIPTION
#' @param tissue PARAM_DESCRIPTION
#' @param targetOpt PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SetupMirListData
#' @export 
SetupMirListData <- function(mirs, orgType, idType, tissue, targetOpt=NULL){

    dataSet$listData <- TRUE;
    data.org <<- dataSet$org <- orgType;
    dataSet$idType <- dataSet$mirnaType <- idType;
    dataSet$targetOpt <- targetOpt;
    dataSet$tissue <- tissue;
    current.msg <<- NULL;

    mir.mat <- .parseListData(mirs);
    mir.vec <- mir.mat[,1];

    if(data.type=="xeno.mir" && idType == "mir_id"){
        mir.vec <- gsub("mir", "miR", mir.vec);;
    }
    if(idType == "mir_id"){
      rownames(mir.mat) <-  tolower(as.vector(mir.vec));
    }else{
      rownames(mir.mat) <-  mir.vec;
    }
    
    mir.mat <- mir.mat[,-1, drop=F];
    dataSet$mir.orig <- mir.mat;

    dataSet<<- dataSet;
    return (nrow(mir.mat));
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param listInput PARAM_DESCRIPTION
#' @param orgType PARAM_DESCRIPTION
#' @param inputType PARAM_DESCRIPTION
#' @param idType PARAM_DESCRIPTION
#' @param tissue PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SetupIndListData
#' @export 
SetupIndListData <- function(listInput, orgType, inputType, idType, tissue){

  data.org <<- dataSet$org <- orgType;
  dataSet$tissue <- tissue;
  current.msg <<- NULL;
  dataSet$listData <- TRUE;

  in.mat <- .parseListData(listInput);
  in.vec <- in.mat[,1];
  rownames(in.mat) <-  in.vec;
  in.mat <- in.mat[,-1, drop=F];

  if(is.null(dataSet$data)){
        dataSet$data <- dataSet$id.types <- list();
  }

  dataSet$data[[inputType]] <- in.mat;
  dataSet$id.types[[inputType]] <- idType;
  dataSet <<- dataSet;
  return (nrow(in.mat));
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param orgType PARAM_DESCRIPTION, Default: 'hsa'
#' @param tissue PARAM_DESCRIPTION
#' @param idType PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SetupItemFromPickList
#' @export 
SetupItemFromPickList <- function(orgType="hsa", tissue, idType){
    if(!exists("picklist.vec")){
        print("Could not find user entered disease list!");
        return(0);
    }
    dataSet$org <- orgType;
    dataSet$tissue <- tissue;
    mir.mat <- .parsePickListItems(picklist.vec);

    if(anal.type == "multilist"){
        dataSet$data[[idType]] <- mir.mat;
    }else{
        dataSet$mir.orig <- mir.mat;
        dataSet$idType <- idType;
    }
    picklist.vec <<- NULL;
    dataSet <<- dataSet;
    return(nrow(mir.mat));
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SetupMirExpressData
#' @export 
SetupMirExpressData <- function(){
    idType <- dataSet$id.current;
    mydata <- data.matrix(dataSet$sig.mat[,"max.logFC",drop=FALSE]);
    if(idType == "mir_id"){ # note, in the mirnet database, all mir ids are lower case! miR=>mir
        mir.vec <- rownames(mydata);
        rownames(mydata) <- tolower(mir.vec);
    }
    dataSet$idType <- idType;
    dataSet$mir.orig <- mydata;
    dataSet<<- dataSet;
    return (nrow(mydata));
}

# "ID", "Accession","Gene", "PMID"
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param netType PARAM_DESCRIPTION
#' @param colInx PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetMirResCol
#' @export 
GetMirResCol <- function(netType, colInx){
  if (anal.type == "multilist"  || anal.type == "snp2mir" ) {
  if(netType == "gene2mir"){
    res <- dataSet$gene2mir[, colInx];
  } else if(netType == "mir2gene"){
    res <- dataSet$mir2gene[, colInx];
  } else if(netType == "mir2mol"){
    res <- dataSet$mir2mol[, colInx];
  } else if(netType == "mol2mir"){
    res <- dataSet$mol2mir[, colInx];
  } else if(netType == "mir2lnc"){
    res <- dataSet$mir2lnc[, colInx];
  } else if(netType == "lnc2mir"){
    res <- dataSet$lnc2mir[, colInx];
  } else if(netType == "mir2tf"){
    res <- dataSet$mir2tf[, colInx];
  } else if(netType == "tf2mir"){
    res <- dataSet$tf2mir[, colInx];
  } else if(netType == "mir2epi"){
    res <- dataSet$mir2epi[, colInx];
  } else if(netType == "epi2mir"){
    res <- dataSet$epi2mir[, colInx];
  } else if(netType == "mir2dis"){
    res <- dataSet$mir2dis[, colInx];
  } else if(netType == "dis2mir"){
    res <- dataSet$dis2mir[, colInx];
  } else if(netType == "snp2mir"){
    res <- dataSet$snp2mir[, colInx];
    }
  } else{
    res <- dataSet$mir.res[, colInx];
  }
    hit.inx <- is.na(res) | res == ""; # note, must use | for element-wise operation
    res[hit.inx] <- "N/A";
    return(res);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param netType PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetMirResRowNames
#' @export 
GetMirResRowNames <- function(netType){
  if (anal.type == "multilist"  || anal.type == "snp2mir" ) {
    if (netType == "gene2mir") {
      resTable <- dataSet$gene2mir;
    } else if (netType == "mir2gene") {
      resTable <- dataSet$mir2gene;
    } else if (netType == "mir2mol") {
      resTable <- dataSet$mir2mol;
    } else if (netType == "mol2mir") {
      resTable <- dataSet$mol2mir;
    } else if (netType == "mir2lnc") {
      resTable <- dataSet$mir2lnc;
    } else if (netType == "lnc2mir") {
      resTable <- dataSet$lnc2mir;
    } else if (netType == "mir2tf") {
      resTable <- dataSet$mir2tf;
    } else if (netType == "tf2mir") {
      resTable <- dataSet$tf2mir;
    } else if (netType == "mir2epi") {
      resTable <- dataSet$mir2epi;
    } else if (netType == "epi2mir") {
      resTable <- dataSet$epi2mir;
    } else if (netType == "mir2dis") {
      resTable <- dataSet$mir2dis;
    } else if (netType == "dis2mir") {
      resTable <- dataSet$dis2mir;
    } else if(netType == "snp2mir"){
      resTable <- dataSet$snp2mir;
    }
  }  else{
    resTable <- dataSet$mir.res;
  }
  if(nrow(resTable) > 1000){
    resTable <- resTable[1:1000, ];
    current.msg <<- "Due to computational constraints, only the top 1000 rows will be displayed.";
  }
  rownames(resTable);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mir.id PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname RemoveMirEntry
#' @export 
RemoveMirEntry <- function(mir.id) {
    inx <- which(rownames(dataSet$mir.res) == mir.id);
    if(length(inx) > 0){
        dataSet$mir.res <<- dataSet$mir.res[-inx,];
    }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param node.id PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname RemoveMirNode
#' @export 
RemoveMirNode <- function(node.id) {
    # node ID is name of either gene or miRNA
    row.ids <- NULL;
    # first try mir names
    hits <- dataSet$mir.res[, 1] == node.id;
    if(sum(hits) > 0){
        row.ids <- rownames(dataSet$mir.res)[hits];
    }else{
        hits <- dataSet$mir.res[, 3] == node.id;
        if(sum(hits) > 0){
            row.ids <- rownames(dataSet$mir.res)[hits];
        }
    }
    if(is.null(row.ids)){
        return("NA");
    }else{
        dataSet$mir.res <<- dataSet$mir.res[!hits,];
        return(row.ids);
    }
}

# batch remove based on
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param col.id PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION
#' @param value PARAM_DESCRIPTION
#' @param action PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname UpdateMirEntries
#' @export 
UpdateMirEntries <- function(col.id, method, value, action) {

    if(col.id == "evidence"){
        col <- dataSet$mir.res$Experiment;
    }else if(col.id == "mir"){
        col <- dataSet$mir.res[,1];
    }else if(col.id == "target"){
        col <- dataSet$mir.res[,3];
    }else if(col.id == "tissue"){
        col <- dataSet$mir.res$Tissue;
    }else if(col.id == "literature"){
        col <- dataSet$mir.res$Literature;
    } else {
        print(paste("unknown column:", col.id));
    }

    if(method == "contain"){
        hits <- grepl(value, col, ignore.case = TRUE);
    }else if(method == "match"){
        hits <- tolower(col) %in% tolower(value);
    }else{ # at least
        if(dataSet$org!= "sma" && dataSet$org!= "gga" && dataSet$org!= "bta" && dataSet$org!= "ssc" && col.id != "evidence"){
            print("This is only for Prediction Score at this moment");
            reuturn("NA");
        } else {
            col.val <- as.numeric(gsub("Predicted miRanda Score:", "", col));
            # note NA will be introduced for non-predicted ones
            na.inx <- is.na(col.val);
            col.val[na.inx] <- max(col.val[!na.inx]);
            hits <- col.val > as.numeric(value);
        }

    }

    if(action == "keep"){
        hits = !hits;
    }

    if(sum(hits) > 0){
        row.ids <- rownames(dataSet$mir.res)[hits];
        dataSet$mir.res <<- dataSet$mir.res[!hits,];
        write.csv(dataSet$mir.res, file="mirnet_mir_target.csv", row.names=FALSE);
        return(row.ids);
    }else{
        return("NA");
    }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetUniqueDiseaseNames
#' @export 
GetUniqueDiseaseNames <- function(){
    db.path <- paste(sqlite.path, "mir2disease", sep="");
    statement <- "SELECT disease FROM disease";
    return(GetUniqueEntries(db.path, statement));
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param orgType PARAM_DESCRIPTION, Default: 'hsa'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetUniqueMoleculeNames
#' @export 
GetUniqueMoleculeNames <- function(orgType="hsa"){
    db.path <- paste(sqlite.path, "mir2molecule", sep="");
    statement <- paste("SELECT molecule FROM ",orgType, sep="");
    return(GetUniqueEntries(db.path, statement));
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param orgType PARAM_DESCRIPTION, Default: 'hsa'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetUniqueEpigeneNames
#' @export 
GetUniqueEpigeneNames <- function(orgType="hsa"){
    db.path <- paste(sqlite.path, "mir2epi", sep="");
    statement <- paste("SELECT epi_regulator FROM ",orgType, sep="");
    return(GetUniqueEntries(db.path, statement));
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SetCurrentDataMulti
#' @export 
SetCurrentDataMulti <- function(){
  dataSet$type <- nms.vec;
  dataSet <<- dataSet;
  return(1);
}
