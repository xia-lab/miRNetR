##################################################
## R script for miRNet
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################
.on.public.web <- FALSE; # only TRUE when on mirnet web server

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
    mir.nmsu <- vector();
    mir.nmsu <<- mir.nmsu;
    dataSet <- list(data=list());
    dataSet$ppiOpts<-list()
    dataSet$ppiOpts$db.name<-"innate"
    dataSet$ppiOpts$require.exp<-T
    dataSet$ppiOpts$min.score<-"900"
    dataSet <<- dataSet
    net.info<- list();
    net.info <<- net.info
    data.type <<- dataType; # mir or xeno.mir
    anal.type <<- analType;
    current.msg <<- "";
    msg.vec <<- vector(mode="character");
    module.count <<- 0;
    if(.on.public.web){
      lib.path <<- "../../data/libs/";
    }else{
      lib.path <<- "https://www.mirnet.ca/resources/data/libs/";
    }

    if(file.exists("/home/glassfish/sqlite/")){ #public server
        sqlite.path <<- "/home/glassfish/sqlite/";
    }else if(file.exists("/Users/xia/Dropbox/sqlite/")){# xia local
        sqlite.path <<- "/Users/xia/Dropbox/sqlite/";
    }else if(file.exists("/home/le/sqlite/gene-id-mapping/")){# le local
        sqlite.path <<- "/home/le/sqlite/mirnet/";
    }else if(file.exists("/home/soufanom/Database/")){# Othman local
        sqlite.path <<- "/home/soufanom/Database/sqlite/";
    }else if(file.exists("/home/zzggyy")){# zgy local
        sqlite.path <<- "/home/zzggyy/Downloads/netsqlite/";
    }else{
      sqlite.path <<- "https://github.com/xia-lab/miRNetR/raw/master/inst/database/sqlite/mirnet/";
    }

    # preload some general package
    library("RSQLite");
    library('Cairo');
    CairoFonts("Arial:style=Regular","Arial:style=Bold","Arial:style=Italic","Helvetica","Symbol")
    print("miRNet init done!");
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
  if (anal.type == "multilist"  || anal.type == "snp2mir" || anal.type == "tf2genemir" || anal.type == "gene2tfmir") {
    res <- dataSet[netType][[1]][, colInx];
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
  if (anal.type == "multilist"  || anal.type == "snp2mir" || anal.type == "tf2genemir" || anal.type == "gene2tfmir") {
   resTable <- dataSet[netType][[1]]
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
#' @param tblnm PARAM_DESCRIPTION
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
RemoveMirEntry <- function(tblnm, mir.id) {
    id <<- mir.id
    inx <- which(rownames(dataSet[tblnm][[1]]) == mir.id);
    if(length(inx) > 0){
        dataSet[tblnm][[1]] <- dataSet[tblnm][[1]][-inx,];
    }
    dataSet<<-dataSet
    return(1)
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
        if(dataSet$org %in% c("sma","gga","bta","ssc") || col.id != "evidence"){
            col.val <- as.numeric(gsub("Predicted miRanda Score:", "", col));
            # note NA will be introduced for non-predicted ones
            na.inx <- is.na(col.val);
            col.val[na.inx] <- max(col.val[!na.inx]);
            hits <- col.val > as.numeric(value);
        } else {
            print("This is only for Prediction Score at this moment");
            return("NA");
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
