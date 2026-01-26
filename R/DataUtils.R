##################################################
## R script for miRNet
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

.onAttach <- function (libname, pkgname){
  .on.public.web <<- FALSE;
  k1 <- paste("miRNetR",
              utils::packageVersion( "miRNetR"),
              "initialized Successfully !")
  k0 <- "\n";
  packageStartupMessage(c(k1,k0));
}

# init resources for analysis
#' Initiate Data
#' @export
Init.Data<-function(dataType, analType, onWeb=T){
  globalConfig <- list();
  .on.public.web <<- onWeb;
  globalConfig$anal.mode <- "web";
  globalConfig <<- globalConfig;
  mir.nmsu <- vector();
  mir.nmsu <<- mir.nmsu;
  dataSet <- list(data=list());
  dataSet$ppiOpts<-list()
  dataSet$ppiOpts$db.name<-"innate"
  dataSet$ppiOpts$require.exp<-T;
  dataSet$ppiOpts$min.score<-"900";
  dataSet$report.format <- "html";
  
  dataSet <<- dataSet
  net.info<- list();
  net.info <<- net.info
  data.type <<- dataType; # mir or xeno.mir
  anal.type <<- analType;
  current.msg <<- "";
  msg.vec <<- vector(mode="character");
  module.count <<- 0;
  
  #api.base <<- "132.216.38.6:8987"
  api.base <<- "http://api.xialab.ca"
  
  if(.on.public.web){
    lib.path <<- "../../data/libs/";
  }else{
    lib.path <<- "https://www.mirnet.ca/resources/data/libs/";     
  }
  
  if(file.exists("/home/glassfish/sqlite/")){ #public server
    sqlite.path <<- "/home/glassfish/sqlite/";
  }else if(file.exists("/Users/xialab/Dropbox/sqlite/")){# xia local
    sqlite.path <<- "/Users/xialab/Dropbox/sqlite/";
  }else if(file.exists("/Users/jeffxia/Dropbox/sqlite/")){# xia local2
    sqlite.path <<- "/Users/jeffxia/Dropbox/sqlite/";
  }else if(file.exists("/home/zzggyy")){# zgy local
    sqlite.path <<- "/media/zzggyy/disk/sqlite/";
  }else if(file.exists("/home/zgy/sqlite")){# zgy local2
    sqlite.path <<- "/home/zgy/sqlite/";
  }else if(file.exists("/home/fiona")){# fiona local
    sqlite.path <<- "/home/fiona/sqlite/";
  }else{
    sqlite.path <<- "https://www.xialab.ca/resources/sqlite/";
  }
  
  # preload some general package
  library("RSQLite");
  library('Cairo');
  CairoFonts("Arial:style=Regular","Arial:style=Bold","Arial:style=Italic","Helvetica","Symbol")
  paramSet <- list(objName="paramSet", jsonNms=list());
  msgSet <- list(objName="msgSet");
  cmdSet <<- list(objName="cmdSet");
  imgSet <- list(objName="imgSet", enrTables=list());
  infoSet <- list();
  
  infoSet$objName <- "infoSet";
  infoSet$paramSet <- paramSet;
  infoSet$cmdSet <- cmdSet;
  infoSet$msgSet <- msgSet;
  infoSet$imgSet <- imgSet;
  saveSet(infoSet);
  
  data.org <<- "NA";
  print("miRNet init done!");
}


# "ID", "Accession","Gene", "PMID"
#' Get Result Column
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

#' Get Result Row
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

#' Remove miRNA Entry
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
#' Update miRNA Entries
#' @export
UpdateMirEntries <- function(col.id, method, value, action, tblnm="") {
  #save.image("updateentries.RData");
  use.mir.res <- T;
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
  if(is.null(col)){
    use.mir.res <- F;
    
    if(col.id == "evidence"){
      col <- dataSet[tblnm][[1]]$Experiment;
    }else if(col.id == "mir"){
      col <- dataSet[tblnm][[1]][,1];
    }else if(col.id == "target"){
      col <- dataSet[tblnm][[1]][,3];
    }else if(col.id == "tissue"){
      col <- dataSet[tblnm][[1]]$Tissue;
    }else if(col.id == "literature"){
      col <- dataSet[tblnm][[1]]$Literature;
    } else {
      print(paste("unknown column:", col.id));
    }
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
  
  print(paste(sum(hits)));
  if(sum(hits) > 0){
    if(use.mir.res){
      row.ids <- rownames(dataSet$mir.res)[hits];
      dataSet$mir.res <- dataSet$mir.res[!hits,];
    }else{
      if(tblnm != "NA"){
        row.ids <- rownames(dataSet[tblnm][[1]])[hits];
        dataSet[tblnm][[1]] <- dataSet[tblnm][[1]][!hits,];
        if("tableType" %in% colnames(dataSet$mir.res)){
          rowsToRemove <- with(dataSet$mir.res, tableType == tblnm & originalRowId %in% row.ids)
          rowsToKeep <- !rowsToRemove
          dataSet$mir.res <- dataSet$mir.res[rowsToKeep, ]
        }
      }
    }
    dataSet<<-dataSet;
    fast.write.csv(dataSet$mir.res, file="mirnet_mir_target.csv", row.names=FALSE);
    return(row.ids);
  }else{
    return("NA");
  }
}

#' Prepare JSON File
#' @export
PrepareJsonFromR <- function(fileNm, type, jsonString, dataSetString){
  library(RJSONIO)
  dataSet <- fromJSON(dataSetString);
  dataSet <<- dataSet
  sink(fileNm);
  cat(jsonString);
  sink();
  return(1)
}


saveSet <- function(obj=NA, set="", output=1){

  if(globalConfig$anal.mode == "api"){
    qs:::qsave(obj, paste0(set, ".qs"));
    # CRITICAL: Prevent race condition - allow file system to sync before Java reads
    Sys.sleep(0.15);
  }else{
    if(set == ""){
      set <- obj$objName;
    }
    if(set == "dataSet"){
      dataSet <<- obj;
    }else if(set == "analSet"){
      analSet <<- obj;
    }else if(set == "imgSet"){
      imgSet <<- obj;
    }else if(set == "paramSet"){
      paramSet <<- obj;
    }else if(set == "msgSet"){
      msgSet <<- obj;
    }else if(set == "cmdSet"){
      cmdSet <<- obj;
    }else if(set == "infoSet"){
      infoSet <<- obj;
    }
    
  }
  return(output);
  
}

readSet <- function(obj=NA, set=""){
  if(globalConfig$anal.mode == "api"){
    path <- "";
    if(exists('user.path')){
      path <- user.path;
    }
    
    if(path != ""){
      obj <- load_qs(paste0(path, set, ".qs"));
    }else{
      obj <- qs:::qread(paste0(set, ".qs"));
    }
  }
  return(obj);
}


#'Record R Commands
#'@param cmd Commands 
#'@export
RecordRCommand <- function(cmd){
  infoSet <- readSet(infoSet, "infoSet"); 
  infoSet$cmdSet$cmdVec <- c(infoSet$cmdSet$cmdVec, cmd);
  saveSet(infoSet, "infoSet");
  return(1);
}

SaveRCommands <- function(){
  infoSet <- readSet(infoSet, "infoSet"); 
  cmds <- paste(infoSet$cmdSet$cmdVec, collapse="\n");
  pid.info <- paste0("# PID of current job: ", Sys.getpid());
  cmds <- c(pid.info, cmds);
  write(cmds, file = "Rhistory.R", append = FALSE);
}

#'Export R Command History
#'@export
GetRCommandHistory <- function(){
  infoSet <- readSet(infoSet, "infoSet"); 
  if(length(infoSet$cmdSet$cmdVec) == 0){
    return("No commands found");
  }
  return(infoSet$cmdSet$cmdVec);
}

ClearRCommandHistory <- function(){
  infoSet <- readSet(infoSet, "infoSet"); 
  infoSet$cmdSet$cmdVec <- c();
}