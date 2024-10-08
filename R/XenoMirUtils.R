#' Setup Xeno miRNA List
#' @param mirs miRNAs.
#' @param orgType Organism type.
#' @param idType ID type.
#' @export
SetupXenoMirListData <- function(mirs, orgType, idType){

    dataSet$listData <- TRUE;
    data.org <<- dataSet$org <- orgType;
    dataSet$idType <- idType;

    current.msg <<- NULL;

    lines <- strsplit(mirs, "\r|\n|\r\n")[[1]];

    mir.lists <- strsplit(lines, "\\s+");
    mir.mat <- do.call(rbind, mir.lists);

    if(dim(mir.mat)[2] == 1){ # add *
        mir.mat <- cbind(mir.mat, rep("*", nrow(mir.mat)));
    }else if(dim(mir.mat)[2] > 2){
        mir.mat <- mir.mat[,1:2];
        current.msg <<- "More than two columns found in the list. Only first two columns will be used.";
    }
    if(idType == "mir_id"){
        mir.vec <- gsub("mir", "miR", mir.vec);
    }
    mir.vec <- mir.mat[,1];
    mir.mat <- mir.mat[,-1, drop=F];
    rownames(mir.mat) <- mir.vec;
    dataSet$mir.orig <- mir.mat;

    dataSet <<- dataSet;
    return (nrow(mir.mat));
}

#' Get miRNA Name
#' @export
GetSeqnm <- function(rowid){
  inx <- which(rownames(dataSet$mir.res) == rowid);
  mirid <- dataSet$mir.res[inx, "miRNA"];
  return(mirid);
}

#' Get miRNA Accession
#' @export
GetSeq <- function(mir.id){
  inx <- which(rownames(dataSet$mir.res) == mir.id);
  seq <- dataSet$mir.res[inx, "Accession"];
  return(seq);
}

# batch remove based on
#' Update Xeno miRNA Entries
#' @export
UpdateXenoMirEntries <- function(col.id, method, value, action) {

    if (col.id == "source"){
        col <- dataSet$mir.res$Source;
    } else if (col.id == "mir"){
        col <- dataSet$mir.res$miRNA;
    } else if (col.id == "target"){
        col <- dataSet$mir.res$Gene;
    } else if (col.id == "miranda"){
        col <- dataSet$mir.res$miRanda;
    } else if (col.id == "tarpmir"){
        col <- dataSet$mir.res$TarPmiR;
    } else if (col.id == "exp"){
        col <- dataSet$mir.res$Expression;
    } else {
        print(paste("unknown column:", col.id));
    }

    if (method == "contain"){
        hits <- grepl(value, col, ignore.case = TRUE);
    } else if (method == "match"){
        hits <- tolower(col) %in% tolower(value);
    } else {
        hits <- col >= as.numeric(value);
    }

    if(action == "keep"){
        hits = !hits;
    }

    if(sum(hits) > 0){
        row.ids <- rownames(dataSet$mir.res)[hits];
        dataSet$mir.res <<- dataSet$mir.res[!hits,];
        fast.write.csv(dataSet$mir.res, file="xeno_mirnet_target.csv", row.names=FALSE);
        return(row.ids);
    }else{
        return("NA");
    }
}

#' Get Unique Source Names
#' @export
GetUniqueSourceNames <- function(orgType){
    db.path <- paste("../../data/libs/xenomir.hosts.json");
    library("RJSONIO");
    browse <- fromJSON(db.path);
    browse <<- browse;
    res <- sort(names(browse[[orgType]]));
    return(res);
}

#' Get Unique Class Names
#' @export
GetUniqueClassNames <- function(orgType, source){
    res <- sort(names(browse[[orgType]][[source]]));
    return(res);
}

#' Get Unique Species Names
#' @export
GetUniqueSpeciesNames <- function(orgType){
    res <- unique(unlist(browse[[orgType]][[1]]));
    return(res);
}

#' Get Updated Class Names
#' @export
GetUpdateClassNames <- function(orgType, source){
    res <- sort(names(browse[[orgType]][[source]]));
    return(res);
}

#' Get Updated Species Names
#' @export
GetUpdateSpeNames <- function(orgType, source){
    res <- unlist(browse[[orgType]][[source]]);
    res <- sort(res);
    return(res);
}

#' Set up items
#' @export
SetupItemFromList <- function(orgType, species){
    mir.mat <- cbind(species, rep("*", length(species)));
    rownames(mir.mat) <- mir.mat[,1];
    mir.mat <- mir.mat[,-1, drop=F];
    dataSet$org <- orgType;
    dataSet$mir.orig <- dataSet$mir.mat <- mir.mat;
    dataSet$det <- species;
    dataSet$pre <- species;
    dataSet <<- dataSet;
    return(nrow(mir.mat));
}

#' Set up source
#' @export
SetupSourceFromList <- function(source){
    if(source == ""){
        print("Please at least choose one source!");
        return(0);
    }
    dataSet$source <- source;

    dataSet <<- dataSet;
    #return(length(source));
}

#' Species Mapping
#' @export
PerformSpeciesMapping <- function(status){
    orgType <- dataSet$org;
    source <- dataSet$source;

    idType <- "exo_species";
    det.vec <- dataSet$det;
    pre.vec <- dataSet$pre;
    det.dic <- QueryXenoMirSQLite(paste(sqlite.path, "xenomirnet", sep=""), det.vec, orgType, idType, source);
    det.res <- det.dic[, c("source", "exo_species", "exo_mirna", "mir_acc", "symbol", "entrez", "data", "exp","miranda", "tarpmir")];

    if (status == "true"){
        source.pre <- "predicted";
        pre.dic <- QueryXenoMirSQLite(paste(sqlite.path, "xenomirnet", sep=""), pre.vec, paste(orgType, "_pre", sep=""), idType, source.pre);
        pre.res <- pre.dic[, c("source", "exo_species", "exo_mirna", "mir_acc", "symbol", "entrez", "data", "exp", "miranda", "tarpmir")];
        res <- rbind(det.res, pre.res);
    } else {
        res <- det.res;
    }

    hit.num <- nrow(res);
    if (hit.num == 0) {
       current.msg <<- "No hits found in the database.";
       print(current.msg);
       return(0);
    } else {
       current.msg <<- paste("A total of unqiue ", hit.num, " pairs of miRNA-gene targets were identified!");
       res$tarpmir <- round(res$tarpmir, 3);
       colnames(res) <- c("Source", "Xeno.species", "miRNA", "Accession", "Gene", "Entrez", "Reference", "Expression", "miRanda", "TarPmiR");
       mir.nms <- res[, "miRNA"];
       gene.nms <- res[,"Gene"];
       net.info$gene.nms <- gene.nms;
       net.info <<-net.info;
       dataSet$seeds <- mir.nms;
       fast.write.csv(res, file="xeno_mirnet_target.csv", row.names=FALSE);

       #control if res is too big for view (table display and network visualization)
       # based on overall score and detected
       trimmed <- FALSE;
       if(hit.num > 10000){
            miranda.score <- (res$miRanda-140)/(200-140);
            tarpmir.score <- (res$TarPmiR-0.5)/(1-0.5);
            sum.score <- miranda.score + tarpmir.score + res$Expression;
            gd.rks <- rank(-sum.score) < 10000; # select top 10000
            res <- res[gd.rks,];
            trimmed <- TRUE;
       }

       mir.list <-  as.list(unique(res$miRNA));
       mir.mat <- do.call(rbind, mir.list);
       mir.mapped <- cbind(mir.mat, rep("*", nrow(mir.mat)));
       rownames(mir.mapped) <- mir.mapped[, 1];
       mir.mapped <- mir.mapped[,-1, drop=F];
       dataSet$mir.mapped <- mir.mapped;
       dataSet$mir.res <- res;
       dataSet$mirtarget <- "gene";
       dataSet$mir2gene <- res
       dataSet$mirtable <- "mir2gene"
       mir.nmsu = mir.mapped;
       dataSet <<- dataSet;
       if(trimmed){
            return(2);
       }else{
            return(1);
       }
    }
}

#' Xeno-miRNA Gene Target Mapping
#' @export
PerformXenoMirGeneMapping <- function(status){

    orgType <- dataSet$org;
    idType <- dataSet$idType;

    mir.mat <- dataSet$mir.orig;
    idVec <- rownames(mir.mat);

    source.vec <- "";

    det.dic <- QueryXenoMirSQLite(paste(sqlite.path, "xenomirnet", sep=""), idVec, orgType, idType, source.vec);
    det.res <- det.dic[, c("source", "exo_species", "exo_mirna", "mir_acc", "symbol", "entrez", "data", "exp", "miranda", "tarpmir")];

    if(status == "true"){
        pre.dic <- QueryXenoMirSQLite(paste(sqlite.path, "xenomirnet", sep=""), idVec, paste(orgType, "_pre", sep=""), idType, source.vec);
        pre.res <- pre.dic[, c("source", "exo_species", "exo_mirna", "mir_acc", "symbol", "entrez", "data", "exp", "miranda", "tarpmir")];
        res <- rbind(det.res, pre.res);
    } else {
        res <- det.res;
    }
    res$tarpmir <- round(res$tarpmir, 3);

    hit.num <- nrow(res)
    if(hit.num == 0){
        current.msg <<- "No hits found in the database. Please check your input";
        print(current.msg);
        return(0);
    }else{
        current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-gene targets were identified!");

        # record the mapped queries and change to same IDs used in network
        uniq.mat <- unique(res[, c("exo_mirna", "symbol", dataSet$idType)]);
        hit.inx <- match(rownames(mir.mat), uniq.mat[, dataSet$idType]);
        if(dataSet$idType %in% c("exo_mirna", "mir_acc")){
            rownames(mir.mat) <- uniq.mat[hit.inx, "exo_mirna"];
        }else{
            rownames(mir.mat) <- uniq.mat[hit.inx, "symbol"];
        }

        dataSet$mir.mapped <- mir.mat;

        # update col names
        colnames(res) <- c("Source", "Xeno.species", "miRNA", "Accession", "Gene", "Entrez", "Literature", "Expression", "miRanda", "TarPmiR");
        gene.nms <- res[,"Gene"];
        mir.nms <- res[, "miRNA"];
        if (dataSet$idType %in% c("exo_mirna", "mir_acc")) {
          dataSet$seeds <- mir.nms
        } else{
          dataSet$seeds <- gene.nms
        }
        net.info$gene.nms <- gene.nms;
        net.info <<-net.info;
        fast.write.csv(res, file="xeno_mirnet_target.csv", row.names=FALSE);

        #control if res is too big for view (table display and network visualization)
        # based on overall score and detected
        trimmed <- FALSE;
       if(hit.num > 10000){
            miranda.score <- (res$miRanda-140)/(200-140);
            tarpmir.score <- (res$TarPmiR-0.5)/(1-0.5);
            sum.score <- miranda.score + tarpmir.score + res$Expression;
            gd.rks <- rank(-sum.score) < 10000; # select top 10000
            res <- res[gd.rks,];
            trimmed <- TRUE;
       }
        dataSet$mir.res <- res;
        dataSet$mir2gene <- res
        mir.nmsu = mir.mat;
        dataSet$mirtable <- "mir2gene"
        dataSet$mirtarget <- "gene";
        dataSet <<- dataSet;
        if(trimmed){
            return(2);
        }else{
          if(.on.public.web){
            return(1);
          }else{
            return(current.msg)
          }
        }
    }
}

##################################################
## R functuions adapted from package miRBaseConverter
## Description: miRBase ID conversion between different version
## Author: Xu, Taosheng taosheng.x@gmail.com
###################################################

#' miRNA Name Conversion From Precursor to Mature
#' @export
miRNA_PrecursorToMature<-function(miRNANames,version=NULL){
  if(is.null(version)){
    c_version <- checkMiRNAVersion(miRNANames,verbose = FALSE)
    cat(paste0("The input miRNA version information: miRBase ",c_version))
  } else {
    c_version=version
  }
  ver_index=match(tolower(c_version),VER)
  if (is.na(ver_index))
    stop("It is a wrong version name, Please check it")
  VMAP <- miRNA_data[[ver_index]][,c(2,6,9)]
  VMAP[,1] <- SYM[VMAP[,1]]
  VMAP[,2] <- SYM[VMAP[,2]]
  VMAP[,3] <- SYM[VMAP[,3]]

  miRNANames <- as.character(miRNANames)
  miRNANames <- gsub(" ","",miRNANames)##Remove the possible space
  uid <- unique(as.vector(miRNANames))
  uid <- na.omit(uid)

  ind <- apply(VMAP,2,function(x){match(uid,x)})
  ind[which(is.na(ind[,1])),1] <- ind[which(is.na(ind[,1])),2]
  ind[which(is.na(ind[,1])),1] <- ind[which(is.na(ind[,1])),3]

  target <- data.frame(
    OriginalName = uid,
    Mature1 = VMAP[ind[,1],2],
    Mature2 = VMAP[ind[,1],3],
    row.names = NULL)
  target = target[match(miRNANames, target$OriginalName),]
  rownames(target)=NULL
  target
}

#' Check miRNA Version
#' @export
checkMiRNAVersion <- function(miRNANames,verbose=TRUE){
  miRNANames <- as.character(miRNANames)
  miRNANames <- gsub(" ","",miRNANames)##Remove the possible space

  uid <- unique(as.vector(miRNANames))
  SYM_ID <- match(uid,SYM)
  SYM_ID <- na.omit(SYM_ID)
  #result=data.frame(matrix(vector(),length(VER), 3,dimnames=list(c(), c("Version","Proportion","Recommend") )))
  #result$Version <- VER
  result <- data.frame("Version"=VER,"Proportion"=NA,"Recommend"="")

  for(i in 1:length(VER)){
    allSYM <- miRNA_data[[i]]
    allSYM <- c(allSYM[, 2],allSYM[, 6],allSYM[, 9])
    allSYM <- na.omit(allSYM)
    num <- length(intersect(allSYM, SYM_ID))
    result$Proportion[i] <- round((num/length(uid))*100, 2)
  }
  ind <- which(result$Proportion == max(result$Proportion))
  ind <- rev(ind)[1]
  if(verbose)
  {
    result[result$Proportion == max(result$Proportion), "Recommend"] <- " ***BEST Matched***"
    result$Proportion=paste(result$Proportion,"%",sep="")
    print.data.frame(result)
  }
  VER[ind]
}

#' Query Xeno-miRNA
#' @export
QueryXenoMirSQLite <- function(db.path, q.vec, table.nm, col.nm, source){
    
  db.path <- paste0(db.path, ".sqlite");
  if(.on.public.web){
    mir.db <- dbConnect(SQLite(), db.path);
  }else{
    msg <- paste("Downloading", db.path);
    db.name <- gsub(sqlite.path, "", db.path);
    if(!file.exists(db.name)){
      print(msg);
      download.file(db.path, db.name);
    }
    mir.db <- dbConnect(SQLite(), db.name);
  }
    query <- paste(shQuote(q.vec), collapse=",");
    if(source == "immune_organ"){
       source <- "immuneorgan";
    }
    if (col.nm == "exo_species" ){
       source <- paste("'", source, "'", sep="");
       statement <- paste("SELECT * FROM ",table.nm, " WHERE (",col.nm," IN (",query,") AND source IN (", source, "))", sep="");
    } else {
       statement <- paste("SELECT * FROM ",table.nm, " WHERE ", col.nm," IN (", query, ")", sep="");
    }
    mir.dic <- .query.sqlite(mir.db, statement);

    #Check the matched miRNA from upload list to database
    if (col.nm == "exo_mirna"){
        mir.lib <- as.vector(unique(mir.dic$exo_mirna));
        notMatch <- setdiff(q.vec, mir.lib);
        if (length(notMatch) > 0){ # Converting miRBase version and mature id.
          load("../../data/libs/mbcdata.rda");
          miRNANames <- gsub(" ","", as.character(notMatch));
          targetVersion <- "v21";
          ver_index <- match(tolower(targetVersion), VER)
          if (!is.na(ver_index)) {
                MiRNAs <- as.matrix(miRNA_data[[ver_index]])
                MiRNAs <- rbind(MiRNAs[, c(1,2,4)], MiRNAs[, c(5,6,7)], MiRNAs[, c(8,9,10)])
                colnames(MiRNAs) <- c("ACC","SYM","SEQ")

                ##check the rows with all NA
                ind <- apply(MiRNAs, 1, function(x) all(is.na(x)))
                VMAP <- as.data.frame(unique(MiRNAs[-ind,]))[, c("ACC", "SYM")];

                uid <- unique(as.vector(miRNANames))
                SYM_ID <- match(uid, SYM)
                df <- data.frame(uid = uid, SYM = SYM_ID)
                df <- merge(df, ACC_SYM)[, c("uid", "ACC")]
                df <- unique( merge(df, VMAP, by="ACC") )

                target <- data.frame(
                    OriginalName = df$uid,
                    TargetName = SYM[df$SYM],
                    Accession = ACC[df$ACC]
                );

                idx <- (target$OriginalName == target$TargetName) | (!target$OriginalName %in% target$TargetName)
                target <- target[idx, , drop=FALSE]

                ## collapse 1:many maps
                splitpaste <- function(x, f) {
                    result <- vapply(split(x, f), paste, character(1), collapse="&")
                    result[!nzchar(result)] <- NA
                    result
                }
                f <- factor(target$OriginalName, levels=uid)
                target <- data.frame(
                    OriginalName = uid,
                    TargetName = splitpaste(target$TargetName, f),
                    Accession = splitpaste(target$Accession, f),
                    row.names=NULL);

                target <- target[match(miRNANames, target$OriginalName),];
                mir.vec <- as.vector(target$TargetName);
                mir.vec <- unique(unlist(strsplit(mir.vec, split="&")));
                if (length(mir.vec) > 0){
                    query2 <- paste(shQuote(mir.vec), collapse=",");
                    statement2 <- paste("SELECT * FROM ", table.nm, " WHERE ", col.nm," IN (", query2, ")", sep="");
                    if(.on.public.web){
                      mir.db <- dbConnect(SQLite(), db.path);
                    }else{
                      msg <- paste("Downloading", db.path);
                      db.name <- gsub(sqlite.path, "", db.path);
                      if(!file.exists(db.name)){
                        print(msg);
                        download.file(db.path, db.name);
                      }
                      mir.db <- dbConnect(SQLite(), db.name);
                    }
                    mir.dic2 <- .query.sqlite(mir.db, statement2);
                    mir.lib2 <- as.vector(unique(mir.dic2$exo_mirna));
                    notMatch2 <- setdiff(mir.vec, mir.lib2);

                    if (length(notMatch2) > 0){
                        notMatch2 <- tolower(notMatch2);
                        res <- miRNA_PrecursorToMature(notMatch2);
                        mir.vec <- unique(c(res2$Mature1, res2$Mature2));
                        query3 <- paste(shQuote(mir.vec), collapse=",");
                        statement3 <- paste("SELECT * FROM ", table.nm, " WHERE ", col.nm," IN (", query3, ")", sep="");
                        if(.on.public.web){
                          mir.db <- dbConnect(SQLite(), db.path);
                        }else{
                          msg <- paste("Downloading", db.path);
                          db.name <- gsub(sqlite.path, "", db.path);
                          if(!file.exists(db.name)){
                            print(msg);
                            download.file(db.path, db.name);
                          }
                          mir.db <- dbConnect(SQLite(), db.name);
                        }
                        mir.dic3 <- .query.sqlite(mir.db, statement3);
                        mir.dic2 <- rbind(mir.dic2, mir.dic3);
                    }
                    mir.dic <- rbind(mir.dic, mir.dic2);
                }
            }
        }
    }
    dbDisconnect(mir.db); CleanMemory();
    return(mir.dic);
}
