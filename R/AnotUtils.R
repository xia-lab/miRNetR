##################################################
## R script for miRNet
## Description: Gene/Compound Annotation
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

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
#' @rdname PerformMirGeneMapping
#' @export
PerformMirGeneMapping <- function(){
    mir.mat <- dataSet$mir.orig;
    idVec <- rownames(mir.mat);

    mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2gene", sep=""), idVec, dataSet$org, dataSet$idType);

    hit.num <- nrow(mir.dic)
    if (hit.num == 0 && dataSet$tissue == "na") {
        current.msg <<- "No hits found in the database. Please check your input.";
        print(current.msg);
        return(0);
    } else if (hit.num == 0 && dataSet$tissue != "na") {
        current.msg <<- "No hits found in the database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
        print(current.msg);
        return(2);
    } else {
        res <- mir.dic[ , c("mir_id", "mir_acc", "symbol", "entrez", "experiment", "pmid", "tissue")];
        rownames(res) <- mir.dic$mirnet;
        current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-gene targets were identified!");

        # record the mapped queries and change to same IDs used in network
        uniq.mat <- unique(mir.dic[, c("mir_id", "symbol", dataSet$idType)]);
        hit.inx <- match(rownames(mir.mat), uniq.mat[, dataSet$idType]);
        if(dataSet$idType %in% c("mir_id", "mir_acc")){
            rownames(mir.mat) <- uniq.mat[hit.inx,"mir_id"];
        }else{
            rownames(mir.mat) <- uniq.mat[hit.inx,"symbol"];
        }
        dataSet$mir.mapped <- mir.mat;

        # update col names
        colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue");
        gene.nms <- res[,"Gene"];
        mir.nms <- res[, "ID"];
        if(dataSet$idType %in% c("mir_id", "mir_acc")){
          dataSet$seeds <- mir.nms;
        }else{
          dataSet$seeds <- gene.nms;
        }
        write.csv(res, file="mirnet_mir_target.csv", row.names=FALSE);
        dataSet$mir.res <- res;
        dataSet$mirtable <- "gene2mir"
        dataSet$gene2mir <- res
        net.info$gene.nms <- gene.nms;
        net.info <<-net.info
        dataSet$mirtarget <- "gene";
        dataSet <<- dataSet;
        return(current.msg);
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
#' @rdname PerformSNPMirGeneMapping
#' @export
PerformSNPMirGeneMapping <- function(){
  snp.mat <- dataSet$mir.orig;
  snpidVec <- rownames(snp.mat);
  idType <- dataSet$idType;

  # first try to match snp2mir, if still have unmatched, do snp2gene
  snp.dic <- Query.miRNetDB(paste(sqlite.path, "snp2mir", sep=""), snpidVec, dataSet$org, dataSet$idType);
  hit.inx <- match(snpidVec, snp.dic$rsid);
  if(NA %in% hit.inx){
    unmatched.snp <- snpidVec[is.na(hit.inx)];
    snp.dic2 <- Query.miRNetDB(paste(sqlite.path, "snp2gene", sep=""), unmatched.snp, dataSet$org, dataSet$idType);
    snp.num <- rbind(snp.dic[,1:2], snp.dic2[,1:2])
  }
  snp.num <- snp.dic;
  hit.num <- nrow(snp.num)
  if (hit.num == 0) {
    current.msg <<- "No hits found in the database. The SNP list has not been annotated to miRNA or miRNA-binding sites. Please check your input.";
    print(current.msg);
    return(0);
  } else {
    if(NA %in% hit.inx){
      # snp2mir2gene
      snp <- na.omit(data.frame(name1 = snp.dic[, idType], id1 = snp.dic[, "rsid"], name2 = snp.dic[, "Mature_Name"], id2 =  snp.dic[, "Mature_Acc"], stringsAsFactors = FALSE));
      hit.num <- nrow(snp);
      idVec <- as.vector(unique(snp.dic[, c("MIRNA_Name")]));
      mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2gene", sep=""), idVec, dataSet$org, "mir_id");

      # for network
      snp.edge <- na.omit(data.frame(Name1=snp.dic[,"Mature_Name"],ID1=snp.dic[,"Mature_Acc"],Name2=snp.dic[,"rsid"],ID2=snp.dic[,"rsid"],stringsAsFactors = FALSE));
      mir.edge <- na.omit(data.frame(Name1=mir.dic[,"mir_id"],ID1=mir.dic[,"mir_acc"],Name2=mir.dic[,"symbol"],ID2=mir.dic[,"entrez"],stringsAsFactors = FALSE));

      # table results
      snp.res <- na.omit(snp.dic[ , c("chr_pos", "rsid", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain")]);
      mir.res <- mir.dic[ , c("mir_id", "mir_acc", "symbol", "entrez", "experiment", "pmid", "tissue")];
      snp.res$Database <- rep("30302893", nrow(snp.res));
      colnames(snp.res) <- c("CHR_POS", "rsID", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain", "DataBase");
      colnames(mir.res) <- c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue");

      # snp2gene2mir
      snp2 <- na.omit(data.frame(name1 = snp.dic2[, idType], id1 = snp.dic2[, "rsid"], name2 = snp.dic2[, "symbol"], id2 =  snp.dic2[, "entrez"], stringsAsFactors = FALSE));
      hit.num2 <- nrow(snp2);
      current.msg <<- paste("A total of", hit.num, "SNPs were mapped to miRNAs &", hit.num2, "SNPs were mapped to miRNA-binding sites!");

      idVec2 <- as.vector(unique(snp.dic2[, c("entrez")]));
      mir.dic2 <- Query.miRNetDB(paste(sqlite.path, "mir2gene", sep=""), idVec2, dataSet$org, "entrez");

      # for network
      snp.edge2 <- na.omit(data.frame(Name1=snp.dic2[,"symbol"],ID1=snp.dic2[,"entrez"],Name2=snp.dic2[,"rsid"],ID2=snp.dic2[,"rsid"],stringsAsFactors = FALSE));
      mir.edge2 <- na.omit(data.frame(Name1=mir.dic2[,"symbol"],ID1=mir.dic2[,"entrez"],Name2=mir.dic2[,"mir_id"],ID2=mir.dic2[,"mir_acc"],stringsAsFactors = FALSE));

      # table results
      snp.res2 <- na.omit(snp.dic2[ , c("chr_pos", "rsid", "transcript_id", "entrez", "symbol")]);
      mir.res2 <- mir.dic2[ , c("mir_id", "mir_acc", "symbol", "entrez", "experiment", "pmid", "tissue")];
      snp.res2$Literature <- rep("24163105", nrow(snp.res2));
      snp.res2$Database <- rep("PolymiRTS_3.0", nrow(snp.res2));
      colnames(snp.res2) <- c("CHR_POS", "rsID", "Transcript_ID", "Entrez", "Gene", "Literature", "Database");
      colnames(mir.res2) <- c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue");

      # rbind snp2mir2gene and snp2gene2mir for network


      # update the data
      gd.inx <- rownames(snp.mat) %in% unique(c(snp.edge[, "Name2"], snp.edge2[, "Name2"]));
      dataSet$mir.mapped <- snp.mat[gd.inx,,drop=F];

      dataSet$seeds <- c(snp.edge[, "Name2"], snp.edge2[, "Name2"]);

      dataSet$mirtarget <- c("gene");
      dataSet$mirtable <- c("snp2mir", "mir2gene", "snp2gene", "gene2mir");

      tf.dic <- Query.miRNetDB(paste(sqlite.path, "snp2tf", sep=""), snpidVec, dataSet$org, dataSet$idType);
      res <- na.omit(tf.dic[ , c("chr_pos", "rsid", "entrez", "symbol", "name")]);
      res$Literature <- rep("NA", nrow(res));
      res$Database <- rep("NA", nrow(res));
      colnames(res) <- c("CHR_POS", "rsID", "Entrez", "Symbol", "Name");
      tf.vec = res[,"Entrez"]
      tf.res = res

      tf.dic <- Query.miRNetDB(paste(sqlite.path, "mir2gene", sep=""), tf.vec, dataSet$org, "entrez");
      res <- tf.dic[ , c("mir_id", "mir_acc", "symbol", "entrez", "experiment", "pmid", "tissue")];
      res$Literature <- rep("NA", nrow(res));
      res$Database <- rep("NA", nrow(res));
      colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue");
      tf.res2 = res

      tf.edge <- na.omit(data.frame(Name1=tf.res[,4],ID1=tf.res[,3],Name2=tf.res[,2],ID2=tf.res[,1],stringsAsFactors = FALSE));
      tf.edge2 <- na.omit(data.frame(Name1=tf.res2[,1],ID1=tf.res2[,2],Name2=tf.res2[,3],ID2=tf.res2[,4],stringsAsFactors = FALSE));
     merge.edge <- rbind(snp.edge, mir.edge, snp.edge2, mir.edge2, tf.edge, tf.edge2);
      rownames(merge.edge) <- 1:nrow(merge.edge);
    dataSet$mir.res <- merge.edge; # for network
      write.csv(merge.edge, file="mirnet_snp_mir_target.csv", row.names=FALSE);

      dataSet$snp2mir <- snp.res;
      dataSet$mir2gene <- mir.res;
      dataSet$snp2gene <- snp.res2;
      dataSet$gene2mir <- mir.res2;
if(nrow(tf.res)>0){
      dataSet$snp2tf <- tf.res;
      dataSet$tf2mir <- tf.res2;
      dataSet$mirtable <- c(dataSet$mirtable, "snp2tf", "tf2mir");
}
    net.info$tf.nms <<- unique(tf.edge$Name1)
    net.info$gene.nms <<- unique(c(snp.edge2$Name1, mir.edge2$Name1, mir.edge$Name2))
    if(nrow(tf.res2)>0){
    mir.nmsu <<- unique(c(mir.edge2$Name2, snp.edge$Name1, mir.edge$Name1, tf.res2[,"ID"]))
    }else{
    mir.nmsu <<- unique(c(mir.edge2$Name2, snp.edge$Name1, mir.edge$Name1))
    }
      dataSet$nodeNumbers <- nrow(merge.edge)
      dataSet <<- dataSet;
      return(1);
    }else{
      # snp2mir2gene
      snp <- na.omit(data.frame(name1 = snp.dic[, idType], id1 = snp.dic[, "rsid"], name2 = snp.dic[, "Mature_Name"], id2 =  snp.dic[, "Mature_Acc"], stringsAsFactors = FALSE));
      hit.num <- nrow(snp);
      current.msg <<- paste("A total of", hit.num, "SNPs were mapped to miRNAs!");

      idVec <- as.vector(unique(snp.dic[, c("MIRNA_Name")]));
      mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2gene", sep=""), idVec, dataSet$org, "mir_id");

      # for network
      snp.edge <- na.omit(data.frame(Name1=snp.dic[,"Mature_Name"],ID1=snp.dic[,"Mature_Acc"],Name2=snp.dic[,"rsid"],ID2=snp.dic[,"rsid"],stringsAsFactors = FALSE));
      mir.edge <- na.omit(data.frame(Name1=mir.dic[,"mir_id"],ID1=mir.dic[,"mir_acc"],Name2=mir.dic[,"symbol"],ID2=mir.dic[,"entrez"],stringsAsFactors = FALSE));

      # table results
      snp.res <- na.omit(snp.dic[ , c("chr_pos", "rsid", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain")]);
      mir.res <- mir.dic[ , c("mir_id", "mir_acc", "symbol", "entrez", "experiment", "pmid", "tissue")];
      snp.res$Database <- rep("30302893", nrow(snp.res));
      colnames(snp.res) <- c("CHR_POS", "rsID", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain", "DataBase");
      colnames(mir.res) <- c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue");

      # rbind  for network
      merge.edge <- rbind(snp.edge, mir.edge);
      rownames(merge.edge) <- 1:nrow(merge.edge);

      write.csv(merge.edge, file="mirnet_snp_mir_target.csv", row.names=FALSE);


      # update the data
      gd.inx <- rownames(snp.mat) %in% unique(c(snp.edge[, "Name2"]));
      dataSet$mir.mapped <- snp.mat[gd.inx,,drop=F];

      dataSet$seeds <- c(snp.edge[, "Name2"]);
      dataSet$mir.res <- merge.edge; # for network
      dataSet$mirtarget <- c("gene");
      dataSet$mirtable <- c("snp2mir", "mir2gene");
      dataSet$snp2mir <- snp.res;
      dataSet$mir2gene <- mir.res;
      dataSet$nodeNumbers <- nrow(merge.edge)
      dataSet <<- dataSet;
      return(1);
    }
  }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param org PARAM_DESCRIPTION
#' @param idType PARAM_DESCRIPTION
#' @param tissue PARAM_DESCRIPTION
#' @param lvlOpt PARAM_DESCRIPTION
#' @param matchMin PARAM_DESCRIPTION, Default: 0.5
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformDataAnnot
#' @export
PerformDataAnnot<-function(org, idType, tissue, lvlOpt, matchMin=0.5){
    data.org <<- dataSet$org <- org;
    dataSet$id.orig <- dataSet$id.current <- idType;
    dataSet$annotation <- NULL;
    dataSet$annotated <- F;
    dataSet$tissue <- tissue;

    # should not contain duplicates, however sanity check
    data.proc <- readRDS("data.proc");
    dataSet$data.anot <- RemoveDuplicates(data.proc, "mean", quiet=T);
    feature.vec <- rownames(data.proc);

    if(idType == "mir_id" | idType=="mir_acc"){
        #do nothing
        matched.len <- length(feature.vec);
    }else{ # genes
        if(tolower(org) != 'na' & tolower(idType) != 'na'){
            anot.id <- doMirGeneAnnotation(feature.vec, idType);
            hit.inx <- !is.na(anot.id);
            matched.len <- sum(hit.inx);
            if(matched.len < length(feature.vec)*0.25){
                perct <- round(matched.len/length(feature.vec),3)*100;
                current.msg <<- paste('Only ', perct, '% ID were matched. Please choose the correct ID type or use default.', sep="");
                return(0);
            }else{
                current.msg <<- paste("ID annotation: ", "Total [", length(anot.id),
                    "] Matched [", matched.len, "] Unmatched [", sum(!hit.inx),"]", collapse="\n");

                if(lvlOpt != 'NA' | idType == "entrez"){
                    # do actual summarization to gene level
                    matched.entrez <- anot.id[hit.inx];
                    data.anot <- data.proc[hit.inx,];
                    rownames(data.anot) <- matched.entrez;
                    current.msg <<- paste(current.msg, "Data is now transformed to gene-level (Entrez) expression.");
                    dataSet$data.anot <- RemoveDuplicates(data.anot, lvlOpt, quiet=F);
                    #dataSet$id.current <- "entrez";
                    dataSet$id.current <- "symbol";
                    dataSet$annotated <- T;
                }else{
                    # record the annotation
                    dataSet$annotation <- anot.id; # this need to be updated to gether with data from now on
                    current.msg <<- paste(current.msg, "No gene level summarization was performed.");
                }
            }
        }else{ # no conversion will be performed
            matched.len <- 9; # dummies
            minLvl <- 1;
            current.msg <<- paste("No annotation was performed. Make sure organism and gene ID are specified correctly!");
        }
    }
    dataSet$data.norm <- dataSet$data.anot; # before
    write.csv(dataSet$data.anot, file="cleaned_annot.csv");
    dataSet <<- dataSet;
    return(matched.len);
}

### convert to gene symbols!!! not entrez
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param id.vec PARAM_DESCRIPTION
#' @param idType PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doMirGeneAnnotation
#' @export
doMirGeneAnnotation <- function(id.vec, idType){
     feature.vec <- id.vec;
     if(idType %in% c("entrez", "symbol", "refseq", "genbank", "emblgene", "embltranscript", "orfid","mir_id","mir_acc")){
         anot.id <- doGeneIDMapping(feature.vec, idType);
     }else{
         anot.id <- doProbeMapping(feature.vec, idType);
     }
     # convert all entrez to symbol
     anot.id <- doEntrez2SymbolMapping(anot.id);
     names(anot.id) <- id.vec;
     return(anot.id);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param id.vec PARAM_DESCRIPTION
#' @param idType PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doAnnotation
#' @export
doAnnotation <- function(id.vec, idType){
     feature.vec <- id.vec;
     if(idType %in% c("entrez", "symbol", "refseq", "genbank", "emblgene", "embltranscript", "orfid","mir_id","mir_acc")){
         anot.id <- doGeneIDMapping(feature.vec, idType);
     }else{
         anot.id <- doProbeMapping(feature.vec, idType);
     }
     names(anot.id) <- id.vec;
     return(anot.id);
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
#' @rdname PerformGeneAnnotation
#' @export
PerformGeneAnnotation <- function(){
    if(!exists("entrez.vec")){
        print("Could not find Entrez ID list!");
        return(0);
    }

    gene.map <-  queryGeneDB("entrez", data.org);
    gene.map[] <- lapply(gene.map, as.character)

    hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
    dat <- cbind(query=entrez.vec, gene.map[hit.inx, c("symbol","name")]);
    write.csv(dat, file="EntrezID2Gene.csv", row.names=F);
    rm(entrez.vec, envir = .GlobalEnv);
    return(1);
}

# probe based on built-in
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param platform PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformProbeAnnotation
#' @export
PerformProbeAnnotation <- function(platform){
    if(!exists("probe.vec")){
        print("Could not find the Probe ID list!");
        return(0);
    }
    entrez.vec <- doProbeMapping(probe.vec, platform);
    dat <- cbind(query=probe.vec, entrez=entrez.vec);
    write.csv(dat, file="Probe2Entrez.csv", row.names=F);
    rm(probe.vec, envir = .GlobalEnv);
    return(1);
}

# from probe ID to entrez ID
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param probe.vec PARAM_DESCRIPTION
#' @param platform PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doProbeMapping
#' @export
doProbeMapping <- function(probe.vec, platform){
    platform.path <- paste(lib.path,  data.org, "/", platform, ".csv", sep="");
    download.file(platform.path);
    probe.map <- read.csv(platform.path, header=T, as.is=T);
    if(is.null(probe.vec)){
        entrez <- probe.map[, "entrez"];
    }else{
        hit.inx <- match(probe.vec, probe.map[, "probe"]);
        entrez <- probe.map[hit.inx, "entrez"];
    }
    rm(probe.map);
    return(entrez);
}


# mapping between genebank, refseq and entrez
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param q.vec PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doGeneIDMapping
#' @export
doGeneIDMapping <- function(q.vec, type){
    require('RSQLite');
    mir.db <- dbConnect(SQLite(), paste(sqlite.geneid.path, data.org, "_genes.sqlite", sep=""));
    query <- paste (shQuote(q.vec),collapse=",");
    if(is.null(q.vec)){
        type.query <- paste("entrez");
        statement <- paste("SELECT * FROM entrez");
        db.map <- dbGetQuery(mir.db, statement);
        q.vec <- db.map[, "gene_id"];
        type = "entrez";
    }
    if(type == "symbol"){
        type.query <- paste("entrez");
        statement <- paste("SELECT * FROM ", type.query, " WHERE symbol IN (",query,")", sep="");
        mir.dic <-.query.sqlite(mir.db, statement);
        hit.inx <- match(q.vec, mir.dic[, "symbol"])
    }
    else if(type == "entrez"){
        type.query <- paste("entrez");
        statement <- paste("SELECT * FROM ", type.query, " WHERE gene_id IN (",query,")", sep="");
        mir.dic <- .query.sqlite(mir.db, statement);
        hit.inx <- match(q.vec, mir.dic[, "gene_id"])
    }
    else{
        # note, some ID can have version number which is not in the database
        # need to strip it off NM_001402.5 => NM_001402
        q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
        q.vec <- q.mat[,1];

        if(type == "genbank"){
            type.query <- paste("entrez_gb");
        }else if(type == "refseq"){
            type.query <- paste("entrez_refseq");
        }else if(type == "emblgene"){
            type.query <- paste("entrez_embl_gene");
        }else if(type == "embltranscript"){
            type.query <- paste("entrez_embl_transcript");
        }else if(type == "orfid"){ # only for yeast
            type.query <- paste("entrez_orf");
        }else{
            print("Unknown data type");
            return(0);
        }
        statement <- paste("SELECT * FROM ", type.query, " WHERE accession IN (",query,")", sep="");
        mir.dic <- .query.sqlite(mir.db, statement);
        hit.inx <- match(q.vec, mir.dic[, "accession"])
    }
    entrezs=mir.dic[hit.inx, "gene_id"];
    mode(entrezs) <- "character";
    rm(mir.dic, q.vec); gc();
    return(entrezs);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param entrez.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doEntrez2SymbolMapping
#' @export
doEntrez2SymbolMapping <- function(entrez.vec){
    gene.map <-  queryGeneDB("entrez", data.org);
    gene.map[] <- lapply(gene.map, as.character)

    hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
    symbols <- gene.map[hit.inx, "symbol"];

    # if not gene symbol, use id by itself
    na.inx <- is.na(symbols);
    symbols[na.inx] <- entrez.vec[na.inx];
    return(symbols);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param entrez.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doSymbol2EntrezMapping
#' @export
doSymbol2EntrezMapping <- function(entrez.vec){
    gene.map <-  queryGeneDB("entrez", data.org);
    gene.map[] <- lapply(gene.map, as.character)

    hit.inx <- match(entrez.vec, gene.map[,"symbol"]);
    symbols <- gene.map[hit.inx, "gene_id"];

    # if not gene symbol, use id by itself
    na.inx <- is.na(symbols);
    symbols[na.inx] <- entrez.vec[na.inx];
    return(symbols);
}

# given a data with duplicates, dups is the one with duplicates
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data PARAM_DESCRIPTION
#' @param lvlOpt PARAM_DESCRIPTION
#' @param quiet PARAM_DESCRIPTION, Default: T
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname RemoveDuplicates
#' @export
RemoveDuplicates <- function(data, lvlOpt, quiet=T){

    all.nms <- rownames(data);
    colnms <- colnames(data);
    dup.inx <- duplicated(all.nms);
    dim.orig  <- dim(data);
    data <- apply(data, 2, as.numeric); # force to be all numeric
    dim(data) <- dim.orig; # keep dimension (will lost when only one item)
    rownames(data) <- all.nms;
    colnames(data) <- colnms;
    if(sum(dup.inx) > 0){
        uniq.nms <- all.nms[!dup.inx];
        uniq.data <- data[!dup.inx,,drop=F];

        dup.nms <- all.nms[dup.inx];
        uniq.dupnms <- unique(dup.nms);
        uniq.duplen <- length(uniq.dupnms);

        for(i in 1:uniq.duplen){
            nm <- uniq.dupnms[i];
            hit.inx.all <- which(all.nms == nm);
            hit.inx.uniq <- which(uniq.nms == nm);

            # average the whole sub matrix
            if(lvlOpt == "mean"){
                uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, mean, na.rm=T);
            }else if(lvlOpt == "median"){
                uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, median, na.rm=T);
            }else if(lvlOpt == "max"){
                uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, max, na.rm=T);
            }else{ # sum
                uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, sum, na.rm=T);
            }
        }
        if(!quiet){
            current.msg <<- paste(current.msg, paste("A total of ", sum(dup.inx), " of duplicates were replaced by their ", lvlOpt, ".", sep=""), collapse="\n");
        }
        return(uniq.data);
    }else{
        if(!quiet){
            current.msg <<- paste(current.msg, "All IDs are unique.", collapse="\n");
        }
        return(data);
    }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param table.nm PARAM_DESCRIPTION
#' @param data.org PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname queryGeneDB
#' @export
queryGeneDB <- function(table.nm, data.org){
    require('RSQLite')

    conv.db <- dbConnect(SQLite(), paste(sqlite.geneid.path, data.org, "_genes.sqlite", sep=""));
    db.map <- dbReadTable(conv.db, table.nm)
    dbDisconnect(conv.db); cleanMem();

    return(db.map)
}
