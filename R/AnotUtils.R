##################################################
## R script for miRNet
## Description: Gene/Compound Annotation
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

PerformMirGeneMapping <- function(input.type="none"){
    if(input.type %in% c("mir2gene_mirtarbase", "mir2gene_tarbase", "mir2gene_mirecords", "mir2gene_miranda")){
      db.type <- gsub("mir2gene_", "", input.type);
    }else{
      db.type <- "mirtarbase";
    }
    mir.mat <- dataSet$mir.orig;

    mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2gene", sep=""), rownames(mir.mat), dataSet$org, dataSet$idType, db.type);

    hit.num <- nrow(mir.dic);
    if (hit.num == 0){
        if(dataSet$tissue == "na") {
            current.msg <<- "No hits found in the database. Please check your input.";
        }else{
            current.msg <<- "No hits found in the database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
        }
        print(current.msg);
        return(0);
    } else {
        current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-gene targets were identified!");
        fast.write.csv(mir.dic, file="mirnet_mir_gene.csv", row.names=FALSE); # this is just for mir2gene results table to show different db source with 0/1 to indicate present
        res <- mir.dic[ , c("mir_id", "mir_acc", "symbol", "entrez", "experiment", "pmid", "tissue")];
        rownames(res) <- mir.dic$mirnet;
        colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue");
        fast.write.csv(res, file="mirnet_mir_target.csv", row.names=FALSE);
        dataSet$mir.res <- res;

        gene.nms <- res[,"Gene"];
        net.info$gene.nms <- gene.nms;
        net.info <<-net.info;

        # record the mapped queries and change to same IDs used in network
        uniq.mat <- unique(mir.dic[, c("mir_id", "symbol", dataSet$idType)]);
        hit.inx <- match(rownames(mir.mat), uniq.mat[, dataSet$idType]);
        if(dataSet$idType %in% c("mir_id", "mir_acc")){
            rownames(mir.mat) <- uniq.mat[hit.inx,"mir_id"];
            dataSet$seeds <- res[, "ID"];
        }else{
            rownames(mir.mat) <- uniq.mat[hit.inx,"symbol"];
            dataSet$seeds <- gene.nms;
        }
        dataSet$mir.mapped <- mir.mat;
        dataSet$mirtable <- "mir2gene"
        dataSet$mir2gene <- res
        dataSet$mirtarget <- "gene";
        dataSet <<- dataSet;
        if(.on.public.web){
          return(1);
        }else{
          return(current.msg);
        }
    }
}

PerformMolMapping <- function(){
  mir.mat <- dataSet$mir.orig;

  mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2molecule", sep=""), rownames(mir.mat), dataSet$org, dataSet$idType);

  hit.num <- nrow(mir.dic);
  if (hit.num == 0){
    if(dataSet$tissue == "na") {
      current.msg <<- "No hits found in the database. Please check your input.";
    }else{
      current.msg <<- "No hits found in the database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
    }
    print(current.msg);
    return(0);
  } else {
    current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-molecule interactions were identified!");

    res <- mir.dic[ , c("mir_id","mir_acc","molecule", "pubchem_id", "method", "pmid", "tissue")];
    rownames(res) <- mir.dic$mirnet;
    colnames(res) <- c("ID","Accession","Molecule", "Pubchem_ID", "Experiment", "Literature", "Tissue");
    fast.write.csv(res, file="mirnet_mir_target.csv", row.names=FALSE);
    dataSet$mir.res <- res;

    mol.nms <- res[,"Molecule"];
    net.info$mol.nms <- mol.nms;
    net.info <<-net.info;

    # record the mapped queries and change to same IDs used in network
    uniq.mat <- unique(mir.dic[, c("mir_id", "molecule", dataSet$idType)]);
    hit.inx <- match(rownames(mir.mat), uniq.mat[, dataSet$idType]);
    if(dataSet$idType %in% c("mir_id", "mir_acc")){
      rownames(mir.mat) <- uniq.mat[hit.inx,"mir_id"];
      dataSet$seeds <- res[, "ID"];
    }else{
      rownames(mir.mat) <- uniq.mat[hit.inx,"molecule"];
      dataSet$seeds <- mol.nms;
    }
    dataSet$mir.mapped <- mir.mat;
    dataSet$mirtable <- "mir2mol"
    dataSet$mir2mol <- res
    dataSet$mirtarget <- "molecule";
    dataSet <<- dataSet;
    if(.on.public.web){
      return(1);
    }else{
      return(current.msg);
    }
  }
}

PerformLncRNAMapping <- function(){
  orgType <- dataSet$org;
  if(orgType != "hsa" ){
    curent.msg <<- "Only human supports lncRNA network."
    print(current.msg);
    return(0);
  }

  mir.mat <- dataSet$mir.orig;
  idType <- dataSet$idType;
  mir.vec <- rownames(mir.mat);
  mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2lncRNA", sep=""), mir.vec, orgType, idType);

  hit.num <- nrow(mir.dic)
  if (hit.num == 0 && dataSet$tissue == "na") {
    current.msg <<- "No hits found in the database. Please check your input. ";
    print(current.msg);
    return(0);
  } else if (hit.num == 0 && dataSet$tissue != "na") {
    current.msg <<- "No hits found in the database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
    print(current.msg);
    return(2);
  } else {
    res <- mir.dic[ , c("mir_id","mir_acc","symbol","entrez", "tissue")];
    rownames(res) <- mir.dic$mirnet;
    current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-lncRNA targets were identified!");

    # update the data
    gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
    dataSet$mir.mapped <- mir.mat[gd.inx,,drop=F];

    colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Tissue");
    res$Experiment <- rep("CLIP-Seq", nrow(res));
    res$Literature <- rep("24297251", nrow(res));
    res <- res[, c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue")];
    lnc.nms <- res[,"Gene"];
    mir.nms <- res[, "ID"];
    if(dataSet$idType %in% c("mir_id", "mir_acc")){
      dataSet$seeds <- mir.nms;
    }else{
      dataSet$seeds <- lnc.nms;
    }
    lnc.nms <- res[,"Gene"];
    net.info$lnc.nms <- lnc.nms;
    net.info <<-net.info;
    fast.write.csv(res, file="mirnet_mir_target.csv", row.names=FALSE);
    dataSet$mir.res <- res;
    dataSet$mirtarget <- "lncrna";
    dataSet$mirtable <- "mir2lnc"
    dataSet$mir2lnc <- res
    dataSet <<- dataSet;
    if(.on.public.web){
      return(1);
    }else{
      return(current.msg);
    }
  }
}

PerformCircRNAMapping <- function(){
  orgType <- dataSet$org;
  if(orgType != "hsa" ){
    curent.msg <<- "Only human supports circRNA network."
    print(current.msg);
    return(0);
  }

  mir.mat <- dataSet$mir.orig;
  idType <- dataSet$idType;
  mir.vec <- rownames(mir.mat);
  mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2circRNA", sep=""), mir.vec, orgType, idType);

  hit.num <- nrow(mir.dic)
  if (hit.num == 0 && dataSet$tissue == "na") {
    current.msg <<- "No hits found in the database. Please check your input. ";
    print(current.msg);
    return(0);
  } else if (hit.num == 0 && dataSet$tissue != "na") {
    current.msg <<- "No hits found in the database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
    print(current.msg);
    return(2);
  } else {
    res <- mir.dic[ , c("mir_id","mir_acc","symbol","entrez", "tissue")];
    rownames(res) <- mir.dic$mirnet;
    current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-circRNA targets were identified!");

    # update the data
    gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
    dataSet$mir.mapped <- mir.mat[gd.inx,,drop=F];

    colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Tissue");
    res$Experiment <- rep("CLIP-Seq", nrow(res));
    res$Literature <- rep("24297251", nrow(res));
    res <- res[, c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue")];
    circ.nms <- res[,"Gene"];
    mir.nms <- res[, "ID"];
    if(dataSet$idType %in% c("mir_id", "mir_acc")){
      dataSet$seeds <- mir.nms;
    }else{
      dataSet$seeds <- circ.nms;
    }
    circ.nms <- res[,"Gene"];
    net.info$circ.nms <- circ.nms;
    net.info <<-net.info;
    fast.write.csv(res, file="mirnet_mir_target.csv", row.names=FALSE);
    dataSet$mir.res <- res;
    dataSet$mirtarget <- "circrna";
    dataSet$mirtable <- "mir2circ"
    dataSet$mir2circ <- res
    dataSet <<- dataSet;
    if(.on.public.web){
      return(1);
    }else{
      return(current.msg);
    }
  }
}

PerformPseudoMapping <- function(){
  orgType <- dataSet$org;
  if(orgType != "hsa" ){
    curent.msg <<- "Only human supports pseudogene network."
    print(current.msg);
    return(0);
  }

  mir.mat <- dataSet$mir.orig;
  idType <- dataSet$idType;
  mir.vec <- rownames(mir.mat);
  mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2pseudogene", sep=""), mir.vec, orgType, idType);

  hit.num <- nrow(mir.dic)
  if (hit.num == 0 && dataSet$tissue == "na") {
    current.msg <<- "No hits found in the database. Please check your input. ";
    print(current.msg);
    return(0);
  } else if (hit.num == 0 && dataSet$tissue != "na") {
    current.msg <<- "No hits found in the database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
    print(current.msg);
    return(2);
  } else {
    res <- mir.dic[ , c("mir_id","mir_acc","symbol","entrez", "tissue")];
    rownames(res) <- mir.dic$mirnet;
    current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-pseudogene targets were identified!");

    # update the data
    gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
    dataSet$mir.mapped <- mir.mat[gd.inx,,drop=F];

    colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Tissue");
    res$Experiment <- rep("CLIP-Seq", nrow(res));
    res$Literature <- rep("24297251", nrow(res));
    res <- res[, c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue")];
    pseudo.nms <- res[,"Gene"];
    mir.nms <- res[, "ID"];
    if(dataSet$idType %in% c("mir_id", "mir_acc")){
      dataSet$seeds <- mir.nms;
    }else{
      dataSet$seeds <- pseudo.nms;
    }
    pseudo.nms <- res[,"Gene"];
    net.info$pseudo.nms <- pseudo.nms;
    net.info <<-net.info;
    fast.write.csv(res, file="mirnet_mir_target.csv", row.names=FALSE);
    dataSet$mir.res <- res;
    dataSet$mirtarget <- "pseudogene";
    dataSet$mirtable <- "mir2pseudo"
    dataSet$mir2pseudo <- res
    dataSet <<- dataSet;
    if(.on.public.web){
      return(1);
    }else{
      return(current.msg);
    }
  }
}

PerformSncRNAMapping <- function(){
  orgType <- dataSet$org;
  if(orgType != "hsa" ){
    curent.msg <<- "Only human supports sncRNA network."
    print(current.msg);
    return(0);
  }

  mir.mat <- dataSet$mir.orig;
  idType <- dataSet$idType;
  mir.vec <- rownames(mir.mat);
  mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2sncRNA", sep=""), mir.vec, orgType, idType);

  hit.num <- nrow(mir.dic)
  if (hit.num == 0 && dataSet$tissue == "na") {
    current.msg <<- "No hits found in the database. Please check your input. ";
    print(current.msg);
    return(0);
  } else if (hit.num == 0 && dataSet$tissue != "na") {
    current.msg <<- "No hits found in the database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
    print(current.msg);
    return(2);
  } else {
    res <- mir.dic[ , c("mir_id","mir_acc","symbol","entrez", "tissue")];
    rownames(res) <- mir.dic$mirnet;
    current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-sncRNA targets were identified!");

    # update the data
    gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
    dataSet$mir.mapped <- mir.mat[gd.inx,,drop=F];

    colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Tissue");
    res$Experiment <- rep("CLIP-Seq", nrow(res));
    res$Literature <- rep("24297251", nrow(res));
    res <- res[, c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue")];
    snc.nms <- res[,"Gene"];
    mir.nms <- res[, "ID"];
    if(dataSet$idType %in% c("mir_id", "mir_acc")){
      dataSet$seeds <- mir.nms;
    }else{
      dataSet$seeds <- snc.nms;
    }
    snc.nms <- res[,"Gene"];
    net.info$snc.nms <- snc.nms;
    net.info <<-net.info;
    fast.write.csv(res, file="mirnet_mir_target.csv", row.names=FALSE);
    dataSet$mir.res <- res;
    dataSet$mirtarget <- "sncrna";
    dataSet$mirtable <- "mir2snc"
    dataSet$mir2snc <- res
    dataSet <<- dataSet;
    if(.on.public.web){
      return(1);
    }else{
      return(current.msg);
    }
  }
}

PerformMir2EpiMapping <- function(){
  orgType <- dataSet$org;
  if(orgType %in% c("bta", "dme","gga","sma", "cel","dre","rno", "ssc") ){
    curent.msg <<- "Only huamn and mouse support the epigene network."
    print(current.msg);
    return(0);
  }

  mir.mat <- dataSet$mir.orig;

  idType <- dataSet$idType;
  mir.vec <- rownames(mir.mat);
  print("Perform epi2mir");
  print(dataSet$tissue);
  mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2epi", sep=""), mir.vec, orgType, idType);

  hit.num <- nrow(mir.dic)
  if (hit.num == 0 && dataSet$tissue == "na") {
    current.msg <<- "No hits found in the database. Please check your input. ";
    print(current.msg);
    return(0);
  } else if (hit.num == 0 && dataSet$tissue != "na") {
    current.msg <<- "No hits found in the database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
    print(current.msg);
    return(2);
  } else {
    res <- mir.dic[ , c("mir_id", "mir_acc", "epi_regulator", "experiment", "condition", "pmid", "tissue")];
    rownames(res) <- mir.dic$mirnet;
    current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-epigene targets were identified!");

    # update the data
    gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
    dataSet$mir.mapped <- mir.mat[gd.inx,,drop=F];

    colnames(res) <- c("ID","Accession","Epigenetics","Experiment", "Condition","Literature", "Tissue");
    epi.nms <- res[,"Epigenetics"];
    mir.nms <- res[, "ID"];
    if(dataSet$idType %in% c("mir_id", "mir_acc")){
      dataSet$seeds <- mir.nms;
    }else{
      dataSet$seeds <- epi.nms;
    }
    epi.nms <- res[,"Epigenetics"];
    net.info$epi.nms <- epi.nms;
    net.info <<-net.info;
    fast.write.csv(res, file="mirnet_mir_target.csv", row.names=FALSE);
    dataSet$mir.res <- res;
    dataSet$mirtarget <- "epigenetics";
    dataSet$mirtable <- "mir2epi"
    dataSet$mir2epi <- res
    dataSet <<- dataSet;
    if(.on.public.web){
      return(1);
    }else{
      return(current.msg);
    }
  }
}

PerformDisMapping <- function(){
  if(dataSet$org != "hsa" ){
    curent.msg <<- "Only huamn support the disease network."
    print(current.msg);
    return(0);
  }

  mir.mat <- dataSet$mir.orig;
  idType <- dataSet$idType;
  mir.vec <- rownames(mir.mat);
  mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2disease", sep=""), mir.vec, "disease", idType);

  hit.num <- nrow(mir.dic)
  if (hit.num == 0 && dataSet$tissue == "na") {
    current.msg <<- "No hits found in the database. Please check your input. ";
    print(current.msg);
    return(0);
  } else if (hit.num == 0 && dataSet$tissue != "na") {
    current.msg <<- "No hits found in the database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
    print(current.msg);
    return(2);
  } else{
    res <- mir.dic[ , c("mir_id", "mir_acc", "disease", "method", "database", "pmid", "tissue")];
    rownames(res) <- mir.dic$mirnet;
    current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-disease associations were identified!");

    # update the data
    gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
    dataSet$mir.mapped <- mir.mat[gd.inx,,drop=F];

    colnames(res) <- c("ID","Accession","Disease","Experiment", "Database", "Literature", "Tissue");
    dis.nms <- res[,"Disease"];
    mir.nms <- res[, "ID"];
    if(dataSet$idType %in% c("mir_id", "mir_acc")){
      dataSet$seeds <- mir.nms;
    }else{
      dataSet$seeds <- dis.nms;
    }
    dis.nms <- res[,"Disease"];
    net.info$dis.nms <- dis.nms;
    net.info <<-net.info;
    fast.write.csv(res, file="mirnet_mir_target.csv", row.names=FALSE);
    dataSet$mir.res <- res;
    dataSet$mirtarget <- "disease";
    dataSet$mirtable <- "mir2dis"
    dataSet$mir2dis <- res
    dataSet <<- dataSet;
    if(.on.public.web){
      return(1);
    }else{
      return(current.msg);
    }
  }
}

PerformTFMapping <- function(){
  orgType <- dataSet$org;
  if(orgType %in% c("bta", "ssc","gga","dme", "sma") ){
    curent.msg <<- "This organism is not supported for transcription factors network research."
    print(current.msg);
    return(0);
  }

  mir.mat <- dataSet$mir.orig;
  idType <- dataSet$idType;
  mir.vec <- rownames(mir.mat);
  mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2tf", sep=""), mir.vec, orgType, idType);

  hit.num <- nrow(mir.dic)
  if (hit.num == 0 && dataSet$tissue == "na") {
    current.msg <<- "No hits found in the database. Please check your input. ";
    print(current.msg);
    return(0);
  } else if (hit.num == 0 && dataSet$tissue != "na") {
    current.msg <<- "No hits found in the database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
    print(current.msg);
    return(2);
  } else {
    res <- mir.dic[ , c("mir_id","mir_acc","symbol","entrez", "pmid", "tissue")];
    rownames(res) <- mir.dic$mirnet;
    current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-TF targets were identified!");

    # update the data
    gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
    dataSet$mir.mapped <- mir.mat[gd.inx,,drop=F];

    colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Literature", "Tissue");
    res$Experiment <- rep("ChIP-seq", nrow(res));
    res <- res[, c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue")];
    tf.nms <- res[,"Gene"];
    mir.nms <- res[, "ID"];
    if(dataSet$idType %in% c("mir_id", "mir_acc")){
      dataSet$seeds <- mir.nms;
    }else{
      dataSet$seeds <- tf.nms;
    }
    tf.nms <- res[,"Gene"];
    net.info$tf.nms <- tf.nms;
    net.info <<-net.info;
    fast.write.csv(res, file="mirnet_mir_target.csv", row.names=FALSE);
    dataSet$mir.res <- res;
    dataSet$mirtarget <- "tf";
    dataSet$nodeNumbers <- nrow(res);
    dataSet$mirtable <- "mir2tf"
    dataSet$mir2tf <- res
    dataSet <<- dataSet;
    if(.on.public.web){
      return(1);
    }else{
      return(current.msg);
    }
  }
}

PerformSNPMirGeneMapping <- function(){
     if(!exists("my.snp.mir.mapping")){ # public web on same user dir
        compiler::loadcmp("../../rscripts/miRNetR/R/utils_mir_snp.Rc");    
    }
    return(my.snp.mir.mapping());
}


### convert to gene symbols!!! not entrez
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

PerformGeneAnnotation <- function(){
    if(!exists("entrez.vec")){
        print("Could not find Entrez ID list!");
        return(0);
    }

    gene.map <-  queryGeneDB("entrez", data.org);
    gene.map[] <- lapply(gene.map, as.character)

    hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
    dat <- cbind(query=entrez.vec, gene.map[hit.inx, c("symbol","name")]);
    fast.write.csv(dat, file="EntrezID2Gene.csv", row.names=F);
    rm(entrez.vec, envir = .GlobalEnv);
    return(1);
}

# from probe ID to entrez ID
doProbeMapping <- function(probe.vec, platform){
    platform.path <- paste(lib.path,  data.org, "/", platform, ".csv", sep="");
    if(.on.public.web){
      probe.map <- read.csv(platform.path, header=T, as.is=T);
    }else{
      destfile <- paste(platform, ".csv", sep="");
      download.file(platform.path, destfile);
      probe.map <- read.csv(destfile, header=T, as.is=T);
    }
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
doGeneIDMapping <- function(q.vec, type){
    require('RSQLite');
    db.path <- paste(sqlite.path, data.org, "_genes.sqlite", sep="");
    if(.on.public.web){
        mir.db <- dbConnect(SQLite(), db.path);
    }else{
        msg <- paste("Downloading", db.path);
        db.name <- gsub(sqlite.path, "", db.path);
        if(!file.exists(db.name)){
          print(msg);
          download.file(db.path, db.name, mode = "wb");
        }
        mir.db <- dbConnect(SQLite(), db.name);
    }
  
    #mir.db <- dbConnect(SQLite(), paste(sqlite.path, data.org, "_genes.sqlite", sep=""));
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

queryGeneDB <- function(table.nm, data.org){
    require('RSQLite')
    
    db.path <- paste(sqlite.path, data.org, "_genes.sqlite", sep="")
    if(.on.public.web){
      conv.db <- dbConnect(SQLite(), db.path);
    }else{
      msg <- paste("Downloading", db.path);
      db.name <- gsub(sqlite.path, "", db.path);
      if(!file.exists(db.name)){
        print(msg);
        download.file(db.path, db.name, mode = "wb");
      }
      conv.db <- dbConnect(SQLite(), db.name);
    }
    #conv.db <- dbConnect(SQLite(), paste(sqlite.path, data.org, "_genes.sqlite", sep=""));
    db.map <- dbReadTable(conv.db, table.nm)
    dbDisconnect(conv.db); cleanMem();

    return(db.map)
}
