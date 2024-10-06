##################################################
## R script for miRNet
## Description: List data I/O and processing
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

#' Setup miRNA List Data
#' @param mirs A list of miRNAs.
#' @param orgType Organism type.
#' @param idType miRNA ID type.
#' @param tissue Tissue type.
#' @param targetOpt Target options (optional).
#' @return miRNA list data initialized.
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
        mir.vec <- gsub("mir", "miR", mir.vec);
    }
    if(idType == "mir_id"){
      rownames(mir.mat) <-  tolower(as.vector(mir.vec));
    }else{
      rownames(mir.mat) <-  mir.vec;
    }

    mir.mat <- mir.mat[,-1, drop=F];
    dataSet$mir.orig <- mir.mat;
    dataSet$data[["miRNA"]] <- mir.mat;
    dataSet$id.types[["miRNA"]] <- idType;
    dataSet$target.types[["miRNA"]] <- targetOpt;
    dataSet<<- dataSet;
    if(.on.public.web){
      return (nrow(mir.mat));
    }else{
      return (paste("A total of",  nrow(mir.mat), "unique items were entered."))
    }
}

#' Setup List Data
#' @param listInput List to query.
#' @param orgType Organism type.
#' @param inputType Input type.
#' @param idType ID type.
#' @param tissue Tissue type.
#' @param target Target.
#' @return List data initialized.
#' @export
SetupIndListData <- function(listInput, orgType, inputType, idType, tissue, target){

print(paste0("================",  orgType));
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
  dataSet$target.types[[inputType]] <- target;
  dataSet <<- dataSet;
  if(.on.public.web){
print(paste0("================123",  orgType));
    return (nrow(in.mat));
  }else{
    return (paste("A total of",  nrow(in.mat), "unique items were entered."))
  }
}

#' Setup Item From Picklist
#' @param orgType Organism type.
#' @param tissue Tissue type.
#' @param idType ID type.
#' @export
SetupItemFromPickList <- function(orgType="hsa", tissue, idType){
    if(!exists("picklist.vec")){
        print("Could not find user entered disease list!");
        return(0);
    }
    dataSet$org <- orgType;
    dataSet$tissue <- tissue;
    mir.mat <- .parsePickListItems(picklist.vec);

    if(anal.type == "multilist" || anal.type == "dis2mir" ){
        dataSet$data[[idType]] <- mir.mat;
    }else{
        dataSet$data[[idType]] <- mir.mat;
        dataSet$mir.orig <- mir.mat;
        dataSet$idType <- idType;
    }
    picklist.vec <<- NULL;
    dataSet <<- dataSet;
    return(nrow(mir.mat));
}

#' Get Unique Disease Names
#' @export
GetUniqueDiseaseNames <- function(){
    db.path <- paste(sqlite.path, "mir2disease.sqlite", sep="");
    statement <- "SELECT disease FROM disease";
    return(GetUniqueEntries(db.path, statement));
}

#' Get Unique Molecule Names
#' @param orgType Organism type.
#' @export
GetUniqueMoleculeNames <- function(orgType="hsa"){
    db.path <- paste(sqlite.path, "mir2molecule.sqlite", sep="");
    statement <- paste("SELECT molecule FROM ",orgType, sep="");
    return(GetUniqueEntries(db.path, statement));
}

#' Get Unique Epigene Names
#' @param orgType Organism type.
#' @export
GetUniqueEpigeneNames <- function(orgType="hsa"){
    db.path <- paste(sqlite.path, "mir2epi.sqlite", sep="");
    statement <- paste("SELECT epi_regulator FROM ",orgType, sep="");
    return(GetUniqueEntries(db.path, statement));
}

#' Set Current Data Multi
#' @export
SetCurrentDataMulti <- function(){
  dataSet$type <- nms.vec;
  dataSet <<- dataSet;
  if(.on.public.web){
    return(1);
  }else{
    return (paste("Targets were entered!"))
  }
}

.init.multilist <- function(){
  anal.type <<- "multilist"
  net.info <<- list();
  mir.mappedu <<- matrix();
  mir.resu <<- data.frame();
  mirtargetu <<- vector();
  mirtableu <<- vector();
  seedsu <<- vector();
  dataSet$directionInx <-vector()
  dataSet$regDirection <-vector()
  dataSet$tfTargetType <-vector()
  edgeNumU <<- vector();
  edgeu.res <<- data.frame();
  nodeu.ids <<- vector();
  nodeu.nms <<- vector();
  mir.nmsu <<- vector();
  tf.nms <<- vector();
  gene.nms <<- vector()
  dataSet<<- dataSet

}

#' Query Multi List
#' @export
QueryMultiList <- function(){
  .init.multilist();
  for(i in 1:length(dataSet$type)){
    if (dataSet$type[i] == "mirna") {
      input.type=paste("mir2", dataSet$targetOpt, sep="");
    } else{
      if(dataSet$type[i] == "protein2protein"){
        input.type= "protein2protein"
      }else if(grepl("2", dataSet$type[i])){
        input.type=dataSet$type[i];
      }else{
        if(dataSet$type[i] != "tf"){
            input.type=paste(dataSet$type[i], "2mir", sep="");
        }
      }
    }
    res <- SearchMultiNet(input.type);
    if(res != 1) return(0);
    net.info <<- .set.net.names(input.type);
  }

  typesu <<- dataSet$type;

  if(length(nodeu.ids) == 0){
    current.msg <<- paste("No interactions have been detected the given interaction types");
    return(c(0,0,1));
  }

  node.res <- data.frame(Id=nodeu.ids, Label=nodeu.nms);
  un.inx= !duplicated(node.res$Id)
  node.res <- node.res[un.inx,];

  edgeu.res[edgeu.res== "<NA>"] = 'NA'
  edgeu.res = na.omit(edgeu.res);

  fast.write.csv(node.res, file="orig_node_list.csv", row.names=FALSE);
  fast.write.csv(mir.resu, file="mirnet_mir_target.csv", row.names=FALSE);

  dataSet$mir.mapped <- na.omit(mir.mappedu);
  dataSet$mir.res <- mir.resu;
  dataSet$mirtarget <- mirtargetu;
  dataSet$mirtable <- unique(mirtableu);
  dataSet$seeds <- seedsu;
  dataSet$nodeNumbers <-edgeNumU
  dataSet <<- dataSet;
  net.info <<- net.info;
  multi.net <<- list(
    node.data = node.res,
    edge.data = edgeu.res
  );
  if(.on.public.web){
    return(1);
  }else{
    return(current.msg);
  }
}

.set.net.names <- function(input.type){

      if (grepl("gene", input.type)) {
        net.info$gene.nms = gene.nms
      }
      if (grepl("lnc", input.type)) {
        net.info$lnc.nms = lnc.nms
      }
      if (grepl("tf", input.type)) {
        net.info$tf.nms = tf.nms
      }
      if (grepl("dis", input.type)) {
        net.info$dis.nms = dis.nms
      }
      if (grepl("mol", input.type)) {
        net.info$mol.nms = mol.nms
      }
      if (grepl("epi", input.type)) {
        net.info$epi.nms = epi.nms
      }
      if (grepl("circ", input.type)) {
        net.info$circ.nms = circ.nms
      }
      if (grepl("pseudo", input.type)) {
        net.info$pseudo.nms = pseudo.nms
      }
      if (grepl("snc", input.type)) {
        net.info$snc.nms = snc.nms
      }
      if (grepl("snp", input.type)) {
        net.info$snp.nms = snp.nms
      }
  return(net.info)
}

#' Query Multi List miRNA
#' @export
QueryMultiListMir <- function(){
  .init.multilist();
  for(i in 1:length(dataSet$type)){
    if(dataSet$type[i] == "protein2protein"){
      input.type=dataSet$type[i]
    }else{
      input.type=paste("mir2", dataSet$type[i], sep="");
    }
    res <- SearchMultiNet(input.type);
    if(res != 1) return(0);
    net.info <<- .set.net.names(dataSet$type[i]);
  }
  typesu <<- dataSet$type;

  if(length(nodeu.ids) == 0){
    current.msg <<- paste("No interactions have been detected the given interaction types");
    return(c(0,0,1));
  }

  node.res <- data.frame(Id=nodeu.ids, Label=nodeu.nms);
  un.inx= !duplicated(node.res$Id)
  node.res <- node.res[un.inx,];

  edgeu.res[edgeu.res== "<NA>"] = 'NA'
  edgeu.res = na.omit(edgeu.res);

  fast.write.csv(node.res, file="orig_node_list.csv", row.names=FALSE);
  fast.write.csv(mir.resu, file="mirnet_mir_target.csv", row.names=FALSE);

  dataSet$mir.mapped <- na.omit(mir.mappedu);
  dataSet$mir.res <- mir.resu;
  dataSet$mirtarget <- mirtargetu;
  dataSet$mirtable <- mirtableu;
  dataSet$seeds <- seedsu;
  dataSet$nodeNumbers <-edgeNumU
  dataSet <<- dataSet;
  net.info <<- net.info;
  multi.net <<- list(
    node.data = node.res,
    edge.data = edgeu.res
  );
  if(.on.public.web){
    return(1);
  }else{
    return(current.msg);
  }
}

.searchMultiNet_mir2gene<-function(input.type){
    if (input.type == "mir2gene"){
      idType <- dataSet$mirnaType;
      mir.mat <- dataSet$mir.orig;
      idVec <- rownames(mir.mat);
      db.type <- "mirtarbase";
      orgType <- dataSet$org;
    }else if(input.type %in% c("mir2gene_mirtarbase", "mir2gene_tarbase", "mir2gene_mirecords", "mir2gene_miranda")){
      idType <- dataSet$mirnaType;
      mir.mat <- dataSet$mir.orig;
      idVec <- rownames(mir.mat);
      db.type <- gsub("mir2gene_", "", input.type);
      orgType <- dataSet$org;
    }else if(input.type=="snpmir2gene"){
      idType <- "mir_id";
      mir.mat <- data.frame(unique(mir.nms))
      rownames(mir.mat) = unique(mir.nms)
      mir.mat[,1] = 0
      idVec <- rownames(mir.mat);
      db.type <- "";
      orgType <- dataSet$org;
    }else{
      if(length(dataSet$data[["gene"]])!= 0){
        idType <- dataSet$id.types[["gene"]];
        mir.mat <- dataSet$data[["gene"]];
        idVec <- rownames(mir.mat);
        db.type <- "";
        orgType <- dataSet$org;
      }else{
        idType <- dataSet$mirnaType;
        mir.mat <- dataSet$mir.orig;
        idVec <- rownames(mir.mat);
        db.type <- "";
        orgType <- dataSet$org;
      }
    }
    mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2gene", sep=""), idVec, orgType, idType, db.type);

    hit.num <- nrow(mir.dic)
    if (hit.num == 0 && dataSet$tissue == "na") {
      current.msg <<- paste("No hits found in the", db.type,"miRNA-gene database. Please check your input.", sep = " ");
      print(current.msg);
      return(0);
    } else if (hit.num == 0 && dataSet$tissue != "na") {
      current.msg <<- paste("No hits found in the", db.type,"miRNA-gene database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.", sep = " ");
      print(current.msg);
      return(2);
    } else {
      fast.write.csv(mir.dic, file="mirnet_mir_gene.csv", row.names=FALSE); # this is just for mir2gene results table to show different db source with 0/1 to indicate present
      res <- mir.dic[ , c("mir_id", "mir_acc", "symbol", "entrez", "experiment", "pmid", "tissue")];
      rownames(res) <- mir.dic$mirnet;
      current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-gene targets were identified from", db.type, sep = " ");

      # record the mapped queries and change to same IDs used in network
      uniq.mat <- unique(mir.dic[, c("mir_id", "symbol", idType)]);
      hit.inx <- match(rownames(mir.mat), uniq.mat[, idType]);
      if(idType %in% c("mir_id", "mir_acc")){
        rownames(mir.mat) <- uniq.mat[hit.inx,"mir_id"];
      }else{
        rownames(mir.mat) <- uniq.mat[hit.inx,"symbol"];
      }

      mir.mappedu <<- rbind(mir.mappedu, mir.mat);

      # update col names
      colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue");
      display.res <- res;
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Entrez"]);    # IDs
      if(nrow(res)!=0){
        row.names(edge.res) <- 1:nrow(res);
      }

      node.ids <- c(ID1=res[,"Accession"], ID2=res[,"Entrez"]);
      node.nms <- c(Name1=res[,"ID"], Name2=res[,"Gene"]);

      edgeu.res <<- rbind(edgeu.res, edge.res); #edgeu.res is an empty dataframe defined in QueryNet
      edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ];
      nodeu.ids <<- c(nodeu.ids, node.ids);
      edgeNumU <<- c(edgeNumU, nrow(edge.res))
      nodeu.nms <<- c(nodeu.nms, node.nms);

      gene.nms <<- c(gene.nms, as.vector(res[,"Gene"]));
      mir.nms <<- res[, "ID"];
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"]);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "gene");

      if(input.type %in% c("mir2gene", "mir2gene_mirtarbase", "mir2gene_tarbase", "mir2gene_mirecords", "mir2gene_miranda")){
        dataSet$mir2gene <<- mir.resu;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(gene.nms)));
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2gene");
      } else{
        dataSet$gene2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(gene.nms)),Mapped=length(unique(mir.nms)));
        seedsu <<- c(seedsu, gene.nms);
        mirtableu <<- c(mirtableu, "gene2mir");
      }
    }
    return(1);
}

.searchMultiNet_molecule2mir <- function(input.type){
    orgType <- dataSet$org;
    if(orgType %in% c("bta", "dme", "gga","sma", "cel", "ssc")){
      curent.msg <<- "This organism is not supported for molecule network research."
      print(current.msg);
      return(0);
    }

    if (input.type == "mir2molecule"){
      idType <- dataSet$mirnaType;
      mir.mat <- dataSet$mir.orig;
      mir.vec <- rownames(mir.mat);
    }else{
      idType <- "molecule";
      mir.mat <- dataSet$data[[idType]];
      mir.vec <- rownames(mir.mat);
    }

    mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2molecule", sep=""), mir.vec, orgType, idType);

    hit.num <- nrow(mir.dic)
    if (hit.num == 0 && dataSet$tissue == "na") {
      current.msg <<- "No hits found in the miRNA-small compound database. Please check your input. ";
      print(current.msg);
      return(0);
    } else if (hit.num == 0 && dataSet$tissue != "na") {
      current.msg <<- "No hits found in the miRNA-small compound database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
      print(current.msg);
      return(2);
    } else {
      res <- mir.dic[, c("mir_id","mir_acc","molecule", "pubchem_id", "method", "pmid", "tissue")];
      rownames(res) <- mir.dic$mirnet;
      current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-gene targets were identified!");

      # update the data
      gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
      mir.mat <- mir.mat[gd.inx,,drop=F];
      mir.mappedu <<- rbind(mir.mappedu, mir.mat);

      colnames(res) <- c("ID","Accession","Molecule", "Pubchem_ID", "Experiment", "Literature", "Tissue");
      display.res <- res;
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Molecule"]);    # IDs
      if(nrow(res)!=0){
        row.names(edge.res) <- 1:nrow(res);
      }

      node.ids <- c(ID1=res[,"Accession"], ID2=res[,"Molecule"]);
      node.nms <- c(Name1=res[,"ID"], Name2=res[,"Molecule"]);

      edgeu.res <<- rbind(edgeu.res, edge.res); #edgeu.res is an empty dataframe defined in QueryNet
      edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ];
      nodeu.ids <<- c(nodeu.ids, node.ids);
      edgeNumU <<- c(edgeNumU, nrow(edge.res))
      nodeu.nms <<- c(nodeu.nms, node.nms);

      mol.nms <<- res[,"Molecule"];
      mir.nms <<- res[, "ID"];
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Molecule"], TargetID=res[,"Pubchem_ID"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"]);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "molecule");
      if(input.type == "mir2molecule"){
        dataSet$mir2mol <<- res;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(mol.nms)));
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2mol");
      } else{
        dataSet$mol2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(mol.nms)),Mapped=length(unique(mir.nms)));
        seedsu <<- c(seedsu, mol.nms);
        mirtableu <<- c(mirtableu, "mol2mir");
      }
    }
    return(1);
}

.searchMultiNet_lncrna2mir <- function(input.type){
    orgType <- dataSet$org;
    if(orgType != "hsa" ){
      curent.msg <<- "Only human supports lncRNA network."
      print(current.msg);
      return(0);
    }

    if (input.type == "mir2lncrna"){
      idType <- dataSet$mirnaType;
      mir.mat <- dataSet$mir.orig;
      mir.vec <- rownames(mir.mat);
    }else{
      idType <- dataSet$id.types[["lncrna"]];
      mir.mat <- dataSet$data[["lncrna"]];
      mir.vec <- rownames(mir.mat);
    }

    mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2lncRNA", sep=""), mir.vec, orgType, idType);

    hit.num <- nrow(mir.dic)
    if (hit.num == 0 && dataSet$tissue == "na") {
      current.msg <<- "No hits found in the miRNA-lncRNA database. Please check your input. ";
      print(current.msg);
      return(0);
    } else if (hit.num == 0 && dataSet$tissue != "na") {
      current.msg <<- "No hits found in the miRNA-lncRNA database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
      print(current.msg);
      return(2);
    } else {
      res <- mir.dic[ , c("mir_id","mir_acc","symbol","entrez", "tissue")];
      rownames(res) <- mir.dic$mirnet;
      current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-lncRNA targets were identified!");

      # update the data
      gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
      mir.mat <- mir.mat[gd.inx,,drop=F];
      mir.mappedu <<- rbind(mir.mappedu, mir.mat);

      colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Tissue");
      res$Experiment <- rep("CLIP-Seq", nrow(res));
      res$Literature <- rep("24297251", nrow(res));
      res <- res[, c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue")];
      display.res <- res;
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Entrez"]);    # IDs
      if(nrow(res)!=0){
        row.names(edge.res) <- 1:nrow(res);
      }

      node.ids <- c(ID1=res[,"Accession"], ID2=res[,"Entrez"]);
      node.nms <- c(Name1=res[,"ID"], Name2=res[,"Gene"]);

      edgeu.res <<- rbind(edgeu.res, edge.res); #edgeu.res is an empty dataframe defined in QueryNet
      edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ];
      nodeu.ids <<- c(nodeu.ids, node.ids);
      edgeNumU <<- c(edgeNumU, nrow(edge.res))
      nodeu.nms <<- c(nodeu.nms, node.nms);

      lnc.nms <<- res[,"Gene"];
      mir.nms <<- res[, "ID"];
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"]);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "lncrna");
      if(input.type == "mir2lncrna"){
        dataSet$mir2lnc <<- res;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(lnc.nms)));
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2lnc");
      } else{
        dataSet$lnc2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(lnc.nms)),Mapped=length(unique(mir.nms)));
        seedsu <<- c(seedsu, lnc.nms);
        mirtableu <<- c(mirtableu, "lnc2mir");
      }
    }
    return(1);
  }


.searchMultiNet_cir2mir <- function(input.type){

    orgType <- dataSet$org;
    if(orgType != "hsa" ){
      curent.msg <<- "Only human supports circRNA network."
      print(current.msg);
      return(0);
    }

    if (input.type == "mir2circrna"){
      idType <- dataSet$mirnaType;
      mir.mat <- dataSet$mir.orig;
      mir.vec <- rownames(mir.mat);
    }else{
      idType <- dataSet$id.types[["circrna"]];
      mir.mat <- dataSet$data[["circrna"]];
      mir.vec <- rownames(mir.mat);
    }

    mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2circRNA", sep=""), mir.vec, orgType, idType);

    hit.num <- nrow(mir.dic)
    if (hit.num == 0 && dataSet$tissue == "na") {
      current.msg <<- "No hits found in the miRNA-circRNA database. Please check your input. ";
      print(current.msg);
      return(0);
    } else if (hit.num == 0 && dataSet$tissue != "na") {
      current.msg <<- "No hits found in the miRNA-circRNA database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
      print(current.msg);
      return(2);
    } else {
      res <- mir.dic[ , c("mir_id","mir_acc","symbol","entrez", "tissue")];
      rownames(res) <- mir.dic$mirnet;
      current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-circRNA targets were identified!");

      # update the data
      gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
      mir.mat <- mir.mat[gd.inx,,drop=F];
      mir.mappedu <<- rbind(mir.mappedu, mir.mat);

      colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Tissue");
      res$Experiment <- rep("CLIP-Seq", nrow(res));
      res$Literature <- rep("24297251", nrow(res));
      res <- res[, c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue")];
      display.res <- res;
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Entrez"]);    # IDs
      if(nrow(res)!=0){
        row.names(edge.res) <- 1:nrow(res);
      }

      node.ids <- c(ID1=res[,"Accession"], ID2=res[,"Entrez"]);
      node.nms <- c(Name1=res[,"ID"], Name2=res[,"Gene"]);

      edgeu.res <<- rbind(edgeu.res, edge.res); #edgeu.res is an empty dataframe defined in QueryNet
      edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ];
      nodeu.ids <<- c(nodeu.ids, node.ids);
      edgeNumU <<- c(edgeNumU, nrow(edge.res))
      nodeu.nms <<- c(nodeu.nms, node.nms);

      circ.nms <<- res[,"Gene"];
      mir.nms <<- res[, "ID"];
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"]);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "circrna");
      if(input.type == "mir2circrna"){
        dataSet$mir2circ <<- res;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(circ.nms)));
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2circ");
      } else{
        dataSet$circ2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(circ.nms)),Mapped=length(unique(mir.nms)));
        seedsu <<- c(seedsu, circ.nms);
        mirtableu <<- c(mirtableu, "circ2mir");
      }
    }
    return(1);
}

.searchMultiNet_pseudogene2mir <- function(input.type){

    orgType <- dataSet$org;
    if(orgType != "hsa" ){
      curent.msg <<- "Only human supports pseudogene network."
      print(current.msg);
      return(0);
    }

    if (input.type == "mir2pseudogene"){
      idType <- dataSet$mirnaType;
      mir.mat <- dataSet$mir.orig;
      mir.vec <- rownames(mir.mat);
    }else{
      idType <- dataSet$id.types[["pseudogene"]];
      mir.mat <- dataSet$data[["pseudogene"]];
      mir.vec <- rownames(mir.mat);
    }

    mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2pseudogene", sep=""), mir.vec, orgType, idType);

    hit.num <- nrow(mir.dic)
    if (hit.num == 0 && dataSet$tissue == "na") {
      current.msg <<- "No hits found in the miRNA-pseudogene database. Please check your input. ";
      print(current.msg);
      return(0);
    } else if (hit.num == 0 && dataSet$tissue != "na") {
      current.msg <<- "No hits found in the miRNA-pseudogene database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
      print(current.msg);
      return(2);
    } else {
      res <- mir.dic[ , c("mir_id","mir_acc","symbol","entrez", "tissue")];
      rownames(res) <- mir.dic$mirnet;
      current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-pseudogene targets were identified!");

      # update the data
      gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
      mir.mat <- mir.mat[gd.inx,,drop=F];
      mir.mappedu <<- rbind(mir.mappedu, mir.mat);

      colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Tissue");
      res$Experiment <- rep("CLIP-Seq", nrow(res));
      res$Literature <- rep("24297251", nrow(res));
      res <- res[, c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue")];
      display.res <- res;
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Entrez"]);    # IDs
      if(nrow(res)!=0){
        row.names(edge.res) <- 1:nrow(res);
      }

      node.ids <- c(ID1=res[,"Accession"], ID2=res[,"Entrez"]);
      node.nms <- c(Name1=res[,"ID"], Name2=res[,"Gene"]);

      edgeu.res <<- rbind(edgeu.res, edge.res); #edgeu.res is an empty dataframe defined in QueryNet
      edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ];
      nodeu.ids <<- c(nodeu.ids, node.ids);
      edgeNumU <<- c(edgeNumU, nrow(edge.res))
      nodeu.nms <<- c(nodeu.nms, node.nms);

      pseudo.nms <<- res[,"Gene"];
      mir.nms <<- res[, "ID"];
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"]);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "pseudogene");
      if(input.type == "mir2pseudogene"){
        dataSet$mir2pseudo <<- res;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(pseudo.nms)));
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2pseudo");
      } else{
        dataSet$pseudo2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(pseudo.nms)),Mapped=length(unique(mir.nms)));
        seedsu <<- c(seedsu, pseudo.nms);
        mirtableu <<- c(mirtableu, "pseudo2mir");
      }
    }
    return(1);
}

.searchMultiNet_sncrna2mir <- function(input.type){
    orgType <- dataSet$org;
    if(orgType != "hsa" ){
      curent.msg <<- "Only human supports sncRNA network."
      print(current.msg);
      return(0);
    }

    if (input.type == "mir2sncrna"){
      idType <- dataSet$mirnaType;
      mir.mat <- dataSet$mir.orig;
      mir.vec <- rownames(mir.mat);
    }else{
      idType <- dataSet$id.types[["sncrna"]];
      mir.mat <- dataSet$data[["sncrna"]];
      mir.vec <- rownames(mir.mat);
    }

    mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2sncRNA", sep=""), mir.vec, orgType, idType);

    hit.num <- nrow(mir.dic)
    if (hit.num == 0 && dataSet$tissue == "na") {
      current.msg <<- "No hits found in the miRNA-sncRNA database. Please check your input. ";
      print(current.msg);
      return(0);
    } else if (hit.num == 0 && dataSet$tissue != "na") {
      current.msg <<- "No hits found in the miRNA-sncRNA database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
      print(current.msg);
      return(2);
    } else {
      res <- mir.dic[ , c("mir_id","mir_acc","symbol","entrez", "tissue")];
      rownames(res) <- mir.dic$mirnet;
      current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-sncRNA targets were identified!");

      # update the data
      gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
      mir.mat <- mir.mat[gd.inx,,drop=F];
      mir.mappedu <<- rbind(mir.mappedu, mir.mat);

      colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Tissue");
      res$Experiment <- rep("CLIP-Seq", nrow(res));
      res$Literature <- rep("24297251", nrow(res));
      res <- res[, c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue")];
      display.res <- res;
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Entrez"]);    # IDs
      if(nrow(res)!=0){
        row.names(edge.res) <- 1:nrow(res);
      }

      node.ids <- c(ID1=res[,"Accession"], ID2=res[,"Entrez"]);
      node.nms <- c(Name1=res[,"ID"], Name2=res[,"Gene"]);

      edgeu.res <<- rbind(edgeu.res, edge.res); #edgeu.res is an empty dataframe defined in QueryNet
      edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ];
      nodeu.ids <<- c(nodeu.ids, node.ids);
      edgeNumU <<- c(edgeNumU, nrow(edge.res))
      nodeu.nms <<- c(nodeu.nms, node.nms);

      snc.nms <<- res[,"Gene"];
      mir.nms <<- res[, "ID"];
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"]);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "sncrna");
      if(input.type == "mir2sncrna"){
        dataSet$mir2snc <<- res;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(snc.nms)));
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2snc");
      } else{
        dataSet$snc2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(snc.nms)),Mapped=length(unique(mir.nms)));
        seedsu <<- c(seedsu, snc.nms);
        mirtableu <<- c(mirtableu, "snc2mir");
      }
    }
    return(1);
}

.searchMultiNet_tf2mir <- function(input.type){
    orgType <- dataSet$org;
    if(orgType %in% c("bta", "ssc","gga","dme", "sma") ){
      curent.msg <<- "This organism is not supported for transcription factors network research."
      print(current.msg);
      return(0);
    }

    if (input.type == "mir2tf"){
      idType <- dataSet$mirnaType;
      mir.mat <- dataSet$mir.orig;
      mir.vec <- rownames(mir.mat);
    }else if(input.type == "gene2tf"){
      idType <- dataSet$id.types[["gene"]];
      mir.mat <- dataSet$data[["gene"]];
      mir.vec <- rownames(mir.mat);
    }else if(input.type == "tf2gene"){
      idType <- dataSet$id.types[["tf"]];
      if(idType=="symbol"){
        idType <- "tfname"; # when tf as input, use "tfname" column
      }else{
        idType <- "tfid";
      }
      mir.mat <- dataSet$data[["tf"]];
      mir.vec <- rownames(mir.mat);
    }else{
      idType <- dataSet$id.types[["tf"]];
      mir.mat <- dataSet$data[["tf"]];
      mir.vec <- rownames(mir.mat);
    }
    if(!input.type %in% c("tf2gene", "gene2tf")){
      mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2tf", sep=""), mir.vec, orgType, idType);
    }else{
      if(input.type == "gene2tf"){
        targetType <- dataSet$target.types[["gene"]];
      }else{
        targetType <- dataSet$target.types[["tf"]];
      }
      dataSet$tfTargetType <<- c(dataSet$tfTargetType, targetType)
      table.nm <- paste(dataSet$org,targetType, sep="_");
      mir.dic <- QueryTFSQLite(table.nm, mir.vec, idType);
    }


    hit.num <- nrow(mir.dic)
    if (hit.num == 0 && dataSet$tissue == "na") {
      current.msg <<- "No hits found in the miRNA-TF database. Please check your input. ";
      print(current.msg);
      return(0);
    } else if (hit.num == 0 && dataSet$tissue != "na") {
      current.msg <<- "No hits found in the miRNA-TF database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
      print(current.msg);
      return(2);
    } else {
      if(!input.type %in% c("tf2gene", "gene2tf")){
        res <- mir.dic[ , c("mir_id","mir_acc","symbol","entrez", "action_type", "pmid", "tissue")];
        rownames(res) <- mir.dic$mirnet;
        current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-TF targets were identified!");

        # update the data
        gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
        mir.mat <- mir.mat[gd.inx,,drop=F];
        mir.mappedu <<- rbind(mir.mappedu, mir.mat);

        colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Direction", "Literature", "Tissue");
        res$Experiment <- rep("ChIP-seq", nrow(res));
        res <- res[, c("ID", "Accession", "Gene", "Entrez", "Direction", "Experiment", "Literature", "Tissue")];
        direction <- res$Direction;
        display.res <- res;
        edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Entrez"]);    # IDs
        if(nrow(res)!=0){
          row.names(edge.res) <- 1:nrow(res);
        }
        node.ids <- c(ID1=res[,"Accession"], ID2=res[,"Entrez"]);
        node.nms <- c(Name1=res[,"ID"], Name2=res[,"Gene"]);
      }else{
        res <- mir.dic[ ,c("tfname","tfid","symbol","entrez", "Direction", "Literature")]; # "tfname" is TF, "symbol" is target gene
        direction <- res$Direction
        literature <- res$Literature
        rownames(res) <- mir.dic$mirnet;
        current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-TF targets were identified!");

        # update the data
        gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
        mir.mat <- mir.mat[gd.inx,,drop=F];
        mir.mappedu <<- rbind(mir.mappedu, mir.mat);

        colnames(res)[c(1:4)] <- c("ID", "Accession", "Gene", "Entrez");
        res$Experiment <- rep("ChIP-seq", nrow(res));
        res <- res[, c("ID", "Accession", "Gene", "Entrez")];
        display.res <- res;
        edge.res <- data.frame(Source=res[,3],Target=res[,1]);    # IDs
        if(nrow(res)!=0){
          row.names(edge.res) <- 1:nrow(res);
        }
        node.ids <- c(ID1=res[,"ID"], ID2=res[,"Gene"]);
        node.nms <- c(Name1=res[,"Accession"], Name2=res[,"Entrez"]);
      }
      edgeu.res <<- rbind(edgeu.res, edge.res); #edgeu.res is an empty dataframe defined in QueryNet
      edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ];
      nodeu.ids <<- c(nodeu.ids, node.ids);
      edgeNumU <<- c(edgeNumU, nrow(edge.res))
      nodeu.nms <<- c(nodeu.nms, node.nms);
      if(input.type %in% c("tf2mir", "mir2tf")){
        tf.nms <<- unique(c(tf.nms, as.vector(res[,"Gene"])));
        mir.nms <<- res[,"ID"];
        direction <- gsub("Activation[(]feedback[)]", "+", direction);
        direction <- gsub("Repression[(]feedback[)]", "-", direction);
        direction <- gsub("Activation", "+", direction);
        direction <- gsub("Repression", "-", direction);
        res$Direction <- direction;
        res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Direction"], Literature=res[,"Literature"], Tissue=res[,"Tissue"]);
        l = nrow(mir.resu);
        dataSet$directionInx =c(dataSet$directionInx, paste0(res[,1],res[,3]));
        dataSet$regDirection =c(dataSet$regDirection, res$Experiment);
        dataSet<<-dataSet;
      }else{
        gene.nms <<- unique(c(gene.nms, as.vector(res[,"Gene"])));
        tf.nms <<- unique(c(tf.nms, as.vector(res[,"ID"])));
        mir.nms <<- "";
          direction <- gsub("Activation", "+", direction)
          direction <- gsub("Repression", "-", direction)
        if(targetType == "encode"){
          res$Literature <- gsub("ENCODE", "29126249", literature);
          res$Experiment <- rep("ChIP-seq", nrow(res));

          res$Direction <- direction
        }else if(targetType == "jaspar"){
          res$Literature <- gsub("JASPAR", "31701148", literature);
          res$Experiment <- rep("Computational", nrow(res));
          res$Direction <- direction
        }else if(targetType == "chea"){
          res$Literature <- gsub("ChEA", "20709693", literature);
          res$Experiment <- rep("ChIP-X", nrow(res));
          res$Direction <- direction
        }else if(targetType == "trrust"){
          res$Experiment <- rep("TRRUST", nrow(res));
          res$Literature <- literature
          res$Direction <- direction
        }else if(targetType == "regnetwork"){
          res$Experiment <- rep("RegNetwork", nrow(res));
          res$Literature <- literature
          res$Direction <- direction
        }
        if(targetType %in% c("trrust", "regnetwork")){
             res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Direction"], Literature=res[,"Literature"],Tissue=rep("Not Applicable", nrow(res)));
             l = nrow(mir.resu)
             dataSet$directionInx =c(dataSet$directionInx, paste0(res[,1],res[,3]))
             dataSet$regDirection =c(dataSet$regDirection, res$Experiment)
             dataSet<<-dataSet
            }else{
             res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"],Tissue=rep("Not Applicable", nrow(res)));
            }
    }
      mir.resu <<- rbind(mir.resu, res);

      mirtargetu <<- c(mirtargetu, "tf");
      if(input.type == "mir2tf"){
        dataSet$mir2tf <<- res;         # save this for network builder and table view
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2tf");
      } else if(input.type == "tf2gene"){
        dataSet$tf2gene <<- res;
        seedsu <<- c(seedsu, tf.nms);
        mirtableu <<- c(mirtableu, "tf2gene");
      } else if(input.type == "gene2tf"){
        dataSet$gene2tf <<- res;
        seedsu <<- c(seedsu, tf.nms);
        mirtableu <<- c(mirtableu, "gene2tf");
      } else{
        dataSet$tf2mir <<- res;
        seedsu <<- c(seedsu, tf.nms);
        mirtableu <<- c(mirtableu, "tf2mir");
      }
    }
    return(1);
}

.searchMultiNet_tf2gene2mir <- function(input.type){
  orgType <- dataSet$org;
  if(orgType %in% c("bta", "ssc","gga","dme", "sma") ){
    curent.msg <<- "This organism is not supported for transcription factors network research."
    print(current.msg);
    return(0);
  }
    idType <- dataSet$id.types[["tf"]];
    if(idType=="symbol"){
      idType <- "tfname";
    }else{
      idType <- "tfid";
    }
    mir.mat <- dataSet$data[["tf"]];
    mir.vec <- rownames(mir.mat);
    targetType <- dataSet$target.types[["tf"]];
    dataSet$tfTargetType <<- c(dataSet$tfTargetType, targetType)
    table.nm <- paste(dataSet$org,targetType, sep="_");
    mir.dic <- QueryTFSQLite(table.nm, mir.vec, idType);
  hit.num <- nrow(mir.dic)
  if (hit.num == 0 && dataSet$tissue == "na") {
    current.msg <<- "No hits found in the miRNA-TF database. Please check your input. ";
    print(current.msg);
    return(0);
  } else if (hit.num == 0 && dataSet$tissue != "na") {
    current.msg <<- "No hits found in the miRNA-TF database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
    print(current.msg);
    return(2);
  } else {
    # map gene2mir
    idVec <- as.vector(unique(mir.dic[, c("entrez")]));
    mir.dic2 <- Query.miRNetDB(paste(sqlite.path, "mir2gene", sep=""), idVec, dataSet$org, "entrez");

      # tf2gene
      res <- mir.dic[ ,c("tfname","tfid","symbol","entrez", "Direction", "Literature")]; # "tfname" is TF, "symbol" is target gene
      direction <- res$Direction
      literature <- res$Literature
      rownames(res) <- mir.dic$mirnet;
      # gene2mir
      res2 <- mir.dic2[ , c("mir_id","mir_acc","symbol","entrez", "experiment", "pmid", "tissue")];
      rownames(res2) <- mir.dic2$mirnet;
      current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-TF targets were identified!");

      # update the data
      gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
      mir.mat <- mir.mat[gd.inx,,drop=F];
      mir.mappedu <<- rbind(mir.mappedu, mir.mat);
      #tf2gene
      colnames(res)[c(1:4)] <- c("ID", "Accession", "Gene", "Entrez");
      res$Experiment <- rep("ChIP-seq", nrow(res));
      res <- res[, c("ID", "Accession", "Gene", "Entrez")];
      display.res <- res;
      edge.res <- data.frame(Source=res[,3],Target=res[,1]);    # IDs
      if(nrow(res)!=0){
        row.names(edge.res) <- 1:nrow(res);
      }
      node.ids <- c(ID1=res[,"ID"], ID2=res[,"Gene"]);
      node.nms <- c(Name1=res[,"Accession"], Name2=res[,"Entrez"]);
      #gene2mir
      # update col names
      colnames(res2) <- c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue");

      edge.res2 <- data.frame(Source=res2[,"Accession"],Target=res2[,"Entrez"]);    # IDs
      if(nrow(res2)!=0){
        row.names(edge.res2) <- 1:nrow(res2);
      }
      node.ids2 <- c(ID1=res2[,"Accession"], ID2=res2[,"Entrez"]);
      node.nms2 <- c(Name1=res2[,"ID"], Name2=res2[,"Gene"]);
      # rbind tf2gene and gene2mir network
      edge.res3 <- rbind(edge.res,edge.res2);
      node.ids3 <- c(node.ids,node.ids2);
      node.nms3 <- c(node.nms,node.nms2);

    edgeu.res <<- rbind(edgeu.res, edge.res3); #edgeu.res is an empty dataframe defined in QueryNet
    edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ];
    nodeu.ids <<- c(nodeu.ids, node.ids3);
    edgeNumU <<- c(edgeNumU, nrow(edge.res3))
    nodeu.nms <<- c(nodeu.nms, node.nms3);

      gene.nms <<- unique(c(gene.nms, as.vector(res[,"Gene"])));
      tf.nms <<- unique(c(tf.nms, as.vector(res[,"ID"])));
      mir.nms <<- res2[,"ID"];
      direction <- gsub("Activation", "+", direction)
      direction <- gsub("Repression", "-", direction)
      if(targetType == "encode"){
        res$Literature <- gsub("ENCODE", "29126249", literature);
        res$Experiment <- rep("ChIP-seq", nrow(res));

        res$Direction <- direction
      }else if(targetType == "jaspar"){
        res$Literature <- gsub("JASPAR", "31701148", literature);
        res$Experiment <- rep("Computational", nrow(res));
        res$Direction <- direction
      }else if(targetType == "chea"){
        res$Literature <- gsub("ChEA", "20709693", literature);
        res$Experiment <- rep("ChIP-X", nrow(res));
        res$Direction <- direction
      }else if(targetType == "trrust"){
        res$Experiment <- rep("TRRUST", nrow(res));
        res$Literature <- literature
        res$Direction <- direction
      }else if(targetType == "regnetwork"){
        res$Experiment <- rep("RegNetwork", nrow(res));
        res$Literature <- literature
        res$Direction <- direction
      }
      if(targetType %in% c("trrust", "regnetwork")){
        res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Direction"], Literature=res[,"Literature"],Tissue=rep("Not Applicable", nrow(res)));
        l = nrow(mir.resu)
        dataSet$directionInx =c(dataSet$directionInx, paste0(res[,1],res[,3]))
        dataSet$regDirection =c(dataSet$regDirection, res$Experiment)
        dataSet<<-dataSet
      }else{
        res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"],Tissue=rep("Not Applicable", nrow(res)));
      }
      res2 <- data.frame(ID=res2[,"ID"], Accession=res2[,"Accession"], Target=res2[,"Gene"], TargetID=res2[,"Entrez"], Experiment=res2[,"Experiment"], Literature=res2[,"Literature"], Tissue=res2[,"Tissue"]);
        mir.resu <<- rbind(mir.resu, res, res2);
        mirtargetu <<- c(mirtargetu, "tf");
      dataSet$tf2gene <<- res;
      dataSet$gene2mir <<- res2;
      seedsu <<- c(seedsu, tf.nms);
      mirtableu <<- c(mirtableu, "tf2gene");
      mirtableu <<- c(mirtableu, "gene2mir");
  }
  return(1);
}

.searchMultiNet_gene2tf2mir <- function(input.type){
  orgType <- dataSet$org;
  if(orgType %notin% c("hsa", "mmu") ){
    curent.msg <<- "This organism is not supported for transcription factors network research."
    print(current.msg);
    return(0);
  }
  idType <- dataSet$id.types[["gene"]];
  mir.mat <- dataSet$data[["gene"]];
  mir.vec <- rownames(mir.mat);
  targetType <- dataSet$target.types[["gene"]];
  dataSet$tfTargetType <<- c(dataSet$tfTargetType, targetType)
  table.nm <- paste(dataSet$org,targetType, sep="_");
  mir.dic <- QueryTFSQLite(table.nm, mir.vec, idType);
  hit.num <- nrow(mir.dic)
  if (hit.num == 0 && dataSet$tissue == "na") {
    current.msg <<- "No hits found in the miRNA-TF database. Please check your input. ";
    print(current.msg);
    return(0);
  } else if (hit.num == 0 && dataSet$tissue != "na") {
    current.msg <<- "No hits found in the miRNA-TF database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
    print(current.msg);
    return(2);
  } else {
    # map tf2mir
    idVec <- as.vector(unique(mir.dic[, c("tfid")]));
    mir.dic2 <- Query.miRNetDB(paste(sqlite.path, "mir2tf", sep=""), idVec, dataSet$org, "entrez");

    # gene2tf
    res <- mir.dic[ ,c("tfname","tfid","symbol","entrez", "Direction", "Literature")]; # "tfname" is TF, "symbol" is target gene
    direction <- res$Direction
    literature <- res$Literature
    rownames(res) <- mir.dic$mirnet;
    # tf2mir
    res2 <- mir.dic2[ , c("mir_id","mir_acc","symbol","entrez", "action_type", "pmid", "tissue")];
    rownames(res2) <- mir.dic2$mirnet;
    current.msg <<- paste("A total of unqiue", hit.num, "pairs of gene-TF targets were identified!");

    # update the data
    gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
    mir.mat <- mir.mat[gd.inx,,drop=F];
    mir.mappedu <<- rbind(mir.mappedu, mir.mat);
    #gene2tf
    colnames(res)[c(1:4)] <- c("ID", "Accession", "Gene", "Entrez");
    res$Experiment <- rep("ChIP-seq", nrow(res));
    res <- res[, c("ID", "Accession", "Gene", "Entrez")];
    display.res <- res;
    edge.res <- data.frame(Source=res[,3],Target=res[,1]);    # IDs
    if(nrow(res)!=0){
      row.names(edge.res) <- 1:nrow(res);
    }
    node.ids <- c(ID1=res[,"ID"], ID2=res[,"Gene"]);
    node.nms <- c(Name1=res[,"Accession"], Name2=res[,"Entrez"]);
    #tf2mir
    # update col names
    colnames(res2) <- c("ID", "Accession", "Gene", "Entrez", "Direction", "Literature", "Tissue");
    res2$Experiment <- rep("ChIP-seq", nrow(res2));
    res2 <- res2[, c("ID", "Accession", "Gene", "Entrez", "Direction", "Experiment", "Literature", "Tissue")];
    direction2 <- res2$Direction;
    edge.res2 <- data.frame(Source=res2[,"Accession"],Target=res2[,"Entrez"]);    # IDs
    if(nrow(res2)!=0){
      row.names(edge.res2) <- 1:nrow(res2);
    }
    node.ids2 <- c(ID1=res2[,"Accession"], ID2=res2[,"Entrez"]);
    node.nms2 <- c(Name1=res2[,"ID"], Name2=res2[,"Gene"]);
    # rbind gene2tf and tf2mir network
    edge.res3 <- rbind(edge.res,edge.res2);
    node.ids3 <- c(node.ids,node.ids2);
    node.nms3 <- c(node.nms,node.nms2);

    edgeu.res <<- rbind(edgeu.res, edge.res3); #edgeu.res is an empty dataframe defined in QueryNet
    edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ];
    nodeu.ids <<- c(nodeu.ids, node.ids3);
    edgeNumU <<- c(edgeNumU, nrow(edge.res3))
    nodeu.nms <<- c(nodeu.nms, node.nms3);

    gene.nms <<- unique(c(gene.nms, as.vector(res[,"Gene"])));
    tf.nms <<- unique(c(tf.nms, as.vector(res[,"ID"])));
    mir.nms <<- res2[,"ID"];
    direction <- gsub("Activation", "+", direction)
    direction <- gsub("Repression", "-", direction)
    if(targetType == "encode"){
      res$Literature <- gsub("ENCODE", "29126249", literature);
      res$Experiment <- rep("ChIP-seq", nrow(res));
      res$Direction <- direction
    }else if(targetType == "jaspar"){
      res$Literature <- gsub("JASPAR", "31701148", literature);
      res$Experiment <- rep("Computational", nrow(res));
      res$Direction <- direction
    }else if(targetType == "chea"){
      res$Literature <- gsub("ChEA", "20709693", literature);
      res$Experiment <- rep("ChIP-X", nrow(res));
      res$Direction <- direction
    }else if(targetType == "trrust"){
      res$Experiment <- rep("TRRUST", nrow(res));
      res$Literature <- literature
      res$Direction <- direction
    }else if(targetType == "regnetwork"){
      res$Experiment <- rep("RegNetwork", nrow(res));
      res$Literature <- literature
      res$Direction <- direction
    }
    if(targetType %in% c("trrust", "regnetwork")){
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Direction"], Literature=res[,"Literature"],Tissue=rep("Not Applicable", nrow(res)));
      l = nrow(mir.resu)
      dataSet$directionInx =c(dataSet$directionInx, paste0(res[,1],res[,3]))
      dataSet$regDirection =c(dataSet$regDirection, res$Experiment)
      dataSet<<-dataSet
    }else{
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"],Tissue=rep("Not Applicable", nrow(res)));
    }
    direction2 <- gsub("Activation[(]feedback[)]", "+", direction2);
    direction2 <- gsub("Repression[(]feedback[)]", "-", direction2);
    direction2 <- gsub("Activation", "+", direction2);
    direction2 <- gsub("Repression", "-", direction2);
    res2$Direction <- direction2;
    res2 <- data.frame(ID=res2[,"ID"], Accession=res2[,"Accession"], Target=res2[,"Gene"], TargetID=res2[,"Entrez"], Experiment=res2[,"Direction"], Literature=res2[,"Literature"], Tissue=res2[,"Tissue"]);
    mir.resu <<- rbind(mir.resu, res, res2);
    mirtargetu <<- c(mirtargetu, "tf");
    dataSet$gene2tf <<- res;
    dataSet$tf2mir <<- res2;
    seedsu <<- c(seedsu, gene.nms);
    mirtableu <<- c(mirtableu, "gene2tf");
    mirtableu <<- c(mirtableu, "tf2mir");
  }
   return(1);
}

.searchMultiNet_mir2tf2gene <- function(input.type){
  orgType <- dataSet$org;
  if(orgType %notin% c("hsa", "mmu") ){
    curent.msg <<- "This organism is not supported for transcription factors network research."
    print(current.msg);
    return(0);
  }
  idType <- dataSet$mirnaType;
  mir.mat <- dataSet$mir.orig;
  mir.vec <- rownames(mir.mat);
  mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2tf", sep=""), mir.vec, orgType, idType);
  hit.num <- nrow(mir.dic)
  if (hit.num == 0 && dataSet$tissue == "na") {
    current.msg <<- "No hits found in the miRNA-TF database. Please check your input. ";
    print(current.msg);
    return(0);
  } else if (hit.num == 0 && dataSet$tissue != "na") {
    current.msg <<- "No hits found in the miRNA-TF database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
    print(current.msg);
    return(2);
  } else {
    # map tf2gene
    idVec <- as.vector(unique(mir.dic[, c("entrez")]));
    table.nm <- paste(dataSet$org,"trrust", sep="_");
    mir.dic2 <- QueryTFSQLite(table.nm, idVec, "tfid");

    # mir2tf
    res <- mir.dic[ , c("mir_id","mir_acc","symbol","entrez", "action_type", "pmid", "tissue")];
    rownames(res) <- mir.dic$mirnet;
    current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-TF interactions were identified!");
    # tf2gene
    res2 <- mir.dic2[ ,c("tfname","tfid","symbol","entrez", "Direction", "Literature")]; # "tfname" is TF, "symbol" is target gene
    direction2 <- res2$Direction
    literature2 <- res2$Literature
    rownames(res2) <- mir.dic2$mirnet;

    # update the data
    gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
    mir.mat <- mir.mat[gd.inx,,drop=F];
    mir.mappedu <<- rbind(mir.mappedu, mir.mat);
    #mir2tf
    colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Direction", "Literature", "Tissue");
    res$Experiment <- rep("ChIP-seq", nrow(res));
    res <- res[, c("ID", "Accession", "Gene", "Entrez", "Direction", "Experiment", "Literature", "Tissue")];
    direction <- res$Direction;
    display.res <- res;
    edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Entrez"]);    # IDs
    if(nrow(res)!=0){
      row.names(edge.res) <- 1:nrow(res);
    }
    node.ids <- c(ID1=res[,"Accession"], ID2=res[,"Entrez"]);
    node.nms <- c(Name1=res[,"ID"], Name2=res[,"Gene"]);
    #tf2gene
    # update col names
    colnames(res2)[c(1:4)] <- c("ID", "Accession", "Gene", "Entrez");
    res2$Experiment <- rep("ChIP-seq", nrow(res2));
    res2 <- res2[, c("ID", "Accession", "Gene", "Entrez")];
    edge.res2 <- data.frame(Source=res2[,3],Target=res2[,1]);    # IDs
    if(nrow(res2)!=0){
      row.names(edge.res2) <- 1:nrow(res2);
    }
    node.ids2 <- c(ID1=res2[,"ID"], ID2=res2[,"Gene"]);
    node.nms2 <- c(Name1=res2[,"Accession"], Name2=res2[,"Entrez"]);
    # rbind mir2tf and tf2gene network
    edge.res3 <- rbind(edge.res,edge.res2);
    node.ids3 <- c(node.ids,node.ids2);
    node.nms3 <- c(node.nms,node.nms2);

    edgeu.res <<- rbind(edgeu.res, edge.res3); #edgeu.res is an empty dataframe defined in QueryNet
    edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ];
    nodeu.ids <<- c(nodeu.ids, node.ids3);
    edgeNumU <<- c(edgeNumU, nrow(edge.res3))
    nodeu.nms <<- c(nodeu.nms, node.nms3);

    gene.nms <<- unique(c(gene.nms, as.vector(res2[,"Gene"])));
    tf.nms <<- unique(c(tf.nms, as.vector(res2[,"ID"])));
    mir.nms <<- res[,"ID"];
    # mir2tf
    direction <- gsub("Activation[(]feedback[)]", "+", direction);
    direction <- gsub("Repression[(]feedback[)]", "-", direction);
    direction <- gsub("Activation", "+", direction);
    direction <- gsub("Repression", "-", direction);
    res$Direction <- direction;
    res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Direction"], Literature=res[,"Literature"], Tissue=res[,"Tissue"]);
    dataSet$directionInx =c(dataSet$directionInx, paste0(res[,1],res[,3]));
    dataSet$regDirection =c(dataSet$regDirection, res$Experiment);
    dataSet<<-dataSet;
    # tf2gene
    res2$Experiment <- rep("TRRUST", nrow(res2));
      res2$Literature <- literature2
      res2$Direction <- direction2
      res2 <- data.frame(ID=res2[,"ID"], Accession=res2[,"Accession"], Target=res2[,"Gene"], TargetID=res2[,"Entrez"], Experiment=res2[,"Direction"], Literature=res2[,"Literature"],Tissue=rep("Not Applicable", nrow(res2)));
      dataSet$directionInx =c(dataSet$directionInx, paste0(res2[,1],res2[,3]))
      dataSet$regDirection =c(dataSet$regDirection, res2$Experiment)
      dataSet<<-dataSet
    mir.resu <<- rbind(mir.resu, res, res2);
    mirtargetu <<- c(mirtargetu, "tf");
    dataSet$tf2gene <<- res2;
    seedsu <<- c(seedsu, mir.nms);
    mirtableu <<- c(mirtableu, "tf2gene");
  }
  return(1);
}

.searchMultiNet_epi2mir <- function(input.type){
    orgType <- dataSet$org;
    if(orgType %in% c("bta", "dme","gga","sma", "cel","dre","rno", "ssc") ){
      curent.msg <<- "Only human and mouse are supported for epigene network."
      print(current.msg);
      return(0);
    }

    if (input.type == "mir2epigene"){
      idType <- dataSet$mirnaType;
      mir.mat <- dataSet$mir.orig;
      mir.vec <- rownames(mir.mat);
    }else{
      idType <- "epi_regulator";
      mir.mat <- dataSet$data[[idType]]
      mir.vec <- rownames(mir.mat);
    }

    mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2epi", sep=""), mir.vec, orgType, idType);

    hit.num <- nrow(mir.dic)
    if (hit.num == 0 && dataSet$tissue == "na") {
      current.msg <<- "No hits found in the miRNA-epigenetic modifier database. Please check your input. ";
      print(current.msg);
      return(0);
    } else if (hit.num == 0 && dataSet$tissue != "na") {
      current.msg <<- "No hits found in the miRNA-epigenetic modifier database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
      print(current.msg);
      return(2);
    } else {
      res <- mir.dic[ , c("mir_id", "mir_acc", "epi_regulator", "experiment", "condition", "pmid", "tissue")];
      rownames(res) <- mir.dic$mirnet;
      current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-epigene targets were identified!");

      # update the data
      gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
      mir.mat <- mir.mat[gd.inx,,drop=F];
      mir.mappedu <<- rbind(mir.mappedu, mir.mat);

      colnames(res) <- c("ID","Accession","Epigenetics","Experiment", "Condition","Literature", "Tissue");
      display.res <- res;
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Epigenetics"]);    # IDs
      if(nrow(res)!=0){
        row.names(edge.res) <- 1:nrow(res);
      }

      node.ids <- c(ID1=res[,"Accession"], ID2=res[,"Epigenetics"]);
      node.nms <- c(Name1=res[,"ID"], Name2=res[,"Epigenetics"]);

      edgeu.res <<- rbind(edgeu.res, edge.res); #edgeu.res is an empty dataframe defined in QueryNet
      edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ];
      nodeu.ids <<- c(nodeu.ids, node.ids);
      edgeNumU <<- c(edgeNumU, nrow(edge.res))
      nodeu.nms <<- c(nodeu.nms, node.nms);

      epi.nms <<- res[,"Epigenetics"];
      mir.nms <<- res[, "ID"];
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Epigenetics"], TargetID=res[,"Epigenetics"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"]);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "epigenetics");
      if(input.type == "mir2epigene"){
        dataSet$mir2epi <<- res;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(epi.nms)));
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2epi");
      } else{
        dataSet$epi2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(epi.nms)),Mapped=length(unique(mir.nms)));
        seedsu <<- c(seedsu, epi.nms);
        mirtableu <<- c(mirtableu, "epi2mir");
      }
    }
    return(1);
}

.searchMultiNet_disease2mir <-function(input.type){
    if(dataSet$org != "hsa" ){
      curent.msg <<- "Only human is supported for disease network."
      print(current.msg);
      return(0);
    }

    if (input.type == "mir2disease"){
      idType <- dataSet$mirnaType;
      mir.mat <- dataSet$mir.orig;
      mir.vec <- rownames(mir.mat);
    }else{
      idType <- "disease";
      mir.mat <- dataSet$data[[idType]];
      mir.vec <- rownames(mir.mat);
    }

    mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2disease", sep=""), mir.vec, "disease", idType);

    hit.num <- nrow(mir.dic)
    if (hit.num == 0 && dataSet$tissue == "na") {
      current.msg <<- "No hits found in the miRNA-disease database. Please check your input. ";
      print(current.msg);
      return(0);
    } else if (hit.num == 0 && dataSet$tissue != "na") {
      current.msg <<- "No hits found in the miRNA-disease database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
      print(current.msg);
      return(2);
    } else{
      res <- mir.dic[ , c("mir_id", "mir_acc", "disease", "method", "database", "pmid", "tissue")];
      rownames(res) <- mir.dic$mirnet;
      current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-disease associations were identified!");

      # update the data
      gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
      mir.mat <- mir.mat[gd.inx,,drop=F];
      mir.mappedu <<- rbind(mir.mappedu, mir.mat);

      colnames(res) <- c("ID","Accession","Disease","Experiment", "Database", "Literature", "Tissue");
      display.res <- res;
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Disease"]);    # IDs
      if(nrow(res)!=0){
        row.names(edge.res) <- 1:nrow(res);
      }

      node.ids <- c(ID1=res[,"Accession"], ID2=res[,"Disease"]);
      node.nms <- c(Name1=res[,"ID"], Name2=res[,"Disease"]);

      edgeu.res <<- rbind(edgeu.res, edge.res); #edgeu.res is an empty dataframe defined in QueryNet
      edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ];
      nodeu.ids <<- c(nodeu.ids, node.ids);
      edgeNumU <<- c(edgeNumU, nrow(edge.res))
      nodeu.nms <<- c(nodeu.nms, node.nms);

      dis.nms <<- res[,"Disease"];
      mir.nms <<- res[, "ID"];
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Disease"], TargetID=res[,"Disease"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"]);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "disease");
      if(input.type == "mir2disease"){
        dataSet$mir2dis <<- res;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(dis.nms)));
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2dis");
      } else{
        dataSet$dis2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(dis.nms)),Mapped=length(unique(mir.nms)));
        seedsu <<- c(seedsu, dis.nms);
        mirtableu <<- c(mirtableu, "dis2mir");
      }
    }
    return(1);
}

.searchMultiNet_protein2protein <-function(input.type){
    mir.vec <- unique(unname(nodeu.ids))
    table.nm = paste(data.org, dataSet$ppiOpts$db.name, sep="_")
    require.exp = dataSet$ppiOpts$require.exp
    min.score = as.numeric(dataSet$ppiOpts$min.score)

    mir.dic <- QueryPpiSQLiteZero(table.nm, mir.vec, require.exp, min.score);

    hit.num <- nrow(mir.dic)
    if (hit.num == 0 && dataSet$tissue == "na") {
      current.msg <<- "No hits found in the PPI database. Please check your input. ";
      print(current.msg);
      return(0);
    } else if (hit.num == 0 && dataSet$tissue != "na") {
      current.msg <<- "No hits found in the PPI database. The gene list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
      print(current.msg);
      return(2);
    } else{
      res <- mir.dic[ , c("id1", "id2", "name1", "name2")];
      rownames(res) <- mir.dic$mirnet;
      current.msg <<- paste("A total of unqiue", hit.num, "pairs of  protein-protein interactions were identified!");

      # update the data

      colnames(res) <- c("ID1","ID2","Symbol1","Symbol2");
      display.res <- res;
      edge.res <- data.frame(Source=res[,"Symbol1"],Target=res[,"Symbol2"]);    # IDs

      if(nrow(res)!=0){
        row.names(edge.res) <- 1:nrow(res);
      }

      node.ids <- c(ID1=res[,"Symbol1"], ID2=res[,"Symbol2"]);
      node.nms <- c(Name1=res[,"ID1"], Name2=res[,"ID2"]);
      edgeu.res <<- rbind(edgeu.res, edge.res); #edgeu.res is an empty dataframe defined in QueryNet
      edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ];
      nodeu.ids <<- c(nodeu.ids, node.ids);
      edgeNumU <<- c(edgeNumU, nrow(edge.res))
      nodeu.nms <<- c(nodeu.nms, node.nms);
      if(table.nm == "hsa_string"){
        res$Literature <- rep("25352553", nrow(res));
      }else if(table.nm == "hsa_innate"){
        res$Literature <- rep("23180781", nrow(res));
      }else if(table.nm == "hsa_rolland"){
        res$Literature <- rep("25416956", nrow(res));
      }else if(table.nm == "hsa_huri"){
        res$Literature <- rep("32296183", nrow(res));
      }

      na.res = rep("Not Applicable", nrow(res))
      res <- data.frame(ID=res[,"Symbol1"], Accession=res[,"ID1"], Target=res[,"Symbol2"], TargetID=res[,"ID2"], Experiment=na.res, Literature=res[,"Literature"], Tissue=na.res);

      mir.resu <<- rbind(mir.resu, res);
      if(input.type == "protein2protein"){
        dataSet$protein2protein <<- res;         # save this for network builder and table view
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "protein2protein");
      }
    }
    return(1);
}

#' Search Multi-Network Data
#' @param input.type miRNA and target, e.g., mir2gene, mir2molecule, mir2disease.
#' @export
SearchMultiNet <- function(input.type){

  node.ids <- vector();
  res <- 0;
  if (input.type %in% c("gene2mir","mir2gene","snpmir2gene", "mir2gene_mirtarbase", "mir2gene_tarbase", "mir2gene_mirecords", "mir2gene_miranda")){
    res <- .searchMultiNet_mir2gene(input.type);
  }else if (input.type %in% c("molecule2mir","mir2molecule","mol2mir","mir2mol")){
    res <- .searchMultiNet_molecule2mir(input.type)
  }else if (input.type  %in% c("lncrna2mir" ,"mir2lncrna" ,"lnc2mir" ,"mir2lnc")){
    res <- .searchMultiNet_lncrna2mir(input.type)
  }else if (input.type %in% c("pseudogene2mir", "mir2pseudogene", "pseudo2mir","mir2pseudo")){
    res <- .searchMultiNet_pseudogene2mir(input.type)
  }else if (input.type %in% c("sncrna2mir" ,"mir2sncrna", "snc2mir","mir2snc")){
    res <- .searchMultiNet_sncrna2mir(input.type)
  }else if (input.type %in% c("tf2mir", "mir2tf", "tf2gene", "gene2tf")){
    res <- .searchMultiNet_tf2mir(input.type)
  }else if (input.type %in% c("tf2gene2mir")){
    res <- .searchMultiNet_tf2gene2mir(input.type)
  }else if (input.type %in% c("gene2tf2mir")){
    res <- .searchMultiNet_gene2tf2mir(input.type)
  }else if (input.type %in% c("mir2tf2gene")){
    res <- .searchMultiNet_mir2tf2gene(input.type)
  }else if(input.type %in% c("epi2mir", "mir2epi")){
    res <- .searchMultiNet_epi2mir(input.type)
  }else if (input.type %in% c("disease2mir", "mir2disease", "dis2mir", "mir2dis")){
    res <- .searchMultiNet_disease2mir(input.type)
  }else if(input.type == "protein2protein"){
    res <- .searchMultiNet_protein2protein(input.type);
  }else if(input.type  %in% c("circrna2mir" ,"mir2circrna", "circ2mir","mir2circ")){
    res <- .searchMultiNet_cir2mir(input.type);
  }

  mir.nmsu <<- unique(c(mir.nmsu, mir.nms));
  return(res);
}

#' Set Protein-Protein Interaction Db
#' @param db Database name.
#' @param req Logical.
#' @param conf Min score.
#' @export
SetPpiDb  <- function(db, req, conf){
  dataSet$ppiOpts$db.name=db;
  if(req=="true"){
    dataSet$ppiOpts$require.exp=TRUE;
  }else{
    dataSet$ppiOpts$require.exp=FALSE;
  }
  dataSet$ppiOpts$min.score=as.numeric(conf)
  dataSet<<-dataSet
}
