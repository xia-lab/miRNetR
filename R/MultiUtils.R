##################################################
## R script for miRNet
## Description: Gene/Compound Annotation
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

.init.multilist <- function(){
  anal.type <<- "multilist"
  net.info <<- list();
  mir.mappedu <<- matrix();
  mir.resu <<- data.frame();
  mirtargetu <<- vector();
  mirtableu <<- vector();
  seedsu <<- vector();
  
  edgeNumU <<- vector();
  edgeu.res <<- data.frame();
  nodeu.ids <<- vector();
  nodeu.nms <<- vector();
  mir.nmsu <<- vector();
  tf.nms <<- vector();
  gene.nms <<- vector()

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
#' @rdname QueryMultiList
#' @export 
QueryMultiList <- function(){
  .init.multilist();
print(dataSet$type)
  for(i in 1:length(dataSet$type)){
    if (dataSet$type[i] == "mirna") {
      input.type=paste("mir2", dataSet$targetOpt, sep="");
      SearchMultiNet(input.type);
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
      SearchMultiNet(input.type);
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
    }
  }
  if("tf" %in% dataSet$type[i]){
    #SearchMultiNet("tf2gene");
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

  write.csv(node.res, file="orig_node_list.csv", row.names=FALSE);
  write.csv(mir.resu, file="mirnet_mir_target.csv", row.names=FALSE);

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
  return(1);
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
#' @rdname QueryMultiListMir
#' @export 
QueryMultiListMir <- function(){
  .init.multilist();

  for(i in 1:length(dataSet$type)){
    if(dataSet$type[i] == "protein2protein"){
      input.type=dataSet$type[i]
    }else{
      input.type=paste("mir2", dataSet$type[i], sep="");
    }
    SearchMultiNet(input.type);
    if (dataSet$type[i] == "gene") {
      net.info$gene.nms = gene.nms
    } else if (dataSet$type[i] == "lncrna") {
      net.info$lnc.nms = lnc.nms
    } else if (dataSet$type[i] == "circrna") {
      net.info$circ.nms = circ.nms
    } else if (dataSet$type[i] == "pseudogene") {
      net.info$pseudo.nms = pseudo.nms
    } else if (dataSet$type[i] == "sncrna") {
      net.info$snc.nms = snc.nms
    } else if (dataSet$type[i] == "tf") {
      net.info$tf.nms = tf.nms
    } else if (dataSet$type[i] == "disease") {
      net.info$dis.nms = dis.nms
    } else if (dataSet$type[i] == "molecule") {
      net.info$mol.nms = mol.nms
    } else if (dataSet$type[i] == "epig") {
      net.info$epi.nms = epi.nms
    }
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

  write.csv(node.res, file="orig_node_list.csv", row.names=FALSE);
  write.csv(mir.resu, file="mirnet_mir_target.csv", row.names=FALSE);

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

  return(1);
}

.searchMultiNet_mir2gene<-function(input.type){

    if (input.type == "mir2gene"){
      idType <- dataSet$mirnaType;
      mir.mat <- dataSet$mir.orig;
      idVec <- rownames(mir.mat);
    }else if(input.type=="snpmir2gene"){
      idType <- "mir_id";
      mir.mat <- data.frame(unique(mir.nms))
      rownames(mir.mat) = unique(mir.nms)
      mir.mat[,1] = 0
      idVec <- rownames(mir.mat);
    }else{
      if(length(dataSet$data[["gene"]])!= 0){
        idType <- dataSet$id.types[["gene"]];
        mir.mat <- dataSet$data[["gene"]];
        idVec <- rownames(mir.mat);
      }else{
        idType <- dataSet$mirnaType;
        mir.mat <- dataSet$mir.orig;
        idVec <- rownames(mir.mat);
      }
    }

    mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2gene", sep=""), idVec, dataSet$org, idType);

    hit.num <- nrow(mir.dic)
    if (hit.num == 0 && dataSet$tissue == "na") {
      current.msg <<- "No hits found in the miRNA-gene database. Please check your input.";
      print(current.msg);
      return(0);
    } else if (hit.num == 0 && dataSet$tissue != "na") {
      current.msg <<- "No hits found in the miRNA-gene database. The miRNA list has not been annotated by this tissue type. Please try NOT to specify the tissue.";
      print(current.msg);
      return(2);
    } else {
      res <- mir.dic[ , c("mir_id", "mir_acc", "symbol", "entrez", "experiment", "pmid", "tissue")];
      rownames(res) <- mir.dic$mirnet;
      current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-gene targets were identified!");

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
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Entrez"],stringsAsFactors = FALSE);    # IDs
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
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"], stringsAsFactors = FALSE);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "gene");

      if(input.type == "mir2gene"){
        dataSet$mir2gene <<- res;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(gene.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2gene");
      } else{
        dataSet$gene2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(gene.nms)),Mapped=length(unique(mir.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, gene.nms);
        mirtableu <<- c(mirtableu, "gene2mir");
      }
    }
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
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Molecule"],stringsAsFactors = FALSE);    # IDs
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
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Molecule"], TargetID=res[,"Pubchem_ID"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"], stringsAsFactors = FALSE);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "molecule");
      if(input.type == "mir2molecule"){
        dataSet$mir2mol <<- res;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(mol.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2mol");
      } else{
        dataSet$mol2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(mol.nms)),Mapped=length(unique(mir.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, mol.nms);
        mirtableu <<- c(mirtableu, "mol2mir");
      }
    }
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
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Entrez"],stringsAsFactors = FALSE);    # IDs
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
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"], stringsAsFactors = FALSE);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "lncrna");
      if(input.type == "mir2lncrna"){
        dataSet$mir2lnc <<- res;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(lnc.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2lnc");
      } else{
        dataSet$lnc2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(lnc.nms)),Mapped=length(unique(mir.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, lnc.nms);
        mirtableu <<- c(mirtableu, "lnc2mir");
      }
    }
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
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Entrez"],stringsAsFactors = FALSE);    # IDs
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
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"], stringsAsFactors = FALSE);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "circrna");
      if(input.type == "mir2circrna"){
        dataSet$mir2circ <<- res;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(circ.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2circ");
      } else{
        dataSet$circ2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(circ.nms)),Mapped=length(unique(mir.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, circ.nms);
        mirtableu <<- c(mirtableu, "circ2mir");
      }
    }
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
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Entrez"],stringsAsFactors = FALSE);    # IDs
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
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"], stringsAsFactors = FALSE);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "pseudogene");
      if(input.type == "mir2pseudogene"){
        dataSet$mir2pseudo <<- res;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(pseudo.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2pseudo");
      } else{
        dataSet$pseudo2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(pseudo.nms)),Mapped=length(unique(mir.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, pseudo.nms);
        mirtableu <<- c(mirtableu, "pseudo2mir");
      }
    }
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
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Entrez"],stringsAsFactors = FALSE);    # IDs
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
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"], stringsAsFactors = FALSE);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "sncrna");
      if(input.type == "mir2sncrna"){
        dataSet$mir2snc <<- res;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(snc.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2snc");
      } else{
        dataSet$snc2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(snc.nms)),Mapped=length(unique(mir.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, snc.nms);
        mirtableu <<- c(mirtableu, "snc2mir");
      }
    }
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
      if(targetType == "encode"){
        table.nm <- paste(targetType, dataSet$org, sep="_");
      }else{
        table.nm <- toupper(targetType);
      }
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
        res <- mir.dic[ , c("mir_id","mir_acc","symbol","entrez", "pmid", "tissue")];
        rownames(res) <- mir.dic$mirnet;
        current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-TF targets were identified!");

        # update the data
        gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
        mir.mat <- mir.mat[gd.inx,,drop=F];
        mir.mappedu <<- rbind(mir.mappedu, mir.mat);

        colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Literature", "Tissue");
        res$Experiment <- rep("ChIP-seq", nrow(res));
        res <- res[, c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue")];
        display.res <- res;
        edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Entrez"],stringsAsFactors = FALSE);    # IDs
        if(nrow(res)!=0){
          row.names(edge.res) <- 1:nrow(res);
        }
        node.ids <- c(ID1=res[,"Accession"], ID2=res[,"Entrez"]);
        node.nms <- c(Name1=res[,"ID"], Name2=res[,"Gene"]);
      }else{
        res <- mir.dic[ , c("symbol","entrez","tfname","tfid")];
        rownames(res) <- mir.dic$mirnet;
        current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-TF targets were identified!");

        # update the data
        gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
        mir.mat <- mir.mat[gd.inx,,drop=F];
        mir.mappedu <<- rbind(mir.mappedu, mir.mat);

        colnames(res) <- c("ID", "Accession", "Gene", "Entrez");
        res$Experiment <- rep("ChIP-seq", nrow(res));
        res <- res[, c("ID", "Accession", "Gene", "Entrez")];
        display.res <- res;
        edge.res <- data.frame(Source=res[,3],Target=res[,1],stringsAsFactors = FALSE);    # IDs
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
        res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"], stringsAsFactors = FALSE);
      }else{
        gene.nms <<- unique(c(gene.nms, as.vector(res[,"Gene"])));
        tf.nms <<- unique(c(tf.nms, as.vector(res[,"ID"])));
        mir.nms <<- "";
        if(targetType == "encode"){
          res$Literature <- rep("29126249", nrow(res));
          res$Experiment <- rep("ChIP-seq", nrow(res));
        }else if(targetType == "jaspar"){
          res$Literature <- rep("31701148", nrow(res));
          res$Experiment <- rep("JASPAR", nrow(res));
        }else if(targetType == "chea"){
          res$Literature <- rep("20709693", nrow(res));
          res$Experiment <- rep("ChIP-X", nrow(res));
        }
          res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"],Tissue=rep("Not Applicable", nrow(res)),stringsAsFactors = FALSE);
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
}

.searchMultiNet_epi2mir <- function(input.type){
orgType <- dataSet$org;
    if(orgType %in% c("bta", "dme","gga","sma", "cel","dre","rno", "ssc") ){
      curent.msg <<- "Only huamn and mouse support the epigene network."
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
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Epigenetics"],stringsAsFactors = FALSE);    # IDs
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
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Epigenetics"], TargetID=res[,"Epigenetics"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"],stringsAsFactors = FALSE);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "epigenetics");
      if(input.type == "mir2epigene"){
        dataSet$mir2epi <<- res;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(epi.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2epi");
      } else{
        dataSet$epi2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(epi.nms)),Mapped=length(unique(mir.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, epi.nms);
        mirtableu <<- c(mirtableu, "epi2mir");
      }
    }
}

.searchMultiNet_disease2mir <-function(input.type){
    if(dataSet$org != "hsa" ){
      curent.msg <<- "Only huamn support the disease network."
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
      edge.res <- data.frame(Source=res[,"Accession"],Target=res[,"Disease"],stringsAsFactors = FALSE);    # IDs
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
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Disease"], TargetID=res[,"Disease"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"], stringsAsFactors = FALSE);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "disease");
      if(input.type == "mir2disease"){
        dataSet$mir2dis <<- res;         # save this for network builder and table view
        dataSet$tableStats <<- data.frame(Query=length(unique(mir.nms)),Mapped=length(unique(dis.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2dis");
      } else{
        dataSet$dis2mir <<- res;
        dataSet$tableStats <<- data.frame(Query=length(unique(dis.nms)),Mapped=length(unique(mir.nms)),stringsAsFactors = FALSE);
        seedsu <<- c(seedsu, dis.nms);
        mirtableu <<- c(mirtableu, "dis2mir");
      }
    }
}

.searchMultiNet_protein2protein <-function(input.type){
mir.vec <- unique(unname(nodeu.ids))
    table.nm = paste("ppi", data.org, dataSet$ppiOpts$db.name, sep="_")
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
      edge.res <- data.frame(Source=res[,"Symbol1"],Target=res[,"Symbol2"],stringsAsFactors = FALSE);    # IDs

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
      if(table.nm == "ppi_hsa_string"){
        res$Literature <- rep("25352553", nrow(res));
      }else if(table.nm == "ppi_hsa_innate"){
        res$Literature <- rep("23180781", nrow(res));
      }else if(table.nm == "ppi_hsa_rolland"){
        res$Literature <- rep("25416956", nrow(res));
      }

      na.res = rep("Not Applicable", nrow(res))
      res <- data.frame(ID=res[,"Symbol1"], Accession=res[,"ID1"], Target=res[,"Symbol2"], TargetID=res[,"ID2"], Experiment=na.res, Literature=res[,"Literature"], Tissue=na.res,stringsAsFactors = FALSE);

      mir.resu <<- rbind(mir.resu, res);
      if(input.type == "protein2protein"){
        dataSet$protein2protein <<- res;         # save this for network builder and table view
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "protein2protein");
      }

    }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param input.type PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SearchMultiNet
#' @export 
SearchMultiNet <- function(input.type){

  node.ids <- vector()
  if (input.type %in% c("gene2mir","mir2gene","snpmir2gene")){
    .searchMultiNet_mir2gene(input.type);
  }else if (input.type %in% c("molecule2mir","mir2molecule","mol2mir","mir2mol")){
    .searchMultiNet_molecule2mir(input.type)
  }else if (input.type  %in% c("lncrna2mir" ,"mir2lncrna" ,"lnc2mir" ,"mir2lnc")){
    .searchMultiNet_lncrna2mir(input.type)
  }else if (input.type %in% c("pseudogene2mir", "mir2pseudogene", "pseudo2mir","mir2pseudo")){
    .searchMultiNet_pseudogene2mir(input.type)
  }else if (input.type %in% c("sncrna2mir" ,"mir2sncrna", "snc2mir","mir2snc")){
    .searchMultiNet_sncrna2mir(input.type)
  }else if (input.type %in% c("tf2mir", "mir2tf", "tf2gene", "gene2tf")){
    .searchMultiNet_tf2mir(input.type)
  }else if(input.type %in% c("epi2mir", "mir2epi")){
    .searchMultiNet_epi2mir(input.type)
  }else if (input.type %in% c("disease2mir", "mir2disease", "dis2mir", "mir2dis")){
    .searchMultiNet_disease2mir(input.type)
  }else if(input.type == "protein2protein"){
    .searchMultiNet_protein2protein(input.type);
  }else if(input.type  %in% c("circrna2mir" ,"mir2circrna", "circ2mir","mir2circ")){
    .searchMultiNet_cir2mir(input.type);
}
  mir.nmsu <<- unique(c(mir.nmsu, mir.nms));
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param db PARAM_DESCRIPTION
#' @param req PARAM_DESCRIPTION
#' @param conf PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SetPpiDb
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
