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
#' @rdname QueryMultiList
#' @export 
QueryMultiList <- function(){

  net.info <<- list();
  mir.mappedu <<- matrix();
  mir.resu <<- data.frame();
  mirtargetu <<- vector();
  mirtableu <<- vector();
  seedsu <<- vector();
  edgeu.res <<- data.frame();
  nodeu.ids <<- vector();
  nodeu.nms <<- vector();

    for(i in 1:length(dataSet$type)){
      if (dataSet$type[i] == "mirna") {
        input.type=paste("mir2", dataSet$targetOpt, sep="");
        SearchMultiNet(input.type);
        if (dataSet$targetOpt == "gene") {
          net.info$gene.nms = gene.nms
        } else if (dataSet$targetOpt == "lncrna") {
          net.info$lnc.nms = lnc.nms
        } else if (dataSet$targetOpt == "tf") {
          net.info$tf.nms = tf.nms
        } else if (dataSet$targetOpt == "disease") {
          net.info$dis.nms = dis.nms
        } else if (dataSet$targetOpt == "molecule") {
          net.info$mol.nms = mol.nms
        } else if (dataSet$targetOpt == "epigene") {
          net.info$epi.nms = epi.nms
        }
      } else{
        input.type=paste(dataSet$type[i], "2mir", sep="");
        SearchMultiNet(input.type);
        if (dataSet$type[i] == "gene") {
          net.info$gene.nms = gene.nms
        } else if (dataSet$type[i] == "lncrna") {
          net.info$lnc.nms = lnc.nms
        } else if (dataSet$type[i] == "tf") {
          net.info$tf.nms = tf.nms
        } else if (dataSet$type[i] == "disease") {
          net.info$dis.nms = dis.nms
        } else if (dataSet$type[i] == "molecule") {
          net.info$mol.nms = mol.nms
        } else if (dataSet$type[i] == "epigene") {
          net.info$epi.nms = epi.nms
        }
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

  if (input.type == "gene2mir" || input.type == "mir2gene"){
    if (input.type == "mir2gene"){
      idType <- dataSet$mirnaType;
      mir.mat <- dataSet$mir.orig;
      idVec <- rownames(mir.mat);
    }else{
      idType <- dataSet$id.types[["gene"]];
      mir.mat <- dataSet$data[["gene"]];
      idVec <- rownames(mir.mat);
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
      nodeu.nms <<- c(nodeu.nms, node.nms);

      gene.nms <<- res[,"Gene"];
      mir.nms <<- res[, "ID"];
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"], stringsAsFactors = FALSE);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "gene");
      
      if(input.type == "mir2gene"){
        dataSet$mir2gene <<- display.res;         # save this for network builder and table view
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2gene");
      } else{
        dataSet$gene2mir <<- display.res;
        seedsu <<- c(seedsu, gene.nms);
        mirtableu <<- c(mirtableu, "gene2mir");
      }
    }
  }else if (input.type == "molecule2mir" || input.type == "mir2molecule"){
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
      nodeu.nms <<- c(nodeu.nms, node.nms);

      mol.nms <<- res[,"Molecule"];
      mir.nms <<- res[, "ID"];
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Molecule"], TargetID=res[,"Pubchem_ID"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"], stringsAsFactors = FALSE);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "molecule");
      if(input.type == "mir2molecule"){
        dataSet$mir2mol <<- display.res;         # save this for network builder and table view
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2mol");
      } else{
        dataSet$mol2mir <<- display.res;
        seedsu <<- c(seedsu, mol.nms);
        mirtableu <<- c(mirtableu, "mol2mir");
      }
    }
  }else if (input.type == "lncrna2mir" || input.type == "mir2lncrna"){
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
      nodeu.nms <<- c(nodeu.nms, node.nms);

      lnc.nms <<- res[,"Gene"];
      mir.nms <<- res[, "ID"];
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"], stringsAsFactors = FALSE);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "lncrna");
      if(input.type == "mir2lncrna"){
        dataSet$mir2lnc <<- display.res;         # save this for network builder and table view
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2lnc");
      } else{
        dataSet$lnc2mir <<- display.res;
        seedsu <<- c(seedsu, lnc.nms);
        mirtableu <<- c(mirtableu, "lnc2mir");
      }
    }
  }else if (input.type == "tf2mir" || input.type == "mir2tf"){
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
    }else{
      idType <- dataSet$id.types[["tf"]];
      mir.mat <- dataSet$data[["tf"]];
      mir.vec <- rownames(mir.mat);
    }

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

      edgeu.res <<- rbind(edgeu.res, edge.res); #edgeu.res is an empty dataframe defined in QueryNet
      edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ];
      nodeu.ids <<- c(nodeu.ids, node.ids);
      nodeu.nms <<- c(nodeu.nms, node.nms);

      tf.nms <<- res[,"Gene"];
      mir.nms <<- res[, "ID"];
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Gene"], TargetID=res[,"Entrez"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"], stringsAsFactors = FALSE);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "tf");
      if(input.type == "mir2tf"){
        dataSet$mir2tf <<- display.res;         # save this for network builder and table view
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2tf");
      } else{
        dataSet$tf2mir <<- display.res;
        seedsu <<- c(seedsu, tf.nms);
        mirtableu <<- c(mirtableu, "tf2mir");
      }
    }
  }else if(input.type == "epigene2mir" || input.type == "mir2epigene"){
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
      nodeu.nms <<- c(nodeu.nms, node.nms);

      epi.nms <<- res[,"Epigenetics"];
      mir.nms <<- res[, "ID"];
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Epigenetics"], TargetID=res[,"Epigenetics"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"],stringsAsFactors = FALSE);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "epigenetics");
      if(input.type == "mir2epigene"){
        dataSet$mir2epi <<- display.res;         # save this for network builder and table view
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2epi");
      } else{
        dataSet$epi2mir <<- display.res;
        seedsu <<- c(seedsu, epi.nms);
        mirtableu <<- c(mirtableu, "epi2mir");
      }
    }
  }else if (input.type == "disease2mir" || input.type == "mir2disease"){
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
      nodeu.nms <<- c(nodeu.nms, node.nms);

      dis.nms <<- res[,"Disease"];
      mir.nms <<- res[, "ID"];
      res <- data.frame(ID=res[,"ID"], Accession=res[,"Accession"], Target=res[,"Disease"], TargetID=res[,"Disease"], Experiment=res[,"Experiment"], Literature=res[,"Literature"], Tissue=res[,"Tissue"], stringsAsFactors = FALSE);
      mir.resu <<- rbind(mir.resu, res);
      mirtargetu <<- c(mirtargetu, "disease");
      if(input.type == "mir2disease"){
        dataSet$mir2dis <<- display.res;         # save this for network builder and table view
        seedsu <<- c(seedsu, mir.nms);
        mirtableu <<- c(mirtableu, "mir2dis");
      } else{
        dataSet$dis2mir <<- display.res;
        seedsu <<- c(seedsu, dis.nms);
        mirtableu <<- c(mirtableu, "dis2mir");
      }
    }
  }
}

