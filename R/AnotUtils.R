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
        dataSet$mirtarget <- "gene";
        dataSet <<- dataSet;
        return(1);
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
#' @rdname PerformMolMapping
#' @export 
PerformMolMapping <- function(){
    orgType <- dataSet$org;
    if(orgType %in% c("bta", "dme", "gga","sma", "cel", "ssc")){
       curent.msg <<- "This organism is not supported for molecule network research."
       print(current.msg);
       return(0);
    }

    mir.mat <- dataSet$mir.orig;
    idType <- dataSet$idType;
    mir.vec <- rownames(mir.mat);
    mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2molecule", sep=""), mir.vec, orgType, idType);

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
        res <- mir.dic[, c("mir_id","mir_acc","molecule", "pubchem_id", "method", "pmid", "tissue")];
        rownames(res) <- mir.dic$mirnet;
        current.msg <<- paste("A total of unqiue", hit.num, "pairs of miRNA-gene targets were identified!");

        # update the data
        gd.inx <- rownames(mir.mat) %in% unique(res[, idType]);
        dataSet$mir.mapped <- mir.mat[gd.inx,,drop=F];

        colnames(res) <- c("ID","Accession","Molecule", "Pubchem_ID", "Experiment", "Literature", "Tissue");
        mol.nms <- res[,"Molecule"];
        mir.nms <- res[, "ID"];
        if(dataSet$idType %in% c("mir_id", "mir_acc")){
          dataSet$seeds <- mir.nms;
        }else{
          dataSet$seeds <- mol.nms;
        }
        write.csv(res, file="mirnet_mir_target.csv", row.names=FALSE);
        dataSet$mir.res <- res;
        dataSet$mirtarget <- "molecule";
        dataSet <<- dataSet;
        return(1);
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
#' @rdname PerformDisMapping
#' @export 
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
        write.csv(res, file="mirnet_mir_target.csv", row.names=FALSE);
        dataSet$mir.res <- res;
        dataSet$mirtarget <- "disease";
        dataSet <<- dataSet;
        return(1);
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
#' @rdname PerformMir2EpiMapping
#' @export 
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
        write.csv(res, file="mirnet_mir_target.csv", row.names=FALSE);
        dataSet$mir.res <- res;
        dataSet$mirtarget <- "epigenetics";
        dataSet <<- dataSet;
        return(1);
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
#' @rdname PerformLncRNAMapping
#' @export 
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
        write.csv(res, file="mirnet_mir_target.csv", row.names=FALSE);
        dataSet$mir.res <- res;
        dataSet$mirtarget <- "lncrna";
        dataSet <<- dataSet;
        return(1);
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
#' @rdname PerformTFMapping
#' @export 
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
            write.csv(res, file="mirnet_mir_target.csv", row.names=FALSE);
            dataSet$mir.res <- res;
            dataSet$mirtarget <- "tf";
            dataSet <<- dataSet;
            return(1);
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
    snp.dic <- Query.miRNetDB(paste(sqlite.path, "snp2mir", sep=""), snpidVec, dataSet$org, dataSet$idType);

    hit.num <- nrow(snp.dic)
    if (hit.num == 0) {
        current.msg <<- "No hits found in the database. The SNP list has not been annotated to miRNA. Please check your input.";
        print(current.msg);
        return(0);
    } else {
        snp <- na.omit(data.frame(name1 = snp.dic[, idType], id1 = snp.dic[, "rsid"], name2 = snp.dic[, "Mature_Name"], id2 =  snp.dic[, "Mature_Acc"], stringsAsFactors = FALSE));
        hit.num <- nrow(snp)
        idVec <- as.vector(unique(snp.dic[, c("MIRNA_Name")]));
        mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2gene", sep=""), idVec, dataSet$org, "mir_id");

        snp.res <- na.omit(snp.dic[ , c("chr_pos", "rsid", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain")]);
        mir.res <- mir.dic[ , c("mir_id", "mir_acc", "symbol", "entrez", "experiment", "pmid", "tissue")];
        merge.res <- merge(snp.res, mir.res, by.x = "Mature_Name", by.y = "mir_id", all.x = T);
        merge.res <- na.omit(merge.res[, c("chr_pos", "rsid", "MIRNA_Name", "MIRNA_Acc", "MIRNA_Domain", "Mature_Name", "Mature_Acc", "symbol", "entrez", "experiment", "pmid", "tissue")])
        rownames(merge.res) <- 1:nrow(merge.res);
        current.msg <<- paste("A total of", hit.num, "SNPs were mapped to miRNAs!");

        # update the data
        gd.inx <- rownames(snp.mat) %in% unique(merge.res[, idType]);
        dataSet$mir.mapped <- snp.mat[gd.inx,,drop=F];

        # update col names
        colnames(merge.res) <- c("CHR_POS", "rsID", "miRNA_Name", "miRNA_Accession", "miRNA_Domain", "Mature_miRNA", "Mature_Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue");
        write.csv(merge.res, file="mirnet_snp_mir_target.csv", row.names=FALSE);
        snp.res$Database <- rep("30302893", nrow(snp.res));
        colnames(snp.res) <- c("CHR_POS", "rsID", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain", "DataBase");
        colnames(mir.res) <- c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue");
        dataSet$seeds <- snp.res[, "rsID"]; 
        dataSet$mir.res <- merge.res;
        dataSet$mirtarget <- "gene";
        dataSet$mirtable <- c("snp2mir", "mir2gene");
        dataSet$snp2mir <- snp.res;
        dataSet$mir2gene <- mir.res;
        dataSet <<- dataSet;
        return(1);
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
#' @rdname PerformSNPMirDisMapping
#' @export 
PerformSNPMirDisMapping <- function(){
    snp.mat <- dataSet$mir.orig;
    snpidVec <- rownames(snp.mat);
    idType <- dataSet$idType;
    snp.dic <- Query.miRNetDB(paste(sqlite.path, "snp2mir", sep=""), snpidVec, dataSet$org, dataSet$idType);

    hit.num <- nrow(snp.dic)
    if (hit.num == 0) {
        current.msg <<- "No hits found in the database. The SNP list has not been annotated to miRNA. Please check your input.";
        print(current.msg);
        return(0);
    } else {
        snp <- na.omit(data.frame(name1 = snp.dic[, idType], id1 = snp.dic[, "rsid"], name2 = snp.dic[, "Mature_Name"], id2 =  snp.dic[, "Mature_Acc"], stringsAsFactors = FALSE));
        hit.num <- nrow(snp)
        idVec <- as.vector(unique(snp.dic[, c("MIRNA_Name")]));
        mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2disease", sep=""), idVec, "disease", "mir_id");

        snp.res <- na.omit(snp.dic[ , c("chr_pos", "rsid", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain")]);
        mir.res <- mir.dic[ , c("mir_id", "mir_acc", "disease", "method", "database", "pmid", "tissue")];
        merge.res <- merge(snp.res, mir.res, by.x = "MIRNA_Name", by.y = "mir_id", all.x = T);
        merge.res <- na.omit(merge.res[, c("chr_pos", "rsid", "MIRNA_Name", "MIRNA_Acc", "MIRNA_Domain", "Mature_Name", "Mature_Acc", "disease", "method", "database", "pmid", "tissue")])
        rownames(merge.res) <- 1:nrow(merge.res);
        current.msg <<- paste("A total of", hit.num, "SNPs were mapped to miRNAs!");

        # update the data
        gd.inx <- rownames(snp.mat) %in% unique(merge.res[, idType]);
        dataSet$mir.mapped <- snp.mat[gd.inx,,drop=F];

        # update col names
        colnames(merge.res) <- c("CHR_POS", "rsID", "miRNA_Name", "miRNA_Accession", "miRNA_Domain", "Mature_miRNA", "Mature_Accession", "Disease", "Experiment", "Database", "Literature", "Tissue");
        write.csv(merge.res, file="mirnet_snp_mir_target.csv", row.names=FALSE);
        snp.res$Database <- rep("30302893", nrow(snp.res));
        colnames(snp.res) <- c("CHR_POS", "rsID", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain", "DataBase");
        colnames(mir.res) <- c("ID","Accession","Disease","Experiment", "Database", "Literature", "Tissue");
        dataSet$seeds <- snp.res[, "rsID"]; 
        dataSet$mir.res <- merge.res;
        dataSet$mirtarget <- "disease";
        dataSet$mirtable <- c("snp2mir", "mir2dis");
        dataSet$snp2mir <- snp.res;
        dataSet$mir2dis <- mir.res;
        dataSet <<- dataSet;
        return(1);
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
#' @rdname PerformSNPMirMolMapping
#' @export 
PerformSNPMirMolMapping <- function(){
    snp.mat <- dataSet$mir.orig;
    snpidVec <- rownames(snp.mat);
    idType <- dataSet$idType;
    snp.dic <- Query.miRNetDB(paste(sqlite.path, "snp2mir", sep=""), snpidVec, dataSet$org, dataSet$idType);

    hit.num <- nrow(snp.dic)
    if (hit.num == 0) {
        current.msg <<- "No hits found in the database. The SNP list has not been annotated to miRNA. Please check your input.";
        print(current.msg);
        return(0);
    } else {
        snp <- na.omit(data.frame(name1 = snp.dic[, idType], id1 = snp.dic[, "rsid"], name2 = snp.dic[, "Mature_Name"], id2 =  snp.dic[, "Mature_Acc"], stringsAsFactors = FALSE));
        hit.num <- nrow(snp)
        idVec <- as.vector(unique(snp.dic[, c("MIRNA_Name")]));
        mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2molecule", sep=""), idVec, "hsa", "mir_id");

        snp.res <- na.omit(snp.dic[ , c("chr_pos", "rsid", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain")]);
        mir.res <- mir.dic[ , c("mir_id", "mir_acc", "molecule", "pubchem_id", "method", "pmid", "tissue")];
        merge.res <- merge(snp.res, mir.res, by.x = "MIRNA_Name", by.y = "mir_id", all.x = T);
        merge.res <- na.omit(merge.res[, c("chr_pos", "rsid", "MIRNA_Name", "MIRNA_Acc", "MIRNA_Domain", "Mature_Name", "Mature_Acc", "molecule", "pubchem_id", "method", "pmid", "tissue")])
        rownames(merge.res) <- 1:nrow(merge.res);
        current.msg <<- paste("A total of", hit.num, "SNPs were mapped to miRNAs!");

        # update the data
        gd.inx <- rownames(snp.mat) %in% unique(merge.res[, idType]);
        dataSet$mir.mapped <- snp.mat[gd.inx,,drop=F];

        # update col names
        colnames(merge.res) <- c("CHR_POS", "rsID", "miRNA_Name", "miRNA_Accession", "miRNA_Domain", "Mature_miRNA", "Mature_Accession", "Molecule", "Pubchem_ID", "Experiment", "Literature", "Tissue");
        write.csv(merge.res, file="mirnet_snp_mir_target.csv", row.names=FALSE);
        snp.res$Database <- rep("30302893", nrow(snp.res));
        colnames(snp.res) <- c("CHR_POS", "rsID", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain", "DataBase");
        colnames(mir.res) <- c("ID","Accession","Molecule", "Pubchem_ID", "Experiment", "Literature", "Tissue");
        dataSet$seeds <- snp.res[, "rsID"]; 
        dataSet$mir.res <- merge.res;
        dataSet$mirtarget <- "molecule";
        dataSet$mirtable <- c("snp2mir", "mir2mol");
        dataSet$snp2mir <- snp.res;
        dataSet$mir2mol <- mir.res;
        dataSet <<- dataSet;
        return(1);
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
#' @rdname PerformSNPMirLncMapping
#' @export 
PerformSNPMirLncMapping <- function(){
  snp.mat <- dataSet$mir.orig;
  snpidVec <- rownames(snp.mat);
  idType <- dataSet$idType;
  snp.dic <- Query.miRNetDB(paste(sqlite.path, "snp2mir", sep=""), snpidVec, dataSet$org, dataSet$idType);

  hit.num <- nrow(snp.dic)
  if (hit.num == 0) {
    current.msg <<- "No hits found in the database. The SNP list has not been annotated to miRNA. Please check your input.";
    print(current.msg);
    return(0);
  } else {
    snp <- na.omit(data.frame(name1 = snp.dic[, idType], id1 = snp.dic[, "rsid"], name2 = snp.dic[, "Mature_Name"], id2 =  snp.dic[, "Mature_Acc"], stringsAsFactors = FALSE));
    hit.num <- nrow(snp)
    idVec <- as.vector(unique(snp.dic[, c("MIRNA_Name")]));
    mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2lncRNA", sep=""), idVec, dataSet$org, "mir_id");

    snp.res <- na.omit(snp.dic[ , c("chr_pos", "rsid", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain")]);
    mir.res <- mir.dic[ , c("mir_id", "mir_acc", "symbol", "entrez", "tissue")];

    merge.res <- merge(snp.res, mir.res, by.x = "MIRNA_Name", by.y = "mir_id", all.x = T);
    merge.res <- na.omit(merge.res[, c("chr_pos", "rsid", "MIRNA_Name", "MIRNA_Acc", "MIRNA_Domain", "Mature_Name", "Mature_Acc", "symbol", "entrez", "tissue")]);
    rownames(merge.res) <- 1:nrow(merge.res);
    current.msg <<- paste("A total of", hit.num, "SNPs were mapped to miRNAs!");

    # update the data
    gd.inx <- rownames(snp.mat) %in% unique(merge.res[, idType]);
    dataSet$mir.mapped <- snp.mat[gd.inx,,drop=F];

    # update col names
    colnames(merge.res) <- c("CHR_POS", "rsID", "miRNA_Name", "miRNA_Accession", "miRNA_Domain", "Mature_miRNA", "Mature_Accession", "Gene", "Entrez", "Tissue");
    merge.res$Experiment <- rep("CLIP-Seq", nrow(merge.res));
    merge.res$Literature <- rep("24297251", nrow(merge.res));
    merge.res <- merge.res[, c("CHR_POS", "rsID", "miRNA_Name", "miRNA_Accession", "miRNA_Domain", "Mature_miRNA", "Mature_Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue")];
    write.csv(merge.res, file="mirnet_snp_mir_target.csv", row.names=FALSE);
    snp.res$Database <- rep("30302893", nrow(snp.res));
    colnames(snp.res) <- c("CHR_POS", "rsID", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain", "DataBase");
    mir.res <- merge.res[, c("miRNA_Name", "miRNA_Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue")];
    colnames(mir.res) <- c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue");
    dataSet$seeds <- snp.res[, "rsID"]; 
    dataSet$mir.res <- merge.res;
    dataSet$mirtarget <- "lncrna";
    dataSet$mirtable <- c("snp2mir", "mir2lnc");
    dataSet$snp2mir <- snp.res;
    dataSet$mir2lnc <- mir.res;
    dataSet <<- dataSet;
    return(1);
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
#' @rdname PerformSNPMirTFMapping
#' @export 
PerformSNPMirTFMapping <- function(){
  snp.mat <- dataSet$mir.orig;
  snpidVec <- rownames(snp.mat);
  idType <- dataSet$idType;
  snp.dic <- Query.miRNetDB(paste(sqlite.path, "snp2mir", sep=""), snpidVec, dataSet$org, dataSet$idType);

  hit.num <- nrow(snp.dic)
  if (hit.num == 0) {
    current.msg <<- "No hits found in the database. The SNP list has not been annotated to miRNA. Please check your input.";
    print(current.msg);
    return(0);
  } else {
    snp <- na.omit(data.frame(name1 = snp.dic[, idType], id1 = snp.dic[, "rsid"], name2 = snp.dic[, "Mature_Name"], id2 =  snp.dic[, "Mature_Acc"], stringsAsFactors = FALSE));
    hit.num <- nrow(snp)
    idVec <- as.vector(unique(snp.dic[, c("MIRNA_Name")]));
    mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2tf", sep=""), idVec, dataSet$org, "mir_id");

    snp.res <- na.omit(snp.dic[ , c("chr_pos", "rsid", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain")]);
    mir.res <- mir.dic[ , c("mir_id", "mir_acc", "symbol","entrez", "pmid", "tissue")];

    merge.res <- merge(snp.res, mir.res, by.x = "MIRNA_Name", by.y = "mir_id", all.x = T);
    merge.res <- na.omit(merge.res[, c("chr_pos", "rsid", "MIRNA_Name", "MIRNA_Acc", "MIRNA_Domain", "Mature_Name", "Mature_Acc", "symbol","entrez", "pmid", "tissue")]);
    rownames(merge.res) <- 1:nrow(merge.res);
    current.msg <<- paste("A total of", hit.num, "SNPs were mapped to miRNAs!");

    # update the data
    gd.inx <- rownames(snp.mat) %in% unique(merge.res[, idType]);
    dataSet$mir.mapped <- snp.mat[gd.inx,,drop=F];

    # update col names
    colnames(merge.res) <- c("CHR_POS", "rsID", "miRNA_Name", "miRNA_Accession", "miRNA_Domain", "Mature_miRNA", "Mature_Accession", "Gene", "Entrez", "Literature", "Tissue");
    merge.res$Experiment <- rep("ChIP-seq", nrow(merge.res));
    merge.res <- merge.res[, c("CHR_POS", "rsID", "miRNA_Name", "miRNA_Accession", "miRNA_Domain", "Mature_miRNA", "Mature_Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue")];
    write.csv(merge.res, file="mirnet_snp_mir_target.csv", row.names=FALSE);
    snp.res$Database <- rep("30302893", nrow(snp.res));
    colnames(snp.res) <- c("CHR_POS", "rsID", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain", "DataBase");
    mir.res <- merge.res[, c("miRNA_Name", "miRNA_Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue")];
    colnames(mir.res) <- c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue");
    dataSet$seeds <- snp.res[, "rsID"]; 
    dataSet$mir.res <- merge.res;
    dataSet$mirtarget <- "tf";
    dataSet$mirtable <- c("snp2mir", "mir2tf");
    dataSet$snp2mir <- snp.res;
    dataSet$mir2tf <- mir.res;
    dataSet <<- dataSet;
    return(1);
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
