my.snp.mir.mapping <- function(){
  snp.mat <- dataSet$mir.orig;
  snpidVec <- rownames(snp.mat);
  idType <- dataSet$idType;

  # first try to match snp2mir, if still have unmatched, do snp2mirbs
  snp.dic <- Query.miRNetDB(paste(sqlite.path, "snp2mir", sep=""), snpidVec, dataSet$org, dataSet$idType);

  hit.inx <- match(snpidVec, snp.dic$rsid);

  na.hits <- is.na(hit.inx);
  if(sum(na.hits) > 0){
    unmatched.snp <- snpidVec[na.hits];
    snp.dic2 <- Query.miRNetDB(paste(sqlite.path, "snp2mirbs", sep=""), unmatched.snp, dataSet$org, dataSet$idType);
    snp.num <- rbind(snp.dic[,1:2], snp.dic2[,1:2])
  }
  snp.num <- snp.dic;
  hit.num <- nrow(snp.num)
  if (hit.num == 0) {
    current.msg <<- "No hits found in the database. The SNP list has not been annotated to miRNA or miRNA-binding sites. Please check your input.";
    print(current.msg);
    return(0);
  } else {
    if(sum(na.hits) > 0){
      # snp2mir2gene
      snp <- na.omit(data.frame(name1 = snp.dic[, idType], id1 = snp.dic[, "rsid"], name2 = snp.dic[, "Mature_Name"], id2 =  snp.dic[, "Mature_Acc"], stringsAsFactors = FALSE));
      hit.num <- nrow(snp);
      idVec <- as.vector(unique(snp.dic[, c("MIRNA_Name")]));
      mir.dic <- Query.miRNetDB(paste(sqlite.path, "mir2gene", sep=""), idVec, dataSet$org, "mir_id");
      mir.dic <- na.omit(mir.dic)

      # for network
        snp.edge <- na.omit(data.frame(
          Name1 = snp.dic[,"Mature_Name"],
          ID1 = snp.dic[,"Mature_Acc"],
          Name2 = snp.dic[,"rsid"],
          ID2 = snp.dic[,"rsid"],
          tableType = "snp2mir", # Data table type
          originalRowId = rownames(snp.dic), # Assuming row order is the original ID
          stringsAsFactors = FALSE
        ));
        mir.edge <- na.omit(data.frame(
          Name1 = mir.dic[,"mir_id"],
          ID1 = mir.dic[,"mir_acc"],
          Name2 = mir.dic[,"symbol"],
          ID2 = mir.dic[,"entrez"],
          tableType = "mir2gene", # Data table type
          originalRowId = rownames(mir.dic), # Assuming row order is the original ID
          stringsAsFactors = FALSE
        ));

      # table results
      snp.res <- na.omit(snp.dic[ , c("chr_pos", "rsid", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain")]);
      mir.res <- mir.dic[ , c("mir_id", "mir_acc", "symbol", "entrez", "experiment", "pmid", "tissue")];
      snp.res$Database <- rep("30302893", nrow(snp.res));
      colnames(snp.res) <- c("CHR_POS", "rsID", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain", "DataBase");
      colnames(mir.res) <- c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue");

      # snp2mirbs2mir
      snp2 <- na.omit(data.frame(name1 = snp.dic2[, idType], id1 = snp.dic2[, "rsid"], name2 = snp.dic2[, "symbol"], id2 =  snp.dic2[, "entrez"], stringsAsFactors = FALSE));
      hit.num2 <- nrow(snp2);
      current.msg <<- paste("A total of", hit.num, "SNPs were mapped to miRNAs &", hit.num2, "SNPs were mapped to miRNA-binding sites!");

      idVec2 <- as.vector(unique(snp.dic2[, c("entrez")]));
      mir.dic2 <- Query.miRNetDB(paste(sqlite.path, "mir2gene", sep=""), idVec2, dataSet$org, "entrez");
      mir.dic2 <- na.omit(mir.dic2);

      # for network
        snp.edge2 <- na.omit(data.frame(
          Name1 = snp.dic2[,"symbol"], 
          ID1 = snp.dic2[,"entrez"], 
          Name2 = snp.dic2[,"rsid"],
          ID2 = snp.dic2[,"rsid"],
          tableType = "snp2mirbs", # Data table type
          originalRowId = rownames(snp.dic2), # Assuming row order is the original ID
          stringsAsFactors = FALSE
        ));
        # Assuming additional queries might yield mir.dic2 or similar
        mir.edge2 <- na.omit(data.frame(
          Name1 = mir.dic2[,"symbol"], 
          ID1 = mir.dic2[,"entrez"], 
          Name2 = mir.dic2[,"mir_id"],
          ID2 = mir.dic2[,"mir_acc"],
          tableType = "mir2gene", # Specify the actual query type
          originalRowId = rownames(mir.dic2), # Assuming row order is the original ID
          stringsAsFactors = FALSE
        ));

      # table results
      snp.res2 <- na.omit(snp.dic2[ , c("chr_pos", "rsid", "transcript_id", "entrez", "symbol")]);
      mir.res2 <- mir.dic2[ , c("mir_id", "mir_acc", "symbol", "entrez", "experiment", "pmid", "tissue")];
      snp.res2$Literature <- rep("24163105", nrow(snp.res2));
      snp.res2$Database <- rep("PolymiRTS_3.0", nrow(snp.res2));
      colnames(snp.res2) <- c("CHR_POS", "rsID", "Transcript_ID", "Entrez", "Gene", "Literature", "Database");
      colnames(mir.res2) <- c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue");

      # rbind snp2mir2gene and snp2mirbs2mir for network

      # update the data
      gd.inx <- rownames(snp.mat) %in% unique(c(snp.edge[, "Name2"], snp.edge2[, "Name2"]));
      dataSet$mir.mapped <- snp.mat[gd.inx,,drop=F];

      dataSet$seeds <- c(snp.edge[, "Name2"], snp.edge2[, "Name2"]);

      dataSet$mirtarget <- c("gene");
      dataSet$mirtable <- c("snp2mir", "mir2gene", "snp2mirbs", "gene2mir");

        tf.dic <- Query.miRNetDB(paste(sqlite.path, "snp2tfbs", sep=""), snpidVec, dataSet$org, dataSet$idType);
        tf.dic <- na.omit(tf.dic)
        res <- na.omit(tf.dic[, c("chr_pos", "rsid", "entrez", "symbol", "name")]);
        res$Literature <- rep("27899579", nrow(res));
        res$Database <- rep("SNP2TFBS", nrow(res));
        res$tableType <- "snp2tfbs"; # Add tableType
        res$originalRowId <- 1:nrow(res); # Assuming row order as original ID
        colnames(res) <- c("CHR_POS", "rsID", "Entrez", "Symbol", "Name", "Literature", "Database", "TableType", "OriginalRowId");
        tf.vec = res["Entrez"]
        tf.res = res

        tf.dic <- Query.miRNetDB(paste(sqlite.path, "mir2gene", sep=""), tf.vec, dataSet$org, "entrez");
        tf.dic <- na.omit(tf.dic)

        res <- tf.dic[, c("mir_id", "mir_acc", "symbol", "entrez", "experiment", "pmid", "tissue")];
        res$Literature <- rep("NA", nrow(res));
        res$Database <- rep("NA", nrow(res));
        res$tableType <- "mir2gene"; # Add tableType
        res$originalRowId <- 1:nrow(res); # Assuming row order as original ID
        colnames(res) <- c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue", "TableType", "OriginalRowId");
        tf.res2 = res

        # Adjust based on the structure of tf.res and tf.res2 to include tableType and originalRowId
        tf.edge <- na.omit(data.frame(
          Name1=tf.res[, "Symbol"],
          ID1=tf.res[, "Entrez"],
          Name2=tf.res[, "Name"], # Adjust if "Name" is not the correct column
          ID2=tf.res[, "rsID"],
          tableType = "snp2tfbs",
          originalRowId = tf.res[, "OriginalRowId"],
          stringsAsFactors = FALSE
        ));

        tf.edge2 <- na.omit(data.frame(
          Name1=tf.res2[, "Gene"], # Adjust if "Gene" is not the correct column
          ID1=tf.res2[, "Entrez"],
          Name2=tf.res2[, "ID"],
          ID2=tf.res2[, "Accession"],
          tableType = "mir2gene",
          originalRowId = tf.res2[, "OriginalRowId"],
          stringsAsFactors = FALSE
        ));
      merge.edge <- rbind(snp.edge, mir.edge, snp.edge2, mir.edge2, tf.edge, tf.edge2);
      rownames(merge.edge) <- 1:nrow(merge.edge);
      dataSet$mir.res <- merge.edge; # for network
      fast.write.csv(merge.edge, file="mirnet_snp_mir_target.csv", row.names=FALSE);

      dataSet$snp2mir <- snp.res;
      dataSet$mir2gene <- mir.res;
      dataSet$snp2mirbs <- snp.res2;
      dataSet$gene2mir <- mir.res2;
      if(nrow(tf.res)>0){
        dataSet$snp2tfbs <- tf.res;
        dataSet$tf2mir <- tf.res2;
        dataSet$mirtable <- c(dataSet$mirtable, "snp2tfbs", "tf2mir");
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
      mir.dic <- na.omit(mir.dic)

      # for network
      snp.edge <- na.omit(data.frame(
          Name1 = snp.dic[,"Mature_Name"],
          ID1 = snp.dic[,"Mature_Acc"],
          Name2 = snp.dic[,"rsid"],
          ID2 = snp.dic[,"rsid"],
          tableType = "snp2mir", # Data table type
          originalRowId = rownames(snp.dic), # Assuming row order is the original ID
          stringsAsFactors = FALSE
        ));        
      mir.edge <- na.omit(data.frame(
          Name1 = mir.dic[,"mir_id"],
          ID1 = mir.dic[,"mir_acc"],
          Name2 = mir.dic[,"symbol"],
          ID2 = mir.dic[,"entrez"],
          tableType = "mir2gene", # Data table type
          originalRowId = rownames(mir.dic), # Assuming row order is the original ID
          stringsAsFactors = FALSE
        ));
      # table results
      snp.res <- na.omit(snp.dic[ , c("chr_pos", "rsid", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain")]);
      mir.res <- mir.dic[ , c("mir_id", "mir_acc", "symbol", "entrez", "experiment", "pmid", "tissue")];
      snp.res$Database <- rep("30302893", nrow(snp.res));
      colnames(snp.res) <- c("CHR_POS", "rsID", "Mature_Name", "Mature_Acc", "MIRNA_Name","MIRNA_Acc", "MIRNA_Domain", "DataBase");
      colnames(mir.res) <- c("ID", "Accession", "Gene", "Entrez", "Experiment", "Literature", "Tissue");

      # rbind  for network
      merge.edge <- rbind(snp.edge, mir.edge);
      rownames(merge.edge) <- 1:nrow(merge.edge);

      fast.write.csv(merge.edge, file="mirnet_snp_mir_target.csv", row.names=FALSE);

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
