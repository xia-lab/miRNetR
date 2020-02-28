##################################################
## R script for miRNet
## Description: GO/Pathway ORA 
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
#' @rdname LoadKEGGLib
#' @export 
LoadKEGGLib<-function(){
    kegg.rda <- paste(lib.path, dataSet$org, "/kegg_", dataSet$org, ".rda", sep="");
    
    print(paste("adding library:", kegg.rda));
    load(kegg.rda);

    current.setlink <- kegg$link;
    current.mset <- kegg$sets;
    set.ids<- names(current.mset); 
    names(set.ids) <- names(current.mset) <- kegg$term;

    current.setlink <<- current.setlink;
    current.setids <<- set.ids;
    current.geneset <<- current.mset;
    current.universe <<- unique(unlist(current.geneset));
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
#' @rdname LoadREACTOMELib
#' @export 
LoadREACTOMELib<-function(){

    reactome.rda <- paste(lib.path, dataSet$org, "/reactome_", dataSet$org, ".rda", sep="");

    print(paste("adding library:", reactome.rda));
    load(reactome.rda);

    current.mset <- reactome$sets;
    set.ids<- names(current.mset); 
    names(set.ids) <- names(current.mset) <- reactome$term;
    current.setlink <<- reactome$link;
    current.setids <<- set.ids;
    current.geneset <<- current.mset;
    current.universe <<- unique(unlist(current.geneset));
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param onto PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname LoadGOLib
#' @export 
LoadGOLib<-function(onto){

    go.rda <- paste(lib.path, dataSet$org, "/go_", tolower(onto), ".rda", sep="");

    print(paste("adding library:", go.rda));
    load(go.rda);

    if(tolower(onto) == "bp"){
        current.link <- go_bp$link;
        current.mset <- go_bp$sets;
        set.ids<- names(current.mset); 
        names(set.ids) <- names(current.mset) <- go_bp$term;
    }else if(tolower(onto) == "mf"){
        current.link <- go_mf$link;
        current.mset <- go_mf$sets;
        set.ids<- names(current.mset); 
        names(set.ids) <- names(current.mset) <- go_mf$term;
    }else{
        current.link <- go_cc$link;
        current.mset <- go_cc$sets;
        set.ids<- names(current.mset); 
        names(set.ids) <- names(current.mset) <- go_cc$term;
    }
    
    current.setlink <<- current.link;
    current.setids <<- set.ids;
    current.geneset <<- current.mset;
    current.universe <<- unique(unlist(current.geneset));
}

# loading miRNA tissue annotation library (human only)
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
#' @rdname LoadTissueLib
#' @export 
LoadTissueLib <- function(){
  tissue.rda <- paste(lib.path, "tissue.rda", sep="");
  load(tissue.rda);
  print(paste("adding library: ", tissue.rda));
  current.mset <- tissue;
  set.ids <- names(current.mset);
  names(set.ids) <- names(current.mset);
  current.setlink <<- "http://bioeng.swjtu.edu.cn/TSmiR";
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}

# loading miRNA functional annotation library Tam 2.0 (human only)
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
#' @rdname LoadFuncLib
#' @export 
LoadFuncLib <- function(){
  func.rda <- paste(lib.path, dataSet$org, "/tam_func.rda", sep="");
  load(func.rda);
  print(paste("adding library: ", func.rda));
  current.mset <- tam_func$sets;
  set.ids <- names(current.mset); 
  names(set.ids) <- names(current.mset) <- tam_func$term;
  current.setlink <<- "http://www.lirmed.com/tam2/";
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}

# loading miRNA hmdd disease annotation library Tam 2.0 (human only)
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
#' @rdname LoadHMDDLib
#' @export 
LoadHMDDLib <- function(){
  hmdd.rda <- paste(lib.path, dataSet$org, "/tam_hmdd.rda", sep="");
  load(hmdd.rda);
  print(paste("adding library: ", hmdd.rda));
  current.mset <- tam_hmdd$sets;
  set.ids <- names(current.mset); 
  names(set.ids) <- names(current.mset) <- tam_hmdd$term;
  current.setlink <<- "http://www.lirmed.com/tam2/";
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}

# loading miRNA cluster annotation library Tam 2.0 (human only)
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
#' @rdname LoadClusterLib
#' @export 
LoadClusterLib <- function(){
  cluster.rda <- paste(lib.path, dataSet$org, "/tam_cluster.rda", sep="");
  load(cluster.rda);
  print(paste("adding library: ", cluster.rda));
  current.mset <- tam_cluster$sets;
  set.ids <- names(current.mset); 
  names(set.ids) <- names(current.mset) <- tam_cluster$term;
  current.setlink <<- "http://www.lirmed.com/tam2/";
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}

# loading miRNA TF annotation library Tam 2.0 (human only)
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
#' @rdname LoadTFLib
#' @export 
LoadTFLib <- function(){
  tf.rda <- paste(lib.path, dataSet$org, "/tam_tf.rda", sep="");
  load(tf.rda);
  print(paste("adding library: ", tf.rda));
  current.mset <- tam_tf$sets;
  set.ids <- names(current.mset); 
  names(set.ids) <- names(current.mset) <- tam_tf$term;
  current.setlink <<- "http://www.lirmed.com/tam2/";
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}

# loading mirfamily library accroding to the species. The names for set.ids are the same as set.ids.
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
#' @rdname LoadmiRFamLib
#' @export 
LoadmiRFamLib <- function(){
  mirfamily.rda <- paste(lib.path, "mirfamily.rda", sep="");
  load(mirfamily.rda);
  print(paste("adding library: ", mirfamily.rda));
  
  if(dataSet$org == "hsa"){
    current.mset <- mirfam$hsa;
    set.ids <- names(current.mset);
    names(set.ids) <- names(current.mset);
    current.setlink <<- "http://www.mirbase.org";
    current.setids <<- set.ids;
    current.geneset <<- current.mset;
    current.universe <<- unique(unlist(current.geneset));
  
    } else if(dataSet$org == "mmu"){
    current.mset <- mirfam$mmu;
    set.ids <- names(current.mset);
    names(set.ids) <- names(current.mset);
    current.setlink <<- "http://www.mirbase.org";
    current.setids <<- set.ids;
    current.geneset <<- current.mset;
    current.universe <<- unique(unlist(current.geneset));
    
    }else if(dataSet$org == "rno"){
    current.mset <- mirfam$rno;
    set.ids <- names(current.mset);
    names(set.ids) <- names(current.mset);
    current.setlink <<- "http://www.mirbase.org";
    current.setids <<- set.ids;
    current.geneset <<- current.mset;
    current.universe <<- unique(unlist(current.geneset));
    
    }else if(dataSet$org == "gga"){
    current.mset <- mirfam$gga;
    set.ids <- names(current.mset);
    names(set.ids) <- names(current.mset);
    current.setlink <<- "http://www.mirbase.org";
    current.setids <<- set.ids;
    current.geneset <<- current.mset;
    current.universe <<- unique(unlist(current.geneset));
    
    }else if(dataSet$org == "bta"){
    current.mset <- mirfam$bta;
    set.ids <- names(current.mset);
    names(set.ids) <- names(current.mset);
    current.setlink <<- "http://www.mirbase.org";
    current.setids <<- set.ids;
    current.geneset <<- current.mset;
    current.universe <<- unique(unlist(current.geneset));
    
    }else if(dataSet$org == "dre"){
    current.mset <- mirfam$dre;
    set.ids <- names(current.mset);
    names(set.ids) <- names(current.mset);
    current.setlink <<- "http://www.mirbase.org";
    current.setids <<- set.ids;
    current.geneset <<- current.mset;
    current.universe <<- unique(unlist(current.geneset));
    
    }else if(dataSet$org == "dme"){
    current.mset <- mirfam$dme;
    set.ids <- names(current.mset);
    names(set.ids) <- names(current.mset);
    current.setlink <<- "http://www.mirbase.org";
    current.setids <<- set.ids;
    current.geneset <<- current.mset;
    current.universe <<- unique(unlist(current.geneset));
    
    }else {
    current.mset <- mirfam$cel;
    set.ids <- names(current.mset);
    names(set.ids) <- names(current.mset);
    current.setlink <<- "http://www.mirbase.org";
    current.setids <<- set.ids;
    current.geneset <<- current.mset;
    current.universe <<- unique(unlist(current.geneset));
    
    }
}

# note: hit.query, resTable must synchronize
# Due to bias in mirRNA targets, for mir2gene, the p value will be tested 
# using 1000 permutation. This should NOT be used for gene2mir
# If users click more than once using the empirical approach 
# We can take advantage of the accumulating the permutation result to increase the significance 

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param fun.type PARAM_DESCRIPTION
#' @param file.nm PARAM_DESCRIPTION
#' @param IDs PARAM_DESCRIPTION
#' @param algo PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformMirTargetEnrichAnalysis
#' @export 
PerformMirTargetEnrichAnalysis <- function(fun.type, file.nm, IDs, algo){

    perm.num <- 1000;

    # prepare lib
    if(tolower(fun.type) == 'kegg'){ 
        LoadKEGGLib();
    }else if(tolower(fun.type) == 'reactome'){ 
        LoadREACTOMELib();
    }else if(tolower(fun.type) == 'mirfamily'){ # when user choose to perform miRNA family enrichment analysis.
        LoadmiRFamLib();
    }else if(tolower(fun.type) == 'tissue'){
        LoadTissueLib();
    }else if(tolower(fun.type) == 'func'){
        LoadFuncLib();
    }else if(tolower(fun.type) == 'hmdd'){
        LoadHMDDLib();
    }else if(tolower(fun.type) == 'cluster'){
        LoadClusterLib();
    }else if(tolower(fun.type) == 'tf'){
        LoadTFLib();
    }else{ # GO
        LoadGOLib(fun.type);
    }

    mirnet.type <- dataSet$mirnet;

    # prepare query, current.mirnet may be subset of all networks
    nodeList <- get.data.frame(current.mirnet, "vertices");
    
    if(data.type == "xeno.mir"){
        if(tolower(fun.type) == 'mirfamily'){
            hit.inx <- dataSet$mir.filtered$miRNA %in% nodeList[,1];
  
            mir.query <- unique(dataSet$mir.filtered$miRNA[hit.inx]);
            my.data <- unique(dataSet$mir.filtered[hit.inx,c("Accession", "miRNA")]); # The original dataset contains miRNA Accession number, you can consider it as entrez id.
  
            ora.vec <- my.data[hit.inx,"Accession"];
            sybls <- my.data[hit.inx,"miRNA"];
            names(ora.vec) <- sybls;
        } else{
            hit.inx <- dataSet$mir.filtered$Gene %in% nodeList[,1];

            mir.query <- unique(dataSet$mir.filtered$miRNA[hit.inx]);
            my.data <- unique(dataSet$mir.filtered[hit.inx,c("Entrez", "Gene")]);
       
            ora.vec <- my.data[hit.inx,"Entrez"];
            sybls <- my.data[hit.inx,"Gene"];
            names(ora.vec) <- sybls;
        }
    }else{
        if (tolower(fun.type) == 'mirfamily'){
            hit.inx <- dataSet$mir.filtered$ID %in% nodeList[,1];
            mir.query <- unique(dataSet$mir.filtered$ID[hit.inx]);
            my.data <- unique(dataSet$mir.filtered[hit.inx, c("Accession", "ID")]); # The original dataset contains miRNA Accession number, you can consider it as entrez id.
  
            ora.vec <- my.data[hit.inx, "Accession"];
            sybls <- my.data[hit.inx, "ID"];
            names(ora.vec) <- sybls;
        } else if (tolower(fun.type) == 'tissue'){
            hit.inx <- dataSet$mir.filtered$ID %in% nodeList[,1];
            mir.query <- unique(dataSet$mir.filtered$ID[hit.inx]);
            my.data <- unique(dataSet$mir.filtered[hit.inx, c("Accession", "ID")]);# The original dataset contains miRNA Accession number, you can consider it as entrez id.
       
            ora.vec <- my.data[hit.inx, "Accession"];
            sybls <- my.data[hit.inx, "ID"];
            names(ora.vec) <- sybls;
        } else if (tolower(fun.type) == 'func'){
          hit.inx <- dataSet$mir.filtered$ID %in% nodeList[,1];
          mir.query <- unique(dataSet$mir.filtered$ID[hit.inx]);
          my.data <- unique(dataSet$mir.filtered[hit.inx, c("Accession", "ID")]);# The original dataset contains miRNA Accession number, you can consider it as entrez id.
          ora.vec <- my.data[hit.inx, "Accession"];
          sybls <- my.data[hit.inx, "ID"];
          names(ora.vec) <- sybls;
        } else if (tolower(fun.type) == 'hmdd'){
          hit.inx <- dataSet$mir.filtered$ID %in% nodeList[,1];
          mir.query <- unique(dataSet$mir.filtered$ID[hit.inx]);
          my.data <- unique(dataSet$mir.filtered[hit.inx, c("Accession", "ID")]);# The original dataset contains miRNA Accession number, you can consider it as entrez id.
          ora.vec <- my.data[hit.inx, "Accession"];
          sybls <- my.data[hit.inx, "ID"];
          names(ora.vec) <- sybls;
        } else if (tolower(fun.type) == 'cluster'){
          hit.inx <- dataSet$mir.filtered$ID %in% nodeList[,1];
          mir.query <- unique(dataSet$mir.filtered$ID[hit.inx]);
          my.data <- unique(dataSet$mir.filtered[hit.inx, c("Accession", "ID")]);# The original dataset contains miRNA Accession number, you can consider it as entrez id.
          ora.vec <- my.data[hit.inx, "Accession"];
          sybls <- my.data[hit.inx, "ID"];
          names(ora.vec) <- sybls;
        } else if (tolower(fun.type) == 'tf'){
          hit.inx <- dataSet$mir.filtered$ID %in% nodeList[,1];
          mir.query <- unique(dataSet$mir.filtered$ID[hit.inx]);
          my.data <- unique(dataSet$mir.filtered[hit.inx, c("Accession", "ID")]);# The original dataset contains miRNA Accession number, you can consider it as entrez id.
          ora.vec <- my.data[hit.inx, "Accession"];
          sybls <- my.data[hit.inx, "ID"];
          names(ora.vec) <- sybls;
        }else {
            hit.inx <- dataSet$mir.filtered$Gene %in% nodeList[, 1];

            mir.query <- unique(dataSet$mir.filtered$ID[hit.inx]);
            my.data <- unique(dataSet$mir.filtered[hit.inx, c("Entrez", "Gene")]);
            ora.vec <- my.data[hit.inx, "Entrez"];
            sybls <- my.data[hit.inx, "Gene"];
            names(ora.vec) <- sybls;
        }
    }
    q.vec <-  unlist(strsplit(IDs, "; "));
    ora.vec <- ora.vec[q.vec];
    
    ora.vec <- ora.vec[!is.na(ora.vec)];
    ora.nms <- names(ora.vec);

    # prepare for the result table
    set.size <- length(current.geneset);

    res.mat<-matrix(0, nrow=set.size, ncol=4);
    colnames(res.mat) <- c("Total", "Expected", "Hits", "Pval");
    rownames(res.mat) <- names(current.geneset);

    # not all query genes can be used, need to cut query to only the universe covered 
    hits.inx <- ora.vec %in% current.universe;
    ora.vec <- ora.vec[hits.inx];
    ora.nms <- ora.nms[hits.inx];

    q.size<-length(ora.vec);

    # get the matched query for each pathway
    hits.query <- lapply(current.geneset, 
        function(x) {
            ora.nms[ora.vec %in% x];
        }
    );
    hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=FALSE);
    names(hits.query) <- names(current.geneset);

    # total unique gene number
    uniq.count <- length(current.universe);
    
    # unique gene count in each pathway
    set.size <- unlist(lapply(current.geneset, length));

    res.mat[,1]<-set.size;
    res.mat[,2]<-q.size*(set.size/uniq.count);
    res.mat[,3]<-hit.num;

    if(algo == 'emp'){
        # empirical sampling
        # do stepped permutations 
        # assume 1000 ==> 200
        # 300, 500 to remove those genesets that are already >20% , to save computing time
        library(fastmatch); # ~20% faster
        perm.out <- matrix(0, nrow=length(current.geneset), ncol = perm.num);
        if(data.type == "xeno.mir"){
            myRandQs <- GetRandomXenoMirTargetGenes(length(mir.query), perm.num);
        }else{
            myRandQs <- GetRandomMirTargetGenes(length(mir.query), perm.num);
        }
        for(i in 1:perm.num){
            perm.out[, i]<- sapply(current.geneset, function(x) {sum(fmatch(myRandQs[[i]], x, nomatch = 0L) > 0L)}); 
        }
    
        #empirical p from permutation - percentage of number large than original
        hmat <- perm.out - hit.num > 0;
        # now, see if we can combine with previous previous permutation with the same query size
        if(is.null(dataSet$perm.res)){ 
            dataSet$perm.res <<- list();
        }
        perm.nm <- paste("Q", length(mir.query), sep="");
        if(perm.nm %in% names(dataSet$perm.res)){ # same query size and we can combine
            hmat <- cbind(dataSet$perm.res[[perm.nm]], hmat);
        }
        dataSet$perm.res[[perm.nm]] <<- hmat;
        perm.pvals <- apply(hmat, 1, sum)/ncol(hmat);

        res.mat[,4] <- perm.pvals; 
    }else{
        # standard hypergeometric tests use lower.tail = F for P(X>x)
        raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
        #res.mat[,4]<- raw.pvals;
    
        fdr.pvals <- p.adjust(raw.pvals, "fdr");
        res.mat[,4] <- fdr.pvals;
    }

    res.mat <- res.mat[hit.num>0,,drop = F];
    hits.query <- hits.query[hit.num>0];

    if(nrow(res.mat)> 1){
        # order by p value
        ord.inx<-order(res.mat[,4]);
        res.mat <- signif(res.mat[ord.inx,],3);
        hits.query <- hits.query[ord.inx];
    }

    #get gene symbols
    resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
    if(nrow(resTable) == 0){
        current.msg <<- "No hits found for your query!";
        print(current.msg);
        return(0);
    }
    current.msg <<- "Functional enrichment analysis was completed";

    # write json
    fun.anot <- hits.query; 

    hit.num <- resTable[,4]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
    fun.ids <- as.vector(current.setids[names(fun.anot)]); 
    if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };

    pval <- resTable[,5]; if(length(pval) ==1) { pval <- matrix(pval) };
    if(algo == "emp"){
        hit.inx <- pval == 0;
        pval[hit.inx] <- paste("<", 1/perm.num);
    }
    json.res <- list(
                    fun.anot = fun.anot,
                    fun.ids = fun.ids,
                    pval = pval,
                    hit.num = hit.num
        );
     json.mat <- toJSON(json.res, .na='null');
     json.nm <- paste(file.nm, ".json", sep="");
     
     sink(json.nm)
     cat(json.mat);
     sink();

    # write csv
    # csv.nm <- paste(file.nm, ".csv", sep="");
    write.csv(resTable, file="mirnet_enrichment.csv", row.names=F);
    return(1);
}

# return a list of gene targets from the same size but randomly selected mirs
# qSize is the query mir vec size
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param qSize PARAM_DESCRIPTION
#' @param perm.num PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetRandomMirTargetGenes
#' @export 
GetRandomMirTargetGenes <- function(qSize, perm.num){

    mir.db <- dbConnect(SQLite(), paste(sqlite.path, "mir2gene", sep=""));
    statement <- paste("SELECT mir_id,entrez FROM ",dataSet$org, sep="");
    mir.dic <- .query.sqlite(mir.db, statement); 

    # now get the unique mirs
    mirs <- unique(mir.dic[,1]);
    
    # sample with replacement
    res <- vector(length=perm.num, mode="list");
    for(i in 1:perm.num){
        q.mir <- sample(mirs, qSize, replace = TRUE);
        hit.inx <- mir.dic[,1] %fin% q.mir;
        ora.vec <- unique(mir.dic[hit.inx, 2]);

        # filter ora.vec based on current universe
        res[[i]] <- ora.vec[ora.vec %fin% current.universe];
    }
    return(res);
}

# return a list of gene targets from the same size but randomly selected mirs
# qSize is the query mir vec size
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param qSize PARAM_DESCRIPTION
#' @param perm.num PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetRandomXenoMirTargetGenes
#' @export 
GetRandomXenoMirTargetGenes <- function(qSize, perm.num){
    mir.db <- dbConnect(SQLite(), paste(lib.path, "xenomirnet", sep=""));
    statement <- paste("SELECT exo_mirna,entrez FROM ",dataSet$org, sep="");
    mir.dic <- .query.sqlite(mir.db, statement);

    # now get the unique mirs
    mirs <- unique(mir.dic[,1]);
    
    # sample without replacement
    res <- vector(length=perm.num, mode="list");
    for(i in 1:perm.num){
        q.mir <- sample(mirs, qSize, replace = FALSE);
        hit.inx <- mir.dic[,1] %fin% q.mir;
        ora.vec <- unique(mir.dic[hit.inx, 2]);

        # filter ora.vec based on current universe
        res[[i]] <- ora.vec[ora.vec %fin% current.universe];
    }
    return(res);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nms PARAM_DESCRIPTION
#' @param operation PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname CalculateMirTargetSet
#' @export 
CalculateMirTargetSet <- function(nms, operation){
    nms <- strsplit(nms, ";")[[1]];

    valid.inx <- nms %in% dataSet$mir.filtered$ID;
    nms <- nms[valid.inx];
    if(length(nms) == 0){
        print("No valid mir found!");
        return("error||No valid miRNA IDs were selected!");
    }

    hit.inx <- dataSet$mir.filtered$ID == nms[1];
    targets <- dataSet$mir.filtered[hit.inx, 3]; # targets at third column
    for(i in 2:length(nms)){
        hit.inx <- dataSet$mir.filtered$ID == nms[i];
        if(operation == "intersect"){
            targets <- intersect(targets, dataSet$mir.filtered[hit.inx,3]);
        }else if(operation == "union"){
            targets <- union(targets, dataSet$mir.filtered[hit.inx,3]);
        }
    }

    # include original queries
    return(paste(unique(c(nms, targets)), collapse="||"));
}
