##################################################
## R script for miRNet
## Description: expression data I/O, normalization and differential expression analysis 
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# read tab delimited file
# stored in dataSet list object
# can have many classes, stored in meta.info
# type: array, count, qpcr
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataName PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ReadTabExpressData
#' @export 
ReadTabExpressData <- function(dataName) {

    dataSet <- ReadTabData(dataName);

    # rename data to data.orig
    int.mat <- dataSet$data;
    dataSet$cls <- dataSet$meta.info[,1];
    dataSet$data <- NULL;
    dataSet$listData <- FALSE;

    msg <- paste("a total of ", ncol(int.mat), " samples and ", nrow(int.mat), " features were found. ");

    # remove NA, null
    row.nas <- apply(is.na(int.mat)|is.null(int.mat), 1, sum);
    good.inx<- row.nas/ncol(int.mat) < 0.5;
    if(sum(!good.inx) > 0){
        int.mat <- int.mat[good.inx,];
        msg <- c(msg, paste("removed ", sum(!good.inx), " features with over 50% missing values"));
    }
    # remove constant values
    filter.val <- apply(int.mat, 1, IQR, na.rm=T);
    good.inx2 <- filter.val > 0;
    if(sum(!good.inx2) > 0){
        int.mat <- int.mat[good.inx2,];
        msg <- c(msg, paste("removed ", sum(!good.inx2), " features with constant values"));
    }

    if(nrow(int.mat) > 2000){
        filter.val <- filter.val[good.inx2];
        rk <- rank(-filter.val, ties.method='random');

        var.num <- nrow(int.mat);
        kept.num <- 0.95*var.num;
        int.mat <- int.mat[rk < kept.num, ];
        # msg <- c(msg, paste("removed 5% features with near-constant values"));
    }

    minVal <- min(int.mat, na.rm=T);
    na.inx <- is.na(int.mat);
    if(sum(na.inx) > 0){
        int.mat[na.inx] <- minVal/2;
        # msg <- c(msg, "the remaining", sum(na.inx), "missing variables were replaced with data min");
    }
    current.msg <<- paste(msg, collapse="; ");
    saveRDS(int.mat, file="data.proc");
    dataSet <<- dataSet;
    if(.on.public.web){
      return (1);  
    }else{
      return (current.msg);
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
#' @rdname SetupMirExpressData
#' @export 
SetupMirExpressData <- function(){
    idType <- dataSet$id.current;
    mydata <- data.matrix(dataSet$sig.mat[,"max.logFC",drop=FALSE]);
    if(idType == "mir_id"){ # note, in the mirnet database, all mir ids are lower case! miR=>mir
        mir.vec <- rownames(mydata);
        rownames(mydata) <- tolower(mir.vec);
    }
    dataSet$idType <- idType;
    dataSet$mir.orig <- mydata;
    dataSet$data[["gene"]] <- mydata;
    dataSet$id.types[["gene"]] <- idType;
    dataSet<<- dataSet;
    if(.on.public.web){
      return (nrow(mydata));
    }else{
      return("Setup complete!")
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
#' @rdname GetClassInfo
#' @export 
GetClassInfo <- function(){
    return(levels(dataSet$cls));
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
#' @rdname GetAnotNames
#' @export 
GetAnotNames<-function(){
    return(rownames(dataSet$data.anot));
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
    if(.on.public.web){
      return(matched.len);
    }else{
      return("Annotation was performed successfully!")
    }
}

# update result based on new cutoff
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param p.lvl PARAM_DESCRIPTION
#' @param fc.lvl PARAM_DESCRIPTION
#' @param direction PARAM_DESCRIPTION
#' @param update PARAM_DESCRIPTION, Default: T
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetSigGenes
#' @export 
GetSigGenes<-function(p.lvl, fc.lvl, direction, update=T){

    resTable <- readRDS("resTable");
    if(update){
        current.msg <<- "";
    }
    # select based on p
    hit.inx.p <- resTable$adj.P.Val <= p.lvl; 

    resTable<-resTable[hit.inx.p,,drop=F];
    if(nrow(resTable) == 0){
        current.msg <<- paste(current.msg, "No significant genes were identified using the given design and cutoff."); 
        return();
    }

    # now rank by logFC, note, the logFC for each comparisons 
    # are returned in resTable before the AveExpr columns 
    # for two-class, only one column, multiple columns can be involved
    # for > comparisons - in this case, use the largest logFC among all comparisons
    fc.lvl <- abs(fc.lvl);
    #if(fc.lvl > 0){ # further filter by logFC
        if(direction == "both"){
            hit.inx.fc <- abs(resTable$max.logFC) >= fc.lvl;
        }else if(direction == "up"){
            hit.inx.fc <- resTable$max.logFC >= fc.lvl;
        }else{
            hit.inx.fc <- resTable$max.logFC <= -fc.lvl;
        }
        resTable<-resTable[hit.inx.fc,,drop=F];
        if(nrow(resTable) == 0){
            current.msg <<- paste(current.msg, "No significant genes were identified using the given design and cutoff."); 
            return();
        }
    #}

    write.csv(signif(resTable,5), file="mirnet_siggenes.csv");
        
    de.Num <- nrow(resTable);
    current.msg <<- paste(current.msg, "A total of", de.Num, "significant genes were identified!"); 
 
    gene.anot <- NULL;
    # display at most 5000 genes for the server (two main reasons)
    # 1) should not have more 22% (human: 23000) DE of all genes (biological)
    # 2) IE canvas can display no more than 6800 pixels (computational)
    if(nrow(resTable) > 5000){
        resTable <- resTable[1:5000,];
        gene.anot <- gene.anot[1:5000,];
        current.msg <<- paste(current.msg, " Due to computational constraints, only top 5000 genes will be used. ", collapse="\n");
    }

    # may need to update data, class and meta.info
    data <- dataSet$data.norm;
    cls <- dataSet$cls; 
    meta.info <- dataSet$meta.info;
    grp.nms <- levels(cls);

    hit.inx <- cls %in% grp.nms;
    if(sum(hit.inx) < length(hit.inx)){
         current.msg <<- paste(current.msg, "Only groups selected for comparisons: ", paste(grp.nms, collapse=", "), "are included.");
         cls <- factor(cls[hit.inx]);
         cls.lvls <- levels(cls);
         data <- data[,hit.inx];
         meta.info <- dataSet$meta.info[hit.inx,];
    }

    dataSet$sig.mat <- resTable;
    dataSet$sig.genes.anot <- gene.anot;
    dataSet$cls.stat <- cls;
    dataSet$meta.stat <- meta.info;

    dataSet <<- dataSet;

    # return DE num
    if(.on.public.web){
      return(de.Num);  
    }else{
      return(current.msg)
    }
}

# note, here also update data type array/count
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param norm.opt PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformArrayDataNormalization
#' @export 
PerformArrayDataNormalization <- function(norm.opt){

    data <- dataSet$data.anot;
    row.nms <- rownames(data);
    col.nms <- colnames(data);

    msg <- NULL;
    if(norm.opt=="log"){
        min.val <- min(data[data>0], na.rm=T)/10;
        data <- log2((data + sqrt(data^2 + min.val^2))/2);
        msg <- "Log2 transformation";
    }else if(norm.opt=="quant"){
        library('preprocessCore');
        data <- normalize.quantiles(data, copy=FALSE);
        msg <- "Quantile normalization";
    }else if(norm.opt == "combine"){
        min.val <- min(data[data>0], na.rm=T)/10;
        data <- log2((data + sqrt(data^2 + min.val^2))/2);
        library('preprocessCore');
        data <- normalize.quantiles(data, copy=FALSE);
        msg <- "Log2 transformation followed by normalization";
    }else{
        msg <-"No log normalization was performed.";
        print(msg);
    }

    rownames(data) <- row.nms;
    colnames(data) <- col.nms;
    dataSet$data.norm <- data;

    current.msg <<- msg;
    dataSet <<- dataSet;
    if(.on.public.web){
      return(1);
    }else{
      return("Normalization was performed successfully!")
    }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param target.grp PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformLimma
#' @export 
PerformLimma<-function(target.grp){

    myargs <- list();
    cls <- dataSet$cls; 
    design <- model.matrix(~ 0 + cls) # no intercept
    colnames(design) <- levels(cls);

    grp.nms <- strsplit(target.grp, " vs. ")[[1]];
    myargs[[1]] <- paste(grp.nms, collapse="-");
    filename = paste("SigGene_", paste(grp.nms, collapse="_vs_"), sep="");

    library(limma);
    myargs[["levels"]] <- design;
    contrast.matrix <- do.call(makeContrasts, myargs);

    fit = lmFit(dataSet$data.norm, design);

    # sanity check
    if(!is.fullrank(design)){
        current.msg <<- paste("This metadata combination is not full rank! Please use other combination."); 
        return(0);
    }

    df.residual <- fit$df.residual;
    if (all(df.residual == 0)){
        current.msg <<- paste("There is not enough replicates in each group (no residual degrees of freedom)!"); 
        return(0);
    }
    fit2 <- contrasts.fit(fit, contrast.matrix);
    fit2 <- eBayes(fit2);
    topFeatures <- topTable(fit2, number=Inf, adjust.method="fdr");

    # add a common column for FC selection
    # this is to prevent no logFC for multi grps or for counts data

    hit.inx <- which(colnames(topFeatures) == "AveExpr");
    maxFC.inx <- hit.inx - 1; # not sure if this is also true for edgeR
    logfc.mat <- topFeatures[,1:maxFC.inx, drop=F];

    # extract the max FC together with direction
    pos.vec <- apply(abs(logfc.mat), 1, which.max);
    pos.mat <- cbind(1:length(pos.vec), pos.vec);
    max.logFC <- logfc.mat[pos.mat]; 

    topFeatures <- cbind(topFeatures, max.logFC = max.logFC);

    dataSet$filename <- filename;
    dataSet <<- dataSet;
    saveRDS(topFeatures, file="resTable");
    if(.on.public.web){
      return(1);
    }else{
      return("Differential Analysis was performed successfully!")
    }
}

###########################
## for RNAseq data
##########################
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param norm.opt PARAM_DESCRIPTION
#' @param disp.opt PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformCountDataNormalization
#' @export 
PerformCountDataNormalization <- function(norm.opt, disp.opt){

    msg <- NULL;
    cls <- dataSet$cls; 
    design <- model.matrix(~ 0 + cls) # no intercept
    colnames(design) <- levels(cls);

    library(edgeR);
    y <- DGEList(counts=dataSet$data.anot, group=dataSet$cls);
    if(norm.opt=="tmm"){
        y <- calcNormFactors(y)
    }else{
        msg <-"No log normalization was performed.";
        print(msg);
    }

    y <- estimateGLMCommonDisp(y, design, verbose=FALSE);

    if(disp.opt=="tagwise"){
        y <- estimateGLMTrendedDisp(y, design);
        y <- estimateGLMTagwiseDisp(y, design, trend=TRUE);
    }
    saveRDS(y, file="edger.y");
    dataSet$design <- design;
    current.msg <<- msg;
    dataSet <<- dataSet;
    if(.on.public.web){
      return(1);
    }else{
      return("Normalization was performed successfully!")
    }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param target.grp PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformEdgeR
#' @export 
PerformEdgeR<-function(target.grp){

    myargs <- list();
    grp.nms <- strsplit(target.grp, " vs. ")[[1]];
    myargs[[1]] <- paste(grp.nms, collapse="-");
    filename = paste("SigGene_", paste(grp.nms, collapse="_vs_"), sep="");

    myargs[["levels"]] <- dataSet$design;
    contrast.matrix <- do.call(makeContrasts, myargs);

    edger.y <- readRDS("edger.y");
    fit <- glmFit(edger.y, dataSet$design);
    lrt <- glmLRT(fit, contrast=contrast.matrix);
    topFeatures<-topTags(lrt,n=Inf)$table;

    # need to change the FDR to adj.P.Val same as limma
    nms <- colnames(topFeatures);
    nms[which(nms == 'FDR')] <- 'adj.P.Val';
    colnames(topFeatures) <- nms; 

    # add a common column for FC selection
    # this is to prevent no logFC for multi grps or for counts data

    hit.inx <- which(colnames(topFeatures) == "logCPM");
    maxFC.inx <- hit.inx - 1; # not sure if this is also true for edgeR
    logfc.mat <- topFeatures[,1:maxFC.inx, drop=F];

    # extract the max FC together with direction
    pos.vec <- apply(abs(logfc.mat), 1, which.max);
    pos.mat <- cbind(1:length(pos.vec), pos.vec);
    max.logFC <- logfc.mat[pos.mat]; 

    topFeatures <- cbind(topFeatures, max.logFC = max.logFC);

    dataSet$filename <- filename;
    dataSet <<- dataSet;
    saveRDS(topFeatures, file="resTable");
    if(.on.public.web){
      return(1);
    }else{
      return("Differential Analysis was performed successfully!")
    }
}

###########################
## for QPCR data
##########################

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param norm.opt PARAM_DESCRIPTION, Default: 'quantile'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformQpcrDataNormalization
#' @export 
PerformQpcrDataNormalization <- function(norm.opt = "quantile") {
    
    data <- dataSet$data.anot;

    library('HTqPCR');

    # need to create a qPCRset with only info needed 
    obj <- new("qPCRset", exprs=data);
    featureNames(obj) <- rownames(data);
    
    if(norm.opt =="deltaCt"){
        if(!exists('delta.genes')){
            current.msg <<- "Cannot find endogenous control genes!";
            print(current.msg);
            return (0);
        }
        norm.obj <- normalizeCtData(obj, norm="deltaCt", deltaCt.genes=delta.genes)
    }else{
        norm.obj <- normalizeCtData(obj, norm=norm.opt)
    }

    dataSet$data.norm <<- exprs(norm.obj);
    if(.on.public.web){
      return(1);
    }else{
      return("Normalization was performed successfully!")
    }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param target.grp PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformHTqPCR
#' @export 
PerformHTqPCR<-function(target.grp, method){

    if(method == "limma"){
        return(PerformLimma(target.grp));
    }

    # need to subselect groups and data
    grp.nms <- strsplit(target.grp, " vs. ")[[1]];
    my.inx <- dataSet$cls %in% grp.nms;
    
    my.cls <- factor(dataSet$cls[my.inx]);
    my.dat <- dataSet$data.norm[, my.inx];

    nonpar <- TRUE;
    if(method == "ttest"){
        nonpar <- FALSE;
    }
    p.value <- GetTtestP(my.dat, my.cls, grp.nms[1], grp.nms[2], FALSE, TRUE, nonpar);
    fdr.p <- p.adjust(p.value, "fdr");
    fc <-GetFC(my.dat, my.cls, grp.nms[1], grp.nms[2]);

    resTable <- data.frame(p.value=p.value, adj.P.Val = fdr.p, ddCt = fc);
    topFeatures <- cbind(resTable, max.logFC = fc);
    rownames(topFeatures)<-rownames(my.dat);
    
    inx<-order(p.value);
    topFeatures<-topFeatures[inx,];

    filename = paste("SigGene_", paste(grp.nms, collapse="_vs_"), sep="");
    dataSet$filename <- filename;
    dataSet <<- dataSet;
    saveRDS(topFeatures, file="resTable");
    if(.on.public.web){
      return(1);
    }else{
      return("Differential Analysis was performed successfully!")
    }
}

# utility method to get p values
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param my.dat PARAM_DESCRIPTION
#' @param my.cls PARAM_DESCRIPTION
#' @param grp1 PARAM_DESCRIPTION
#' @param grp2 PARAM_DESCRIPTION
#' @param paired PARAM_DESCRIPTION, Default: FALSE
#' @param equal.var PARAM_DESCRIPTION, Default: TRUE
#' @param nonpar PARAM_DESCRIPTION, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetTtestP
#' @export 
GetTtestP <- function(my.dat, my.cls, grp1, grp2, paired=FALSE, equal.var=TRUE, nonpar=F){

    inx1 <- which(my.cls==grp1);
    inx2 <- which(my.cls==grp2);

    if(nonpar){
         p.value <- apply(as.matrix(my.dat), 1, function(x) {
                tmp <- try(wilcox.test(x[inx1], x[inx2], paired = paired));
                if(class(tmp) == "try-error") {
                    return(NA);
                }else{
                    return(tmp$p.value);
                }
        })
    }else{
        if(nrow(my.dat) < 1000){
            p.value <- apply(as.matrix(my.dat), 1, function(x) {
                tmp <- try(t.test(x[inx1], x[inx2], paired = paired, var.equal = equal.var));
                if(class(tmp) == "try-error") {
                    return(NA);
                }else{
                    return(tmp$p.value);
                }
            })
        }else{ # use fast version
            library(genefilter);
            p.value <- try(rowttests(t(as.matrix(my.dat)), my.cls)$p.value);
            if(class(p.value) == "try-error") {
               p.value <- NA;
            }
        }
    }
    return(p.value);
}

# utility method to calculate FC
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param my.dat PARAM_DESCRIPTION
#' @param my.cls PARAM_DESCRIPTION
#' @param grp1 PARAM_DESCRIPTION
#' @param grp2 PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetFC
#' @export 
GetFC <- function(my.dat, my.cls, grp1, grp2){
    m1 <- rowMeans(my.dat[, which(my.cls==grp1)]);
    m2 <- rowMeans(my.dat[, which(my.cls==grp2)]);

    # create a named matrix of sig vars for display
    fc <- signif (m1-m2, 5);    
    return(fc);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param imgNm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PlotDataOverview
#' @export 
PlotDataOverview<-function(imgNm){
    dat <- dataSet$data.anot;
    library('lattice');
    subgene=10000;
    if (nrow(dat)>subgene) {
        set.seed(28051968);
        sg  = sample(nrow(dat), subgene)
        Mss = dat[sg,,drop=FALSE]
    } else {
        Mss = dat
    }

    subsmpl=100;
    if (ncol(Mss)>subsmpl) {
        set.seed(28051968);
        ss  = sample(ncol(Mss), subsmpl)
        Mss = Mss[,ss,drop=FALSE]
    } else {
        Mss = Mss
    }

    sample_id = rep(seq_len(ncol(Mss)), each = nrow(Mss));
    values  = as.numeric(Mss)
    formula = sample_id ~ values

  Cairo(file=imgNm, width=460, height=420, type="png", bg="white");
    box = bwplot(formula, groups = sample_id, layout = c(1,1), as.table = TRUE,
        strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
        horizontal = TRUE,
        pch = "|",  col = "black", do.out = FALSE, box.ratio = 2,
        xlab = "", ylab = "Samples",
        fill = "#1c61b6AA",
        panel = panel.superpose,
        scales = list(x=list(relation="free"), y=list(axs="i")),
        ylim = c(ncol(Mss)+0.7,0.3),
        prepanel = function(x, y) {
          list(xlim = quantile(x, probs = c(0.01, 0.99), na.rm=TRUE))
        },
        panel.groups = function(x, y, ...) {
          panel.bwplot(x, y, ...)
        })
  print(box);
  dev.off();
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
