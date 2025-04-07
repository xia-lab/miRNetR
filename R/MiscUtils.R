##################################################
## R scripts for miRNet
## Various utility methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# new range [a, b]
#' Rescale to New Range
#' @param qvec Vector to be rescaled.
#' @param a Lower bound.
#' @param b Upper bound.
#' @export
rescale2NewRange <- function(qvec, a, b){
    q.min <- min(qvec);
    q.max <- max(qvec);
    if(length(qvec) < 50){
        a <- a*2;
    }
    if(q.max == q.min){
        new.vec <- rep(8, length(qvec));
    }else{
        coef.a <- (b-a)/(q.max-q.min);
        const.b <- b - coef.a*q.max;
        new.vec <- coef.a*qvec + const.b;
    }
    return(new.vec);
}

`%fin%` <- function(x, table) {
  fmatch(x, table, nomatch = 0L) > 0L
}

`%notin%` <- Negate(`%in%`);

# normalize to zero mean and unit variance
#' Auto Normalization
#' @param x Variable to normalize.
#' @export
AutoNorm<-function(x){
    (x - mean(x))/sd(x, na.rm=T);
}

#' Log Normalization
#' @param x Variable to normalize.
#' @param min.val Min value.
#' @export
LogNorm<-function(x, min.val){
    log2((x + sqrt(x^2 + min.val^2))/2)
}

# #FFFFFF to rgb(1, 0, 0)
#' Hex color to RGBA format
#' @param cols Hexadecimal color code (e.g., "#FFFFFF").
#' @export
hex2rgba <- function(cols){
  return(apply(sapply(cols, col2rgb), 2, function(x){paste("rgba(", x[1], ",", x[2], ",", x[3], ",0.8)", sep="")}));
}

# re-arrange one vector elements according to another vector values
# usually src is character vector to be arranged
# target is numberic vector of same length
#' Synchronize Two Vectors
#' @param src.vec Source character vector to be arranged.
#' @param tgt.vec Target numeric vector of the same length.
#' @export
sync2vecs <- function(src.vec, tgt.vec){
    if(length(unique(src.vec)) != length(unique(tgt.vec))){
        print("must be of the same unique length!");
        return();
    }

    tgt.vec[is.na(tgt.vec)] <- min(tgt.vec, na.rm=T)-1;
    names(src.vec) <- sort(unique(tgt.vec));
    tgt.vec <- as.character(tgt.vec);
    as.character(src.vec[tgt.vec]);
}

# col vec is for low high null
#' Colors Based on Expression
#' @param nd.vec Numeric vector containing values to be mapped to colors.
#' @param col.vec Color vector.
#' export
getExpColors <- function(nd.vec, col.vec){
    nvec <- rep("", length(nd.vec));
    m.inx <- is.null(nd.vec) | is.na(nd.vec);
    nvec[m.inx] <- col.vec[3];
    pos.inx <- nd.vec > 0;
    nvec[pos.inx] <- col.vec[2];
    nvec[!pos.inx] <- col.vec[1];
    as.character(nvec);
}

#' Get Color Schema
#' @param my.grps Groups.
#' @export
GetColorSchema <- function(my.grps){
    # test if total group number is over 9
     grp.num <- length(levels(my.grps));

     if(grp.num > 9){
        pal12 = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                    "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
                    "#FFFF99", "#B15928");
        dist.cols <- colorRampPalette(pal12)(grp.num);
        lvs <- levels(my.grps);
        colors <- vector(mode="character", length=length(my.grps));
        for(i in 1:length(lvs)){
            colors[my.grps == lvs[i]] <- dist.cols[i];
        }
     }else{
        colors <- as.numeric(my.grps)+1;
     }
    return (colors);
}

# borrowed from Hmisc
#' Check/Convert Values to Numeric
#' @param x Character vector to check/convert.
#' @param what "test" (check if numeric) or "vector" (convert to numeric).
#' @param extras Additional values to be treated as missing or empty.
#' @export
all.numeric <- function (x, what = c("test", "vector"), extras = c(".", "NA")){
    what <- match.arg(what)
    old <- options(warn = -1)
    on.exit(options(old));
    x <- sub("[[:space:]]+$", "", x);
    x <- sub("^[[:space:]]+", "", x);
    inx <- x %in% c("", extras);
    xs <- x[!inx];
    isnum <- !any(is.na(as.numeric(xs)))
    if (what == "test")
        isnum
    else if (isnum)
        as.numeric(x)
    else x
}

# utils to remove from within, leading and trailing spaces
#' Clear String
#' @param query Query.
#' @export
ClearStrings<-function(query){
    # remove leading and trailing space
    query<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", query, perl=TRUE);

    # kill multiple white space
    query <- gsub(" +",".",query);
    query <- gsub("/", ".", query);
    query <- gsub("-", ".", query);
    return (query);
}

# need to obtain the full path to convert (from imagemagik) for cropping images
#' Get Bash Full Path
#' @export
GetBashFullPath<-function(){
    path <- system("which bash", intern=TRUE);
    if((length(path) == 0) && (typeof(path) == "character")){
        print("Could not find bash in the PATH!");
        return("NA");
    }
    return(path);
}

# overwrite ave, => na.rm=T
#' Average
#' @param x Numeric vector to calculate mean.
#' @param ... Optional grouping factors.
#' @export
myave <- function (x, ...) {
    n <- length(list(...))
    if (n) {
        g <- interaction(...)
        split(x, g) <- lapply(split(x, g), mean, na.rm=T)
    }
    else x[] <- FUN(x, na.rm=T)
    x
}

# log scale ratio
#' Calculate Pairwise Difference
#' @param mat Numeric matrix.
#' @export
CalculatePairwiseDiff <- function(mat){
    f <- function(i, mat) {
       z <- mat[, i-1] - mat[, i:ncol(mat), drop = FALSE]
       colnames(z) <- paste(colnames(mat)[i-1], colnames(z), sep = "/")
       z
    }
    res <- do.call("cbind", sapply(2:ncol(mat), f, mat));
    round(res,5);
}

#' Get Extend Range
#' @param vec Vector input.
#' @param unit Unit to extend the range.
#' @export
GetExtendRange<-function(vec, unit=10){
    var.max <- max(vec);
    var.min <- min(vec);
    exts <- (var.max - var.min)/unit;
    c(var.min-exts, var.max+exts);
}

# perform scale on row (scale is on column)
#' Row Scale
#' @param x Numeric matrix to standardize.
#' @export
RowScale <- function(x){
    x <- sweep(x, 1L, rowMeans(x, na.rm = T), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = T)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE);
    x
}

# utils to remove from within, leading and trailing spaces
#' Clear Factor Strings
#' @param cls.nm Character string to be prefixed.
#' @param query Character vector of string to be cleaned.
#' @export
ClearFactorStrings<-function(cls.nm, query){
    # remove leading and trailing space
    query<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", query, perl=TRUE);

    # kill multiple white space
    query <- gsub(" +","_",query);
    # remove non alphabets and non numbers
    query <- gsub("[^[:alnum:] ]", "_", query);

    # test all numbers (i.e. Time points)
    chars <- substr(query, 0, 1);
    num.inx<- chars >= '0' & chars <= '9';
    if(all(num.inx)){
        query = as.numeric(query);
        nquery <- paste(cls.nm, query, sep="_");
        query <- factor(nquery, levels=paste(cls.nm, sort(unique(query)), sep="_"));
    }else{
        query[num.inx] <- paste(cls.nm, query[num.inx], sep="_");
        query <- factor(query);
    }
    return (query);
}

# parse two-column list as a string input from text area in web page
.parseListData <- function(my.input){
  if(grepl("\n", my.input[1])) {
    lines <- strsplit(my.input, "\r|\n|\r\n")[[1]];
  }else{
    lines <- my.input;
  }
    my.lists <- strsplit(lines, "\\s+");
    my.mat <- do.call(rbind, my.lists);
    if(dim(my.mat)[2] == 1){ # add *
        my.mat <- cbind(my.mat, rep("*", nrow(my.mat)));
    }else if(dim(my.mat)[2] > 2){
        my.mat <- my.mat[,1:2];
        current.msg <- "More than two columns found in the list. Only first two columns will be used.";
        print(currret.msg);
    }
    return(my.mat);
}

.parsePickListItems <- function(my.vec){
    my.mat <- cbind(my.vec, rep("*", length(my.vec)));
    rownames(my.mat) <- my.mat[,1];
    my.mat <- my.mat[,-1, drop=F];
    return(my.mat);
}

# read tab delimited file
# stored in dataSet list object
# can have many classes, stored in meta.info
#' Read Tab Delimited Data
#' @param dataName Data name.
#' @export
ReadTabData <- function(dataName) {

    msg <- NULL;
    # using the powerful fread function, 10 times faster, note: default return data.table, turn off
    dat1 <- .readDataTable(dataName);

    # look for #CLASS, could have more than 1 class labels, store in a list
    meta.info <- list();
    cls.inx <- grep("^#CLASS", dat1[,1]);
    if(length(cls.inx) > 0){
        for(i in 1:length(cls.inx)){
            inx <- cls.inx[i];
            cls.nm <- substring(dat1[inx, 1],2); # discard the first char #
            if(nchar(cls.nm) > 6){
                cls.nm <- substring(cls.nm, 7); # remove class
            }
            cls.lbls <- dat1[inx, -1];
            # test NA
            na.inx <- is.na(cls.lbls);
            cls.lbls[na.inx] <- "NA";
            cls.lbls <- ClearFactorStrings(cls.nm, cls.lbls);

            meta.info[[cls.nm]] <- cls.lbls;
        }
    }else{
        current.msg <<- "No metadata labels #CLASS found in your data!";
        return("F");
    }

    meta.info <- data.frame(meta.info);
    dat1 <- .to.numeric.mat(dat1);

    list(
        name = basename(dataName),
        data=dat1,
        meta.info=meta.info
    );
}

.readDataTable <- function(fileName){
    if(length(grep('\\.zip$',fileName,perl=TRUE))>0){
       fileName <- unzip(fileName);
       if(length(fileName) > 1){
            # test if "__MACOSX" or ".DS_Store"
            osInx <- grep('MACOSX',fileName,perl=TRUE);
            if(length(osInx) > 0){
                fileName <- fileName[-osInx];
            }
            dsInx <- grep('DS_Store',fileName,perl=TRUE);
            if(length(dsInx) > 0){
                fileName <- fileName[-dsInx];
            }
            dat.inx <- grep(".[Tt][Xx][Tt]$", fileName);
            if(length(dat.inx) != 1){
                current.msg <<- "More than one text files (.txt) found in the zip file.";
                return(0);
            }
       }
    }
    dat <- tryCatch(
            data.table::fread(fileName, header=TRUE, check.names=FALSE, blank.lines.skip=TRUE, data.table=FALSE),
            error=function(e){
                print(e);
                return(.my.slowreaders(fileName));
            },
            warning=function(w){
                print(w);
                return(.my.slowreaders(fileName));
            });

    if(any(dim(dat) == 0)){
        dat <- .my.slowreaders(fileName);
    }
    return(dat);
}

.my.slowreaders <- function(fileName){
  print("Using slower file reader ...");
  formatStr <- substr(fileName, nchar(fileName)-2, nchar(fileName))
  if(formatStr == "txt"){
    dat <- try(read.table(fileName, header=TRUE, comment.char = "", check.names=F, as.is=T));
  }else{ # note, read.csv is more than read.table with sep=","
    dat <- try(read.csv(fileName, header=TRUE, comment.char = "", check.names=F, as.is=T));
  }
  return(dat);
}

#' Query miRNet DB
#' @param db.path Path to database.
#' @param q.vec Vector to query.
#' @param table.nm Table name.
#' @param col.nm Column names.
#' @param db.nm Name of the database.
#' @export
Query.miRNetDB <- function(db.path, q.vec, table.nm, col.nm, db.nm = "mirtarbase"){

  db.path <- paste0(db.path, ".sqlite");
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
  query <- paste (shQuote(q.vec),collapse=",");
  if(grepl("mir2gene", db.path) && col.nm %in% c("mir_id", "mir_acc") ){
    statement <- paste("SELECT * FROM ", table.nm, " WHERE ", col.nm," IN (", query, ")"," AND ", db.nm," == 1 ", sep="");
  }else{
    statement <- paste("SELECT * FROM ", table.nm, " WHERE ", col.nm," IN (", query, ")", sep="");
  }
print(statement);
  mir.dic <- .query.sqlite(mir.db, statement);

  #Check the matched miRNA from upload list to database
  if (col.nm == "mir_id"){
    mir.lib <- as.vector(unique(mir.dic$mir_id));
    notMatch <- setdiff(q.vec, mir.lib);

    if (length(notMatch) > 0){ # Converting miRBase version and mature id.
      print("Converting ids for different miRBase versions ....");
      notMatch <- gsub("mir", "miR", notMatch);
      if(.on.public.web){
        load("../../data/libs/mbcdata.rda");
      }else{
        mbcdata.rda <- paste(lib.path, "/mbcdata.rda", sep="");
        destfile <- paste("mbcdata.rda");
        download.file(mbcdata.rda, destfile, mode = "wb");
        load(destfile);
      }
      
      miRNANames <- gsub(" ","", as.character(notMatch));
      targetVersion <- "v22";

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

        # map to mature form
        VMAP <-miRNA_data[[ver_index]][,c(2,6,9)]
        # [1] "Precursor" "Mature1"   "Mature2"
        VMAP[,1]=SYM[VMAP[,1]]
        VMAP[,2]=SYM[VMAP[,2]] # mature 1
        VMAP[,3]=SYM[VMAP[,3]] # mature 2
        miRNANames <- gsub("miR", "mir", miRNANames);
        miRNANames=as.character(miRNANames)
        miRNANames=gsub(" ","",miRNANames)##Remove the possible space
        uid = unique(as.vector(miRNANames))

        uid=na.omit(uid)
        ind=apply(VMAP,2,function(x){match(uid,x)})
        if(length(miRNANames) == 1){
          ind <- rbind(ind, rep(NA, 3));
        }else if(length(miRNANames) > 1){
          ind <- ind;
        }
        ind[which(is.na(ind[,1])),1]=ind[which(is.na(ind[,1])),2]
        ind[which(is.na(ind[,1])),1]=ind[which(is.na(ind[,1])),3]
        target2 <- data.frame(
          OriginalName = uid,
          Mature1 = VMAP[ind[,1],2],
          Mature2 = VMAP[ind[,1],3],
          row.names=NULL)
        target2=target2[match(miRNANames, target2$OriginalName),]

        #merge
        mir.vec <- as.vector(target$TargetName);
        mir.vec2 <- c(as.vector(target2$Mature1),as.vector(target2$Mature2));
        mir.vec3 <- paste(target$OriginalName,"-3p",sep="");
        mir.vec4 <- paste(target$OriginalName,"-5p",sep="");
        mir.vec <- na.omit(c(mir.vec, mir.vec2,mir.vec3,mir.vec4))
        mir.vec <- unique(unlist(strsplit(mir.vec, split="&")));
        query2 <- paste(shQuote(mir.vec), collapse=",");
        if(grepl("mir2gene", db.path) && col.nm %in% c("mir_id", "mir_acc")){
          statement2 <- paste("SELECT * FROM ", table.nm, " WHERE ", col.nm," IN (", query2, ")", " AND ", db.nm," == 1 ", sep="");
        }else{
          statement2 <- paste("SELECT * FROM ", table.nm, " WHERE ", col.nm," IN (", query2, ")", sep="");
        }
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
        mir.dic2 <- .query.sqlite(mir.db, statement2);

        # now add back to the main data
        mir.dic <- rbind(mir.dic, mir.dic2);

        # remove duplicates
        dup.inx <- duplicated(mir.dic$mirnet);
        mir.dic <- mir.dic[!dup.inx, ];
      }
    }
  }

  if(col.nm == "mir_acc"){
    # when use Accession number, it can match both new and old version, use the new one (old one is ranked later)
    dup.inx <- duplicated(mir.dic[, c("mir_acc", "symbol")]);
    mir.dic <- mir.dic[!dup.inx, ];
  }

  # Perform tissue annotation if specified
  if (nrow(mir.dic) > 0){
    tissue <- dataSet$tissue;
    if(dataSet$org != "hsa" || tissue == "na"){
      mir.dic[, "tissue"] <- "Unspecified";
    }else{
      path <- paste(lib.path, "hsa/mir_tissue.csv", sep="");
      mir_tissue <- read.csv(file=path);
      if (tissue == "Others"){
        ts_df <- mir_tissue[!(mir_tissue$tissue %in% ts_count$tissue), ];
        ts_ano <- aggregate(tissue ~ mir_acc, data = ts_df, paste, collapse = "//");
        mir.dic <- merge(mir.dic, ts_ano, by="mir_acc");
      } else if (tissue != "na" && tissue != "Others"){
        tissue.acc <- unique(mir_tissue[which(mir_tissue$tissue == tissue), "mir_acc"]);
        mir.dic <- mir.dic[which(mir.dic$mir_acc %in% tissue.acc), ];
        if (nrow(mir.dic) > 0){
          mir.dic[, "tissue"] <- tissue;
        }
      } else {
        mir.acc <- unique(mir.dic$mir_acc);
        ts_df <- mir_tissue[which(mir_tissue$mir_acc %in% mir.acc), ];
        if (nrow(ts_df) > 0){
          ts_df <- ts_df[order(ts_df$tissue), ];
          ts_ano <- aggregate(tissue ~ mir_acc, data = ts_df, paste, collapse = "//");
          mir.dic <- merge(mir.dic, ts_ano, by="mir_acc", all.x=T);
          mir.dic[is.na(mir.dic$tissue), "tissue"] <- "Not specified";
        } else {
          mir.dic[, "tissue"] <- "Not specified";
        }
      }
    }
  }

  return(mir.dic);
}

#' Query Transcription Factor
#' @param table.nm Table name.
#' @param q.vec Vector to query.
#' @param col.nm Column names.
#' @export
QueryTFSQLite <- function(table.nm, q.vec, col.nm){
  require('RSQLite');
  db.path <- paste(sqlite.path, "tf2gene.sqlite", sep="");
  if(.on.public.web){
    tf.db <- dbConnect(SQLite(), db.path);
  }else{
    msg <- paste("Downloading", db.path);
    db.name <- gsub(sqlite.path, "", db.path);
    if(!file.exists(db.name)){
      print(msg);
      download.file(db.path, db.name);
    }
    tf.db <- dbConnect(SQLite(), db.name);
  }
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", table.nm, " WHERE ", col.nm," IN (",query,")", sep="");
  return(.query.sqlite(tf.db, statement));
}

#' Clean Memory
#' @export
CleanMemory <- function() { 
    gc(); 
}

#' Get Unique Entries
#' @param db.path Path to database.
#' @param statement Statement.
#' @export
GetUniqueEntries <- function(db.path, statement){
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
    res <- .query.sqlite(mir.db, statement);
    res <- sort(unique(as.character(res[,1])));
    return (res);
}

#' Generate Breaks
#' @export
generate_breaks = function(x, n, center = F){
    if(center){
        m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
        res = seq(-m, m, length.out = n + 1)
    }
    else{
        res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
    }
    return(res)
}

#' Compute Color Gradient
#' @param nd.vec Numeric vector.
#' @param background Background color.
#' @param centered Logical.
#' @export
ComputeColorGradient <- function(nd.vec, background="black", centered){
    library("RColorBrewer");
    if(sum(nd.vec<0, na.rm=TRUE) > 0){
        centered <- T;
    }else{
        centered <- F;
    }
    color <- GetColorGradient(background, centered);
    breaks <- generate_breaks(nd.vec, length(color), center = centered);
    return(scale_vec_colours(nd.vec, col = color, breaks = breaks));
}

#' Get Color Gradient
#' @param backgroun Background.
#' @param center Logical
#' @export
GetColorGradient <- function(background, center){
    if(background == "black"){
        if(center){
            return(c(colorRampPalette(c("#31A231", "#5BC85B", "#90EE90", "#C1FFC1"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)));
        }else{
            return(colorRampPalette(rev(heat.colors(9)))(100));
        }
    }else{ # white background
        if(center){
            return(c(colorRampPalette(c("#137B13", "#31A231", "#5BC85B", "#90EE90"))(50), colorRampPalette(c("#FF7783", "#E32636", "#BD0313", "#96000D"))(50)));
        }else{
            return(colorRampPalette(c("grey", "orange", "red", "darkred"))(100));
        }
    }
}

#' Scale Vector Colours
#' @export
scale_vec_colours = function(x, col = rainbow(10), breaks = NA){
    return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])
}

#' Scale Matrix Colours
#' @export
scale_colours = function(mat, col = rainbow(10), breaks = NA){
    mat = as.matrix(mat)
    return(matrix(scale_vec_colours(as.vector(mat), col = col, breaks = breaks), nrow(mat), ncol(mat), dimnames = list(rownames(mat), colnames(mat))))
}

###########
# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           capture.output(format(utils::object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

# shorthand
#' Show Memory Use
#' @export
ShowMemoryUse <- function(..., n=30) {
    library(pryr);
    sink(); # make sure print to screen
    print(mem_used());
    print(sessionInfo());
    print(.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n));
    print(warnings());
}


# private method for all sqlite queries
.query.sqlite <- function(db.con, statement, offline=TRUE){
  rs <- dbSendQuery(db.con, statement);
  res <- fetch(rs, n=-1); # get all records
  dbClearResult(rs);
  if(offline){
    dbDisconnect(db.con);
  }
  CleanMemory();
  return(res);
}

#' Query Protein-Protein Interaction
#' @export
QueryPpiSQLiteZero <- function(table.nm, q.vec, requireExp, min.score){
    require('RSQLite')
  db.path <- paste(sqlite.path, "ppi.sqlite", sep="");
  if(.on.public.web){
    ppi.db <- dbConnect(SQLite(), db.path);
  }else{
    msg <- paste("Downloading", db.path);
    db.name <- gsub(sqlite.path, "", db.path);
    if(!file.exists(db.name)){
      print(msg);
      download.file(db.path, db.name);
    }
    ppi.db <- dbConnect(SQLite(), db.name);
  }
  query <- paste(shQuote(q.vec),collapse=",");

    if(grepl("string$", table.nm)){
        if(requireExp){
            statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")) OR (name1 IN (", query, "))OR (name2 IN (", query, ")))  AND combined_score >=", min.score, " AND experimental > 0", sep="");
        }else{
            statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")) OR (name1 IN (", query, "))OR (name2 IN (", query, ")))  AND combined_score >=", min.score, sep="");
        }
    }else{
        statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")) OR (name1 IN (", query, "))OR (name2 IN (", query, ")))", sep="");
    }

    ppi.res <- .query.sqlite(ppi.db, statement);
    hit.inx1 <- ppi.res[,1] %in% q.vec
    hit.inx2 <- ppi.res[,2] %in% q.vec
    ppi.res1 <- ppi.res[(hit.inx1 & hit.inx2),]

    hit.inx3 <- ppi.res[,3] %in% q.vec
    hit.inx4 <- ppi.res[,4] %in% q.vec
    ppi.res2 <- ppi.res[(hit.inx3 & hit.inx4),]
    ppi.res = rbind(ppi.res1,ppi.res2)

    return(ppi.res);
}

#' Capitalize First Letter
#' @export
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}

#' Group Color Palette
#' @export
gg_color_hue <- function(grp.num, filenm=NULL) {
    grp.num <- as.numeric(grp.num)
    pal18 <- c("#3cb44b", "#f032e6", "#ffe119", "#e6194B", "#f58231", "#bfef45", "#fabebe", "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#42d4f4","#000075", "#ff4500");
    if(grp.num <= 18){ # update color and respect default
        colArr <- pal18[1:grp.num];
    }else{
        colArr <- colorRampPalette(pal18)(grp.num);
    }
    if(is.null(filenm)){
        return(colArr);
    }else{
        sink(filenm);
        cat(toJSON(colArr));
        sink();
        return(filenm);
    }
}

# obtain a numeric matrix, exclude comments if any
.to.numeric.mat <- function(dat1){
  # now remove all comments in dat1
  # assign rownames after covert to matrix as data.frame does not allow duplicate names
  comments.inx <- grep("^#", dat1[,1]);
  if(sum(comments.inx) > 0){
    row.nms <- dat1[-comments.inx,1];
    dat1 <- dat1[-comments.inx,-1];
  }else{
    row.nms <- dat1[,1];
    dat1 <- dat1[,-1];
  }
  dimensions <- dim(dat1)
  col.nms <- colnames(dat1)
  dat1 <- sapply(dat1, as.numeric);
  dat1 <- matrix(data=dat1, ncol=dimensions[2], nrow=dimensions[1])
  rownames(dat1) <- row.nms;
  colnames(dat1) <- col.nms;
  return(dat1);
}

fast.write.csv <- function(dat, file, row.names=TRUE){
    tryCatch(
        {
           if(is.data.frame(dat)){
                # there is a rare bug in data.table (R 3.6) which kill the R process in some cases 
                data.table::fwrite(dat, file, row.names=row.names);
           }else{
                write.csv(dat, file, row.names=row.names);  
           }
        }, error=function(e){
            print(e);
            fast.write.csv(dat, file, row.names=row.names);   
        }, warning=function(w){
            print(w);
            fast.write.csv(dat, file, row.names=row.names); 
        });
}

makeReadable <- function(str){
    result <- switch(str,
                 pct = "Percent",
                 abs = "Absolute",
                 log = "Log2",
                 rle = "RLE",
                 array = "Microarray",
                 count= "RNA-Seq",
                 hsa = "H. sapiens (human)",
                 mmu = "M. musculus (mouse)",
                 rno = "R. norvegicus (rat)",
                 cel = "C. elegans (roundworm)",
                 dme = "D. melanogaster (fruitfly)",
                 dre = "D. rerio (zebrafish)",
                 sce = "S. cerevisiae (yeast)",
                 eco = "E. coli",
                 ath = "A. thaliana (Arabidopsis)",
                 bta = "B. taurus (cow)",
                 gga = "G. gallus (chicken)",
                 mun = "M. unguiculatus (Mongolian gerbil)",
                 bsu = "B. subtilis",
                 pae = "P. aeruginosa",
                 mtb = "M. tuberculosis",
                 smm = "S. mansoni (schistosomiasis)",
                 tbr = "T. brucei (trypanosoma)",
                 pfa = "P. falciparum (malaria)",
                 cjo = "C. japonica (japanese quail)",
                 xla = "X. laevis (African clawed frog)",
                 ppr = "P. promelas (fathead minnow; custom)",
                 fhm = "P. promelas (fathead minnow; NCBI)",
                 nlf = "L. pipiens (northern leopard frog)",
                 omk = "O. mykiss (rainbow trout)",
                 ham = "H. americanus (American lobster)",
                 cdi = "C. dilutus",
                 dma = "D. magna",
                 rsu = "R. subcapitata",
                 haz = "H. azteca",
                 fcd = "F. candida",
                 "entrez" = "Entrez ID",
                 "refseq" = "RefSeq ID",
                   "gb" = "Genbank ID",
                   "symbol" = "Official Gene Symbol",
                   "embl_gene" = "Ensembl Gene ID",
                   "embl_transcript" = "Ensemble Transcript ID",
                   "embl_protein" = "Ensembl Protein ID",
                   "uniprot" = "Uniprot Accession ID",
                   "hgu95a" = "Affymetrix Human Genome U95 (chip hgu95a)",
                   "hgu95av2" = "Affymetrix Human Genome U95 (chip hgu95av2)",
                   "hgu95b" = "Affymetrix Human Genome U95 (chip hgu95b)",
                   "hgu95c" = "Affymetrix Human Genome U95 (chip hgu95c)",
                   "hgu95d" = "Affymetrix Human Genome U95 (chip hgu95d)",
                   "hgu95e" = "Affymetrix Human Genome U95 (chip hgu95e)",
                   "hgu133a" = "Affymetrix Human Genome U133 (chip hgu133a)",
                   "hgu133b" = "Affymetrix Human Genome U133 (chip hgu133b)",
                   "hgu133plus2" = "Affymetrix Human Genome U133plus2 (hgu133plus2)",
                   "hgu133plus2pm" = "Affymetrix Human Genome U133plus2_PM (hgu133plus2pm)",
                   "lumiht12v3" = "Illumina HumanHT-12 V3 BeadArray",
                   "lumiht12v4" = "Illumina HumanHT-12 V4 BeadArray",
                   "lumiref8v2" = "Illumina HumanRef-8 V2 BeadArray",
                   "lumiref8v3" = "Illumina HumanRef-8 V3 BeadArray",
                   "lumiwg6v2" = "Illumina HumanWG-6 V2 BeadArray",
                   "lumiwg6v3" = "Illumina HumanWG-6 V3 BeadArray",
                   "agi4100a" = "Agilent Human 1 cDNA Microarray (4100A)",
                   "agi4101a" = "Agilent Human 2 cDNA Microarray (4101A)",
                   "agi4110b" = "Agilent Human 1A cDNA Microarray (4110B)",
                   "agi4111a" = "Agilent Human 1B cDNA Microarray (4111A)",
                   "agi4112a" = "Agilent Human Genome Whole Microarray (4x44k/4112)",
                   "agi4845a" = "Agilent Human AMADID 026652 Microarray (4845A)",
                   "lumiwg6v1" = "Illumina MouseWG-6 v1.0 Bead Array",
                   "lumiwg6v11" = "Illumina MouseWG-6 v1.1 Bead Array",
                   "lumiwg6v2" = "Illumina MouseWG-6 v2.0 Bead Array",
                   "lumiref8v1" = "Illumina MouseRef-8 v1.0 Bead Array",
                   "lumiref8v2" = "Illumina MouseRef-8 v2.0 Bead Array",
                   "mgu74a" = "Affymetrix Murine Genome U74v2 (chip mgu74a)",
                   "mgu74av2" = "Affymetrix Murine Genome U74v2 (chip mgu74av2)",
                   "mgu74b" = "Affymetrix Murine Genome U74v2 (chip mgu74b)",
                   "mgu74bv2" = "Affymetrix Murine Genome U74v2 (chip mgu74bv2)",
                   "mgu74c" = "Affymetrix Murine Genome U74v2 (chip mgu74c)",
                   "mgu74cv2" = "Affymetrix Murine Genome U74v2 (chip mgu74cv2)",
                   "moe430a" = "Affymetrix Mouse Expression Set 430 (chip moe430a)",
                   "moe430b" = "Affymetrix Mouse Expression Set 430 (chip moe430b)",
                   "moe430_2" = "Affymetrix GeneChip Mouse Genome 430 2.0",
                   "mgi_st1" = "Affymetrix Mouse Gene 1.0 ST Array",
                   "mgu4101a" = "Agilent Mouse Array (chip mgug4104a)",
                   "mgu4120a" = "Agilent Mouse Array (chip mgug4120a)",
                   "mgu4121a" = "Agilent Mouse Array (chip mgug4121a)",
                   "mgu4122a" = "Agilent Mouse Array (chip mgug4122a)",
                   "kegg" = "KEGG",
                    "reactome" = "Reactome",
                    "go_bp" = "GO:BP",
                    "go_mf" = "GO:MF",
                    "go_cc" = "GO:CC",
                    "panth" = "PANTHER Slim",
                    "motif_set" = "Motif",
                 str)
}


CheckDetailsTablePerformed <-function(type){
  performed <- T;
  if(type == "node"){
    performed <- file.exists("node_table.csv");
  }else if(type %in% c( "network_enr", "regNetwork_enr", "gba_enr", "module_enr", "defaultEnr")){
    clean_type <- gsub("_enr", "", type);
    performed <- !is.null(infoSet$imgSet$enrTables[[clean_type]]);
  }
  print(paste("checkPerformed=", type, "====",performed));

return(performed)
}


GetNodeMat <- function(){
  if(is.null(dataSet$imgSet$node_table)){
    df <- .readDataTable('node_table.csv')
    df[,-c(1:2)] <- lapply(df[,-c(1:2)], function(col) as.numeric(as.character(col)))
    dataSet$imgSet$node_table <<- df;
  }
  return(as.matrix(dataSet$imgSet$node_table[,-c(1:2)]))  # ensure matrix of numerics
}

GetNodeRowNames <- function(){
  if(is.null(dataSet$imgSet$node_table)){
  df <- .readDataTable('node_table.csv')
    dataSet$imgSet$node_table <<- df;

  }
  dataSet$imgSet$node_table$Id;
}

GetNodeGeneSymbols <- function(){
  if(is.null(dataSet$imgSet$node_table)){
  df <- .readDataTable('node_table.csv')
    dataSet$imgSet$node_table <<- df;

  }
  dataSet$imgSet$node_table$Label;
}

GetNodeColNames <- function(){
  if(is.null(dataSet$imgSet$node_table)){
  df <- .readDataTable('node_table.csv')
    dataSet$imgSet$node_table <<- df;

  }
  return(colnames(dataSet$imgSet$node_table[,-c(1:2)]));

}

CleanNumber <-function(bdata){
  if(sum(bdata==Inf)>0){
    inx <- bdata == Inf;
    bdata[inx] <- NA;
    bdata[inx] <- 999999;
  }
  if(sum(bdata==-Inf)>0){
    inx <- bdata == -Inf;
    bdata[inx] <- NA;
    bdata[inx] <- -999999;
  }
  bdata;
}


GetEnrResultMatrix <-function(type){
  infoSet <- readSet(infoSet, "infoSet"); 

  imgSet <- infoSet$imgSet;
  res <- imgSet$enrTables[[type]]$res.mat
  res <- suppressWarnings(apply(res, 2, as.numeric)); # force to be all numeric
  return(signif(as.matrix(res), 5));
}

GetEnrResultColNames<-function(type){
  infoSet <- readSet(infoSet, "infoSet"); 

  imgSet <- infoSet$imgSet;
  res <- imgSet$enrTables[[type]]$res.mat
  colnames(res);
}

GetEnrResSetIDs<-function(type){
  infoSet <- readSet(infoSet, "infoSet"); 

  imgSet <- infoSet$imgSet; 
  res <- imgSet$enrTables[[type]]$table;
  return(res$IDs);
}

GetEnrResSetNames<-function(type){
  infoSet <- readSet(infoSet, "infoSet"); 

  imgSet <- infoSet$imgSet;  res <- imgSet$enrTables[[type]]$table;
  if("Pathway" %in% colnames(res)){
  return(res$Pathway);
  }else if("Name" %in% colnames(res)){
  return(res$Name);
  }else{
    return(res[,1]);
  }

}


PerformDefaultEnrichment <- function(file.nm, fun.type, algo="ora"){
  require("igraph");
  net.nm <- names(mir.nets)[1];
  my.ppi <- mir.nets[[net.nm]];
  IDs <- V(my.ppi)$name;
  names(IDs) <- IDs;
  save.type <- "defaultEnr";
  PerformMirTargetEnrichAnalysis("", fun.type, file.nm, IDs, algo, mode="serial",save.type);

  return(1);
}


GetSetIDLinks <- function(type=""){
  infoSet <- readSet(infoSet, "infoSet"); 

  imgSet <- infoSet$imgSet;
  fun.type <- imgSet$enrTables[[type]]$library;

  ids <- imgSet$enrTables[[type]]$table$IDs
  pathways <- imgSet$enrTables[[type]]$table$Pathway
  print("GetSetIDLinks");
  print(imgSet$enrTables[[type]]$library);

    if(fun.type %in% c("go_bp", "go_mf", "go_cc")){
        annots <- paste("<a href='https://www.ebi.ac.uk/QuickGO/term/", ids, "' target='_blank'>", pathways, "</a>", sep="");
    }else if(fun.type %in% c("go_panthbp", "go_panthmf", "go_panthcc")){
        annots <- paste("<a href='https://www.pantherdb.org/panther/categoryList.do?searchType=basic&fieldName=all&organism=all&fieldValue=", ids, "&listType=5' target='_blank'>", pathways, "</a>", sep="");
    }else if(fun.type == "kegg"){
        annots <- paste("<a href='https://www.genome.jp/dbget-bin/www_bget?pathway+", ids, "' target='_blank'>", pathways, "</a>", sep="");
    }else if(fun.type == "reactome"){
        annots <- paste("<a href='https://reactome.org/content/query?q=", ids, "' target='_blank'>", pathways, "</a>", sep="");
    }else{
        annots <- ids;
    }
  
  print(head(annots));
  return(annots);
}
