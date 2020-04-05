##################################################
## R scripts for miRNet
## Various utility methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# new range [a, b]
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param qvec PARAM_DESCRIPTION
#' @param a PARAM_DESCRIPTION
#' @param b PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname rescale2NewRange
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
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname AutoNorm
#' @export 
AutoNorm<-function(x){
    (x - mean(x))/sd(x, na.rm=T);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param min.val PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname LogNorm
#' @export 
LogNorm<-function(x, min.val){
    log2((x + sqrt(x^2 + min.val^2))/2)
}

# #FFFFFF to rgb(1, 0, 0)
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param cols PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname hex2rgba
#' @export 
hex2rgba <- function(cols){
  return(apply(sapply(cols, col2rgb), 2, function(x){paste("rgba(", x[1], ",", x[2], ",", x[3], ",0.8)", sep="")}));
}

# re-arrange one vector elements according to another vector values
# usually src is character vector to be arranged
# target is numberic vector of same length
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param src.vec PARAM_DESCRIPTION
#' @param tgt.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname sync2vecs
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
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nd.vec PARAM_DESCRIPTION
#' @param col.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname getExpColors
#' @export 
getExpColors <- function(nd.vec, col.vec){
    nvec <- rep("", length(nd.vec));
    m.inx <- is.null(nd.vec) | is.na(nd.vec);
    nvec[m.inx] <- col.vec[3];
    pos.inx <- nd.vec > 0;
    nvec[pos.inx] <- col.vec[2];
    nvec[!pos.inx] <- col.vec[1];
    as.character(nvec);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param my.grps PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetColorSchema
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
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param what PARAM_DESCRIPTION, Default: c("test", "vector")
#' @param extras PARAM_DESCRIPTION, Default: c(".", "NA")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname all.numeric
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

# utils to remove from
# within, leading and trailing spaces
# remove /
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param query PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ClearStrings
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
#' @rdname GetBashFullPath
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
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname myave
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
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mat PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname CalculatePairwiseDiff
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param vec PARAM_DESCRIPTION
#' @param unit PARAM_DESCRIPTION, Default: 10
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetExtendRange
#' @export 
GetExtendRange<-function(vec, unit=10){
    var.max <- max(vec);
    var.min <- min(vec);
    exts <- (var.max - var.min)/unit;
    c(var.min-exts, var.max+exts);
}

# perform scale on row (scale is on column)
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname RowScale
#' @export 
RowScale <- function(x){
    x <- sweep(x, 1L, rowMeans(x, na.rm = T), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = T)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE);
    x
}

# utils to remove from
# within, leading and trailing spaces
# remove /
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param cls.nm PARAM_DESCRIPTION
#' @param query PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ClearFactorStrings
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
    lines <- strsplit(my.input, "\r|\n|\r\n")[[1]];
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
#' @rdname ReadTabData
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
    # now remove all comments in dat1
    # assign rownames after covert to matrix as data.frame does not allow duplicate names
    comments.inx <- grep("^#", dat1[,1]);
    dat1.nms <- dat1[-comments.inx,1];
    dat1<-dat1[-comments.inx,-1];
    dat1 <- data.matrix(dat1);
    rownames(dat1) <- dat1.nms;

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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param db.path PARAM_DESCRIPTION
#' @param q.vec PARAM_DESCRIPTION
#' @param table.nm PARAM_DESCRIPTION
#' @param col.nm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Query.miRNetDB
#' @export 
Query.miRNetDB <- function(db.path, q.vec, table.nm, col.nm){
  db.path <- paste0(db.path, ".sqlite");
  db.url <- paste(sqlite.path, db.path, sep="");
  msg <- paste("Downloading", db.path, "from", db.url);
  print(msg);
  download.file(db.url, db.path);
  mir.db <- dbConnect(SQLite(), db.path);
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", table.nm, " WHERE ", col.nm," IN (", query, ")", sep="");
  mir.dic <- .query.sqlite(mir.db, statement);

  #Check the matched miRNA from upload list to database
  if (col.nm == "mir_id"){
    mir.lib <- as.vector(unique(mir.dic$mir_id));
    notMatch <- setdiff(q.vec, mir.lib);

    if (length(notMatch) > 0){ # Converting miRBase version and mature id.
      print("Converting ids for different miRBase versions ....");
      notMatch <- gsub("mir", "miR", notMatch);
      load("../../data/libs/mbcdata.rda");

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
        df <- data.frame(uid = uid, SYM = SYM_ID, stringsAsFactors=FALSE)
        df <- merge(df, ACC_SYM)[, c("uid", "ACC")]
        df <- unique( merge(df, VMAP, by="ACC") )

        target <- data.frame(
          OriginalName = df$uid,
          TargetName = SYM[df$SYM],
          Accession = ACC[df$ACC],
          stringsAsFactors = FALSE
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
          row.names=NULL, stringsAsFactors = FALSE);

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
          row.names=NULL, stringsAsFactors = FALSE)
        target2=target2[match(miRNANames, target2$OriginalName),]

        #merge
        mir.vec <- tolower(as.vector(target$TargetName));
        mir.vec2 <- c(tolower(as.vector(target2$Mature1)),tolower(as.vector(target2$Mature2)));
        mir.vec3 <- tolower(paste(target$OriginalName,"-3p",sep=""));
        mir.vec4 <- tolower(paste(target$OriginalName,"-5p",sep=""));
        mir.vec <- na.omit(c(mir.vec, mir.vec2,mir.vec3,mir.vec4))
        mir.vec <- unique(unlist(strsplit(mir.vec, split="&")));
        query2 <- paste(shQuote(mir.vec), collapse=",");
        statement2 <- paste("SELECT * FROM ", table.nm, " WHERE ", col.nm," IN (", query2, ")", sep="");
        mir.db <- dbConnect(SQLite(), db.path);
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
      download.file(path,"mir_tissue.csv");
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param table.nm PARAM_DESCRIPTION
#' @param q.vec PARAM_DESCRIPTION
#' @param col.nm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname QueryTFSQLite
#' @export 
QueryTFSQLite<- function(table.nm, q.vec, col.nm){
  require('RSQLite');
  db.path <- paste0("tf2gene.sqlite");
  db.url <- paste(sqlite.path, db.path, sep="");
  msg <- paste("Downloading", db.path, "from", db.url);
  print(msg);
  download.file(db.url, db.path);
  tf.db <- dbConnect(SQLite(), db.path);
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", table.nm, " WHERE ", col.nm," IN (",query,")", sep="");
  return(.query.sqlite(tf.db, statement));
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param n PARAM_DESCRIPTION, Default: 10
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname cleanMem
#' @export 
cleanMem <- function(n=10) { for (i in 1:n) gc() }

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param n PARAM_DESCRIPTION
#' @param center PARAM_DESCRIPTION, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname generate_breaks
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nd.vec PARAM_DESCRIPTION
#' @param background PARAM_DESCRIPTION, Default: 'black'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ComputeColorGradient
#' @export 
ComputeColorGradient <- function(nd.vec, background="black"){
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param background PARAM_DESCRIPTION
#' @param center PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetColorGradient
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param col PARAM_DESCRIPTION, Default: rainbow(10)
#' @param breaks PARAM_DESCRIPTION, Default: NA
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname scale_vec_colours
#' @export 
scale_vec_colours = function(x, col = rainbow(10), breaks = NA){
    return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mat PARAM_DESCRIPTION
#' @param col PARAM_DESCRIPTION, Default: rainbow(10)
#' @param breaks PARAM_DESCRIPTION, Default: NA
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname scale_colours
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
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @param n PARAM_DESCRIPTION, Default: 30
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ShowMemoryUse
#' @export 
ShowMemoryUse <- function(..., n=30) {
    library(pryr);
    sink(); # make sure print to screen
    print(mem_used());
    print(sessionInfo());
    print(.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n));
    print(warnings());
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
#' @rdname CleanMemory
#' @export 
CleanMemory <- function(){
    for (i in 1:10){
        gc(reset = T);
    }
}

# private method for all sqlite queries
.query.sqlite <- function(db.con, statement, offline=TRUE){
  rs <- dbSendQuery(db.con, statement);
  res <- fetch(rs, n=-1); # get all records
  dbClearResult(rs);
  if(offline){
    dbDisconnect(db.con);
  }
  cleanMem();
  return(res);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param table.nm PARAM_DESCRIPTION
#' @param q.vec PARAM_DESCRIPTION
#' @param requireExp PARAM_DESCRIPTION
#' @param min.score PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname QueryPpiSQLiteZero
#' @export 
QueryPpiSQLiteZero <- function(table.nm, q.vec, requireExp, min.score){
    require('RSQLite');
    db.path <- paste0("ppi.sqlite");
    db.url <- paste(sqlite.path, db.path, sep="");
    msg <- paste("Downloading", db.path, "from", db.url);
    print(msg);
    download.file(db.url, db.path);
    ppi.db <- dbConnect(SQLite(), db.path);
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname simpleCap
#' @export 
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param grp.num PARAM_DESCRIPTION
#' @param filenm PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname gg_color_hue
#' @export 
gg_color_hue <- function(grp.num, filenm=NULL) {
    grp.num <- as.numeric(grp.num)
    pal18 <- c( "#911eb4", "#3cb44b", "#4363d8",  "#f032e6", "#ffe119", "#e6194B", "#f58231", "#bfef45", "#fabebe", "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#42d4f4","#000075");
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
