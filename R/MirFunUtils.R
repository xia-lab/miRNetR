##################################################
## R script for miRNet
## Description: GO/Pathway ORA
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# note: hit.query, resTable must synchronize
# Due to bias in mirRNA targets, for mir2gene, the p value will be tested
# using 1000 permutation. This should NOT be used for gene2mir
# If users click more than once using the empirical approach
# We can take advantage of the accumulating the permutation result to increase the significance

PerformMirTargetEnrichAnalysis <- function(adjust.type, fun.type, file.nm, IDs, algo, mode="serial"){
     if(!exists("my.mir.target.enrich")){ # public web on same user dir
        compiler::loadcmp("../../rscripts/miRNetR/R/_utils_mir_target_enrich.Rc");    
    }
    return(my.mir.target.enrich(adjust.type, fun.type, file.nm, IDs, algo, mode));
}

PerformAPIMirTargetEnrichAnalysis <- function(adjust.type, fun.type, file.nm, IDs, algo, mode = "parallel"){
  
  # send all objects needed for API
  toSend <- list(dataSet = dataSet, 
                 adjust.type = adjust.type,
                 fun.type = fun.type,
                 file.nm = file.nm,
                 ids = IDs,
                 algo = algo,
                 current.mirnet = current.mirnet,
                 data.type = data.type)
  
  library(httr)
  base <- api.base
  endpoint <- "/mir_target_enrich"
  call <- paste(base, endpoint, sep="")
  #print(call)
  
  saveRDS(toSend, "tosend.rds")
  request <- httr::POST(url = call, 
                        body = list(rds = upload_file("tosend.rds", "application/octet-stream")))
  
  # check if successful
  if(request$status_code != 200){
    current.msg <<- c("Failed to connect to Xia Lab API Server!")
    return(0)
  }
  
  # now process return
  request <- httr::content(request, "raw")
  request <- unserialize(request)
  
  if(is.null(request$json.nm)){
    current.msg <<- c("Error! Enrichment analysis via api unsuccessful!")
    return(0)
  }else{
    
    # create json
    sink(request$json.nm)
    cat(request$json.mat);
    sink();
    
    # write csv
    fast.write.csv(request$resTable, file="mirnet_enrichment.csv", row.names=F);

    # update dataSet
    dataSet <<- request$dataSet
  }
  
  current.msg <<- "Functional enrichment analysis was completed!"
  return(1)
}


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
