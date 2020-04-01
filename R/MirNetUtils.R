##################################################
## R scripts for miRNet
## Description: miRNA-gene network analysis methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net.type PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname CreateMirNets
#' @export 
CreateMirNets <- function(net.type){
  dataSet$mirnet <<- net.type;
  if(net.type == "multilist" || net.type == "snp2mir2gene" || net.type == "snp2mir2dis" || net.type == "snp2mir2mol" || net.type == "snp2mir2lnc" || net.type == "snp2mir2tf"){
    require('igraph');
    for(i in 1:length(dataSet$mirtable)){
        if(i == 1){
            dataSet$mir.res =  dataSet[dataSet$mirtable[i]][[1]];
        }else{
            dataSet$mir.res = rbind(dataSet$mir.res, dataSet[dataSet$mirtable[i]][[1]]);
        }   
    }
    my.nodes <- dataSet$mir.res[, c(1, 3)];
    nd.nms <- c(dataSet$mir.res[, 1], dataSet$mir.res[, 3]);
    nd.ids <- c(dataSet$mir.res[, 2], dataSet$mir.res[, 4]);
    names(nd.ids) <- nd.nms;
    dups <- duplicated(nd.ids); #note using unique will lose the names attribute
    dataSet$node.anot <<- nd.ids[!dups];
    mir.graph <- simplify(graph.data.frame(my.nodes, directed=FALSE));
  }else{
    if(data.type == "xeno.mir"){
      my.nodes <- dataSet$mir.res[, c("miRNA", "Gene")];
      nd.nms <- c(dataSet$mir.res[, "miRNA"], dataSet$mir.res[, "Gene"]);
      nd.ids <- c(dataSet$mir.res[, "Accession"], dataSet$mir.res[, "Entrez"]);
    }else {
      my.nodes <- dataSet$mir.res[, c(1, 3)];
      nd.nms <- c(dataSet$mir.res[, 1], dataSet$mir.res[, 3]);
      nd.ids <- c(dataSet$mir.res[, 2], dataSet$mir.res[, 4]);
    }
    names(nd.ids) <- nd.nms;
    dups <- duplicated(nd.ids); #note using unique will lose the names attribute
    dataSet$node.anot <<- nd.ids[!dups];
    library(igraph);
    mir.graph <- simplify(graph.data.frame(my.nodes, directed=FALSE));
if(length(dataSet$nodeNumbers) == 0){
    dataSet$nodeNumber = nrow(dataSet$mir.res);
    dataSet<<- dataSet
}
  }
  # add node expression value
  match.index <- match(V(mir.graph)$name, rownames(dataSet$mir.mapped));

  expr.vals <- dataSet$mir.mapped[match.index, 1];
  mir.graph <- set.vertex.attribute(mir.graph, "abundance", index = V(mir.graph), value = expr.vals);
  mir.graph <<- mir.graph;

  substats <- DecomposeMirGraph(net.type, mir.graph, 2);
  if(!is.null(substats)){
    mir.graph <<- mir.graph;
    mir.query <- nrow(dataSet$mir.mapped);
    #mir.query <- nrow(dataSet$mir.orig); #original query
    if(net.type == "multilist" || net.type == "snp2mir2gene" || net.type == "snp2mir2dis" || net.type == "snp2mir2mol" || net.type == "snp2mir2lnc" || net.type == "snp2mir2tf"){
      mir.count <- length(unique(my.nodes[,1]));#matched mir
      tgt.count <- length(unique(my.nodes[,2]));#matched target
    }else{
      mir.count <- length(unique(my.nodes[,1]));#matched mir
      tgt.count <- length(unique(my.nodes[,2]));#matched target
    }
    return(c(mir.query, mir.count, tgt.count, ecount(mir.graph), length(mir.nets), substats));
  }else{
    return(0);
  }
}

# decompose to individual connected subnetworks, discard very small ones (defined by minNodeNum)
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net.type PARAM_DESCRIPTION
#' @param gObj PARAM_DESCRIPTION
#' @param minNodeNum PARAM_DESCRIPTION, Default: 2
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname DecomposeMirGraph
#' @export 
DecomposeMirGraph <- function(net.type, gObj, minNodeNum = 2){
    comps <-decompose.graph(gObj, min.vertices=minNodeNum);

    if(length(comps) == 0){
        current.msg <<- paste("No connected nodes found after this filtering!");
        return(NULL);
    }

    # first get stats
    queries <- unique(dataSet$seeds);
    net.stats <- as.data.frame(matrix(0, ncol = 3, nrow = length(comps)));
    dataSet$query.nums <- vector()
    dataSet$type.nums <- vector()
    for(i in 1:length(comps)){
        g <- comps[[i]];
        if(vcount(g) > 0){
            res=GetQueriesNumByType(g)
            nodeNum = GetNodeNumByType(g)
            dataSet$type.nums = c(dataSet$type.nums, nodeNum)
            dataSet$query.nums = c(dataSet$query.nums, res)
            net.stats[i,] <- c(
                               vcount(g),
                               ecount(g),
                               sum(queries %in% V(g)$name)
                              );
        }
    }

    # now sort graph based on node size and add names
    ord.inx <- order(net.stats[,1], decreasing=TRUE);
    net.stats <- net.stats[ord.inx,];
    comps <- comps[ord.inx];
    names(comps) <- rownames(net.stats) <- paste("mirnet", 1:length(comps), sep="");

    net.stats <- cbind(rownames(net.stats), net.stats);
    colnames(net.stats) <- c("Name", "Node", "Edge", "Query");

    # note, we report stats for all nets (at least 2 nodes);
    # but only contruct at least min node
    hit.inx <- net.stats$Node >= minNodeNum;
    comps <- comps[hit.inx];
    sub.stats <- NULL;
    json.res <- rep(list(list()), length(comps));
    i <- 0;
    for(nm in names(comps)){
        sub.stats <- c(sub.stats, vcount(comps[[nm]]));
    }

    # now save the components
    mir.nets <<- comps;
    net.stats <<- net.stats[,-1];  # remove the first name col

    # update the mir.res edge table
    # both side of the edge must present in all.nodes
    all.nodes <- V(gObj)$name;
    if(net.type == "multilist" || net.type == "snp2mir2gene"  || net.type == "snp2mir2dis" || net.type == "snp2mir2mol" || net.type == "snp2mir2lnc" || net.type == "snp2mir2tf"){
      hit.inx <- (dataSet$mir.res[,1] %in% all.nodes) & (dataSet$mir.res[,3] %in% all.nodes);
    }else{
      if(data.type == "xeno.mir"){
        hit.inx <- (dataSet$mir.res[, "miRNA"] %in% all.nodes) & (dataSet$mir.res[, "Gene"] %in% all.nodes);
      }else{
        hit.inx <- (dataSet$mir.res[,1] %in% all.nodes) & (dataSet$mir.res[,3] %in% all.nodes);
      }
    }
    dataSet$mir.filtered <- dataSet$mir.res[hit.inx, ];
    dataSet<<-dataSet
    return(sub.stats);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net.type PARAM_DESCRIPTION
#' @param nd.type PARAM_DESCRIPTION, Default: 'all'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ReduceEdgeDensity
#' @export 
ReduceEdgeDensity <- function(net.type, nd.type="all"){
    all.nms <- V(mir.graph)$name;
    edge.mat <- get.edgelist(mir.graph);
    dgrs <- degree(mir.graph);
    nodes2rm <- NULL;

    set.seed(8574);
    if(length(all.nms) > 50){
        # only get top 50 with highest density (degree)
        inx <- rank(-dgrs) < 50;
        seed.vec <- all.nms[inx];
    }else{
        seed.vec <- all.nms;
    }
    paths.list <-list();
    # now calculate the shortest paths only between these densely connected nodes
    for(pos in 1:length(seed.vec)){
        paths.list[[pos]] <- get.shortest.paths(mir.graph, seed.vec[pos], seed.vec[-pos])$vpath;
    }
    nds.inxs <- unique(unlist(paths.list));
    nodes2rm <- all.nms[-nds.inxs];

    # keep queries
    if(nd.type == "mir"){ # only apply removing to miRNA nodes
        mir.nms <- unique(edge.mat[,1]);
        nodes2rm <- nodes2rm[nodes2rm %in% mir.nms];
    }else if(nd.type=="other"){
        my.nms <- unique(edge.mat[,2]);
        nodes2rm <- nodes2rm[nodes2rm %in% my.nms];
    }else{
        #nothing to do
    }
    path.list <- NULL; gc();
    nodes2rm <- unique(nodes2rm);
    mir.graph <- simplify(delete.vertices(mir.graph, nodes2rm));
    current.msg <<- paste("A total of", length(nodes2rm) , "was reduced.");
    substats <- DecomposeMirGraph(net.type, mir.graph, 2);
    if(!is.null(substats)){
        mir.graph <<- mir.graph;
#CHECK
        return(c(nrow(dataSet$mir.orig), length(unique(dataSet$mir.filtered[,1])), length(unique(dataSet$mir.filtered[,3])), ecount(mir.graph), length(mir.nets), substats));
    }else{
        return(0);
    }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net.type PARAM_DESCRIPTION
#' @param nd.type PARAM_DESCRIPTION
#' @param min.dgr PARAM_DESCRIPTION
#' @param min.btw PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname FilterMirNet
#' @export 
FilterMirNet <- function(net.type, nd.type, min.dgr, min.btw){
    all.nms <- V(mir.graph)$name;
    edge.mat <- get.edgelist(mir.graph);
    dgrs <- degree(mir.graph);
    nodes2rm.dgr <- nodes2rm.btw <- NULL;

    if(nd.type == "mir"){
        mir.nms <- unique(edge.mat[,1]);
        hit.inx <- all.nms %in% mir.nms;
    }else if(nd.type=="other"){
        my.nms <- unique(edge.mat[,2]);
        hit.inx <- all.nms %in% my.nms;
    }else{ # all
        hit.inx <- rep(TRUE, length(all.nms));
    }

    if(min.dgr > 0){
        rm.inx <- dgrs <= min.dgr & hit.inx;
        nodes2rm.dgr <- V(mir.graph)$name[rm.inx];
    }
    if(min.btw > 0){
        btws <- betweenness(mir.graph);
        rm.inx <- btws <= min.btw & hit.inx;
        nodes2rm.btw <- V(mir.graph)$name[rm.inx];
    }

    nodes2rm <- unique(c(nodes2rm.dgr, nodes2rm.btw));
    mir.graph <- simplify(delete.vertices(mir.graph, nodes2rm));
    current.msg <<- paste("A total of", length(nodes2rm) , "was reduced.");
    substats <- DecomposeMirGraph(net.type, mir.graph, 2);
    if(!is.null(substats)){
        mir.graph <<- mir.graph;
        return(c(nrow(dataSet$mir.orig), length(unique(dataSet$mir.filtered[,1])), length(unique(dataSet$mir.filtered[,3])), ecount(mir.graph), length(mir.nets), substats));
    }else{
        return(0);
    }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net.type PARAM_DESCRIPTION
#' @param ids PARAM_DESCRIPTION
#' @param id.type PARAM_DESCRIPTION
#' @param remove PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname FilterMirNetByList
#' @export 
FilterMirNetByList <- function(net.type, ids, id.type, remove){

    lines <- strsplit(ids, "\r|\n|\r\n")[[1]];
    lines<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", lines, perl=TRUE);
    nms.vec = unique(unlist(dataSet$mir.res[, c(1,2,3,4)]))
    # need to first convert to correct id used in the graph
    hit.inx <- nms.vec %in% lines;
    nodes2rm = nms.vec[hit.inx]
    
    if(remove== "true"){
    nodes2rm <- nodes2rm[nodes2rm %in% V(mir.graph)$name];    # make sure they are in the igraph object
    }else{
    nodes2rm <- V(mir.graph)$name[!(V(mir.graph)$name %in% nodes2rm)];    # make sure they are in the igraph object
    }
    mir.graph <- simplify(delete.vertices(mir.graph, nodes2rm));
    current.msg <<- paste("A total of", length(nodes2rm) , "was reduced.");
    substats <- DecomposeMirGraph(net.type, mir.graph, 2);
    if(!is.null(substats)){
        mir.graph <<- mir.graph;
        return(c(nrow(dataSet$mir.orig), length(unique(dataSet$mir.filtered[,1])), length(unique(dataSet$mir.filtered[,3])), ecount(mir.graph), length(mir.nets), substats));
    }else{
        return(0);
    }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param g PARAM_DESCRIPTION
#' @param filenm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname convertIgraph2JSON
#' @export 
convertIgraph2JSON <- function(g, filenm){
    g <<- g
    nms <- V(g)$name;
    lbls <- as.character(dataSet$node.anot[nms]);

    # get layers
    if(anal.type == "multilist" || anal.type == "snp2mir"){
      my.nodes <- dataSet$mir.res[, c(1, 3)];
    } else if(anal.type == "snp2mir"){
      my.nodes <- dataSet$mir.res[,c(2, 6, 8)];
    } else{
      my.nodes <- dataSet$mir.res[, c(1, 3)];
    }

    m <- as.matrix(my.nodes);
    layers = ceiling(match(V(g)$name, m)/nrow(m));

    # setup shape (mir square, gene circle)
    shapes <- rep("circle", length(nms));

    # get edge data
    edge.mat <- get.edgelist(g);
    edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], type=rep("arrow", nrow(edge.mat)));

    # now get coords
    pos.xy <- PerformLayOut(g, layers, "Default");


    node.btw <- as.numeric(betweenness(g));
    node.dgr <- as.numeric(degree(g));

    if(anal.type %notin% c("array", "rnaseq", "qpcr")){
      node.exp <- as.character(get.vertex.attribute(g, name="abundance", index = V(g)));
    }else{
      node.exp <- as.numeric(get.vertex.attribute(g, name="abundance", index = V(g)));
    }

    if(vcount(g) > 1000){
        minSize = 2;
    }else if(vcount(g) > 300){
        minSize = 2;
    }else{
        minSize = 3;
    }

    # node size to betweenness values
    node.sizes <- as.numeric(rescale2NewRange(node.btw, minSize, 8));

    # update mir node size, shape and color
    if(anal.type == "multilist" || anal.type == "snp2mir"){
    mir.inx <- nms %in% mir.nmsu;
    }else{
    mir.inx <- nms %in% edge.mat[,2];
    }
    shapes[mir.inx] <- "square";

    # update snp node size, shape and color
    snp.inx <- as.vector(sapply(nms, function(x) substr(x, 1, 2) == "rs"));
    shapes[snp.inx] <- "diamond";

    #if(anal.type == "multilist"){
    # update disease node shape
    dis.inx <- nms %in% net.info$dis.nms;
    shapes[dis.inx] <- "equilateral";
    circ.inx <- nms %in% net.info$circ.nms
    pseudo.inx <- nms %in% net.info$pseudo.nms
    lnc.inx <- nms %in% net.info$lnc.nms;
    tf.inx <- nms %in% net.info$tf.nms;
    epi.inx <- nms %in% net.info$epi.nms;
    snc.inx <- nms %in% net.info$snc.nms;
    mol.inx <- nms %in% net.info$mol.nms;
    gene.inx <- nms %in% net.info$gene.nms;
    # extract other molecule index
    #}
    # slightly highlight mir node in general

    if(substring(anal.type, 1,3) == "mir"){ #highlight miRNA if they are the query
      node.sizes[mir.inx] <- node.sizes[mir.inx] + 1;
    }else{
      node.sizes[mir.inx] <- node.sizes[mir.inx] + 0.4;
    }
    node.types <- rep("Gene", length(node.dgr));
    node.types[mir.inx] <- "miRNA";
    node.types[snp.inx] <- "SNP";
    node.types[lnc.inx] <- "lncRNA"; 
    node.types[snc.inx] <- "sncRNA";
    node.types[circ.inx] <- "circRNA";
    node.types[dis.inx] <- "Disease"; 
    node.types[tf.inx] <- "TF";
    node.types[mol.inx] <- "Compound"; 
    node.types[epi.inx] <- "Epigenetic";
    node.types[pseudo.inx] <- "Pseudogene";

    node.cols = rep("#ff4500", length(node.dgr));
    ntype = unique(node.types)
    color.vec = gg_color_hue(length(ntype))
    for(i in 1:length(ntype)){
        if(ntype[i] != "Gene"){
            node.cols[which(node.types ==ntype[i])]=color.vec[i]
        }
    }

    node.cols[mir.inx] <- "#306EFF"; # dark blue
    # update mir node color
    topo.colsw <- node.cols;
    node.cols[mir.inx] <- "#98F5FF";
    topo.colsb <- node.cols;


    freq = table(node.types)

    duplicated.types=node.types
    for(i in 1:length(unique(node.types))){
        duplicated.types[duplicated.types == names(freq[i])]=order(freq)[i]
    }

    if((length(freq) == 3 && "Gene" %in% names(freq) && "TF" %in% names(freq) && "miRNA" %in% names(freq))){
    V(g)$layers = as.numeric(duplicated.types);
    duplicated.types[node.types == "Gene"] = 2
    duplicated.types[node.types == "TF"] = 3
    duplicated.types[node.types == "miRNA"] = 1
    }else if(length(freq) == 3 && "SNP" %in% names(freq) && "Gene" %in% names(freq)&& "miRNA" %in% names(freq)){
      V(g)$layers = as.numeric(duplicated.types);
     duplicated.types[node.types == "miRNA"] = 2
      duplicated.types[node.types == "SNP"] = 1
      duplicated.types[node.types == "Gene"] = 3
    }else if(length(freq) == 3 && "lncRNA" %in% names(freq) && "Gene" %in% names(freq)&& "miRNA" %in% names(freq)){
      duplicated.types[node.types == "miRNA"] = 2
      duplicated.types[node.types == "lncRNA"] = 1
      duplicated.types[node.types == "Gene"] = 3
      V(g)$layers = as.numeric(duplicated.types);
    }else{
    duplicated.types = as.numeric(duplicated.types)
    V(g)$layers = duplicated.types
    }

    V(g)$group = as.numeric(duplicated.types); #concentric circle


    node.cols[mir.inx] <- "#306EFF"; # dark blue
    # update mir node color
    topo.colsw <- node.cols;
    node.cols[mir.inx] <- "#98F5FF";
    topo.colsb <- node.cols;
    # color based on expression
    bad.inx <- is.na(node.exp) | node.exp==0;
    if(!all(bad.inx)){
        exp.val <- node.exp;
        node.colsb.exp <- getExpColors(node.exp, c("#78ff4d", "#FA8072", "#ebebeb"));
        node.colsw.exp <- getExpColors(node.exp, c("#269b00", "#b30000", "#333333"));
        node.colsb.exp[bad.inx] <- "#d3d3d3";
        node.colsw.exp[bad.inx] <- "#c6c6c6";
    }else{
        node.colsb.exp <- rep("#d3d3d3",length(node.exp));
        node.colsw.exp <- rep("#c6c6c6",length(node.exp));
    }

    seed.inx <- nms %in% unique(dataSet$seeds);
    seed_arr <- rep("notSeed",length(node.dgr));
    seed_arr[seed.inx] <- "seed";

    # now create the json object
    nodes <- vector(mode="list");
    for(i in 1:length(node.sizes)){
        nodes[[i]] <- list(
                  id=nms[i],
                  size=node.sizes[i],
                  molType=node.types[i],
                  type=shapes[i],
                  seedArr =seed_arr[i],
                  url=lbls[i],
                  colorb=topo.colsb[i],
                  colorw=topo.colsw[i],
                  x=pos.xy[i,1],
                  y=pos.xy[i,2],
                  attributes=list(
                    expr = node.exp[i],
                    expcolb=node.colsb.exp[i],
                    expcolw=node.colsw.exp[i],
                    degree=node.dgr[i], # actual degree in orginal network
                    between=node.btw[i])
                );
    }

    current.mirnet <<- g
    # save node table
    nd.tbl <- data.frame(Id=nms, Label=nms, Degree=node.dgr, Betweenness=node.btw);
    fileNm <- paste("node_table_", substring(filenm, 0, nchar(filenm)-5), ".csv", sep="")
    write.csv(nd.tbl, file=fileNm, row.names=FALSE);

    # covert to json
    library(RJSONIO);
    netData <- list(mirnet=dataSet$mirnet, mirtarget=dataSet$mirtarget, organism=dataSet$org, nodes=nodes, edges=edge.mat);
    sink(filenm);
    cat(toJSON(netData));
    sink();

    # also save to GraphML
    write.graph(g, file="mirnet.graphml", format="graphml");
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net.nm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PrepareGraphML
#' @export 
PrepareGraphML <- function(net.nm){
    write.graph(mir.nets[[net.nm]], file=paste(net.nm, ".graphml", sep=""), format="graphml");
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param table.nm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PrepareCSV
#' @export 
PrepareCSV <- function(table.nm){
  if(anal.type == "multilist" || anal.type == "snp2mir" || anal.type == "tf2genemir" || anal.type == "gene2tfmir"){
    write.csv(dataSet[[table.nm]], file=paste(table.nm, ".csv", sep=""), row.names = FALSE);
  } else {
    write.csv(dataSet$mir.res, file=paste(dataSet$mirnet, ".csv", sep=""), row.names = FALSE);
  }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mir.nm PARAM_DESCRIPTION
#' @param file.nm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PrepareMirNet
#' @export 
PrepareMirNet <- function(mir.nm, file.nm){
   my.mirnet <- mir.nets[[mir.nm]];
   current.mirnet <<- my.mirnet;
   convertIgraph2JSON(my.mirnet, file.nm);
   return(1);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param g PARAM_DESCRIPTION
#' @param layers PARAM_DESCRIPTION
#' @param algo PARAM_DESCRIPTION
#' @param focus PARAM_DESCRIPTION, Default: ''
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformLayOut
#' @export 
PerformLayOut <- function(g, layers, algo, focus=""){
    vc <- vcount(g);
    if(algo == "Default"){
        if(vc > 1000) {
           # pos.xy <- layout.fruchterman.reingold(g, area=30*vc^2);
            pos.xy <- layout.lgl(g);
        }else if(vc < 100){
            pos.xy <- layout.kamada.kawai(g);
        }else{
            pos.xy <- layout.fruchterman.reingold(g, area=40*vc^2);
        }
    }else if(algo == "FrR"){
        pos.xy <- layout.fruchterman.reingold(g, area=34*vc^2);
    }else if(algo == "circle"){
        pos.xy <- layout.circle(g);
    }else if(algo == "random"){
        pos.xy <- layout.random(g);
    }else if(algo == "lgl"){
        pos.xy <- layout.lgl(g);
    }else if(algo == "gopt"){
        pos.xy <- layout.graphopt(g)
    }else if(algo == "circular_tripartite"){
library(ggforce)
      l <- layout_with_sugiyama(g, layers = V(g)$group*(vc/3) +30)
      layout <- l$layout

  radial <- radial_trans(
    r.range = rev(range(layout[,2])),
    a.range = range(layout[,1]),
    offset = 0
  )
  coords <- radial$transform(layout[,2], layout[,1])
  layout[,1] <- coords$x
  layout[,2] <- coords$y
  pos.xy= layout
    }else if(algo == "tripartite"){
      l <- layout_with_sugiyama(g, layers = V(g)$layers*(vc/4))
      pos.xy <- -l$layout[,2:1]
    }else if(algo == "concentric"){
      library(graphlayouts)
      # the fist element in the list for concentric is the central node.
      if(focus==""){
      inx=1;
      }else{
      inx = which(V(g)$name == focus)
      }
      coords <- layout_with_focus(g,inx)
      pos.xy <- coords$xy
    }else if(algo == "backbone"){
    library(graphlayouts)
    if(length(V(g)$name)<2000){
        coords = layout_with_stress(g)
        pos.xy = coords
    }else{
        coords = layout_with_sparse_stress(g,pivots=100)
        pos.xy = coords
    }
    
    }else if(algo == "mds"){
    library(graphlayouts)
    coords = layout_with_pmds(g,length(V(g)$name)/10)
    pos.xy = coords/100
    rownames(pos.xy) = NULL
    }
    pos.xy;
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param algo PARAM_DESCRIPTION
#' @param filenm PARAM_DESCRIPTION
#' @param focus PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname UpdateNetworkLayout
#' @export 
UpdateNetworkLayout <- function(algo, filenm, focus){
    # get layers
    if(anal.type == "multilist" || anal.type == "snp2mir"){
      my.nodes <- dataSet$mir.res[, c(1, 3)];
    } else if(anal.type == "snp2mir"){
      my.nodes <- dataSet$mir.res[,c(2, 6, 8)];
    } else{
        my.nodes <- dataSet$mir.res[, c(1, 3)];
      }
    m <- as.matrix(my.nodes);
    layers = ceiling(match(V(current.mirnet)$name, m)/nrow(m));

    pos.xy <- PerformLayOut(current.mirnet, layers, algo, focus);
    nms <- V(current.mirnet)$name;
    nodes <- vector(mode="list");
    for(i in 1:length(nms)){
        nodes[[i]] <- list(
                  id=nms[i],
                  x=pos.xy[i,1],
                  y=pos.xy[i,2]
                );
    }
    # now only save the node pos to json
    library(RJSONIO);
    netData <- list(nodes=nodes);
    sink(filenm);
    cat(toJSON(netData));
    sink();
    return(filenm);
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
#' @rdname GetNetNames
#' @export 
GetNetNames <- function(){
    rownames(net.stats);
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
#' @rdname GetTableNames
#' @export 
GetTableNames <- function(){
  if(anal.type == "multilist"|| anal.type == "gene"){
    dataSet$mirtable;
  }else if(anal.type == "gene2tfmir" ){
    dataSet$type
 }else if(anal.type == "snp2mir" || anal.type == "tf2genemir"){
    dataSet$mirtable;
  } else {
    dataSet$mirnet
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
#' @rdname GetSeedsColumn
#' @export 
GetSeedsColumn <- function(){
    tbls = dataSet$mirtable;
    vec = vector();
    for( i in 1:length(tbls)){
        nms = strsplit(tbls[i], "2")[[1]];
        orignms = nms
        nms= gsub("mir", "miRNA",nms)
        nms= gsub("tf", "TF",nms)
        nms= gsub("rna", "RNA",nms)
        nms= gsub("pseudo", "pseudo-gene",nms)
        nms= gsub("epi", "Epigenetic Modif.",nms)
        nms= gsub("dis", "Disease",nms)
        nms= gsub("snp", "SNP",nms)
        nms= gsub("mol", "Compounds",nms)
        nms= gsub("lnc", "lncRNA",nms)
        nms= gsub("snc", "sncRNA",nms)
        nms= gsub("xeno", "xeno-miRNA",nms)
        if(nms[1] == "protein"){
         vec[i]=paste0(nms[1],":" ,length(unique(dataSet[tbls[i]][[1]][,1])) + length(unique(dataSet[tbls[i]][[1]][,3])) )
        }else if(tbls[i] %in% c("tf2gene", "gene2tf") ){
        vec[i]=paste0(nms[1],":" ,length(unique(dataSet[tbls[i]][[1]][,1])),", ",nms[2],": ",length(unique(dataSet[tbls[i]][[1]][,3])))
        }else if(tbls[i] %in% c("tf2mir")){
        vec[i]=paste0(nms[1],":" ,length(unique(dataSet[tbls[i]][[1]][,3])),", ",nms[2],": ",length(unique(dataSet[tbls[i]][[1]][,1])))
        }else if(orignms[1] == "mir"){
        vec[i]=paste0(nms[1],":" ,length(unique(dataSet[tbls[i]][[1]][,1])),", ",nms[2],": ",length(unique(dataSet[tbls[i]][[1]][,3])))
        }else if(orignms[2] == "mir"){
        vec[i]=paste0(nms[1],":" ,length(unique(dataSet[tbls[i]][[1]][,3])),", ",nms[2],": ",length(unique(dataSet[tbls[i]][[1]][,1])))
        }else{
        vec[i]=paste0(nms[1],":" ,length(unique(dataSet[tbls[i]][[1]][,1])),", ",nms[2],": ",length(unique(dataSet[tbls[i]][[1]][,3])))
        }
}
    return(vec)
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
#' @rdname GetNetStats
#' @export 
GetNetStats <- function(){
    as.matrix(net.stats);
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
#' @rdname GetNetsNameString
#' @export 
GetNetsNameString <- function(){
    paste(rownames(net.stats), collapse="||");
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net.type PARAM_DESCRIPTION
#' @param max.len PARAM_DESCRIPTION, Default: 200
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetMinConnectedGraphs
#' @export 
GetMinConnectedGraphs <- function(net.type, max.len = 200){
  set.seed(8574);
  # first get shortest paths for all pair-wise seeds
  my.seeds <- unique(dataSet$seeds);
  sd.len <- length(my.seeds);
  paths.list <-list();

  # first trim mir.graph to remove no-seed nodes of degree 1
  dgrs <- degree(mir.graph);
  keep.inx <- dgrs > 1 | (names(dgrs) %in% my.seeds);
  nodes2rm <- V(mir.graph)$name[!keep.inx];
  mir.graph <-  simplify(delete.vertices(mir.graph, nodes2rm));

  # need to restrict the operation b/c get.shortest.paths is very time consuming
  # for top max.len highest degrees
  if(sd.len > max.len){
    hit.inx <- names(dgrs) %in% my.seeds;
    sd.dgrs <- dgrs[hit.inx];
    sd.dgrs <- rev(sort(sd.dgrs));
    # need to synchronize all (dataSet$seeds) and top seeds (my.seeds)
    dataSet$seeds <- names(sd.dgrs);
    my.seeds <- dataSet$seeds[1:max.len];
    sd.len <- max.len;
    current.msg <<- paste("The minimum connected network was computed using the top", sd.len, "seed nodes in the network based on their degrees.");
  }else{
    current.msg <<- paste("The minimum connected network was computed using all seed nodes in the network.");
  }
  # now calculate the shortest paths for
  # each seed vs. all other seeds (note, to remove pairs already calculated previously)
  for(pos in 1:sd.len){
    paths.list[[pos]] <- get.shortest.paths(mir.graph, my.seeds[pos], dataSet$seeds[-(1:pos)])$vpath;
  }
  nds.inxs <- unique(unlist(paths.list));
  nodes2rm <- V(mir.graph)$name[-nds.inxs];
  g <- simplify(delete.vertices(mir.graph, nodes2rm));

  nodeList <- get.data.frame(g, "vertices");
  nodeList <- as.data.frame(nodeList[, 1]);
  colnames(nodeList) <- c("ID");
  write.csv(nodeList, file="orig_node_list.csv", row.names=F, quote=F);

  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  write.csv(edgeList, file="orig_edge_list.csv", row.names=F, quote=F);

  path.list <- NULL;
  substats <- DecomposeMirGraph(net.type, g);
  if(!is.null(substats)){
    mir.graph <<- g;
    return(c(nrow(dataSet$mir.orig), length(unique(dataSet$mir.filtered[,1])), length(unique(dataSet$mir.filtered[,3])), ecount(mir.graph), length(mir.nets), substats));
  }else{
    return(0);
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
#' @rdname UpdateSubnetStats
#' @export 
UpdateSubnetStats <- function(){
    old.nms <- names(mir.nets);
    net.stats <- ComputeSubnetStats(mir.nets);
    ord.inx <- order(net.stats[,1], decreasing=TRUE);
    net.stats <- net.stats[ord.inx,];
    rownames(net.stats) <- old.nms[ord.inx];
    net.stats <<- net.stats;
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param comps PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ComputeSubnetStats
#' @export 
ComputeSubnetStats <- function(comps){

    net.stats <- as.data.frame(matrix(0, ncol = 3, nrow = length(comps)));
    colnames(net.stats) <- c("Node", "Edge", "Query");
    queries <- rownames(dataSet$mir.mapped);
    for(i in 1:length(comps)){
        g <- comps[[i]];
        net.stats[i,] <- c(vcount(g),ecount(g),sum(queries %in% V(g)$name));
    }
    return(net.stats);
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
#' @rdname GetQueryNum
#' @export 
GetQueryNum <-function(){
    return(dataSet$query.nums)
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
#' @rdname GetTypeNum
#' @export 
GetTypeNum <-function(){
    return(dataSet$type.nums)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param g PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetQueriesNumByType
#' @export 
GetQueriesNumByType <- function(g){
g<<-g
    inx <- unique(rownames(dataSet$mir.mapped)) %in% V(g)$name;
    queries <- unique(rownames(dataSet$mir.mapped))[inx];
if(length(dataSet$mirtable)>1){
vec <- dataSet$mirtable
vec <- vec[vec != "protein2protein"]
    vec = strsplit(vec, "2");
    vec = unlist(vec);
    vec = unique(vec)
    res = vector();
for( i in 1:length(vec)){
        if(vec[i] == "mir"){
        nms = mir.nmsu
        lbl = "miRNA"
        }else if(vec[i] == "circ"){
        nms = net.info$circ.nms
        lbl = "circRNA";
        }else if(vec[i] == "snc"){
        nms = net.info$snc.nms
        lbl = "sncRNA";
        }else if(vec[i] == "lnc"){
        nms = net.info$lnc.nms
        lbl = "lncRNA";
        }else if(vec[i] == "dis"){
        nms = net.info$dis.nms
        lbl = "Disease";
        }else if(vec[i] == "mol"){
        nms = net.info$mol.nms
        lbl = "Compound";
        }else if(vec[i] == "tf"){
        nms = net.info$tf.nms
        lbl = "TF";
        }else if(vec[i] == "gene"){
        nms = net.info$gene.nms
        lbl = "Gene";
        }else if(vec[i] == "epi"){
        nms = net.info$epi.nms
        lbl = "Epigenetic modif.";
        }else if(vec[i] == "snp"){
              nms = unique(c(dataSet$snp2mir$rsID, dataSet$snp2gene$rsID))
        lbl = "SNP";
        }else if(vec[i] == "xeno"){
              nms = unique(c(dataSet$xeno2gene[,3]))
        lbl = "Xeno-miRNA";
        }

    nms = unique(nms)
    if( sum(nms %in% queries)>0){
        res = paste0(lbl,": ", sum(nms %in% queries), "; ", res) 
    }
    }
    return(res)
}else{
    sum(queries %in% V(g)$name)
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
#' @rdname GetUnmappedNum
#' @export 
GetUnmappedNum <- function(){
    totalOrig=0
    allData = vector()
    if(length(dataSet$data)>0){
        for(i in 1:length(names(dataSet$data))){
            nm=names(dataSet$data)[i]
            totalOrig = totalOrig + nrow(dataSet$data[[nm]])
            allData = c(allData, as.vector(rownames(dataSet$data[[nm]])));
        }
    return(totalOrig-sum(unique(allData) %in% V(mir.graph)$name))
    }else{
    queries <- unique(rownames(dataSet$mir.mapped));
    return(length(queries)-sum(queries %in% V(mir.graph)$name))
}
}

# from to should be valid nodeIDs
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param from PARAM_DESCRIPTION
#' @param to PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetShortestPaths
#' @export 
GetShortestPaths <- function(from, to){

    paths <- get.all.shortest.paths(current.mirnet, from, to)$res;
    if(length(paths) == 0){
        return (paste("No connection between the two nodes!"));
    }

    path.vec <- vector(mode="character", length=length(paths));
    for(i in 1:length(paths)){
        path.inx <- paths[[i]];
        path.ids <- V(current.mirnet)$name[path.inx];
        #path.sybls <- V(current.mirnet)$Label[path.inx];
        path.sybls <- path.ids;
        pids <- paste(path.ids, collapse="->");
        psbls <- paste(path.sybls, collapse="->");
        path.vec[i] <- paste(c(pids, psbls), collapse=";")
    }

    if(length(path.vec) > 50){
        path.vec <- path.vec[1:50];
    }

    all.paths <- paste(path.vec, collapse="||");
    return(all.paths);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nodeids PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ExtractMirNetModule
#' @export 
ExtractMirNetModule<- function(nodeids){
    set.seed(8574);
    nodes <- strsplit(nodeids, ";")[[1]];
    g <- current.mirnet;

    # try to see if the nodes themselves are already connected
    hit.inx <- V(g)$name %in% nodes;
    gObj <- induced.subgraph(g, V(g)$name[hit.inx]);

    # now find connected components
    comps <-decompose.graph(gObj, min.vertices=1);

    if(length(comps) == 1){ # nodes are all connected
        g <- comps[[1]];
    }else{
        # extract modules
        paths.list <-list();
        sd.len <- length(nodes);
        for(pos in 1:sd.len){
            paths.list[[pos]] <- get.shortest.paths(g, nodes[pos], nodes[-(1:pos)])$vpath;
        }
        nds.inxs <- unique(unlist(paths.list));
        nodes2rm <- V(g)$name[-nds.inxs];
        g <- simplify(delete.vertices(g, nodes2rm));
    }
    nodeList <- get.data.frame(g, "vertices");
    if(nrow(nodeList) < 3){
        return ("NA");
    }
    module.count <- module.count + 1;
    module.nm <- paste("module", module.count, sep="");

    colnames(nodeList) <- c("Id", "Label");
    ndFileNm = "mirnet_node_list.csv";
    write.csv(nodeList, file=ndFileNm, row.names=F, quote=F);

    edgeList <- get.data.frame(g, "edges");
    edgeList <- cbind(rownames(edgeList), edgeList);
    colnames(edgeList) <- c("Id", "Source", "Target");
    edgFileNm = "mirnet_edge_list.csv";
    write.csv(edgeList, file=edgFileNm, row.names=F, quote=F);

    filenm <- paste(module.nm, ".json", sep="");
    convertIgraph2JSON(g, filenm);

    # record the module
    mir.nets[[module.nm]] <<- g;
    UpdateSubnetStats();
    module.count <<- module.count
    return (filenm);
}


# exclude nodes in current.net (networkview)
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nodeids PARAM_DESCRIPTION
#' @param filenm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ExcludeNodes
#' @export 
ExcludeNodes <- function(nodeids, filenm){
    nodes2rm <- strsplit(nodeids, ";")[[1]];
    current.mirnet <- delete.vertices(current.mirnet, nodes2rm);

    # need to remove all orphan nodes
    bad.vs<-V(current.mirnet)$name[degree(current.mirnet) == 0];
    current.mirnet <<- delete.vertices(current.mirnet, bad.vs);

    # return all those nodes that are removed
    nds2rm <- paste(c(bad.vs, nodes2rm), collapse="||");

    # update topo measures
    node.btw <- as.numeric(betweenness(current.mirnet));
    node.dgr <- as.numeric(degree(current.mirnet));

    if(anal.type %notin% c("array", "rnaseq", "qpcr")){
      node.exp <- as.character(get.vertex.attribute(current.mirnet, name="abundance", index = V(current.mirnet)));
    }else{
      node.exp <- as.numeric(get.vertex.attribute(current.mirnet, name="abundance", index = V(current.mirnet)));
    }
    nms <- V(current.mirnet)$name;
    nodes <- vector(mode="list");
    for(i in 1:length(nms)){
        nodes[[i]] <- list(
                  id=nms[i],
                  expr = node.exp[i],
                  degree=node.dgr[i],
                  between=node.btw[i],
                  expr = node.exp[i]
                );
    }
    # now only save the node pos to json
    library(RJSONIO);
    netData <- list(deletes=nds2rm,nodes=nodes);
    sink(filenm);
    cat(toJSON(netData));
    sink();
    return(filenm);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param g PARAM_DESCRIPTION
#' @param group PARAM_DESCRIPTION, Default: 1
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname layout_in_circles
#' @export 
layout_in_circles <- function(g, group=1) {
    layout <- lapply(split(V(g), group), function(x) {
        layout_in_circle(induced_subgraph(g,x))
    })
    layout <- Map(`*`, layout, seq_along(layout))
    x <- matrix(0, nrow=vcount(g), ncol=2)
    split(x, group) <- layout
    x
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
#' @rdname GetNetStatMatrix
#' @export 
GetNetStatMatrix <-function(){

  return(signif(as.matrix(res), 5));
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
#' @rdname GetNetStatColNames
#' @export 
GetNetStatColNames <-function(){
  return(signif(as.matrix(res), 5));
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
#' @rdname GetNetNms
#' @export 
GetNetNms <-function(){
return(rownames(dataSet$resTable));
}


# use microservice
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
#' @rdname ComputePCSFNet
#' @export 
ComputePCSFNet <- function(){
   .prepare.pcsf();
   .perform.computing();
   .save.pcsf();
}

# in public web, this is done by microservice
.perform.computing <- function(){
    dat.in <- readRDS("dat.in.rds");
    dat.in$my.res <- dat.in$my.fun();
    saveRDS(dat.in, file="dat.in.rds");
}

.prepare.pcsf <- function(){
    my.fun <- function(){
        require('PCSF');
        edg <- get.edgelist(dat.in$data);
        edg <- as.data.frame(edg);
        edg$V3 <- rep(1, nrow(edg));
        colnames(edg) <- c("from", "to", "cost");
        ppi <- construct_interactome(edg);
        g <- PCSF(ppi, dat.in$terminals, w = 5, b = 100, mu = 0.0005);
        return(g);
    }
    expr.vec = as.vector(dataSet$mir.mapped)
    inx = expr.vec == "*"
    expr.vec[inx] = 1;
    expr.vec = as.numeric(expr.vec)
    names(expr.vec) = rownames(dataSet$mir.mapped)
    dat.in <- list(data=mir.graph, terminals = expr.vec, my.fun = my.fun);
    saveRDS(dat.in, file="dat.in.rds");
    return(1)
}

.save.pcsf <- function(){

    dat.in <- readRDS("dat.in.rds");
    g <- dat.in$my.res;

    nodeList <- get.data.frame(g, "vertices");
    colnames(nodeList) <- c("Id", "Label");
    write.csv(nodeList, file="orig_node_list.csv", row.names=F, quote=F);

    edgeList <- get.data.frame(g, "edges");
    edgeList <- cbind(rownames(edgeList), edgeList);
    colnames(edgeList) <- c("Id", "Source", "Target");
    write.csv(edgeList, file="orig_edge_list.csv", row.names=F, quote=F);

    path.list <- NULL;
    substats <- DecomposeMirGraph(dataSet$mirnet, g);
    if(!is.null(substats)){
        mir.graph <<- g;
        return(c(nrow(dataSet$mir.orig), length(unique(dataSet$mir.filtered[,1])), nrow(nodeList), nrow(edgeList), length(mir.nets), substats));
    }else{
        return(0);
    }
}

# support walktrap, infomap and lab propagation
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param method PARAM_DESCRIPTION, Default: 'walktrap'
#' @param use.weight PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname FindCommunities
#' @export 
FindCommunities <- function(method="walktrap", use.weight=FALSE){
    library(igraph)
    # make sure this is the connected
    current.net <- current.mirnet
    g <- current.net;
    if(!is.connected(g)){
        g <- decompose.graph(current.net, min.vertices=2)[[1]];
    }
    total.size <- length(V(g));

    if(method == "walktrap"){
        fc <- walktrap.community(g);
    }else if(method == "infomap"){
        fc <- infomap.community(g);
    }else if(method == "labelprop"){
        fc <- label.propagation.community(g);
    }else{
        print(paste("Unknown method:", method));
        return ("NA||Unknown method!");
    }

    if(length(fc) == 0 || modularity(fc) == 0){
        return ("NA||No communities were detected!");
    }

    # only get communities
    communities <- communities(fc);
    community.vec <- vector(mode="character", length=length(communities));
    gene.community <- NULL;
    qnum.vec <- NULL;
    pval.vec <- NULL;
    rowcount <- 0;
    nms <- V(g)$name;
    #hit.inx <- match(nms, multi.net$node.data[,1]);
    #sybls <- multi.net$node.data[hit.inx,2];
    #names(sybls) <- V(g)$name;
    for(i in 1:length(communities)){
        # update for igraph 1.0.1
        path.ids <- communities[[i]];
        psize <- length(path.ids);
        if(psize < 5){
            next; # ignore very small community
        }
        hits <- dataSet$seeds %in% path.ids;
        qnums <- sum(hits);
        if(qnums == 0){
            next; # ignor community containing no queries
        }

        rowcount <- rowcount + 1;
        pids <- paste(path.ids, collapse="->");
        ##path.sybls <- V(g)$name[path.inx];
        #path.sybls <- sybls[path.ids];
        com.mat <- cbind(path.ids, path.ids, rep(i, length(path.ids)));
        gene.community <- rbind(gene.community, com.mat);
        qnum.vec <- c(qnum.vec, qnums);

        # calculate p values (comparing in- out- degrees)
        #subgraph <- induced.subgraph(g, path.inx);
        subgraph <- induced.subgraph(g, path.ids);
        in.degrees <- degree(subgraph);
        #out.degrees <- degree(g, path.inx) - in.degrees;
        out.degrees <- degree(g, path.ids) - in.degrees;
        ppval <- wilcox.test(in.degrees, out.degrees)$p.value;
        ppval <- signif(ppval, 3);
        pval.vec <- c(pval.vec, ppval);

        # calculate community score
        community.vec[rowcount] <- paste(c(psize, qnums, ppval, pids), collapse=";");
    }

    ord.inx <- order(pval.vec, decreasing=F);
    community.vec <- community.vec[ord.inx];
    qnum.vec <- qnum.vec[ord.inx];
    ord.inx <- order(qnum.vec, decreasing=T);
    community.vec <- community.vec[ord.inx];

    all.communities <- paste(community.vec, collapse="||");
    colnames(gene.community) <- c("Id", "Label", "Module");
    write.csv(gene.community, file="module_table.csv", row.names=F);
    return(all.communities);
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param g PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetNodeNumByType
#' @export 
GetNodeNumByType <- function(g){
    queries <-V(g)$name;
if(length(dataSet$mirtable)>1){
vec <- dataSet$mirtable
vec <- vec[vec != "protein2protein"]
    vec = strsplit(vec, "2");
    vec = unlist(vec);
    vec = unique(vec)
    res = "";
for( i in 1:length(vec)){
        if(vec[i] == "mir"){
        nms = mir.nmsu
        lbl = "miRNA"
        }else if(vec[i] == "circ"){
        nms = net.info$circ.nms
        lbl = "circRNA";
        }else if(vec[i] == "snc"){
        nms = net.info$snc.nms
        lbl = "sncRNA";
        }else if(vec[i] == "lnc"){
        nms = net.info$lnc.nms
        lbl = "lncRNA";
        }else if(vec[i] == "dis"){
        nms = net.info$dis.nms
        lbl = "Disease";
        }else if(vec[i] == "mol"){
        nms = net.info$mol.nms
        lbl = "Compound";
        }else if(vec[i] == "tf"){
        nms = net.info$tf.nms
        lbl = "TF";
        }else if(vec[i] == "gene"){
        nms = net.info$gene.nms
        lbl = "Gene";
        }else if(vec[i] == "epi"){
        nms = net.info$epi.nms
        lbl = "Epigenetic modif.";
        }else if(vec[i] == "snp"){
              nms = unique(c(dataSet$snp2mir$rsID, dataSet$snp2gene$rsID))
        lbl = "SNP";
        }else if(vec[i] == "xeno"){
              nms = unique(c(dataSet$xeno2gene[,3]))
        lbl = "Xeno-miRNA";
        }
        nms = unique(nms)
    if( sum(nms %in% queries)>0 && !grepl(lbl, res)){
        res = paste0(lbl,": ", sum(nms %in% queries), "; ", res) 
    }
    }
    return(res)
}else{
    return(length(V(g)$name))
}
}
