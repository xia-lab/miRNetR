##################################################
## R scripts for miRNet
## Description: Data IO functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

#' Read Graph Data
#' @export
ReadGraphData <- function(fileName, fileType, org) {
  require("igraph");
  data.org<<- org;
  types_arr <<- "";

  fileTypeu <<- fileType;
  current.msg <<- NULL;

  if(endsWith(fileName, ".json")){
    fileType = "json"
  }else if(endsWith(fileName, ".graphml")){
    fileType = "graphml"
  }else if(endsWith(fileName, ".txt")){
    fileType = "txt"
  }else if(endsWith(fileName, ".sif")){
    fileType = "sif"
  }

  if(fileType == "graphml"){
    graphX = tryCatch({
      read_graph(fileName, format = "graphml")
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, finally = {

    })
  }else if(fileType == "sif"){
    graphX = tryCatch({
      read.sif(fileName, format="igraph", directed = FALSE)
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, finally = {

    })
  }else if(fileType == "txt"){
    df <- read.table(fileName, header=FALSE, stringsAsFactors = FALSE)
    df = as.matrix(df)
    graphX = tryCatch({
      graph_from_edgelist(df)
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, finally = {

    })
  }else if(fileType == "json"){
    require("RJSONIO");
    dat = fromJSON(fileName);
    dfn = unlist(dat$elements$nodes);
    conv = data.frame(id1=dfn[which(names(dfn)=='data.id')], name1=dfn[which(names(dfn)=='data.name')]);
    dfe = unlist(dat$elements$edges);
    dffe = data.frame(id1=dfe[which(names(dfe) == "data.source")], id2=dfe[which(names(dfe) == "data.target")]);
    dfint = merge(conv, dffe, by="id1");
    colnames(conv) = c("id2", "name2");
    df = merge(conv, dfint, by="id2");
    df = df[,c("name1", "name2")];
    df=as.matrix(df)

    graphX = tryCatch({
      graph_from_edgelist(df, directed=FALSE)
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, finally = {

    })
  }else{
    current.msg <<- "Unknown format, please make sure that the file is saved in the supported formats!";
    return(0)
  }

  if(!is_igraph(graphX)){
    current.msg <<- "Failed to parse your file, please make sure that the file is formatted correctly";
    return(0)
  }
  current.msg <<- "Sucessfully parsed your graph file!";
  print(current.msg);
  nms <- V(graphX)$name;
  if(length(nms)<1){
    nms <- V(graphX)$id;
    graphX = set_vertex_attr(graphX, "name", value=nms)
  }
  node.data = data.frame(nms, nms);
  seed.proteins <<- nms;
  seed.genes <<- seed.proteins;
  e=as_edgelist(graphX)
  edge.data= data.frame(Source=e[,1], Target=e[,2])
  seed.expr <<- rep(0, length(node.data));
  substats <- DecomposeMirGraph("graphfile",graphX);
  net.nm <- names(mir.nets)[1];
  net.nmu <<- net.nm;
  current.net.nm <<- net.nm;
  g <- mir.nets[[net.nm]];
  # ppi.net <<- list(db.type="abc",
  #                  db.type="ppi",
  #                  order=1,
  #                  seeds=nms,
  #                  table.nm=" ",
  #                  node.data = node.data,
  #                  edge.data = edge.data
  # );

  convertIgraph2JSONFromFile(net.nm, "mirnet_0.json");
  return(1);
}


read.sif <- function (sif.file, format = "graphNEL", directed = FALSE, header = TRUE, sep = "\t", ...) {

  net <- read.csv(file = sif.file, sep = sep, colClasses = "character", header = header, ...)

  # Assume form: node1 linktype node2 side.info..
  if ( ncol(net) > 2 ) {

    # remove NA nodes
    nas <- apply(net, 1, function (x) {any(is.na(x[c(1,3)]))})
    if (any(nas)) {
      net <- net[!nas, ]
      warning("NAs removed from network node list, ", sum(nas), " edges removed.")
    }

    net <- graph.edgelist(as.matrix(net[, -2]), directed = directed)

  } else if ( ncol(net) == 2 ) { # assume form: node1 node2

    # remove NA nodes
    nas <- apply(net, 1, function (x) {any(is.na(x))})
    if (any(nas)) {
      net <- net[!nas, ]
      warning("NAs removed from network node list, ", sum(nas), " edges removed.")
    }

    net <- graph.edgelist(cbind(net[,1],net[,2]), directed = directed)
  }

  if (format == "graphNEL") { net <- igraph.to.graphNEL(net) }
  # if (format == "igraph") { net <- igraph.from.graphNEL(igraph.to.graphNEL(net)) }

  net
}


#' Convert Igraph to JSON
#' @export
convertIgraph2JSONFromFile <- function(net.nm, filenm){
  g <- mir.nets[[net.nm]];

  # annotation
  nms <- V(g)$name;

  # get layers
  m <- as.matrix(nms);
  layers = ceiling(match(V(g)$name, m)/nrow(m));

  # setup shape (gene circle, other squares)
  shapes <- rep("circle", length(nms));

  # get edge data
  edge.mat <- as_edgelist(g);

  edge.mat1 = data.frame(edge.mat)
  edge.mat1$color = "target"
  edge.mat1 = as.matrix(edge.mat1)
  edge.color = rep("#d3d3d3",nrow(edge.mat));

  #edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], color = edge.mat1[,3]);
  edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], color=edge.color);

  # get the note data
  if(length(E(g)$weight)>0){
    E(g)$weight = abs(E(g)$weight)
    weight = E(g)$weight
    mir.nets[[net.nm]] <<- g
  }
  node.btw <- as.numeric(betweenness(g));
  node.clo <- as.numeric(closeness(g));
  #node.adh <- as.numeric(adhesion(g));
  node.eig <- eigen_centrality(g);
  node.eig = as.numeric(node.eig$vector);
  node.tra=transitivity(g,type=c("local"))

  node.dgr <- as.numeric(degree(g));
  node.exp <- as.numeric(vertex_attr(g, name="Expression", index = V(g)));

  if(length(node.exp) == 0){
    node.exp <- rep(0,length(node.dgr));
  }

  # node size to degree values
  if(vcount(g)>500){
    min.size = 1;
  }else if(vcount(g)>200){
    min.size = 2;
  }else{
    min.size = 3;
  }

  minval = min(node.dgr, na.rm=T);
  maxval = max(node.dgr, na.rm=T);
  result = maxval-minval;

  if(result == 0){
    node.sizes <- rep((log(node.dgr))^2, length(nms));
  }else{
    #node.sizes <- as.numeric(rescale2NewRange(node.btw, min.size, 8));
    #node.sizes <- as.numeric(rescale2NewRange((log(node.dgr))^2, min.size, 9))*3 +10;
    node.sizes <- as.numeric(rescale2NewRange((log(node.dgr))^2, min.size, 9)) +1;
  }

  centered = T;
  notcentered = F;
  # update node color based on betweenness
  require("RColorBrewer");
  topo.val <- log(node.btw+1);
  topo.colsb <- ComputeColorGradient(topo.val, "black", notcentered);
  topo.colsw <-  ComputeColorGradient(topo.val, "white", notcentered);

  # color based on expression
  bad.inx <- is.na(node.exp) | node.exp==0;
  if(!all(bad.inx)){
    exp.val <- node.exp;
    node.colsb.exp <- ComputeColorGradient(exp.val, "black", centered);
    node.colsw.exp <- ComputeColorGradient(exp.val, "white", centered);
    node.colsb.exp[bad.inx] <- "#d3d3d3";
    node.colsw.exp[bad.inx] <- "#c6c6c6";
    # node.colsw.exp[bad.inx] <- "#99ddff";
  }else{
    node.colsb.exp <- rep("#d3d3d3",length(node.exp));
    node.colsw.exp <- rep("#c6c6c6",length(node.exp));
  }

  node_attr = list.vertex.attributes(g);

  attr=list();
  for(j in 1:length(node_attr)){
    attr[[node_attr[j]]] = vertex_attr(g, node_attr[j])
  }
  attr_names = names(attr);
  attr_nd = list();
  arr = list()
  for(i in 1:length(node.sizes)){
    for(j in 1:length(attr)){
      attr_nd[node_attr[j]] = as.character(unlist(attr[node_attr[j]])[i])
    }
    arr[[i]] = attr_nd;
  }
  network_prop = list();
  for(i in 1:length(node.sizes)){
    network_prop[[i]]  <- list(
      closeness = node.clo[i],
      eigen = node.eig[i]
    )
  }
  lblsu <<- nms;
  # now get coords
  pos.xy <- PerformLayOut(g, layers, "Default");
  # now create the json object
  nodes <- vector(mode="list");
  for(i in 1:length(node.sizes)){
    nodes[[i]] <- list(
      id=nms[i],
      size=node.sizes[i],
      molType="Unspecified", # need to add details
      type=shapes[i], # need to add details
      seedArr ="seed",
      url=nms[i],
      x = pos.xy[i,1],
      y = pos.xy[i,2],
      colorb=topo.colsb[i],
      colorw=topo.colsw[i],
      #topocolb=topo.colsb[i],
      #topocolw=topo.colsw[i],
      #expcolb=node.colsb.exp[i],
      #expcolw=node.colsw.exp[i],
      #user=network_prop[[i]],
      attributes=list(
        expr = node.exp[i],
        degree=node.dgr[i],
        between=node.btw[i]
      )
    );
  }

  # save node table
  nd.tbl <- data.frame(Id=nms, Degree=node.dgr, Betweenness=round(node.btw,2));
  # order
  ord.inx <- order(nd.tbl[,2], nd.tbl[,3], decreasing =TRUE)
  nd.tbl <- nd.tbl[ord.inx, ];
  fast.write.csv(nd.tbl, file="node_table.csv", row.names=FALSE);
  # covert to json
  require(RJSONIO);
  # dg <- decompose(g)
  # if(length(dg)>1){
  #   modules = "NA"
  # }else{
  #   current.mirnet <- g;
  #   modules = FindCommunities("walktrap", FALSE);
  # }
  netData <- list(mirnet="graphfile", mirtarget="NA", organism=data.org, nodes=nodes, edges=edge.mat);
  sink(filenm);
  cat(toJSON(netData));
  sink();
}

getGraphStatsFromFile <- function(){
  g <- mir.nets[[net.nmu]];
  nms <- V(g)$name;
  edge.mat <- as_edgelist(g);
  return(c(length(nms), nrow(edge.mat)));
}
