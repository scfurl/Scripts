
writeMasterNetworkFile<-function (adjMat, masterFile=NULL, weighted = TRUE, 
                            threshold = 0.02, nodeNames = NULL, moduleNames = NULL, nodeAttr = NULL, 
                            includeColNames = TRUE, writeDEData=FALSE, expressionset=NULL, classvec=NULL, comparison=NULL) 
{
  adjMat = as.matrix(adjMat)
  adjMat[is.na(adjMat)] = 0
  nRow = nrow(adjMat)
  checkAdjMat(adjMat, min = -1, max = 1)
  if (is.null(nodeNames)) 
    nodeNames = dimnames(adjMat)[[1]]
  if (is.null(nodeNames)) 
    stop("Cannot determine node names: nodeNames is NULL and adjMat has no dimnames.")
  rowMat = matrix(c(1:nRow), nRow, nRow, byrow = TRUE)
  colMat = matrix(c(1:nRow), nRow, nRow)
  adjDst = as.dist(adjMat)
  dstRows = as.dist(rowMat)
  dstCols = as.dist(colMat)
  edges = abs(adjDst) > threshold
  nEdges = sum(edges)
  edgeData = data.frame(fromNode = nodeNames[dstRows[edges]], 
                        toNode = nodeNames[dstCols[edges]], weight = if (weighted) 
                          adjDst[edges]
                        else rep(1, nEdges), direction = rep("undirected", nEdges))
  nodesPresent = rep(FALSE, ncol(adjMat))
  nodesPresent[dstRows[edges]] = TRUE
  nodesPresent[dstCols[edges]] = TRUE
  nNodes = sum(nodesPresent)
  if(writeDEData==TRUE){
    #check for expression set and classvec
    if(class(expressionset)!="ExpressionSet"){stop("ExpressionSet not correct")}
    if(class(classvec)!="factor"){stop("Classvec not correct")}
    if(class(comparison)!="numeric"){stop("comparison not correct")}
    mean1<-exprs(data)[nodeNames,classvec %in% levels(classvec)[comparison][1]]
    mean2<-exprs(data)[nodeNames,classvec %in% levels(classvec)[comparison][2]]
    FC<-2^(apply(mean1,1,mean)-apply(mean2,1,mean))
    edgeData$FC<-FC[match(edgeData$fromNode, names(FC))]
  }
  if (!is.null(moduleNames)){
    modLookup<-data.frame(nodes=nodeNames, modules=moduleNames)
    edgeData$modules<-modLookup$modules[match(edgeData$fromNode, modLookup$nodes)]
  }
  if (!is.null(masterFile)) {
    write.table(edgeData, file = masterFile, quote = FALSE, 
                row.names = FALSE, col.names = includeColNames)
  }
}

exportNetwork.SF<-function (adjMat, edgeFile = NULL, nodeFile = NULL, weighted = TRUE, 
          threshold = 0.5, nodeNames = NULL, altNodeNames = NULL, nodeAttr = NULL, 
          includeColNames = TRUE, writeDEData=FALSE, expressionset=NULL, classvec=NULL, 
          comparison=NULL, abbreviated=TRUE, directory=getwd()) 
{
  adjMat = as.matrix(adjMat)
  adjMat[is.na(adjMat)] = 0
  nRow = nrow(adjMat)
  checkAdjMat(adjMat, min = -1, max = 1)
  if (is.null(nodeNames)) 
    nodeNames = dimnames(adjMat)[[1]]
  if (is.null(nodeNames)) 
    stop("Cannot determine node names: nodeNames is NULL and adjMat has no dimnames.")
  rowMat = matrix(c(1:nRow), nRow, nRow, byrow = TRUE)
  colMat = matrix(c(1:nRow), nRow, nRow)
  adjDst = as.dist(adjMat)
  dstRows = as.dist(rowMat)
  dstCols = as.dist(colMat)
  edges = abs(adjDst) > threshold
  nEdges = sum(edges)
  edgeData = data.frame(fromNode = nodeNames[dstRows[edges]], 
                        toNode = nodeNames[dstCols[edges]], weight = if (weighted) 
                          adjDst[edges]
                        else rep(1, nEdges), direction = rep("undirected", nEdges), 
                        fromAltName = if (is.null(altNodeNames)) 
                          rep("NA", nEdges)
                        else altNodeNames[dstRows[edges]], toAltName = if (is.null(altNodeNames)) 
                          rep("NA", nEdges)
                        else altNodeNames[dstCols[edges]])
  nodesPresent = rep(FALSE, ncol(adjMat))
  nodesPresent[dstRows[edges]] = TRUE
  nodesPresent[dstCols[edges]] = TRUE
  nNodes = sum(nodesPresent)
  nodeData = data.frame(nodeName = nodeNames[nodesPresent], 
                        altName = if (is.null(altNodeNames)) 
                          rep("NA", nNodes)
                        else altNodeNames[nodesPresent], nodeAttribute = if (is.null(nodeAttr)) 
                          rep("NA", nNodes)
                        else nodeAttr[nodesPresent])
  if(writeDEData==TRUE){
    #check for expression set and classvec
    if(class(expressionset)!="ExpressionSet"){stop("ExpressionSet not correct")}
    if(class(classvec)!="factor"){stop("Classvec not correct")}
    if(class(comparison)!="numeric"){stop("comparison not correct")}
    mean1<-exprs(data)[nodeNames,classvec %in% levels(classvec)[comparison][1]]
    mean2<-exprs(data)[nodeNames,classvec %in% levels(classvec)[comparison][2]]
    FC<-2^(apply(mean1,1,mean)-apply(mean2,1,mean))
    nodeData$Expression<-FC[match(nodeData$nodeName, names(FC))]
  }
  if(abbreviated==TRUE){
    if (!is.null(edgeFile)) 
      write.table(edgeData[,1:3], file = file.path(directory, edgeFile), quote = FALSE, sep="\t",
                  row.names = FALSE, col.names = includeColNames)
    if (!is.null(nodeFile)) 
      write.table(nodeData[,c(1,4)], file = file.path(directory, nodeFile), quote = FALSE,  sep="\t",
                  row.names = FALSE, col.names = includeColNames)
  }
  else{
  if (!is.null(edgeFile)) 
    write.table(edgeData, file = file.path(directory, edgeFile), quote = FALSE, sep="\t",
                row.names = FALSE, col.names = includeColNames)
  if (!is.null(nodeFile)) 
    write.table(nodeData, file = file.path(directory, nodeFile), quote = FALSE,  sep="\t",
                row.names = FALSE, col.names = includeColNames)}
  return(list(edgeData = edgeData, nodeData = nodeData))
}


breakup<-function(charvec, n){
  return(split(charvec, ceiling(seq_along(charvec)/n)))
}


geneCorPlot<-function(dataNorm, gene1, pvalue, perscreen=c(4,4), cols){
  ####This function takes gene1 and plots its correlation with a set of target genes = "genes" by pvalue threshold = pvalue; perscreen is the vector layout to pass to par to set mfrow (allow 1-5 for each)
  .pardefault <- par(no.readonly = T)
  temp<-apply(dataNorm, 1, function(x) cor.test(dataNorm[gene1,], x))
  p.value<-lapply(temp, function(x) x[3])
  p.value1<-unlist(p.value)
  genes<-rownames(dataNorm)[which(p.value1<pvalue)]
  par(mfrow = perscreen)
  par(cex = (1-(floor(mean(perscreen))/10)))
  par(mar = c(4, 4, 0.5, 0.5))
  numplots<-perscreen[1]*perscreen[2]
  for(i in 1:ceiling(length(genes)/numplots)){
    toplot<-print(breakup(genes, numplots)[[i]])
    for(j in 1:numplots){
      if(is.na(toplot[j])==FALSE){
        plot(x=dataNorm[toplot[j],], y=dataNorm[gene1,], ylab=gene1,xlab=toplot[j], col=cols)
      }
    }
  }
  par(.pardefault)
  return(genes)
}
