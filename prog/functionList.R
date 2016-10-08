#OUTLIER TEST FUNCTION
outlierDendrogram = function(datamatrix,threshold){
# sample network based on squared Euclidean distance, data transposed:
A=adjacency(datamatrix,type="signed") #adjacency assumes columns correspond to genes, but here testing samples
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
#plot(Z.k)
# Designate samples as outlying (See Horvath 2011 Chapter 12)
# if their Z.k value is below the threshold
thresholdZ.k=threshold # often -2.5
# the color vector indicates outlyingness (red)
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
# calculate the cluster tree using flahsClust or hclust
sampleTree = flashClust(as.dist(1-A), method = "average")
plotDendroAndColors(sampleTree,colors=outlierColor,cex.dendroLabels=.5)
}

#Make data.frame numeric
makeNumeric <- function (df) {
rowname.df = rownames(df)
colname.df = colnames(df)
df = as.data.frame(df)
df.num = sapply(df,function(x) as.numeric(as.character(x)))
rownames(df.num) = rowname.df
colnames(df.num) = colname.df
return(df.num)
}

#Rescale and center by gene
center.and.scale.by.feature<-function(dataset){
  expr.mean=apply(dataset,1,mean) #Mean of each column
  expr.sd=apply(dataset,1,sd) #Standard deviation of each column
  expr.norm=matrix(NA,nrow(dataset),ncol(dataset))
  for(i in 1:nrow(dataset)){
    temp=(dataset[i,]-expr.mean[i])/expr.sd[i]
    expr.norm[i,]=temp
  }
  expr.norm=as.data.frame(expr.norm)
  expr.norm=sapply(expr.norm, function(x) as.numeric(as.character(x)))
  rownames(expr.norm)=rownames(dataset) #CHECK THIS IS RIGHT
  colnames(expr.norm)=colnames(dataset) #CHECK THIS IS RIGHT
  return(expr.norm)
}

#Rescale and center by sample
center.and.scale.by.sample<-function(dataset){
  expr.mean=apply(dataset,2,mean) #Mean of each column
  expr.sd=apply(dataset,2,sd) #Standard deviation of each column
  expr.norm=matrix(NA,nrow(dataset),ncol(dataset))
  for(i in 1:ncol(dataset)){
    temp=(dataset[,i]-expr.mean[i])/expr.sd[i]
    expr.norm[,i]=temp
  }
  expr.norm=as.data.frame(expr.norm)
  expr.norm=sapply(expr.norm, function(x) as.numeric(as.character(x)))
  rownames(expr.norm)=rownames(dataset) #CHECK THIS IS RIGHT
  colnames(expr.norm)=colnames(dataset) #CHECK THIS IS RIGHT
  return(expr.norm)
}

###################################################
#WGCNA supporting functions

#This section contains supporting functions for WGCNA

##extractModules: a function which uses the module assignment list as input and writes individual files with the probeset ids for each module
dir=NULL
extractModules<-function(colorh1,datExpr,anno,write=F,file_prefix="",dir=NULL)
{
  module<-list()
  if (!is.null(dir))
  {
    dir.create(dir)
    file_prefix=paste(dir,"/",file_prefix,sep="")
  }
  i<-1
  for (c in unique(colorh1))
  {
    module[[i]]<-(anno[colnames(datExpr)[which(colorh1==c)],1])
    if (write) {write.table(rownames(anno)[which(colorh1==c)],file=paste(file_prefix,"_",c,".txt",sep=""),quote=F,row.names=F,col.names=F)}
    i<-i+1
  }
  names(module)<-unique(colorh1)
  module
}


##EigenGenes : this is used by the plotting function to display close together similar modules based on their eigen values
getEigenGeneValues<-function(datRef,colorh1,datAll)
{
  eigenGenesCoef<-list()
  i<-0
  for (c in unique(colorh1))
  {
    i<-i+1
    eigenGenesCoef[[i]]<-prcomp(scale(datRef[,which(colorh1 == c)]))$rotation[,1]
  }
  names(eigenGenesCoef)<-unique(colorh1)
  values<-NULL
  for( c in unique(colorh1))
  {
    v<-rbind(datAll)[,which(colorh1 == c)] %*%  eigenGenesCoef[[c]]
    values<-cbind(values,sign(mean(v))*v)
  }
  colnames(values)<-unique(colorh1)
  values
}
####plotting function for comparative heatmap
plotC1C2Heatmap<-function(colorh1C1C2,AdjMat1C1,AdjMat1C2, datC1, datC2,ordering=NULL,file=paste("DifferentialPlot",name,".png",sep=""))
{
  if (is.null(ordering))
  {
    h<-hclust(as.dist(1-abs(cor(getEigenGeneValues(datC1[,which(colorh1C1C2!="grey")],colorh1C1C2[which(colorh1C1C2!="grey")],rbind(datC1,datC2)[,which(colorh1C1C2!="grey")])))))
    for (c in h$label[h$order])
    {
      ordering<-c(ordering,which(colorh1C1C2 ==c))
    }
  }
  mat_tmp<-(AdjMat1C1[ordering,ordering])
  mat_tmp[which(row(mat_tmp)>col(mat_tmp))]<-(AdjMat1C2[ordering,ordering][which(row(mat_tmp)>col(mat_tmp))])
  diag(mat_tmp)<-0
  mat_tmp<-sign(mat_tmp)*abs(mat_tmp)^(1/2)
  png(file=file,height=1000,width=1000)
  image(mat_tmp,col=rev(brewer.pal(11,"RdYlBu")),axes=F,asp=1,breaks=seq(-1,1,length.out=12))
  dev.off()
  unique(colorh1C1C2[ordering])
}

##This function plots side by side the color bar of module assignments, and the change in mean expression of the modules between the two conditions.
plotExprChange<-function(datC1,datC2, colorhC1C2,ordering=NULL)
{
  if (is.null(ordering))
  {
    h<-hclust(as.dist(1-abs(cor(getEigenGeneValues(datC1[,which(colorh1C1C2!="grey")],colorh1C1C2[which(colorh1C1C2!="grey")],rbind(datC1,datC2)[,which(colorh1C1C2!="grey")])))))
    for (c in h$label[h$order])
    {
      ordering<-c(ordering,which(colorh1C1C2 ==c))
    }
  }
  mycolors<-colorh1C1C2[ordering]
  plot(x=0:length(which(mycolors!="grey")),y=rep(1,length(which(mycolors!="grey"))+1),col="white",axes=F,xlab="",ylab="",ylim=c(0,1))
  rr=c(244,239,225,215,209,193,181,166,151,130,110)
  gg=c(228,204,174,160,146,117,94,58,44,45,45)
  bb=c(176,140,109,105,102,91,84,74,70,68,66)
  MyColours<-NULL
  for ( i in 1:11)
  {
    MyColours=c(MyColours,rgb(rr[i],gg[i],bb[i],maxColorValue=255)  )
  }
  exprDiff<-NULL
  l<-0
  for (c in setdiff(unique(mycolors),"grey"))
  {
    meanC1<-mean(t(datC1)[colnames(datC1)[which(colorh1C1C2 == c)],])
    meanC2<-mean(t(datC2)[colnames(datC2)[which(colorh1C1C2 == c)],])
    exprDiff<-rbind(exprDiff,c(meanC1,meanC2))
    r<-l+length(which(mycolors==c))
    rect(l,0.85,r,1,col=c,border=F)
    rect(l,0,r,.4,col=MyColours[floor(meanC2*2)-10],border="white",lwd=2)
    rect(l,0.4,r,.8,col=MyColours[floor(meanC1*2)-10],border="white",lwd=2)
    l<-r
  }
  exprDiff
}

#plotMatrix is a function used to make Additional File 2: Figure S1 plot displaying the # permutation results.
plotMatrix<-function(mat)
{
  mat[which(row(mat)>col(mat))]<-1001  
  image(mat,col=c(gray.colors(4),"white"),breaks=c(0,0.1,50,100,1000,1001),xaxt='n',yaxt='n',xlim=c(-0.2,1.2),ylim=c(-0.2,1.2),bty='n',asp=1)

  text(0:(nrow(mat)-1)/(nrow(mat)-1),1.1,rownames(mat),cex=1,col=rownames(mat))

  text(-0.15,0:(ncol(mat)-1)/(ncol(mat)-1),colnames(mat),cex=1,col=colnames(mat))

  text(apply(matrix(0:(nrow(mat)-1)/(nrow(mat)-1)),1,rep,ncol(mat)),rep(0:(ncol(mat)-1)/(ncol(mat)-1),nrow(mat)),as.numeric(t(mat)),col="white",cex=1.5)
}



