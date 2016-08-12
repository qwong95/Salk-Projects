###Dir Brain
#######################
#Load data from aba_tissueexpression.R
#LOAD DEVELOPING HUMAN
load("C:/Users/Quinn Wong/Documents/SalkProjects/abaresults.rda")
#or LOAD ADULT HUMAN
#load("~SalkProjects/abaresults_adult.rda")
#######################
# Load libraries
library(tiff) # for loading in tiffs
library(fields) # for imaging with legends (base R 'image' doesn't show legends)
library(RColorBrewer)

#Initialize Class
setClass(Class="Comp", representation(tissueExp1 = "numeric",tissueExp2 = "table",composite = "matrix",random.matrix = "data.frame"))

#Load private functions
LoadFiles <- function(i, Dir){
  file <- paste(Dir,i,".tif",sep="")
  suppressWarnings(region <- readTIFF(source=file)[,,1])
}
countTiss <- function(g){
  #Purpose: identify the tissues that the gene is expressed in
  a <- rowmeta[rowmeta$gene_symbol == g,"probeset_id"]
  a2 <- colSums(abatissues[match(a,rownames(abatissues)),])
  return(colnames(abatissues[,a2 > 0]))
}
TissueSummary <- function(genes){
  if(!is.null(ncol(na.exclude(abatissuesBygenes[genes,])))){
    tissueswithgenes <- colSums(na.exclude(abatissuesBygenes[genes,]))
  }else{
    tissueswithgenes <- na.exclude(abatissuesBygenes[genes,])
  }
  a <- table(c(as.character(colmeta$structure_acronym),as.character(colmeta$structure_acronym),names(tissueswithgenes)))
  a <- a[a==2]
  a[a==2] <- 0
  tissueswithgenes <- c(a, tissueswithgenes)
  tissueswithgenes <- tissueswithgenes[order(names(tissueswithgenes))]
  return((tissueswithgenes))
}
reColor <- function(i,tissueExp,dim,Abrev,Files,slice){
  region <- matrix(data=unlist(Files[Abrev==i]),nrow=dim[1],ncol=dim[2])
  converted <- as.character(na.exclude(conversion[as.character(conversion[,paste("panel",slices[slice],sep="")]) == as.character(Abrev[Abrev==i]),"acronym"]))
  converted_tissexp <- tissueExp[converted]
  newValue <- as.numeric(ifelse(test=length(converted_tissexp) > 0,yes=sum(unlist(converted_tissexp)),no=0.0001))
  region[region == 1] <- 0
  region[region != 0] <- newValue
  return(as.vector(region))
}
BrainMap <- function(dim,tissueExp,Abrev,Files,slice){
  tmp <- do.call("cbind",lapply(X=Abrev,FUN=reColor, tissueExp = tissueExp,dim=dim, Abrev = Abrev,Files = Files,slice=slice))
  tmp2 <- matrix(apply(X=tmp,MARGIN=1,FUN=sum),nrow=dim[1],ncol=dim[2])
  #tmp2[tmp2 < 1 & tmp2 != 0] <- 0.1
  return(tmp2)
}
RandomTissueSummary <- function(i, genes, samplesize){
  genes2 <- genes[(!is.na(rowmeta[match(genes,as.character(rowmeta$gene_symbol)),"probeset_id"]))]
  genes.1 <-  sample(x=genes2,size=samplesize ,replace=FALSE)
  tissueExp2 <- TissueSummary(genes.1)
  return(tissueExp2)
}
PlotBrain <- function(map,Breaks=9,colorname="Reds"){
  map <- as.matrix(map)
  map[!is.finite(map)] <- 0
  redGradient <- brewer.pal(n=round(Breaks/2),name="Reds")
  blueGradient <- rev(brewer.pal(n=round(Breaks/2),name="Blues"))

  fullGradient <- c(blueGradient,redGradient)
  map[map > 2 & map < max(map)] <- 2#Outlinefill] <- 2
  cutNormal.low <- seq(0.0002,2,length.out = Breaks)[1:(Breaks-1)]
  cutNormal.hi <- seq(0.0002,2,length.out = Breaks)[2:Breaks]
  
  normalized.color <- map
  for ( i in 1:length(cutNormal.low)){
    normalized.color[map >= cutNormal.low[i] & map <= cutNormal.hi[i]] <- fullGradient[i]
  }
  normalized.color[map == 1] <- "white"
  normalized.color[map == 0.0001] <- "#dddcdc"
  normalized.color[map == max(map)] <- "black"#Outlinefill] <- "black"
  grid.newpage()
  plot(1,1)
  grid.newpage()
  grid.raster(as.raster(normalized.color),interpolate=TRUE)
}
GetGenes <- function(genes,regions = 1,region.name = NA){
  normTissueExp <- tissueExp1 / tissueExp2
  normTissueExp <- normTissueExp[order(normTissueExp,decreasing=TRUE)]
  if(is.na(region.name)){
    tmp <- abatissuesBygenes[,names(normTissueExp)[1:regions]]
    if(regions > 1){
      tmp <- apply(tmp,1,mean)}
      tmp <- ifelse(tmp == 0, 0, 1)
      put.genes <- rownames(abatissuesBygenes[tmp ==1,])
      a <- table(c(put.genes,genes))
    return(names(a[a==2]))
  }else{
    Columns <- grep(pattern=region.name,x=colnames(abatissuesBygenes),ignore.case=TRUE)
    if(length(Columns) > 1){
      Test <- apply(abatissuesBygenes[,Columns],MARGIN=1,FUN=sum)
    }else{
      Test <- abatissuesBygenes[,Columns]
    }
    put.genes <- rownames(abatissuesBygenes[Test > 0,])
    a <- table(c(put.genes,genes))
    return(names(a[a==2]))
  }
}
PValue.onetail <- function(regions,tissueExp1,random.matrix){
  a <- na.exclude(random.matrix[regions,])
  b <- tissueExp1[regions]
  if(length(regions) ==1){
    return(sum(a > b) / length(a))
  }else{
    out <- vector()
    for (i in 1:nrow(a)){
      calc <- (b[i] / mean(as.numeric(na.exclude(a[i,]))))
      if ( !is.na(calc) & calc > 1){
        out <- c(out, sum(a[i,] > b[i]) / ncol(a))
      }else{
        out <- c(out, sum(a[i,] < b[i]) / ncol(a))
      }
    }
    out <- as.list(out)
    names(out) <- regions
    return(out)
  }
}
InABA <- function(genes){
  genes2 <- genes[(!is.na(rowmeta[match(genes,as.character(rowmeta$gene_symbol)),"probeset_id"]))]
  
}
SpatialEnrichment <- function(slice = 5, genes,background = NULL, reps = 100, backgroundOff = FALSE){
  #Get info about the slice
  Dir <-  DIR[slice] #pick which slice
  Abrev <- as.vector(unlist(abrev[Dir]))
  Fileloc <- unlist(fileloc[Dir])
  Files <- lapply(X=as.vector(Abrev),FUN=LoadFiles, Dir = Dir)
  dim <- as.vector(unlist(DIM[Dir]))
  #Get info about the probes
  probes <- rowmeta[match(genes,rowmeta$gene_symbol),"probeset_id"]
  surviving <- samplesize <- length(rowmeta[match(probes[!is.na(probes)],rowmeta$probeset_id),"gene_symbol"])
  message(paste("Identifying Tissues where query genes are expressed,",surviving,"of",length(genes),"genes queried are present in ABA dataset"))
  tissueExp1 <- tissueExp <- TissueSummary(genes)
  composite1 <- BrainMap(dim = dim ,tissueExp = tissueExp,  Abrev = Abrev,Files = Files,slice=slice)
  if(backgroundOff == FALSE){
    message(paste("Identifying Tissues where a random selection of",surviving,"genes are expressed using",reps , "iterations..."))
    if(!is.null(background)){
      message("Using user-provided set of background genes")
      starting.b <- length(background)
      probes.b <- rowmeta[match(background,rowmeta$gene_symbol),"probeset_id"]
      background <- rowmeta[match(probes.b[!is.na(probes.b)],rowmeta$probeset_id),"gene_symbol"]
      surviving.b <- length(background)
      message(paste(surviving.b,"of",starting.b,"background genes are present in ABA dataset"))
    }else{
      message("Using all genes in ABA microarray as background")
      background <- unique(as.character(rowmeta$gene_symbol))
    }
    random.matrix <- as.data.frame(sapply(X=c(1:reps),FUN=RandomTissueSummary, genes= background, samplesize = surviving))
    tissueExp2 <-  tissueExp <- as.table(apply(X=random.matrix,MARGIN=1,FUN=mean))
    composite2 <- BrainMap(dim = dim, tissueExp = tissueExp, Abrev = Abrev,Files = Files)
    #message("Calculating bootstrapped p-values...")
    #boot <- PValue.onetail(names(tissueExp1))
    message("Creating Spatially-enriched Brain Image...")
    normalized.1 <- (composite1+1)/(composite2+1)
    normalized <- normalized.1
    normalized[composite1 == 0.0001] <- 0.0001
    Outlinefill <- max(normalized.1[!is.na(normalized.1)]) + sd(normalized.1[!is.na(normalized.1)])
    normalized[outline[[Dir]] < 1] <- Outlinefill
    comp <- new(Class="Comp",
               tissueExp1 = tissueExp1,
               tissueExp2 =  tissueExp2,
               composite = normalized,
               random.matrix = random.matrix
               )
  }else{
    message(paste("Background correction turned off. Showing raw gene enrichment"))
    normalized <- normalized.1 <- composite1 #Scale?
    #normalized[composite1 == 0.0001] <- 0.0001
    Outlinefill <- max(normalized.1[!is.na(normalized.1)]) + sd(normalized.1[!is.na(normalized.1)])
    normalized[outline[[Dir]] < 1] <- Outlinefill
    comp <- new(Class="Comp",
                tissueExp1 = tissueExp1,
                tissueExp2 = table(NA),
                composite = normalized,
                random.matrix = data.frame()
                )
    
  }
  return(comp)
}
Boot <- function(tissueExp1, tissueExp2, random.matrix){
  boot <- PValue.onetail(names(tissueExp1),tissueExp1, random.matrix = random.matrix)
  boot2 <- data.frame(a=as.numeric(boot),b=as.numeric(boot),
                      count.sample = tissueExp1[names(boot)],
                      count.random = tissueExp2[names(boot)],
                      change = tissueExp1[names(boot)] / tissueExp2[names(boot)]
  )
  rownames(boot2) <- names(boot)
  boot2 <- (boot2[order(boot2$change,decreasing=TRUE),])
  return(boot2)
}
#######################
# Get file locations 
LoadDevelopingHuman <- function(){
  #slices <- c(1,6,8,12,19,24,29,32,39,45)
  slices <- c(29)
  ## Get File Locations
  DIR <- c(paste("~/SalkProjects/developingHuman/AllenBrainSlicescopy/",paste("Slice_",slices,"/",sep=""),sep=""))
  fileloc <- c()
  i <- 0
  for (d in DIR){
    i <- i + 1
    fileloc[[i]] <- paste(d,list.files(path=d,pattern="tif"),sep="")
  }
  names(fileloc) <- DIR
  
  ## Get Abreviations
  abrev <- c()
  for (i in 1:length(DIR)){
    tmp <- do.call("rbind",strsplit(x=list.files(path=DIR[i],pattern="tif"),fixed=TRUE,split="."))[,1]
    abrev[[i]] <-  tmp[-c(which(tmp == "Outline"))]
  }
  names(abrev) <- DIR
  
  #conversion of abreviation to the regions that are in this image
  conversion <- as.data.frame(read.table(as.matrix("~/SalkProjects/colmeta_hierarchy_2.txt"),header=TRUE))
  #Dimensions of the image
  DIM <- c()
  options(warn = -1)
  for (d in 1:length(DIR)){
    fordim <- paste(DIR[d],"CP.tif",sep="")
    DIM[[d]] <- dim(readTIFF(source=fordim)[,,1])
  }
  names(DIM) <- DIR
  #Outline image
  outline <- c()
  for (i in 1:length(DIR)){
     tmpfile <- paste(DIR[i],"Outline.tif",sep="")
     outline[[i]] <- readTIFF(source=tmpfile)[,,1]
  }
  options(warn = 0)
  names(outline) <- DIR
 
  fileloc <<- fileloc
  abrev <<- abrev
  conversion <<- conversion
  DIM <<- DIM
  outline <<- outline
  slices <<- slices
  DIR <<- DIR
  
}
LoadAdultHuman <- function(){
  slices <- c(5,17,25,31,51,59,67,77,83,91)#43 missing Outline.tiff
  ## Get File Locations
  DIR <- c(paste("~/Documents/SalkProjects/ME/brainImageR_ABA/brainImageR_rawdata/adultHuman/AllenBrainSlicescopy/",paste("Slice_",slices,"/",sep=""),sep=""))
  fileloc <- c()
  i <- 0
  for (d in DIR){
    i <- i + 1
    fileloc[[i]] <- paste(d,list.files(path=d,pattern="tif"),sep="")
  }
  names(fileloc) <- DIR
  
  ## Get Abreviations
  abrev <- c()
  for (i in 1:length(DIR)){
    tmp <- do.call("rbind",strsplit(x=list.files(path=DIR[i],pattern="tif"),fixed=TRUE,split="."))[,1]
    abrev[[i]] <-  tmp[-c(which(tmp == "Outline"))]
  }
  names(abrev) <- DIR
  
  #conversion of abreviation to the regions that are in this image
  conversion <- as.data.frame(read.table(as.matrix("~/Documents/SalkProjects/ME/brainImageR_ABA/brainImageR_rawdata/adultHuman/AllenBrainSlicescopy/Human34_Hierarchy.txt"),header=TRUE))
  #Dimensions of the image
  DIM <- c()
  options(warn = -1)
  for (d in 1:length(DIR)){
    f <- list.files(DIR[d])
    fordim <- paste(DIR[d],f[1],sep="")
    DIM[[d]] <- dim(readTIFF(source=fordim)[,,1])
  }
  names(DIM) <- DIR
  #Outline image
  outline <- c()
  for (i in 1:length(DIR)){
    tmpfile <- paste(DIR[i],"Outline.tif",sep="")
    outline[[i]] <- readTIFF(source=tmpfile)[,,1]
  }
  options(warn = 0)
  names(outline) <- DIR
  
  fileloc <<- fileloc
  abrev <<- abrev
  conversion <<- conversion
  DIM <<- DIM
  outline <<- outline
  slices <<- slices
  DIR <<- DIR
  
}


#######################
# RUN THE PROGRAM
#######################
# pick your Dir
setwd("C:/Users/Quinn Wong/Documents/SalkProjects")
mydata <- read.table("Genes.txt")
genes <- as.character(mydata$V1)
LoadDevelopingHuman()
#LoadAdultHuman()
######
composite <- SpatialEnrichment(slice = 1, genes,reps = 20)
PlotBrain(composite@composite, Breaks = 9)
### Find the genes in that region
regiongenes <- GetGenes(genes,tissueExp1,region.name="SZ")
boot <- Boot(tissueExp1 = composite@tissueExp1, tissueExp2 = composite@tissueExp2, random.matrix = composite@random.matrix)
