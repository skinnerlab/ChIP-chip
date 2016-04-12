####################################################
## Primary Researcher: Ben and Carlos
##		Washington State University
## 		Skinner labratory
##
## Bioinformatician: Matt Settles
##	Washington State University
##  Normalize by Adjusted GC content
##
##  Code Modified by: Md. Muksitul Haque
################## PREAMBLE ############

base <- "Enter base Directory" # The base directory
experimentName <- "Enter Experiment Name"          # Experiment Name

setwd(base)
rawPath <- file.path(getwd(),"Data")
figurePath <- file.path(getwd(), "Figures")
tablePath <- file.path(getwd(), "Tables")
pairPath <- file.path(getwd(),"Raw_Data_Files")
designPath <- file.path( getwd(),"Design_Information")
probeAnnoPath <- file.path( getwd(),"Probe_Anno")
library(limma)

compChr <- paste("chr",c(1:20,"1_random","2_random","3_random","4_random",
		"5_random","6_random","7_random","8_random",
		"9_random","11_random","12_random","15_random",
		"16_random","19_random","20_random","Un","X","X_random"),sep="") 

design <- c(1,1,1)
groupings <- c(1,2,3)

targetsFile <- "SampleKey.txt"  # SampleKeyFile from Nimblegen
designids <- c("Enter DesignIDs")           # design ID of the Experiment
ext <- ".pair"
################################################################################
### LOAD DATA FILES
#############################
### Load Raw Data ###
#############################

for (i in 1:length(designids)){
	# i <- 1
	design_id <- designids[i] 
	
		hybes  <- readTargets(targetsFile)
		hybes <- hybes[order(hybes$DESIGN_ID,hybes$ORD_ID,hybes$DYE),]
		hybes <- hybes[which(hybes$DESIGN_ID == design_id),]
		
		source(file.path(base,"readNimblegen.R"))
		design_name <- hybes[,"DESIGN_NAME"][1]
		
		narray <- length(unique(hybes[,"CHIP_ID"]))
		cy3FileName <- grep("[Cc][Yy]3", hybes$DYE)
		cy5FileName <- grep("[Cc][Yy]5",  hybes$DYE)
		files <- paste(hybes[c(cy3FileName,cy5FileName),"CHIP_ID",drop=TRUE],"_",rep(c(532,635),each=narray),ext,sep="")
		files <- matrix(files, ncol = 2, byrow = FALSE)
		RG <- readPairFiles(files, verbose=TRUE, path=pairPath)
		RG$Rb <- RG$Gb <- NULL ## background not used
		spottypes <- readSpotTypes("spottypes.txt")
		RG$genes$Status <- controlStatus(spottypes, RG$genes)
		RG$genes <- RG$genes[, match(c("PROBE_ID","Status"),colnames(RG$genes))]
		RG$targets <- cbind(hybes[cy3FileName,c("ORD_ID","CHIP_ID","DESIGN_NAME","DESIGN_ID","SAMPLE_SPECIES")],
				matrix(hybes$SAMPLE_DESCRIPTION[c(cy3FileName,cy5FileName)],ncol=2,byrow=FALSE,dimnames=list(NULL,c("Cy3","Cy5"))))
		## remove duplicated probes, if any
		if (!all(!duplicated(RG$genes$PROBE_ID)))
			RG <- RG[-which(duplicated(RG$genes$PROBE_ID)==TRUE),]
		
		#############################
		### Read ndf and Pos files
		#############################
			designName <- RG$targets$DESIGN_NAME[1]
			NDF <- read.table( file.path(designPath,paste(designName,"ndf",sep=".")) , 
					sep="\t",header=TRUE,as.is=TRUE)
			POS <- read.table( file.path(designPath,paste(designName,"pos",sep=".")) , 
					sep="\t",header=TRUE,as.is=TRUE)
			NDF <- NDF[match(RG$genes$PROBE_ID,NDF$PROBE_ID),]
			NDF$CHROMOSOME <- "RANDOM"
			NDF[match(POS$PROBE_ID,NDF$PROBE_ID),"CHROMOSOME"] <- POS$CHROMOSOME
			NDF$LENGTH <- nchar(NDF$PROBE_SEQUENCE)
			NDF$COUNT <- 0
			NDF[match(POS$PROBE_ID,NDF$PROBE_ID),"COUNT"] <- POS$COUNT
			library(Biostrings)
			NDF$GC <-sapply(gregexpr("[CG]",NDF$PROBE_SEQUENCE) ,length)
			NDF$CG <-sapply(gregexpr("CG",NDF$PROBE_SEQUENCE) ,length)
			NDF <- NDF[order(NDF$CHROMOSOME,NDF$POSITION),]
			POS <- NDF
			rm(NDF)
			save(POS,file=file.path(probeAnnoPath,paste("ndfpos",design_id,"RData",sep=".")))
			RG <- RG[match(POS$PROBE_ID,RG$genes$PROBE_ID),]
			save(RG,file=file.path(rawPath,paste("UnmodifiedRawData",design_id,"RData",sep=".")))
	
	### Remove Saturated Probes
	sat <- 2^(15.5)
	removeA1 <- c(apply(RG$R,2,function(x) which (x > sat)),
			apply(RG$G,2,function(x) which (x > sat)))
	## DUPLICATED PROBES COUNT > 1
	duplicated <- which(POS$COUNT > 1)
	RG$R[as.numeric((unique(c(unlist(removeA1),duplicated)))),] <- NA
	RG$G[as.numeric((unique(c(unlist(removeA1),duplicated)))),] <- NA
	## removed saturated and duplicated probes
	save(RG,file=file.path(rawPath,paste("PreprocRawData",design_id,"RData",sep=".")))
	
	################################################33
	# Log Values
	logR <- log2(RG$R)
	logG <- log2(RG$G)
	## Adjust and GC content
	narray <- ncol(RG$R)

	min <- 2
	GCcount = table(as.factor(POS$GC))
	Rmean = apply(logR,2, function (x) tapply(x,POS$GC, mean,na.rm=TRUE))
	Gmean = apply(logG,2, function (x) tapply(x,POS$GC, mean,na.rm=TRUE))
	Rvar = apply(logR,2, function(x) tapply(x,POS$GC,var,na.rm=TRUE))
	Gvar = apply(logG,2, function(x) tapply(x,POS$GC,var,na.rm=TRUE))
	gccov = matrix(NA,nrow=length(GCcount),ncol=ncol(logR))
	for(i in 1:ncol(logR)) 
		gccov[,i] <- tapply(1:nrow(logR),POS$GC,
				function(x) cov(cbind(logR[x,i],logG[x,i]),use="pairwise.complete.obs")[1,2])
	
	GCsum <- list(GCcount=GCcount,Rmean=Rmean,Gmean=Gmean,Rvar=Rvar,Gvar=Gvar,gccov=gccov)
	
	## no background correction
	object <- MA.RG(RG, bc.method = "none", offset = 0)
	if (is.vector(object$M)) 
		object$M <- as.matrix(object$M)
	if (is.vector(object$A)) 
		object$A <- as.matrix(object$A)
	
	## GC correction 
	source(file.path(codePath,"maplot.R"))
	plotit=TRUE
	pdf(file=file.path(figurePath,paste(experimentName,design_id,"loesscurvesbyGC","pdf",
							sep=".")),height=7,width=7,pointsize=8)
	for (j in 1:narray) {
		if (plotit){
			plot(object$A[,j],object$M[,j],type="n",cex = 0.5,ylim=c(-5,5),xlab="A",ylab="M",main=paste("Loess Curves by GC, Chipset", design_id,", Array:",RG$targets$CHIP_ID[1]))
			abline(0, 0, col = "blue")
		}
		cat("array", j,"\tGC:")
		for (GCp in levels(factor(POS$GC)) ) {
			cat(GCp,":")
			spots <- which(POS$GC == as.numeric(GCp) )
			cat(length(spots),"..")
			if (length(na.exclude(object$M[spots,j])) > 1){
				y <- object$M[spots, j]
				x <- object$A[spots, j]
				#w <- weights[spots, j]
				lfit <- loessFit(y, x, span = 0.3, iterations = 4)
				object$M[spots, j] <- lfit$residuals
				if(plotit){
					o <- order(x)
					A <- x[o]
					M <- lfit$fitted[o]
					o <- which(!duplicated(x))
	#				lines(approx(A[o], M[o]), col = "red", lwd = 1, lty = 1)
				}
			} else {
				object$M[spots,j] <- NA
			}
		}
		cat("\n")
	}
	dev.off()
	
	
	object <- normalizeBetweenArrays(object) # A-quantile , quantile norm the A values
	
	############## MA2C
	Rmean = apply(logR,2, function (x) tapply(x,POS$GC, mean,na.rm=TRUE))
	numer <- (logR - Rmean[as.character(POS$GC),]) - (logG - Gmean[as.character(POS$GC),])
	denom <- sqrt(Gvar + Rvar - 2*gccov)[as.character(POS$GC),]
	MA2C <- numer/denom  ## adjusted further for correlation
	MA2C <- scale(MA2C,center=FALSE)
	MA2CA <-  (logR + logG)/2
	
	save(object,MA2C,MA2CA,GCsum,file=file.path(rawPath,paste("GCNormalized",design_id,"RData",sep=".")))
}



########## compute running medians
source("computeRunningMedians.R")
################################################################################
### LOAD DATA FILES
#############################
### Load Raw Data, one array at a time ###
#############################
## smooth M
for (i in 1:length(designids)){
	design_id <- designids[i] 
	load(file=file.path(rawPath,paste("GCNormalized",design_id,"RData",sep=".")))
	load(file=file.path(probeAnnoPath,paste("ndfpos",design_id,"RData",sep=".")))
	POS <- POS[match(object$genes$PROBE_ID,POS$PROBE_ID),]
	ord <- order(POS$CHROMOSOME,POS$POSITION)
	POS <- POS[ord,]
	object<- object[ord,]
	
	chr <- intersect(unique(POS$CHROMOSOME),compChr)
	
	smooth300.gc <- computeRM(
			object = object ,
			probeAnno = POS[,c("PROBE_ID","CHROMOSOME","POSITION","LENGTH","COUNT")] ,
			design=design ,
			groupings = groupings ,
			allChr = chr ,
			winHalfSize = 300 ,
			min.probes = 4 ,
			quant=0.5 ,
			combineReplicates = FALSE , 
			checkUnique = FALSE ,
			uniqueCodes = c(1) ,
			verbose = TRUE 
	)
	save(smooth300.gc,file=file.path(rawPath,paste("GCNormSmoothed",design_id,"RData",sep="."))) 
}


maxp <- 1e-7
for (i in 1:length(designids)){
	design_id <- designids[i] 
	load(file=file.path(probeAnnoPath,paste("ndfpos",design_id,"RData",sep=".")))
	load(file=file.path(rawPath,paste("GCNormSmoothed",design_id,"RData",sep=".")))

	nas <- which(!duplicated(rownames(smooth300.gc)))
        smooth300.gc <- smooth300.gc[nas,]
 
	rowM <- rowMedians(smooth300.gc$M)
	rowA <- rowMedians(smooth300.gc$A)
	sdA <- apply(smooth300.gc$A,1,sd,na.rm=TRUE)
	meanM <- mean(rowM,na.rm=TRUE);sdM <- sd(rowM,na.rm=TRUE)
	rawp <- 2*pnorm(abs(rowM-meanM),mean=0,sd=sdM,lower.tail=FALSE)
	res <- data.frame(M=rowM,A=rowA,sdA = sdA,t=NA,rawp=rawp,row.names=rownames(smooth300.gc))
	write.table(res,file=file.path(tablePath,paste("tresults",design_id,"csv",sep=".")),sep=",",row.names=TRUE,col.names=TRUE)

	probes <- res[which(res$rawp < maxp),]
	probes <- data.frame(probes,POS[match(rownames(probes),POS$PROBE_ID),c("CHROMOSOME","POSITION")])
	probes$POSITION <- probes$POSITION + ceiling(POS[match(rownames(probes),POS$PROBE_ID),c("LENGTH")]/2)
	
	probes$Dir <-  cut(probes$M, c(-10000,0 ,10000), labels = c("Cont","Vinc"))
	probes$rank= floor(-log10(probes$rawp))
	probes$cluster <- 1
	probes <- probes[order(probes$CHROMOSOME,as.numeric(as.character(probes$POSITION))),]
	for (i in 2:nrow(probes)){
		if ( probes[i,"CHROMOSOME"] == probes[i-1,"CHROMOSOME"]){
			if( (as.numeric(as.character(probes[i,"POSITION"])) - 300) < (as.numeric(as.character(probes[i-1,"POSITION"])) + 300) & (probes[i,"Dir"] == probes[i-1,"Dir"])){
				probes[i,"cluster"] <- probes[i-1,"cluster"]
			} else {
				probes[i,"cluster"] <- probes[i-1,"cluster"] + 1
			}
		} else {
			probes[i,"cluster"] <- probes[(i-1),"cluster"] + 1
		}       
	}
	clusterID <- tapply(1:length(probes$cluster),probes$cluster,function(x) paste(probes[["CHROMOSOME"]][x][1],":",min(probes[["POSITION"]][x])-300,"-",max(probes[["POSITION"]][x]+300),sep=""))
	probes$clusterloc <- clusterID[probes$cluster]
	
	#################### Clusters #####################################
	library(BSgenome)
	library("BSgenome.Rnorvegicus.UCSC.rn4")
	
	sregions <- split(probes,probes$cluster)
	regions <- sapply(sregions, function(x) {
				z <- c(ClusterID=x[1,"clusterloc"],Chromosome = x[1,"CHROMOSOME"],cSTART = min(x$POSITION)-300, cSTOP=max(x$POSITION)+300, meanM = mean(x$M),meanA = mean(x$A),minP = min(x$rawp),nprobes = nrow(x))
				a <- c(Score = max(x$rank))
				t(c(z,a))
			})
	regions <- as.data.frame(t(regions))
	colnames(regions) <- c("ClusterID","Chromosome","cSTART","cSTOP","meanM","meanA","minP","nProbes","Score")
	Seq <- apply(regions[,c("Chromosome","cSTART","cSTOP")],1,function(x) as.character(unmasked(getSeq(Rnorvegicus, x[1], start=as.integer(x[2]), end=as.integer(x[3]),as.character=FALSE))))
	base <- alphabetFrequency(DNAStringSet(Seq),baseOnly=TRUE)
	regions$Length <- nchar(Seq)
	regions$GC <- rowSums(base[,2:3])
	regions$CpG <- dinucleotideFrequency(DNAStringSet(Seq))[,"CG"]	
	regions$Seq <- Seq
	write.table(probes,file=file.path(tablePath,paste("ProbeResults",design_id,"csv",sep=".")),sep=",",row.names=TRUE,col.names=TRUE)
	write.table(regions,file=file.path(tablePath,paste("RegionsResults",design_id,"csv",sep=".")),sep=",",row.names=TRUE,col.names=TRUE,quote=FALSE)
}

for (i in 1:length(designids)){
	design_id <- designids[i] 
	
	tmp <- read.table(file=file.path(tablePath,paste("RegionsResults",design_id,"csv",sep=".")),sep=",")
	if (i == 1)
		cresults <- tmp
	else
		cresults <- rbind(cresults,tmp)
}
write.table(cresults,file=file.path(tablePath,paste("RegionsCombined","csv",sep=".")),sep=",",row.names=TRUE,col.names=TRUE,quote=FALSE)


############ Visualization of the Chip-Chip Data in Barplots #####################################################################

# Look for RG$genes$PROBE_ID and p-value to plot each regions
for (i in 1:length(designids))
{
design_id <- designids[i] 
load(file=file.path(probeAnnoPath,paste("ndfpos",design_id,"RData",sep=".")))
load(file=file.path(rawPath,paste("GCNormSmoothed",design_id,"RData",sep=".")))
# load smooth300.gc list then convert it to RG
# get the POS list to get probe_sequence_Id
RG <- RG.MA(smooth300.gc)
RG$genes$SEQ_ID <- POS[match(POS$PROBE_ID, RG$genes$PROBE_ID),"SEQ_ID"]

temp <- RG$genes$SEQ_ID
g <- factor(temp)
temp <- split(temp, g)
# remove the random SEQ_ID
temp <- temp[1:21554]

# read in the regionscombined files [after annotation]
regions <- read.table("RegionsCombined.csv")
plot_table <-  matrix(0, ncol=2, nrow=nrow(regions), dimnames=list(c(1:nrow(regions)), c("Probe_Sequence_ID", "GeneName")))
regions <- data.frame(Chromosome=regions[,"Chromosome"], cSTART=regions[,"cSTART"], cSTOP=regions[,"cSTOP"], Treatment=regions[,"Treatment"],minp=regions[,"minP"],GeneName=regions[,"GeneName"])

chromosome <- ProbeRegionStart <- ProbeRegionEnd <- Probe_Start <- Probe_End <- Probe_chromosome <- plot_table <- c()
# Find probe_id_sequence
for(i in 1:nrow(regions))
{
chromosome <- as.character(regions[i,"Chromosome"])
ProbeRegionStart <- as.numeric(regions[i,"cSTART"])
ProbeRegionEnd <- as.numeric(regions[i,"cSTOP"])

	for(j in 1:length(temp))
	{
	temp1 <- as.character(names(temp[j]))
	Probe_chromosome <- as.character(strsplit(temp1, ":")[[1]][1])
	Probe_Start <- as.numeric(strsplit(strsplit(temp1, ":")[[1]][2],"-")[[1]][1])
	Probe_End <- as.numeric(strsplit(strsplit(temp1, ":")[[1]][2],"-")[[1]][2])
		if(Probe_chromosome==chromosome)
		{
			if( ((Probe_Start >= ProbeRegionStart) & (ProbeRegionEnd >= Probe_End)) |
			    ((Probe_End >= ProbeRegionStart) & (ProbeRegionEnd >= Probe_End)) |
			    ((Probe_Start >= ProbeRegionStart) & (ProbeRegionEnd >= Probe_Start)) |
			    ((Probe_Start <= ProbeRegionStart) & (ProbeRegionEnd <= Probe_End)) )
				{
				print(regions[i,])
				print(temp1)
				plot_table[i,"Probe_Sequence_ID"] <- temp1
				plot_table[i,"GeneName"] <- as.character(regions[i,"GeneName"])
				}
		}
	}
	print(i)
}

# get the p-value for each probe and use it in the plot
nas <- which(!duplicated(rownames(smooth300.gc)))
smooth300.gc <- smooth300.gc[nas,]
library(Biobase)
rowM <- rowMedians(smooth300.gc$M)
rowA <- rowMedians(smooth300.gc$A)
sdA <- apply(smooth300.gc$A,1,sd,na.rm=TRUE)
meanM <- mean(rowM,na.rm=TRUE);sdM <- sd(rowM,na.rm=TRUE)
rawp <- 2*pnorm(abs(rowM-meanM),mean=0,sd=sdM,lower.tail=FALSE)
res <- data.frame(M=rowM,A=rowA,sdA = sdA,t=NA,rawp=rawp,row.names=rownames(smooth300.gc))
#maxp <- 1e-7
#res <- res[which(res$rawp < maxp),]
RG$genes$rawp <- res[match(RG$genes$PROBE_ID, rownames(res)),"rawp"]

for(i in 1:nrow(plot_table))
{

	 jpeg(file=file.path(figurePath,paste(as.character(plot_table[i,"GeneName"]),".", as.character(gsub(":",".", plot_table[i,"Probe_Sequence_ID"])), "jpeg",sep=".")) , quality=75, width=5000, height=1575,pointsize=24)
	barplot(log2(rowMeans(RG$R[which(RG$genes$SEQ_ID==plot_table[i,"Probe_Sequence_ID"]),]))-log2(rowMeans(RG$G[which(RG$genes$SEQ_ID==plot_table[i,"Probe_Sequence_ID"]),])), 
	legend=RG$genes[which(RG$genes$SEQ_ID==as.character(plot_table[i,"Probe_Sequence_ID"])),"rawp"], lwd=2, col="black", main=paste(plot_table[i,"Probe_Sequence_ID"], sub=paste("ChipSet",design_id)))
	print(i)
	dev.off()
}

