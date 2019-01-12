#!/usr/bin/env Rscript


#----------------------------------
#- Load / install required packages
#----------------------------------
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

if (!require("DESeq2")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("DESeq2", suppressUpdates=TRUE)
  library("DESeq2")
}
if (!require("ReportingTools")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("ReportingTools", suppressUpdates=TRUE)
  library("ReportingTools")
}
if (!require("hwriter")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("hwriter", suppressUpdates=TRUE)
  library("hwriter")
}
if (!require("ini")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("ini", suppressUpdates=TRUE)
  library("ini")
}
if (!require("ggrepel")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("ggrepel", suppressUpdates=TRUE)
  library("ggrepel")
}
if (!require("png")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("png", suppressUpdates=TRUE)
  library("png")
}

#------------
#- Parse args
#------------

args<-read.ini("deseq2.conf")

inputdir = args$deseq2$INPUTDIR
metadata = args$deseq2$METADATA
design = args$deseq2$DESIGN
condition = args$deseq2$CONDITION
treatment = args$deseq2$TREATMENT
control = args$deseq2$CONTROL
pval = as.numeric(args$deseq2$PVAL)
fc = as.numeric(args$deseq2$FC)

#------------
#- Get counts
#------------
countData<-read.table(paste0(inputdir,"/featureCounts/merged_gene_counts.txt"), sep="\t", header=T, check.names=FALSE)
geneID<-countData$ENSEMBL_ID
countData<-select(countData, -ENSEMBL_ID)
rownames(countData)<-geneID

#--------------
#- Get metadata
#--------------
colData <- read.table(metadata, sep="\t",header=T)

#-------------------------------------------------------------------------
#- Make sure that first column has the sample names, that these samples
#-  are all in countData, and that they are in the same order
#-------------------------------------------------------------------------
rownames(colData)<-colData[,1]
colData<-colData[,-c(1)]
if(length(rownames(colData)[!rownames(colData) %in% colnames(countData)]) > 0) {	
	stop("ERROR: the following samples are not in the featureCounts matrix: ", paste(rownames(colData)[!rownames(colData) %in% colnames(countData)], collapse=", "),
	      ". Please ,ake sure that the first column of your metadata file has the sample IDs.")
}else{
  countData = countData[ , rownames(colData) ]
  if(!all(rownames(colData) == colnames(countData))) {	
    stop("Something is wrong, this should never happen [rownames(colData) ne colnames(countData)??]")
  }
}

#-------------------------------
#- Specify design and contrast
#-------------------------------
default="no"
if(design == '-'){
	design = paste0(colnames(colData)[1:ncol(colData)], collapse=" + ")
	default="yes"
}else{
  if(!(condition %in% colnames(colData))) {	
    stop("ERROR: condition '",condition,"' is not a column of the metadata file, please specify a valid condition")
  }
  if(!(treatment %in% levels(colData[colnames(colData) == condition][,1]))) {	
    stop("ERROR: treatment '",treatment,"' is not a level in condition ",condition,". Please specify a valid treatment")
  }
  if(!(control %in% levels(colData[colnames(colData) == condition][,1]))) {	
    stop("ERROR: control '",control,"' is not a level in condition ",condition,". Please specify a valid control")
  }
}

#-------------
#- create DDS
#-------------
dds <- DESeqDataSetFromMatrix(countData = countData,
            colData = colData,
            design = eval(parse(text=paste0("~ ", design))))

#-------------------------------------------
#- remove genes with 0 counts in all samples
#-------------------------------------------
dds <- dds[ rowSums(counts(dds)) > 1, ]

#---------
#- Run DGE
#----------
dds <- DESeq(dds)



#---------
#- PCA
#----------
rld <- rlog(dds, blind=TRUE)
pcaData<- plotPCA(rld, intgroup=colnames(colData), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaplot_name<-"PCAplot.png"
png(file=pcaplot_name)
if(ncol(colData)>1){
  ggplot(pcaData, aes(PC1, PC2, color=colData[,1], shape=colData[,2])) +
    geom_point(size=3) +
    ggtitle("PCA plot") +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    theme_classic() + theme(legend.position = "bottom", legend.title=element_blank())
}else{
  ggplot(pcaData, aes(PC1, PC2, color=colData[,1])) +
    geom_point(size=3) +
    ggtitle("PCA plot") +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    theme_classic() + theme(legend.position = "bottom", legend.title=element_blank())
}
dev.off()


if (default == "no") {
	  
	#-------------------------------------------------------------
	#- Get results using  design and contrast defined by the user
	#-------------------------------------------------------------
	contrast = c(condition,treatment, control)
	       
	#- alpha is fdr threshold for summary display only
	res<-results(dds, contrast = contrast, alpha=0.05)
	resNorm <- lfcShrink(dds, contrast = contrast, res=res, type="normal")
	    
	#----------------------------
	#- File with all the results
	#----------------------------
	write.table(resNorm, file=paste0(condition, "_", treatment, "_vs_", control, "_results.txt"), sep="\t",quote=F,row.names=T, col.names=NA)
	    
	#-----------------------
	#- MA plot
	#-----------------------
	maplot_name<-paste0(condition, "_", treatment, "_vs_", control, "_MAplot.png")
	png(file=maplot_name)
	maplot <- plotMA(resNorm, ylim=c(-5,5))
	dev.off()
	    
	#------------------------
	#- Volcano plot
	#------------------------
	  
	vplot_name<-paste0(condition, "_", treatment, "_vs_", control, "_VolcanoPlot.png")
	
	ha<-data.frame(gene = rownames(resNorm), logFC = resNorm$log2FoldChange, Pval = resNorm$pvalue)
	ha = within(ha, {Col="Other"})
	ha[abs(ha$logFC) > fc, "Col"] <- paste0("|logFC| > ",fc)
	ha[!is.na(ha$Pval) & ha$Pval <pval, "Col"] <- paste0("Pval < ",pval)
	ha[!is.na(ha$Pval) & ha$Pval <pval & abs(ha$logFC)>fc, "Col"] <- paste0("Pval < ",pval, " & |logFC| > ", fc)
	ha$Col<-factor(ha$Col, levels=c(paste0("Pval < ",pval),paste0("|logFC| > ",fc),paste0("Pval < ",pval, " & |logFC| > ", fc), "Other"), labels=c(paste0("Pval < ",pval),paste0("|logFC| > ",fc),paste0("Pval < ",pval, " & |logFC| > ", fc), "Other"))
	
	png(file=vplot_name)
	p<-ggplot(ha, aes(logFC,  -log10(Pval))) +
	  geom_point(aes(col=Col), alpha=0.5) + 
	  ggtitle(paste0(condition, " - ", treatment, " vs ", control)) +
	  geom_hline(yintercept = -log10(pval), linetype = 2, alpha = 0.5) + 
	  geom_vline(xintercept = fc, linetype = 2, alpha = 0.5) +
	  geom_vline(xintercept = -fc, linetype = 2, alpha = 0.5) +
	  scale_color_manual(values=c("indianred1", "gold2", "cornflowerblue", "gray47")) + 
	  theme_classic() + theme(legend.position = "bottom", legend.title=element_blank())
	
	#- add gene labels for gene swith pval anf fc below theresholds
	ha2 <- ha %>%
	  filter(Pval < pval & (logFC >= fc | logFC <= -fc)) %>% 
	  select(gene,logFC,Pval, Col)
	
	p<-p+geom_text_repel(data=ha2,aes(label=gene), box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"), force=10, size=3,segment.size=0.25,segment.alpha=0.5)
	print(p)
	dev.off()
	
	#-----------------------
	#- Report
	#-----------------------
	    
	allplot_name<-paste0(condition, "_", treatment, "_vs_", control, "_AllPlot.png")
	img1 <- readPNG(maplot_name)
	img2 <- readPNG(vplot_name)
	png(file=allplot_name, width=1200, height=600)
	par(mai=rep(0,4)) # no margins
	layout(matrix(1:2, ncol=2, byrow=TRUE))
	for(i in 1:2) {
	  plot(NA,xlim=0:1,ylim=0:1,bty="n",axes=0,xaxs = 'i',yaxs='i')
	  rasterImage(eval(parse(text=paste0("img", i))),0,0,1,1)
	}
	dev.off()

	reportdir="report"
	des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2', title = 'RNA-seq DGE analysis using DESeq2 (\'normal\' shrinkage) with user-defined design',reportDirectory = reportdir)

	himg <- hwriteImage(paste0("figuresRNAseq_analysis_with_DESeq2/",pcaplot_name))
	publish(hwrite(himg, br=TRUE, center=T), des2Report)
	publish(paste0("<center><h3><u>Contrast:</u>  ",condition, " - ", treatment, " vs ", control), des2Report)
	publish("<h4>MA/Volcano plots", des2Report)
	himg <- hwriteImage(paste0("figuresRNAseq_analysis_with_DESeq2/",allplot_name))
	publish(hwrite(himg, br=TRUE, center=T), des2Report)
	publish("<h4>Top 100 differentially expressed genes", des2Report)
	publish(resNorm,des2Report, reportDir=reportdir, pvalueCutoff=1, n=100, DataSet=dds, factor=colData(dds)[[i]])
	
	finish(des2Report)

}else{

	#-------------------------------------------------------------
	#- Get results using default design and all possible contrasts 
	#-------------------------------------------------------------
	reportdirALL="report"
	des2ReportALL <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2', title = 'RNA-seq DGE analysis using DESeq2 (\'normal\' shrinkage) in default (no design) mode',reportDirectory = reportdirALL)
	himg <- hwriteImage(paste0("figuresRNAseq_analysis_with_DESeq2/",pcaplot_name))
	publish(hwrite(himg, br=TRUE,center=TRUE), des2ReportALL)
	n=1
	for (i in 1:ncol(colData)) {
		pairs<-combn(unique(colData[,i]),2)
		for (j in 1:ncol(pairs)) {

			#- alpha is fdr threshold for summary display only
			res<-results(dds, contrast = c(colnames(colData)[i],as.character(pairs[,j])), alpha=0.05)
			resNorm <- lfcShrink(dds, contrast = c(colnames(colData)[i],as.character(pairs[,j])), res=res, type="normal")
	
			#-----------------------
			#- File with all results
			#-----------------------
			write.table(resNorm, file=paste0(colnames(colData)[i], "_", paste0(as.character(pairs[,j]),collapse="_vs_"), "_results.txt"), sep="\t",quote=F,row.names=T, col.names=NA)

			#-----------------------
			#- MA plot
			#-----------------------
			maplot_name<-paste0(colnames(colData)[i], "_", paste0(as.character(pairs[,j]),collapse="_vs_"), "_MAplot.png")
			png(file=maplot_name)
			maplot <- plotMA(resNorm, ylim=c(-5,5))
			dev.off()

			#------------------------
			#- Volcano plot
			#------------------------

			vplot_name<-paste0(colnames(colData)[i], "_", paste0(as.character(pairs[,j]),collapse="_vs_"), "_VolcanoPlot.png")
			
			ha<-data.frame(gene = rownames(resNorm), logFC = resNorm$log2FoldChange, Pval = resNorm$pvalue)
			ha = within(ha, {Col="Other"})
			ha[abs(ha$logFC) > fc, "Col"] <- paste0("|logFC| > ",fc)
			ha[!is.na(ha$Pval) & ha$Pval <pval, "Col"] <- paste0("Pval < ",pval)
			ha[!is.na(ha$Pval) & ha$Pval <pval & abs(ha$logFC)>fc, "Col"] <- paste0("Pval < ",pval, " & |logFC| > ", fc)
			ha$Col<-factor(ha$Col, levels=c(paste0("Pval < ",pval),paste0("|logFC| > ",fc),paste0("Pval < ",pval, " & |logFC| > ", fc), "Other"), labels=c(paste0("Pval < ",pval),paste0("|logFC| > ",fc),paste0("Pval < ",pval, " & |logFC| > ", fc), "Other"))

			png(file=vplot_name)
			p<-ggplot(ha, aes(logFC,  -log10(Pval))) +
			geom_point(aes(col=Col), alpha=0.5) + 
			ggtitle(paste0(colnames(colData)[i], "_", paste0(as.character(pairs[,j]),collapse="_vs_"))) + 
			geom_hline(yintercept = -log10(pval), linetype = 2, alpha = 0.5) + 
			geom_vline(xintercept = fc, linetype = 2, alpha = 0.5) +
			geom_vline(xintercept = -fc, linetype = 2, alpha = 0.5) +
			scale_color_manual(values=c("indianred1", "gold2", "cornflowerblue", "gray47")) + 
			theme_classic() + theme(legend.position = "bottom", legend.title=element_blank())

			#- add gene labels for gene swith pval anf fc below theresholds
			ha2 <- ha %>%
			   filter(Pval < pval & (logFC >= fc | logFC <= -fc)) %>% 
			   select(gene,logFC,Pval, Col)
			
			p<-p+geom_text_repel(data=ha2,aes(label=gene), box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"), force=10, size=3,segment.size=0.25,segment.alpha=0.5)
			print(p)
			dev.off()

			#-----------------------
			#- Report
			#-----------------------
	
			allplot_name<-paste0(colnames(colData)[i], "_", paste0(as.character(pairs[,j]),collapse="_vs_"), "_AllPlot.png")
			img1 <- readPNG(maplot_name)
			img2 <- readPNG(vplot_name)
			png(file=allplot_name, width=1200, height=600)
			par(mai=rep(0,4)) # no margins
			layout(matrix(1:2, ncol=2, byrow=TRUE))
			for(l in 1:2) {
			  plot(NA,xlim=0:1,ylim=0:1,bty="n",axes=0,xaxs = 'i',yaxs='i')
			  rasterImage(eval(parse(text=paste0("img", l))),0,0,1,1)
			}
		        dev.off()

			publish(paste0("<center><h3><u>Contrast ", n, ":</u>  ", colnames(colData)[i], " - ",  paste0(as.character(pairs[,j]),collapse=" vs ")), des2ReportALL)
			publish("<h4>MA/Volcano plots", des2ReportALL)
			himg <- hwriteImage(paste0("figuresRNAseq_analysis_with_DESeq2/",allplot_name))
			publish(hwrite(himg, br=TRUE,center=TRUE), des2ReportALL)
			publish("<h4>Top 100 differentially expressed genes", des2ReportALL)
			publish(resNorm,des2ReportALL, reportDir=reportdirALL, pvalueCutoff=1, n=100, DataSet=dds, factor=colData(dds)[[i]])
      			n=n+1
		}	
	}
	finish(des2ReportALL)
}


