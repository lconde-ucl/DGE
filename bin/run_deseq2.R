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
if (!require("rhdf5")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("rhdf5", suppressUpdates=TRUE)
  library("rhdf5")
}
if (!require("tximport")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("tximport", suppressUpdates=TRUE)
  library("tximport")
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
kallisto = args$deseq2$KALLISTO
tx2gene.file = args$deseq2$TX2GENE


#--------------
#- Get metadata
#--------------
colData <- read.table(metadata, sep="\t",header=T)
samples<-colData[,1]
colnames<-colnames(colData)
colData<-as.data.frame(colData[,-c(1)])
rownames(colData)<-samples
colnames(colData)<-colnames[-1]


#-------------------------------
#- Specify design and contrast
#-------------------------------
default="no"
if(design == 'false'){
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

#-----------------------------------------------------------------------
#- Get counts and create dds
#- Make sure that first column has the sample names, that these samples
#-  are all in countData, and that they are in the same order
#-----------------------------------------------------------------------
counts="featureCounts"
if(kallisto == 'false'){

  countData<-read.table(paste0(inputdir,"/featureCounts/merged_gene_counts.txt"), sep="\t", header=T, check.names=FALSE)
  geneID<-countData$ENSEMBL_ID
  countData<-select(countData, -ENSEMBL_ID)
  rownames(countData)<-geneID 

   if(length(rownames(colData)[!rownames(colData) %in% colnames(countData)]) > 0) {	
    stop("ERROR: the following samples are not in the featureCounts matrix: ", paste(rownames(colData)[!rownames(colData) %in% colnames(countData)], collapse=", "),
         ". Please make sure that the first column of your metadata file has the sample IDs.")
  }else{
    countData = countData[ , rownames(colData) ]
    if(!all(rownames(colData) == colnames(countData))) {	
      stop("Something is wrong, this should never happen [rownames(colData) ne colnames(countData)??]")
    }
  }
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = eval(parse(text=paste0("~ ", design))))
}else{
  counts="kallisto"
  
  files <- file.path(inputdir, "kallisto", rownames(colData), "abundance.h5")
  names(files)<-basename(dirname(files))
  if(length(rownames(colData)[!rownames(colData) %in% names(files)]) > 0) {	
    stop("ERROR: the following samples don't have a kallisto abundance.h5 file: ", paste(rownames(colData)[!rownames(colData) %in% names(files)], collapse=", "),
         ". Please make sure that the first column of your metadata file has the sample IDs.")
  }else{
    tx2gene<-read.table(tx2gene.file, header=T, sep="\t")
    tx2gene_ext<-select(tx2gene, -ENSEMBL_GENE_ID)
    txi.kallisto <- tximport(files, type = "kallisto", tx2gene=tx2gene_ext)
    countData<-txi.kallisto$counts
    countData = countData[ , rownames(colData) ] #- this is not really necessary as this should already be sorted properly
    if(!all(rownames(colData) == colnames(countData))) {	
      stop("Something is wrong, this should never happen [rownames(colData) ne colnames(countData)??]")
    }  
    dds <- DESeqDataSetFromTximport(txi=txi.kallisto,
                                    colData=colData,
                                    design=eval(parse(text=paste0("~ ", design))))
  }
 
}

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
if(ncol(colData) > 50){
	transformation <- rlog(dds, blind=TRUE)
}else{
	trasnformation <- vst(dds, blind=TRUE)
}
pcaData<- plotPCA(transformation, intgroup=colnames(colData), returnData=TRUE)
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

#-----------------------------------------------------------------------
#- Modify one of the functions in ReportingTools so that the norm count 
#- 	plots per gene are in different folders for each contrast
#-	so that plots for a gene that come up in 2 analyses are not overwritten
#- modified from https://rdrr.io/bioc/ReportingTools/src/R/addReportColumns-methods.R
#-----------------------------------------------------------------------
setMethod("modifyReportDF",
          signature = signature(
            object = "DESeqResults"),
          definition = function(df, htmlRep, object, DataSet, factor, 
                             make.plots = TRUE, contrast, ...)
          {
            if("EntrezId" %in% colnames(df)){
              df <- entrezGene.link(df)
            }
            if(make.plots){
              dots <- list(...)
              par.settings <- list()
              if("par.settings" %in% names(dots))
                par.settings <- dots$par.settings
              
              figure.dirname <- paste('figures', htmlRep$shortName,"/",contrast, sep='')  
              figure.directory <- file.path(dirname(path(htmlRep)), 
                                            figure.dirname)
              dir.create(figure.directory, recursive = TRUE)
              
              df <- ReportingTools:::eSetPlot(df, DataSet, factor, figure.directory,
                             figure.dirname, par.settings = par.settings, 
                             ylab.type = "Normalized Counts")
              df
            }
            df
          }
)

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

	reportdir="deseq2_report"
	title=paste0('RNA-seq DGE analysis using DESeq2 with user-defined design and ',counts, ' counts')
	des2Report <- HTMLReport(shortName = 'DESeq2_nextflow_pipeline_results', title = title, reportDirectory = reportdir)

	himg <- hwriteImage(paste0("figuresDESeq2_nextflow_pipeline_results/",pcaplot_name))
	publish(hwrite(himg, br=TRUE, center=T), des2Report)
	publish(paste0("<center><h3><u>Contrast:</u>  ",condition, " - ", treatment, " vs ", control), des2Report)
	publish("<h4>MA/Volcano plots", des2Report)
	himg <- hwriteImage(paste0("figuresDESeq2_nextflow_pipeline_results/",allplot_name))
	publish(hwrite(himg, br=TRUE, center=T), des2Report)
	publish("<h4>Top 100 differentially expressed genes", des2Report)
	publish(resNorm,des2Report, contrast=paste0(condition, "_", treatment, "_vs_", control), pvalueCutoff=1, n=100, DataSet=dds, factor=colData(dds)[[condition]])
	
	finish(des2Report)

}else{

	#-------------------------------------------------------------
	#- Get results using default design and all possible contrasts 
	#-------------------------------------------------------------
	reportdirALL="deseq2_report"
	title=paste0('RNA-seq DGE analysis using DESeq2 with default (no design specified) mode and ',counts, ' counts')
	des2ReportALL <- HTMLReport(shortName = 'DESeq2_nextflow_pipeline_results', title = title, reportDirectory = reportdirALL)

	himg <- hwriteImage(paste0("figuresDESeq2_nextflow_pipeline_results/",pcaplot_name))
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

			#- add gene labels for genes with pval anf fc below theresholds
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
			himg <- hwriteImage(paste0("figuresDESeq2_nextflow_pipeline_results/",allplot_name))
			publish(hwrite(himg, br=TRUE,center=TRUE), des2ReportALL)
			publish("<h4>Top 100 differentially expressed genes", des2ReportALL)
			publish(resNorm,des2ReportALL, contrast = paste0(colnames(colData)[i], "_", paste0(as.character(pairs[,j]),collapse="_vs_")), pvalueCutoff=1, n=100, DataSet=dds, factor=colData(dds)[[i]])
			n=n+1
		}	
	}
	finish(des2ReportALL)
}


