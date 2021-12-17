#!/usr/bin/env Rscript


library(ggplot2) 
library(dplyr)


args = commandArgs(trailingOnly=TRUE)

file <- args[1]
perm <- as.numeric(args[2])


data<-read.table(file, sep="\t", header=T)

#- Get range of NES for ylim
rangeNES=round(max(abs(max(data$NES)), abs(min(data$NES)))+0.5)

#- Replace FDR = 0 to FDR = 0.001 (as this was run on 1000 iterations)
data$FDR_p[data$FDR_p == 0] <- (1/perm)

##- no need to readjust pvalues as there is only one ranked list of genes
#data<-data %>%  mutate(FDR=p.adjust(FDR_p, method="BH"))

#- Add numeric variable for GENESET to be able to draw horizontal lines
data <- transform(data, GENESET0 = as.numeric(GENESET))

#- get plot (based on code from https://www.biostars.org/p/168044/)

png("gsea.results.png", width=800)
p <- ggplot(data, aes(NES, GENESET)) + 
    geom_point(aes(colour=FDR_p, size=RATIO)) +
    geom_segment(mapping = aes(yend=GENESET0, xend = 0), size=0.5, colour="gray50") +
    geom_point(aes(colour=FDR_p, size=RATIO)) +
    scale_color_gradient(limits=c(0, 0.05), low="red", high="white") +
    geom_vline(xintercept=0, size=0.5, colour="gray50") +
    theme(strip.text.x = element_text(size = 8), 
	  panel.background=element_rect(fill="gray95", colour="gray95"),
          panel.grid.major=element_line(size=0.25,linetype='solid', colour="gray90"), 
          panel.grid.minor=element_line(size=0.25,linetype='solid', colour="gray90"),
          axis.title.y=element_blank()) +
    expand_limits(x=c(-rangeNES,rangeNES)) +
    scale_x_continuous(breaks=seq(-rangeNES, rangeNES, 2)) +
    facet_grid(.~RANK)
print(p)
dev.off()


