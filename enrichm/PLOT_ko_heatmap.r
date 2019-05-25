#!/usr/bin/env Rscript
library("optparse")
library("ggplot2")
library("data.table")

option_list = list(
  make_option(c("-i", "--matrix"), type="character",
              help="Input KO matrix", metavar="character"),
    make_option(c("-m", "--metadata"), type="character",
              help="Input metadata", metavar="character"),
    make_option(c("-o", "--output"), type="character",
              help="output file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

matrix 		= fread(opt$matrix)
metadata	= fread(opt$metadata, header=F)
model 		= prcomp(t(matrix[,2:ncol(matrix)]))
model_summary = summary(model)
PC_weight	= round(model_summary$importance[2,1:2] * 100, 1)

PC=data.frame(model$x[,1:2])
PC$genome = rownames(PC)
PC$Group = metadata[match(PC$genome, metadata$V1),]$V2

PCA_PLOT =	ggplot(PC,
                   aes(PC1,
                          PC2,
                          fill = Group)) +
                geom_point(shape = 21,
                           size = 2,
                           stroke = 0.5) +
                geom_hline(yintercept = 0, col = 'gray70', size=0.5) +
                geom_vline(xintercept = 0, col = 'gray70', size=0.5) +
                xlab(paste(PC_weight[1], '% of variance', sep='')) +
                ylab(paste(PC_weight[2], '% of variance', sep='')) +
                theme_bw() +
                theme(text = element_text(size		= 12,
                                          family 	= 'Arial'))

svg(opt$output, height = 3, width=5)
print(PCA_PLOT)
dev.off()
