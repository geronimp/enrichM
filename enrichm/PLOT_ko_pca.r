#!/usr/bin/env Rscript
library("optparse")
library("ggplot2")
library("data.table")
library("gridExtra")

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

matrix			= fread(opt$matrix)
metadata		= fread(opt$metadata, header=F)
matrix			= matrix[rowSums(matrix[,2:ncol(matrix), with = FALSE])>0]
matrix[matrix > 0] = 1
model_all			= prcomp(t(matrix[,2:ncol(matrix), with=FALSE]))
model_all_summary	= summary(model_all)
PC_all_weight		= round(model_all_summary$importance[2,1:2] * 100, 1)

PC_all=data.frame(model_all$x[,1:2])
PC_all$genome = rownames(PC_all)
PC_all$Group = metadata[match(PC_all$genome, metadata$V1),]$V2


model_sub			= prcomp(t(matrix[,metadata$V1, with=FALSE]))
model_sub_summary	= summary(model_sub)
PC_sub_weight		= round(model_sub_summary$importance[2,1:2] * 100, 1)

PC_sub=data.frame(model_sub$x[,1:2])
PC_sub$genome = rownames(PC_sub)
PC_sub$Group = metadata[match(PC_sub$genome, metadata$V1),]$V2

PCA_PLOT_1 = ggplot(PC_all,
                aes(PC1,
                       PC2,
                       fill = Group)) +
             geom_hline(yintercept = 0, col = 'gray70', size=0.5) +
             geom_vline(xintercept = 0, col = 'gray70', size=0.5) +
             geom_point(shape = 21,
                        size = 2,
                        stroke = 0.5) +
             xlab(paste(PC_all_weight[1], '% of variance', sep='')) +
             ylab(paste(PC_all_weight[2], '% of variance', sep='')) +
             theme_bw() +
             theme(text = element_text(family 	= 'Arial',
                                       size = 12))

PCA_PLOT_2 = ggplot(PC_sub,
                aes(PC1,
                       PC2,
                       fill = Group)) +
             geom_hline(yintercept = 0, col = 'gray70', size=0.5) +
             geom_vline(xintercept = 0, col = 'gray70', size=0.5) +
             geom_point(shape = 21,
                        size = 2,
                        stroke = 0.5) +
             xlab(paste(PC_sub_weight[1], '% of variance', sep='')) +
             ylab(paste(PC_sub_weight[2], '% of variance', sep='')) +
             theme_bw() +
             theme(text = element_text(family 	= 'Arial',
                                       size = 12))

svg(opt$output, height = 5, width=15)
print(grid.arrange(PCA_PLOT_1, PCA_PLOT_2, nrow = 1))
dev.off()
