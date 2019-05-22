#!/usr/bin/env Rscript
library("optparse")
library("ggplot2")
library("data.table")
library("gridExtra")
library("grid")

option_list = list(
  make_option(c("-i", "--input"), type="character",
              help="Input KO matrix", metavar="character"),
  make_option(c("-k", "--ko"), type="character",
              help="", metavar="character"),
  make_option(c("-p", "--pvalue"), type="character",
              help="Input KO matrix", metavar="character"),
    make_option(c("-o", "--output"), type="character",
              help="output file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

ko00000 = fread(opt$k, header=F)
input = fread(opt$i)

generate_plots = function(res_1, res_2, ko00000, header, output){

    res_1_ko00000 = ko00000[match(res_1$annotation, ko00000$V4)]
    res_1_match = res_1[match(res_1_ko00000$V4, res_1$annotation)]
    res_1_ko00000$pvalue = res_1_match$pvalue
    res_1_ko00000$ratio = res_1_match$group_1_true/(res_1_match$group_1_true + res_1_match$group_1_false)
    res_1_ko00000_V1 = data.frame(table(res_1_ko00000$V1))
    res_1_ko00000_V2 = data.frame(table(res_1_ko00000$V2))
    res_1_ko00000_V3 = data.frame(table(res_1_ko00000$V3))

    res_2_ko00000 = ko00000[match(res_2$annotation, ko00000$V4)]
    res_2_match = res_2[match(res_2_ko00000$V4, res_2$annotation)]
    res_2_ko00000$pvalue = res_2_match$pvalue
    res_2_ko00000$ratio = res_2_match$group_2_true/(res_2_match$group_2_true + res_2_match$group_2_false)
    res_2_ko00000_V1 = data.frame(table(res_2_ko00000$V1))
    res_2_ko00000_V2 = data.frame(table(res_2_ko00000$V2))
    res_2_ko00000_V3 = data.frame(table(res_2_ko00000$V3))
    print(nrow(res_1))
    print(nrow(res_2))
    h = paste(header[1], header[2], sep='_')
    if (nrow(res_1)>0) {
        p1 = ggplot(res_1_ko00000_V1, aes(Var1, Freq)) +
            geom_bar(fill = 'black', stat='identity') +
            theme_bw(base_family = "Arial",
                     base_size = 12) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  text = element_text(family = 'Arial', size = 12)) +
            xlab("") +
            ggtitle(paste(header[1] ," KEGG heirarchy level 1")) +
            ylab("Number of KOs")

        p2 = ggplot(res_1_ko00000_V2, aes(Var1, Freq)) +
            geom_bar(fill = 'black', stat='identity') +
            theme_bw(base_family = "Arial",
                     base_size = 12) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  text = element_text(family = 'Arial', size = 12)) +
            xlab("") +
            ggtitle(paste(header[1] ," KEGG heirarchy level 2")) +
            ylab("Number of KOs")

        p3 = ggplot(res_1_ko00000_V3, aes(Var1, Freq)) +
            geom_bar(fill = 'black', stat='identity') +
            theme_bw(base_family = "Arial",
                     base_size = 12) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  text = element_text(family = 'Arial', size = 12)) +
            xlab("") +
            ggtitle(paste(header[1] ," KEGG heirarchy level 3")) +
            ylab("Number of KOs")

        p4 = ggplot(res_1_ko00000, aes(V3, pvalue, size = ratio)) +
            geom_point() +
            theme_bw(base_family = "Arial",
                     base_size = 12) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  text = element_text(family = 'Arial', size = 12)) +
            xlab("") +
            scale_y_log10() +
            ggtitle(paste(header[1] ," KEGG heirarchy level 3")) +
            ylab("Corrected p-value")
    } else {
        p1 = ggplot()
        p2 = ggplot()
        p3 = ggplot()
        p4 = ggplot()
    }
    if (nrow(res_2)>0) {
        p5 = ggplot(res_2_ko00000_V1, aes(Var1, Freq)) +
            geom_bar(fill = 'black', stat='identity') +
            theme_bw(base_family = "Arial",
                     base_size = 12) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  text = element_text(family = 'Arial', size = 12)) +
            xlab("") +
            ggtitle(paste(header[2] ," KEGG heirarchy level 1")) +
            ylab("Number of KOs")

        p6 = ggplot(res_2_ko00000_V2, aes(Var1, Freq)) +
            geom_bar(fill = 'black', stat='identity') +
            theme_bw(base_family = "Arial",
                     base_size = 12) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  text = element_text(family = 'Arial', size = 12)) +
            xlab("") +
            ggtitle(paste(header[2] ," KEGG heirarchy level 2")) +
            ylab("Number of KOs")

        p7 = ggplot(res_2_ko00000_V3, aes(Var1, Freq)) +
            geom_bar(fill = 'black', stat='identity') +
            theme_bw(base_family = "Arial",
                     base_size = 12) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  text = element_text(family = 'Arial', size = 12)) +
            xlab("") +
            ggtitle(paste(header[2] ," KEGG heirarchy level 3")) +
            ylab("Number of KOs")

        p8 = ggplot(res_2_ko00000, aes(V3, pvalue, size = ratio)) +
            geom_point() +
            theme_bw(base_family = "Arial",
                     base_size = 12) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  text = element_text(family = 'Arial', size = 12)) +
            xlab("") +
            scale_y_log10() +
            ggtitle(paste(header[2] ," KEGG heirarchy level 3")) +
            ylab("Corrected p-value")
    } else {
        p5 = ggplot()
        p6 = ggplot()
        p7 = ggplot()
        p8 = ggplot()
    }

    svg(paste(output, paste(h, "KEGG_heirarchy_level_1.svg", sep='_'), sep = '/'), height = 18, width = 15)
    grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p5), size="last"))
    dev.off()
    svg(paste(output, paste(h, "KEGG_heirarchy_level_2.svg", sep='_'), sep = '/'), height = 18, width = 15)
    grid.draw(rbind(ggplotGrob(p2), ggplotGrob(p6), size="last"))
    dev.off()
    svg(paste(output, paste(h, "KEGG_heirarchy_level_3.svg", sep='_'), sep = '/'), height = 18, width = 15)
    grid.draw(rbind(ggplotGrob(p3), ggplotGrob(p7), size="last"))
    dev.off()
    svg(paste(output, paste(h, "KEGG_heirarchy_level_3_evals.svg", sep='_'), sep = '/'), height = 18, width = 15)
    grid.draw(rbind(ggplotGrob(p4), ggplotGrob(p8), size="last"))
    dev.off()

}

hq_input = input[input$pvalue<as.numeric(opt$p)]
combos=paste(hq_input$group_1, hq_input$group_2, sep='_')
hq_input$combos = combos

for (i in unique(combos)) {
    i_hq_input = hq_input[hq_input$combos==i,]
    if (grepl("cdf.tsv", opt$i)) {
        g1_res = i_hq_input[i_hq_input$group_1_mean>i_hq_input$group_2_count]
        g2_res = i_hq_input[i_hq_input$group_1_mean<i_hq_input$group_2_count]
    }
    if (grepl("fisher.tsv", opt$i)) {
        g1_res = i_hq_input[i_hq_input$group_1_true>i_hq_input$group_2_true]
        g2_res = i_hq_input[i_hq_input$group_1_true<i_hq_input$group_2_true]
    }
    header = paste(i_hq_input[1,2:3, with = FALSE])
    generate_plots(g1_res, g2_res, ko00000, header, opt$o)
}

