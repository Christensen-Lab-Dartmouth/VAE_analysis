#####################
# Correlation testing between VAE nodes and top 300K CpGs
#
#
#####################

#####################
# Set up the environment
#####################
     require(data.table)
     require(limma)

#####################
# Load the data
#####################
     #dir = '/Users/alexandertitus/Documents/github/DNAm_data_generation/code'
     
     dir = 'C:/Users/atitus/Documents/github/DNAm_data_generation/results'
     BRCA.covFile = '../BRCAtarget_covariates.csv'

     setwd(dir)

     BRCA.covs = data.frame(fread(BRCA.covFile), row.names=1)
     BRCA.covs$sample.typeInt = ifelse(BRCA.covs$sample.type == 'Solid Tissue Normal', 0, 1)
     BRCA.covs = BRCA.covs[!is.na(BRCA.covs$age.Dx), ]
     
     dir = 'C:/Users/atitus/Documents/github/DNAm_data_generation'
     setwd(dir)
     
     ## Betas
     beta.file = 'data/TCGA_BRCA_top300kMAD_cpg.tsv'
     betas = data.frame(fread(beta.file))
     rownames(betas) = betas[,1]
     betas = betas[,2:ncol(betas)]
     betas = betas[rownames(betas) %in% BRCA.covs$Basename, ]
     betas = betas[order(rownames(betas), decreasing=T), ]
     
     
     BRCA.covs = BRCA.covs[order(BRCA.covs$Basename, decreasing=T), ]
     
     ## VAE nodes
     vae.file = 'results/encoded_methyl_onehidden_warmup_batchnorm_300K-100.tsv'
     vae = data.frame(fread(vae.file))
     colnames(vae) = vae[1,]
     rownames(vae) = vae[,1]
     vae = vae[2:nrow(vae), 2:ncol(vae)]
     vae = vae[rownames(vae) %in% BRCA.covs$Basename,]
     
     vae = vae[order(rownames(vae), decreasing=T), ]
     
     ## Check sample concordance
     all(rownames(vae) == rownames(betas))
     all(BRCA.covs$Basename == rownames(betas))
     
     
#####################
# Correlations
#####################  
     node = 1
     vaeNode = vae[, as.character(node)]
     
     cor.func = function(x){return(cor(x, vaeNode))}
     betas = t(betas)
     
     correlations = apply(betas, 2, cor.func)
     correlations = data.frame(correlations)
     correlations = cbind('CpG' = rownames(correlations), correlations)
     
     nodeLabel = paste('VAE', node, sep = '')
     correlations = cbind(nodeLabel = rep(nodeLabel, nrow(correlations)), correlations)
     
     #correlations40 = correlations[order(abs(correlations$correlations), decreasing=T), ]
     #correlations42 = correlations[order(abs(correlations$correlations), decreasing=T), ]
     #correlations93 = correlations[order(abs(correlations$correlations), decreasing=T), ]
     #correlations1 = correlations[order(abs(correlations$correlations), decreasing=T), ]
     
     betas = t(betas)     
     
     
#####################
# Correlation diagram
##################### 
     ## Correlation elbow
     cors = correlations40$correlations
     file.name = paste('correlation_elbow_node40', '.png', sep = '')
     png(file.name, width = 1000, height = 1000, res = 100)
     plot(abs(cors), main=paste('Node ', node, sep = ''), ylab = '|Correlation|', xlab='Index')
     dev.off()
     
     cors = correlations42$correlations
     file.name = paste('correlation_elbow_node42', '.png', sep = '')
     png(file.name, width = 1000, height = 1000, res = 100)
     plot(abs(cors), main=paste('Node ', node, sep = ''), ylab = '|Correlation|', xlab='Index')
     dev.off()
     
     cors = correlations93$correlations
     file.name = paste('correlation_elbow_node93', '.png', sep = '')
     png(file.name, width = 1000, height = 1000, res = 100)
     plot(abs(cors), main=paste('Node ', node, sep = ''), ylab = '|Correlation|', xlab='Index')
     dev.off()
     
     cors = correlations1$correlations
     file.name = paste('correlation_elbow_node1', '.png', sep = '')
     png(file.name, width = 1000, height = 1000, res = 100)
     plot(abs(cors), main=paste('Node ', node, sep = ''), ylab = '|Correlation|', xlab='Index')
     dev.off()
     
     
#####################
# Correlation diagram
##################### 
     ## Correlation with methylation
     cors = correlations1$correlations
     cpg = correlations1$CpG[1]
     cg14789818 = betas[,cpg]
     node = 40
     vaeNode = vae[, as.character(node)]

     df = data.frame(cbind(vaeNode, cg14789818))
     BRCA.covsSub = BRCA.covsSub[order(BRCA.covsSub$Basename, decreasing=T), ]
     all(BRCA.covsSub$Basename == rownames(df))
     df = cbind(df, BRCA.covsSub$sample.typeInt)
     colnames(df) = c('activation', 'beta', 'sample')
     df$sample = factor(df$sample)
     
     file.name = paste('methylation_V_node', node, '.png', sep = '')
     png(file.name, width = 1000, height = 1000, res = 100)
     
     ggscatter(df, x = 'activation', y = 'beta',
               color = 'sample', ylab = 'CpG cg13985135 beta value',
               xlab = 'Activation of node 40',
               add = "reg.line",  # Add regressin line
               add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
               conf.int = TRUE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
               cor.coeff.args = list(method = "pearson", label.x = 0.35, label.y = 0.35, label.sep = "\n"))

     dev.off()
     
     
#####################
# Chord diagram 
#####################
     library(circlize)
     
     threshold = 0.6
     
     correlations = correlations1
     summary(correlations)
     
     cor.sub = correlations[(correlations$correlations >= threshold | 
                                  correlations$correlations < -1*threshold), ]
     
     file.name = paste('VAE', node, '_ChordPlot.png', sep = '')
     
     png(file.name, width = 2000, height = 2000, res = 300)
          
          col_fun = colorRamp2(range(cor.sub$correlations), c("yellow", "blue"), 
                               transparency = 0.5)
          
          chordDiagram(cor.sub, 
                       col = col_fun,
                       link.sort = TRUE, link.decreasing = TRUE,
                       annotationTrack = c("grid"),
                       preAllocateTracks = 1)
          
          circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
               xlim = get.cell.meta.data("xlim")
               ylim = get.cell.meta.data("ylim")
               sector.name = get.cell.meta.data("sector.index")
               circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
               circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
          }, bg.border = NA)
          
          
          title(paste("Correlations (>=", threshold , ") between CpGs and VAE node ",  node, sep = ''), cex = 0.8)
          
     dev.off()
     
     circos.clear()
     
     
#####################
# Genomic context
#####################
     ## https://www.bioconductor.org/help/workflows/methylationArrayAnalysis/
     #biocLite('missMethyl')
     library(missMethyl)
     library(Gviz)
     require(minfi)
     
     threshold = 0.5
     node = 'All'
     #correlations = correlations1
     #correlations = correlations40
     #correlations = correlations42
     #correlations = correlations93

     summary(correlations)
     
     cor.sub = correlations[(correlations$correlations >= threshold | 
                                  correlations$correlations < -1*threshold), ]
     
     
     #cpgs = rownames(cor.sub)
     #cpg_set = cor.sub
     cpg_set = rbind(cpg_set, cor.sub)
     
     cpgs = cpg_set$CpG     

     ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
     head(ann450k)
     
     anno.sub = ann450k[rownames(ann450k) %in% cpgs, ]
     anno.sub = data.frame(anno.sub)
     anno.sub = cbind('NodeCor' = cor.sub$correlations, anno.sub)
     
     file.name = paste('../results/anno450K_node', node, '.csv', sep = '')
     write.csv(anno.sub, file.name)     
     
     
     ## GSA & GO analysis
     # load Broad human curated (C2) gene sets
     # http://bioinf.wehi.edu.au/software/MSigDB/
     load('human_c2_v5p2.rdata')
     
     # analysis
     par(mfrow=c(1,1))
     all = colnames(betas)
     gst <- gometh(sig.cpg=cpgs, all.cpg=all, plot.bias=TRUE)
     gsa <- gsameth(sig.cpg=cpgs, all.cpg=all, collection=Hs.c2)
     
     # Top 10 GO categories
     go = topGO(gst, number = 30); go

     # Top 10 gene sets
     gsa = topGSA(gsa, number=20); gsa
     
     gofile.name = paste('results/go_anno_node', node, '.csv', sep = '')
     write.csv(go, gofile.name) 
     
     gsafile.name = paste('results/gsa_anno_node', node, '.csv', sep = '')
     write.csv(gsa, gsafile.name) 
     
     
#####################
# Enhancer calculations
#####################
     threshold = 0.65
     node = 13
     #correlations = correlations1
     #correlations = correlations40
     #correlations = correlations42
     #correlations = correlations93
     correlations = correlations13
     
     summary(correlations)
     
     cor.sub = correlations[(correlations$correlations >= threshold | 
                                  correlations$correlations < -1*threshold), ]
     
     
     cpgs = rownames(cor.sub)
     #cpg_set = cor.sub
     cpg_set = rbind(cpg_set, cor.sub)
     
     #cpg_set.unique = cpg_set[unique(cpg_set$CpG), ]
     
     ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
     head(ann450k)
     
     anno.sub = ann450k[rownames(ann450k) %in% cpgs, ]
     anno.sub = data.frame(anno.sub)
     anno.sub = cbind('NodeCor' = cor.sub$correlations, anno.sub)
     
     nodeEnhancer = nrow(anno.sub[anno.sub$Enhancer == 'TRUE', ])
     nodeNoEnhancer = nrow(anno.sub[anno.sub$Enhancer == '', ])
     annoEnhancer = nrow(ann450k[ann450k$Enhancer == 'TRUE', ])
     annoNoEnhancer = nrow(ann450k[ann450k$Enhancer == '', ])
     
     enhancers <- matrix(c(nodeEnhancer, annoEnhancer, 
                           nodeNoEnhancer, annoNoEnhancer), nrow = 2,
                           dimnames =
                                list(c("NodeRelated", "NotNodeRelated"),
                                     c("Enhancer", "NotEnhancer")))
     
     node.results = fisher.test(enhancers)
     
     #enhancerResults = data.frame('Node',
                                  'Est', 
                                  'Conf95low', 
                                  'Conf95high', 
                                  'Pvalue')
     
     enhancerResults = rbind(enhancerResults,
                             c(node,
                               node.results$estimate, 
                               node.results$conf.int[1],
                               node.results$conf.int[2],
                               node.results$p.value))
     
     #colnames(enhancerResults) = enhancerResults[1, ]
     #enhancerResults = enhancerResults[2:nrow(enhancerResults), ]
     enhancerResults$Est = as.numeric(enhancerResults$Est)
     enhancerResults$Conf95low = as.numeric(enhancerResults$Conf95low)
     enhancerResults$Conf95high = as.numeric(enhancerResults$Conf95high)
     enhancerResults$Pvalue = as.numeric(enhancerResults$Pvalue)
     
     
#####################
# ggplot theme
#####################   
     theme_Publication <- function(base_size=14, base_family="helvetica") {
          library(grid)
          library(ggthemes)
          (theme_foundation(base_size=base_size, base_family=base_family)
               + theme(plot.title = element_text(face = "bold",
                                                 size = rel(1.2), hjust = 0.5),
                       text = element_text(),
                       panel.background = element_rect(colour = NA),
                       plot.background = element_rect(colour = NA),
                       panel.border = element_rect(colour = NA),
                       axis.title = element_text(face = "bold",size = rel(1)),
                       axis.title.y = element_text(angle=90,vjust =2),
                       axis.title.x = element_text(vjust = -0.2),
                       axis.text = element_text(), 
                       axis.line = element_line(colour="black"),
                       axis.ticks = element_line(),
                       panel.grid.major = element_line(colour="#f0f0f0"),
                       panel.grid.minor = element_blank(),
                       legend.key = element_rect(colour = NA),
                       legend.position = "bottom",
                       legend.direction = "horizontal",
                       legend.key.size= unit(0.2, "cm"),
                       legend.margin = unit(0, "cm"),
                       legend.title = element_text(face="italic"),
                       plot.margin=unit(c(10,5,5,5),"mm"),
                       strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                       strip.text = element_text(face="bold")
               ))
          
     }
     
     scale_fill_Publication <- function(...){
          library(scales)
          discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
          
     }
     
     scale_colour_Publication <- function(...){
          library(scales)
          discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
          
     }
     
     
#####################
# Enhancer OR
##################### 
     fp <- ggplot(data=enhancerResults, aes(x=Est, 
                                            y=as.character(Node), 
                                            xmin=Conf95low, 
                                            xmax=Conf95high)) +
          geom_point(color = 'black') +
          geom_text(aes(label=format(round(Est, 2), nsmall = 2)), hjust=0, vjust=-1) +
          geom_errorbarh(height=.02) +
          ylab('Node') +
          geom_vline(xintercept=1, color='black', linetype='dashed') +
          scale_x_continuous(limits=c(0,max(enhancerResults$Conf95high)), name='Odds ratio w/ 95% CI') +
          ggtitle('OR for enhancer in each node set of CpGs') + theme_Publication()
     
     png('OR_enhancer_by_node.png', width = 2000, height = 2000, res = 300)
     fp
     dev.off()
     
     
#####################
# Visualize genomic context
#####################

     ## Visualization
     #http://zuguang.de/circlize_book/book/high-level-genomic-functions.html#genomic-heatmap
     cpg_set.unique = cpg_set[unique(cpg_set$CpG), ]
     
     anno.sub = ann450k[which(rownames(ann450k) %in% cpg_set.unique$CpG), ]
     anno.sub = anno.sub[order(rownames(anno.sub)), ]
     cpg_set.unique = cpg_set.unique[order(rownames(cpg_set.unique)), ]
     all(rownames(anno.sub) == rownames(cpg_set.unique))
     
     bed = cbind(anno.sub$chr, 
                 anno.sub$pos, 
                 anno.sub$pos + 10,
                 cpg_set.unique$correlations,
                 rownames(anno.sub),
                 cpg_set.unique$CpG,
                 cpg_set.unique$nodeLabel)
     bed = data.frame(bed)
     colnames(bed) = c('chr', 'start', 'end', 'value1', 'AnnoCpG', 'CorrCpG', 'Node')
     all(bed$AnnoCpG == bed$CorrCpG)
     
     bed$start = as.numeric(bed$start)
     bed$end = as.numeric(bed$end)
     bed$value1 = as.numeric(bed$value1)
     bed = bed[order(bed$chr), ]
     
     bed_list = list('VAE1' = bed[bed$Node == 'VAE1', ], 
                     'VAE40' = bed[bed$Node == 'VAE40', ],
                     'VAE42' = bed[bed$Node == 'VAE42', ],
                     'VAE93' = bed[bed$Node == 'VAE93', ])
     
     circlize_plot = function() {
          circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr", 1:22))
          circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("#FF000080", 
                                                                        "#0000FF80", 
                                                                        "orangered1",
                                                                        "lightsteelblue"))
          
          circos.genomicDensity(bed[bed$Node == 'VAE1', ], col = c("#FF000080"), track.height = 0.1)
          circos.genomicDensity(bed[bed$Node == 'VAE40', ], col = c("#0000FF80"), track.height = 0.1)
          circos.genomicDensity(bed[bed$Node == 'VAE42', ], col = c("orangered1"), track.height = 0.1)
          circos.genomicDensity(bed[bed$Node == 'VAE93', ], col = c("lightsteelblue"), track.height = 0.1)
          circos.clear()
     }
     
     
     ## Corr > 0.6
     png('genomic_context_corr0.6.png', width = 2000, height = 2000, res = 300)
     circlize_plot()
     dev.off()
     
     out.data = merge(anno.sub, bed, by.x = 'Name', by.y = 'AnnoCpG')
     write.csv(out.data, file = 'results/genomic_context_anno_corr0.6.csv')
     
     ## Corr > 0.75
     png('results/genomic_context_corr0.75.png', width = 2000, height = 2000, res = 300)
     circlize_plot()
     dev.off()
     
     out.data = merge(anno.sub, bed, by.x = 'Name', by.y = 'AnnoCpG')
     write.csv(out.data, file = 'results/genomic_context_anno_corr0.75.csv')
          
     
#####################
# VAE activations
#####################
     vae.sub = vae[, c(1, 40, 42, 93)]
     vae.sub$SubjectID = rownames(vae.sub)
     VAE.covs = merge(vae.sub, BRCA.covs, by.x = 'SubjectID', by.y = 'Basename')
     temp = melt(VAE.covs[, c(1, 40, 42, 93, 'sample.typeInt', "sample.type", "PAM50.RNAseq")], 
                 id=c('sample.typeInt', 'sample.type', 'PAM50.RNAseq'))
     
     colnames(temp) = c('SampleTypeInt', 'SampleType', 'PAM50orNormal', 'Node', 'Activation')
     temp$SampleTypeInt = factor(temp$SampleTypeInt)
     
     library(ggpubr)
     png('Node_activationsInt.png', width = 2000, height = 2000, res = 300)
     my_comparisons <- list( c("1", "40"), c("40", "42"), c("42", "93"), 
                             c("1", "42"), c("1", "93"))
     ggboxplot(temp, x = "Node", y = "Activation", 
               color = "SampleTypeInt") +  
          stat_compare_means(label.y = 4)  
     dev.off()
     
     temp = temp[temp$PAM50 != '', ]
     png('Node_activationsPAM50.png', width = 2000, height = 2000, res = 300)
     ggboxplot(temp, x = "Node", y = "Activation", 
               color = "PAM50orNormal") +  
          stat_compare_means(label.y = 4)  
     dev.off()
   
     
     

#####################
# VAE activations
#####################