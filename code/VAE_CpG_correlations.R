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

     node_corrs = function(node, betas, vae) {
          vaeNode = vae[, as.character(node)]
          
          cor.func = function(x){return(cor(x, vaeNode))}
          
          correlations = apply(betas, 2, cor.func)
          correlations = data.frame(correlations)
          correlations = cbind('CpG' = rownames(correlations), correlations)
          
          nodeLabel = paste('VAE', node, sep = '')
          correlations = cbind(nodeLabel = rep(nodeLabel, nrow(correlations)), correlations)
          results = correlations[order(abs(correlations$correlations), decreasing=T), ]
          
          return(results)
     }

     subset_cors = function(threshold, correlations){
          cor.sub = correlations[(correlations$correlations >= threshold | 
                                       correlations$correlations < -1*threshold), ]
          return(cor.sub)
     }
     
     threshold = 0.5
     
     correlations1 = node_corrs(1, betas, vae)
     correlations1 = subset_cors(threshold, correlations1)
     
     correlations8 = node_corrs(8, betas, vae)
     correlations8 = subset_cors(threshold, correlations8)
     
     correlations13 = node_corrs(13, betas, vae)
     correlations13 = subset_cors(threshold, correlations13)
     
     correlations16 = node_corrs(16, betas, vae)
     correlations16 = subset_cors(threshold, correlations16)
     
     correlations33 = node_corrs(33, betas, vae)
     correlations33 = subset_cors(threshold, correlations33)
     
     correlations40 = node_corrs(40, betas, vae)
     correlations40 = subset_cors(threshold, correlations40)
     
     correlations48 = node_corrs(48, betas, vae)
     correlations48 = subset_cors(threshold, correlations48)
     
     correlations49 = node_corrs(49, betas, vae)
     correlations49 = subset_cors(threshold, correlations49)
     
     correlations55 = node_corrs(55, betas, vae)
     correlations55 = subset_cors(threshold, correlations55)
     
     correlations56 = node_corrs(56, betas, vae)
     correlations56 = subset_cors(threshold, correlations56)
     
     correlations58 = node_corrs(58, betas, vae)
     correlations58 = subset_cors(threshold, correlations58)
     
     correlations66 = node_corrs(66, betas, vae)
     correlations66 = subset_cors(threshold, correlations66)
     
     correlations68 = node_corrs(68, betas, vae)
     correlations68 = subset_cors(threshold, correlations68)
     
     correlations86 = node_corrs(86, betas, vae)
     correlations86 = subset_cors(threshold, correlations86)
     
     correlations93 = node_corrs(93, betas, vae)
     correlations93 = subset_cors(threshold, correlations93)
     
     
#####################
# Correlation diagram
##################### 
     plot_corr_elbow = function(correlationsN, node, threshold, plot_line = F, line = NA){
          ## Correlation elbow
          cors = correlationsN$correlations
          file.name = paste('correlation_elbow_node', node, '.png', sep = '')
          
          png(file.name, width = 1000, height = 1000, res = 100)
          plot(abs(cors), main=paste('Node ', node, sep = ''), 
               ylab = '|Correlation|', xlab='Index', ylim = c(threshold, 0.9))
          if(plot_line == T){
               abline(h = line)   
               abline(v = 1000, col = "black")
               abline(v = 5000, col = "green")
          }
          dev.off()
     }
     
     plot_corr_elbow(correlations1, 1, threshold, plot_line = T, line = 0.75)
     plot_corr_elbow(correlations8, 8, threshold, plot_line = T, line = 0.62)
     plot_corr_elbow(correlations13, 13, threshold, plot_line = T, line = 0.65)
     plot_corr_elbow(correlations16, 16, threshold, plot_line = T, line = 0.75)
     plot_corr_elbow(correlations33, 33, threshold, plot_line = T, line = 0.6)
     plot_corr_elbow(correlations40, 40, threshold, plot_line = T, line = 0.7)
     plot_corr_elbow(correlations48, 48, threshold, plot_line = T, line = 0.72)
     plot_corr_elbow(correlations49, 49, threshold, plot_line = T, line = 0.8)
     plot_corr_elbow(correlations55, 55, threshold, plot_line = T, line = 0.63)
     plot_corr_elbow(correlations56, 56, threshold, plot_line = T, line = 0.75)
     plot_corr_elbow(correlations58, 58, threshold, plot_line = T, line = 0.8)
     plot_corr_elbow(correlations66, 66, threshold, plot_line = T, line = 0.65)
     plot_corr_elbow(correlations68, 68, threshold, plot_line = T, line = 0.57)
     plot_corr_elbow(correlations86, 86, threshold, plot_line = T, line = 0.77)
     plot_corr_elbow(correlations93, 93, threshold, plot_line = T, line = 0.6)
     
     
#####################
# Correlation diagram
##################### 
     ## Correlation with methylation
     
     # cors = correlations1$correlations
     # cpg = correlations1$CpG[1]
     # CpGprobe = betas[,cpg]
     # node = 40
     # vaeNode = vae[, as.character(node)]
     # 
     # df = data.frame(cbind(vaeNode, cg14789818))
     # BRCA.covsSub = BRCA.covsSub[order(BRCA.covsSub$Basename, decreasing=T), ]
     # all(BRCA.covsSub$Basename == rownames(df))
     # df = cbind(df, BRCA.covsSub$sample.typeInt)
     # colnames(df) = c('activation', 'beta', 'sample')
     # df$sample = factor(df$sample)
     # 
     # file.name = paste('methylation_V_node', node, '.png', sep = '')
     # png(file.name, width = 1000, height = 1000, res = 100)
     # 
     # ggscatter(df, x = 'activation', y = 'beta',
     #           color = 'sample', ylab = 'CpG cg13985135 beta value',
     #           xlab = 'Activation of node 40',
     #           add = "reg.line",  # Add regressin line
     #           add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
     #           conf.int = TRUE, # Add confidence interval
     #           cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
     #           cor.coeff.args = list(method = "pearson", label.x = 0.35, label.y = 0.35, label.sep = "\n"))
     # 
     # dev.off()
     
          # library(circlize)
     # 
     # threshold = 0.6
     # 
     # correlations = correlations1
     # summary(correlations)
     # 
     # cor.sub = correlations[(correlations$correlations >= threshold | 
     #                              correlations$correlations < -1*threshold), ]
     # 
     # file.name = paste('VAE', node, '_ChordPlot.png', sep = '')
     # 
     # png(file.name, width = 2000, height = 2000, res = 300)
     #      
     #      col_fun = colorRamp2(range(cor.sub$correlations), c("yellow", "blue"), 
     #                           transparency = 0.5)
     #      
     #      chordDiagram(cor.sub, 
     #                   col = col_fun,
     #                   link.sort = TRUE, link.decreasing = TRUE,
     #                   annotationTrack = c("grid"),
     #                   preAllocateTracks = 1)
     #      
     #      circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
     #           xlim = get.cell.meta.data("xlim")
     #           ylim = get.cell.meta.data("ylim")
     #           sector.name = get.cell.meta.data("sector.index")
     #           circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
     #           circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
     #      }, bg.border = NA)
     #      
     #      
     #      title(paste("Correlations (>=", threshold , ") between CpGs and VAE node ",  node, sep = ''), cex = 0.8)
     #      
     # dev.off()
     # 
     # circos.clear()
     # 
#####################
# Chord diagram 
#####################

     
#####################
# Genomic context
#####################
     ## https://www.bioconductor.org/help/workflows/methylationArrayAnalysis/
     #install.packages('matrixStats')
     #source("https://bioconductor.org/biocLite.R")
     #biocLite('missMethyl')
     library(missMethyl)
     library(Gviz)
     require(minfi)
     library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
     
     ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
     head(ann450k)
     
     go_pathway_analysis = function(correlationsN, node, threshold, annotations){
          correlations = correlationsN
          
          cor.sub = correlations[(correlations$correlations >= threshold | 
                                       correlations$correlations < -1*threshold), ]
          
          cpgs = rownames(cor.sub)
          cpg_set = cor.sub
          cpg_set = rbind(cpg_set, 1, cor.sub)
          
          cpgs = cpg_set$CpG     
          
          ann450k = annotations
          
          anno.sub = ann450k[rownames(ann450k) %in% cpgs, ]
          anno.sub = data.frame(anno.sub)
          anno.sub = cbind('NodeCor' = cor.sub$correlations, anno.sub)
          
          file.name = paste('anno450K_node', node, '.csv', sep = '')
          write.csv(anno.sub, file.name)     
          
          
          ## GSA & GO analysis
          # load Broad human curated (C2) gene sets
          # http://bioinf.wehi.edu.au/software/MSigDB/
          load('../human_c2_v5p2.rdata')
          
          # analysis
          par(mfrow=c(1,1))
          all = colnames(betas)
          cpgs = as.character(cpgs)
          gst <- gometh(sig.cpg=cpgs, all.cpg=all, plot.bias=TRUE)
          gsa <- gsameth(sig.cpg=cpgs, all.cpg=all, collection=Hs.c2)
          
          # Top 10 GO categories
          go = topGO(gst, number = 30)
          
          # Top 10 gene sets
          gsa = topGSA(gsa, number=20)
          
          gofile.name = paste('go_anno_node', node, '.csv', sep = '')
          write.csv(go, gofile.name) 
          
          gsafile.name = paste('gsa_anno_node', node, '.csv', sep = '')
          write.csv(gsa, gsafile.name) 
     }

     go_pathway_analysis(correlations1, 1, threshold = 0.75, ann450k)
     go_pathway_analysis(correlations8, 8, threshold = 0.62, ann450k)
     go_pathway_analysis(correlations13, 13, threshold = 0.65, ann450k)
     go_pathway_analysis(correlations16, 16, threshold = 0.75, ann450k)
     go_pathway_analysis(correlations33, 33, threshold = 0.6, ann450k)
     go_pathway_analysis(correlations40, 40, threshold = 0.7, ann450k)
     go_pathway_analysis(correlations48, 48, threshold = 0.72, ann450k)
     go_pathway_analysis(correlations49, 49, threshold = 0.8, ann450k)
     go_pathway_analysis(correlations55, 55, threshold = 0.63, ann450k)
     go_pathway_analysis(correlations56, 56, threshold = 0.75, ann450k)
     go_pathway_analysis(correlations58, 58, threshold = 0.8, ann450k)
     go_pathway_analysis(correlations66, 66, threshold = 0.65, ann450k)
     go_pathway_analysis(correlations68, 68, threshold = 0.57, ann450k)
     go_pathway_analysis(correlations86, 86, threshold = 0.77, ann450k)
     go_pathway_analysis(correlations93, 93, threshold = 0.6, ann450k)
         
     
#####################
# Enhancer calculations
#####################

     ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
     
     enhancer_analysis = function(correlations_list, node_list, threshold, annotations){
          
          enhancerResults = c('Node', 'Est', 'Conf95low', 'Conf95high', 'Pvalue')
          
          for( i in 1:length(correlations_list) ){
               correlations = correlations_list[[i]]
               
               nodeName = paste('VAE', node_list[i], sep='')

               cor.sub = correlations[(correlations$correlations >= threshold | 
                                            correlations$correlations < -1*threshold), ]
               
               
               cpgs = rownames(cor.sub)
               cpg_set = cor.sub
               cpg_set.unique = cpg_set[rownames(cpg_set) %in% unique(cpg_set$CpG), ]
               
               ann450k = annotations
               
               anno.sub = ann450k[rownames(ann450k) %in% cpg_set.unique[cpg_set.unique$nodeLabel == nodeName, ]$CpG, ]
               anno.sub = data.frame(anno.sub)
               anno.sub = anno.sub[order(rownames(anno.sub)), ]
               
               cpg_set.unique = cpg_set.unique[order(cpg_set.unique$CpG), ]
               
               temp = cpg_set.unique[cpg_set.unique$nodeLabel == nodeName, ]
               all(temp$CpG == rownames(anno.sub))
               
               anno.sub = cbind('NodeCor' = temp$correlations, anno.sub)
               
               nodeEnhancer = nrow(anno.sub[anno.sub$Enhancer == 'TRUE', ]); nodeEnhancer
               nodeNoEnhancer = nrow(anno.sub[anno.sub$Enhancer == '', ]); nodeNoEnhancer
               annoEnhancer = nrow(ann450k[ann450k$Enhancer == 'TRUE', ]); annoEnhancer
               annoNoEnhancer = nrow(ann450k[ann450k$Enhancer == '', ]); annoNoEnhancer
               
               enhancers <- matrix(c(nodeEnhancer, annoEnhancer, 
                                     nodeNoEnhancer, annoNoEnhancer), nrow = 2,
                                   dimnames =
                                        list(c("NodeRelated", "NotNodeRelated"),
                                             c("Enhancer", "NotEnhancer")))
               
               node.results = fisher.test(enhancers); node.results
               
               
               
               row = c(nodeName, node.results$estimate, 
                       node.results$conf.int[1],
                       node.results$conf.int[2],
                       node.results$p.value)
               
               enhancerResults = rbind(enhancerResults, row)
          }
          
          enhancerResults = data.frame(enhancerResults)
          enhancerResults2 = enhancerResults
          colnames(enhancerResults2) = as.character(unlist(enhancerResults2[1, ]))
          enhancerResults2 = enhancerResults2[2:nrow(enhancerResults2), ]
          enhancerResults2$Est = as.numeric(as.character(enhancerResults2$Est))
          enhancerResults2$Conf95low = as.numeric(as.character(enhancerResults2$Conf95low))
          enhancerResults2$Conf95high = as.numeric(as.character(enhancerResults2$Conf95high))
          enhancerResults2$Pvalue = as.numeric(as.character(enhancerResults2$Pvalue))
          
          return(enhancerResults2)
     }
    
     

     correlations_list = list(correlations1, correlations8, correlations13, correlations16, 
                              correlations33, correlations40, correlations48,
                              correlations49, correlations55, correlations56, correlations58,
                              correlations66, correlations68, correlations86, correlations93)
     
     node_list = c(1, 8, 13, 16, 33, 40, 48, 49, 55, 56, 58, 66, 68, 86, 93) 
     
     enhancer_results = enhancer_analysis(correlations_list, node_list, threshold, ann450k)
     enhancer_results$Node = factor(enhancer_results$Node, 
                                    levels = enhancer_results$Node[order(as.numeric(substr(enhancer_results$Node, 4, 100)))])
     View(enhancer_results)
     
     
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
     library(ggplot2)
     fp <- ggplot(data=enhancer_results, aes(x=Est, 
                                            y=Node, 
                                            xmin=Conf95low, 
                                            xmax=Conf95high)) +
          geom_point(color = 'black') +
          geom_text(aes(label=format(round(Est, 2), nsmall = 2)), hjust=0, vjust=-1) +
          geom_errorbarh(height=.02) +
          ylab('Node') +
          geom_vline(xintercept=1, color='black', linetype='dashed') +
          scale_x_continuous(limits=c(0, max(enhancer_results$Conf95high)), name='Odds ratio w/ 95% CI') +
          ggtitle('OR for enhancer in each node set of CpGs') + theme_Publication()
     
     png('OR_enhancer_by_node.png', width = 2000, height = 2500, res = 300)
     fp
     dev.off()
     
     
#####################
# Visualize genomic context
#####################
     #biocLite('cirlize')
     library(circlize)
     
     ## Visualization
     #http://zuguang.de/circlize_book/book/high-level-genomic-functions.html#genomic-heatmap
     cpg_set.unique = cpg_set.unique = cpg_set[rownames(cpg_set) %in% unique(cpg_set$CpG), ]
     
     anno.sub = ann450k[which(rownames(ann450k) %in% cpg_set.unique$CpG), ]
     anno.sub = data.frame(anno.sub)
     
     anno.sub = anno.sub[order(rownames(anno.sub)), ]
     cpg_set.unique = cpg_set.unique[order(rownames(cpg_set.unique)), ]
     all(rownames(anno.sub) == rownames(cpg_set.unique))
     
     bed = cbind(anno.sub$chr, 
                 anno.sub$pos, 
                 anno.sub$pos + 10,
                 cpg_set.unique$correlations,
                 rownames(anno.sub),
                 rownames(cpg_set.unique),
                 as.character(cpg_set.unique$nodeLabel))
     bed = data.frame(bed)
     colnames(bed) = c('chr', 'start', 'end', 'value1', 'AnnoCpG', 'CorrCpG', 'Node')
     all(bed$AnnoCpG == bed$CorrCpG)
     
     bed$start = as.numeric(as.character(bed$start))
     bed$end = as.numeric(as.character(bed$end))
     bed$value1 = as.numeric(as.character(bed$value1))
     
     bed = bed[with(bed, order(bed$chr)), ]
     
     bed_list = list('VAE1' = bed[bed$Node == 'VAE1', ],
                     'VAE13' = bed[bed$Node == 'VAE13', ],
                     'VAE40' = bed[bed$Node == 'VAE40', ])#,
                     #'VAE42' = bed[bed$Node == 'VAE42', ],
                     #'VAE93' = bed[bed$Node == 'VAE93', ])
     
     circlize_plot = function() {
          circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr", 1:22))
          circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("#FF000080", 
                                                                        "blue",
                                                                        "#0000FF80", 
                                                                        "orangered1",
                                                                        "lightsteelblue"))
          
          circos.genomicDensity(bed[bed$Node == 'VAE1', ], col = c("#FF000080"), track.height = 0.1)
          circos.genomicDensity(bed[bed$Node == 'VAE13', ], col = c("blue"), track.height = 0.1)
          circos.genomicDensity(bed[bed$Node == 'VAE40', ], col = c("#0000FF80"), track.height = 0.1)
          #circos.genomicDensity(bed[bed$Node == 'VAE42', ], col = c("orangered1"), track.height = 0.1)
          #circos.genomicDensity(bed[bed$Node == 'VAE93', ], col = c("lightsteelblue"), track.height = 0.1)
          circos.clear()
     }
     
     
     ## Corr > 0.6
     png('results/genomic_context_corr.png', width = 2000, height = 2000, res = 300)
     circlize_plot()
     dev.off()
     
     out.data = merge(anno.sub, bed, by.x = 'Name', by.y = 'AnnoCpG')
     write.csv(out.data, file = 'results/genomic_context_anno_corr.csv')
          
     
#####################
# VAE activations
#####################
     library(reshape)
     vae.sub = vae[, c(1, 8, 13, 16, 33, 40, 48, 49, 55, 56, 58, 66, 68, 86, 93)  ]
     vae.sub$SubjectID = rownames(vae.sub)
     VAE.covs = merge(vae.sub, BRCA.covs, by.x = 'SubjectID', by.y = 'Basename')
     temp = melt(VAE.covs[, c(1, 8, 13, 16, 33, 40, 48, 49, 55, 56, 58, 66, 68, 86, 93, 
                              'sample.typeInt', "sample.type", "PAM50.RNAseq")], 
                 id=c('sample.typeInt', 'sample.type', 'PAM50.RNAseq'))
     
     colnames(temp) = c('SampleTypeInt', 'SampleType', 'PAM50orNormal', 'Node', 'Activation')
     temp$SampleTypeInt = factor(temp$SampleTypeInt)
     
     #install.packages('ggpubr')
     library(ggpubr)
     png('Node_activationsInt.png', width = 2000, height = 2000, res = 300)
     
     ggplot(temp, aes(x = reorder(Node, Activation, FUN = median), y = Activation, fill = SampleTypeInt)) + 
          geom_boxplot() + facet_grid(SampleTypeInt ~ .) + xlab("Node") + 
          ylab("Latent node activation value")
     
     dev.off()
     
     temp = temp[temp$PAM50 != '', ]
     png('Node_activationsPAM50.png', width = 3000, height = 2000, res = 300)
     ggplot(temp, aes(x = reorder(Node, Activation, FUN = median), y = Activation, fill = PAM50orNormal)) + 
          geom_boxplot() + facet_grid(PAM50orNormal ~ .) + xlab("Node") + 
          ylab("Latent node activation value")
     
     dev.off()
   
     
     
#####################
# multinomial regression
#####################
     
     library("nnet")
     
     set.seed(100)
     temp = VAE.covs[VAE.covs$PAM50 != '', ]
     trainingRows <- sample(1:nrow(temp), 0.7*nrow(temp))
     training <- temp[trainingRows, ]
     test <- temp[-trainingRows, ]
     
     train <- multinom(PAM50.RNAseq ~ `1` + `13` + `16` +  
                                      `40` + `49` + `55` + 
                                      `56` + `58` + `66` + 
                                      `68` + `86` + `93`, data = training)
     summary(train)
     z <- summary(train)$coefficients/summary(train)$standard.errors
     z
     p <- (1 - pnorm(abs(z), 0, 1))*2
     p
     exp(coef(train))
     
     pred <- predict (train, test, "probs") # predict on new data
     pred
     pred_class <- predict(train, test)
     pred_class
     
     table(pred_class, test$PAM50.RNAseq)
     mean(as.character(pred_class) == as.character(test$PAM50.RNAseq))
     
     
     