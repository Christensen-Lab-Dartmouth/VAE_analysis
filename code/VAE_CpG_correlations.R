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
     setwd(dir)
     BRCA.covs = data.frame(fread(BRCA.covFile), row.names=1)
     BRCA.covs$sample.typeInt = ifelse(BRCA.covs$sample.type == 'Solid Tissue Normal', 0, 1)
     
     dir = 'C:/Users/atitus/Documents/github/DNAm_data_generation'
     setwd(dir)
     
     ## Betas
     beta.file = 'data/TCGA_BRCA_top300kMAD_cpg.tsv'
     betas = data.frame(fread(beta.file))
     rownames(betas) = betas[,1]
     betas = betas[,2:ncol(betas)]
     betas = betas[order(rownames(betas), decreasing=T), ]
     
     BRCA.covs = BRCA.covs[order(BRCA.covs$Basename, decreasing=T), ]
     
     ## VAE nodes
     vae.file = 'results/encoded_methyl_onehidden_warmup_batchnorm_300K-100.tsv'
     vae = data.frame(fread(vae.file))
     colnames(vae) = vae[1,]
     rownames(vae) = vae[,1]
     vae = vae[2:nrow(vae), 2:ncol(vae)]
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

     correlations = apply(betas, 2, cor.func)
     correlations = data.frame(correlations)
     correlations = cbind('CpG' = rownames(correlations), correlations)
     
     nodeLabel = paste('VAE', node, sep = '')
     correlations = cbind(nodeLabel = rep(nodeLabel, nrow(correlations)), correlations)
     
     correlations = correlations[order(abs(correlations$correlations), decreasing=T), ]
     
     
#####################
# Correlation diagram
##################### 
     ## Correlation elbow
     cors = correlations$correlations
     file.name = paste('correlation_elbow_node', node, '.png', sep = '')
     png(file.name, width = 1000, height = 1000, res = 100)
     plot(abs(cors), main=paste('Node ', node, sep = ''), ylab = '|Correlation|', xlab='Index')
     dev.off()
     
     ## Correlation with methylation
     cors = correlations$correlations
     cpg = correlations$CpG[1]
     cg14789818 = betas[,cpg]
     vaeNode = vae[, as.character(node)]

     df = data.frame(cbind(vaeNode, cg14789818))
     
     file.name = paste('methylation_V_node', node, '.png', sep = '')
     png(file.name, width = 1000, height = 1000, res = 100)
     
     ggplot(data = df, aes(x = vaeNode, y = cg14789818)) + 
          geom_point(fill = '#C7BBC9', pch = 21, size = 5) + 
          geom_smooth(method ='lm', formula = y~x) + 
          xlab('VAE node activation') + 
          ylab('cg14789818 methylation value')  +
          ggtitle(paste('Node ', node, sep = ''))
     dev.off()
     
     
#####################
# Chord diagram 
#####################
     library(circlize)
     summary(correlations)
     
     threshold = 0.8
     
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
     
     cpgs = rownames(cor.sub)
     
     ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
     head(ann450k)
     
     anno.sub = ann450k[rownames(ann450k) %in% cpgs, ]
     anno.sub = data.frame(anno.sub)
     anno.sub = cbind('NodeCor' = cor.sub$correlations, anno.sub)
     
     file.name = paste('results/anno450K_node', node, '.csv', sep = '')
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
     go = topGO(gst, number = 10); go

     # Top 10 gene sets
     gsa = topGSA(gsa, number=10); gsa
     
     gofile.name = paste('results/go_anno_node', node, '.csv', sep = '')
     write.csv(go, gofile.name) 
     
     gsafile.name = paste('results/gsa_anno_node', node, '.csv', sep = '')
     write.csv(gsa, gsafile.name) 
     
     
#####################
# OmicCircos
#####################
     library(OmicCircos) 
     head(anno.sub)
     