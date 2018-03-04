
node1 = read.csv('/Users/alexandertitus/Documents/github/DNAm_data_generation/results/anno450K_node1.csv')
table(node1$Relation_to_Island)
table(node1$Relation_to_Island, node1$Enhancer)

node13 = read.csv('/Users/alexandertitus/Documents/github/DNAm_data_generation/results/anno450K_node13.csv')
table(node13$Relation_to_Island)
table(node13$Relation_to_Island, node13$Enhancer)

node16 = read.csv('/Users/alexandertitus/Documents/github/DNAm_data_generation/results/anno450K_node16.csv')
table(node16$Relation_to_Island)
table(node16$Relation_to_Island, node16$Enhancer)

node40 = read.csv('/Users/alexandertitus/Documents/github/DNAm_data_generation/results/anno450K_node40.csv')
table(node40$Relation_to_Island)
table(node40$Relation_to_Island, node40$Enhancer)

node42 = read.csv('/Users/alexandertitus/Documents/github/DNAm_data_generation/results/anno450K_node42.csv')
table(node42$Relation_to_Island)
table(node42$Relation_to_Island, node42$Enhancer)

node56 = read.csv('/Users/alexandertitus/Documents/github/DNAm_data_generation/results/anno450K_node56.csv')
table(node56$Relation_to_Island)
table(node56$Relation_to_Island, node56$Enhancer)

node58 = read.csv('/Users/alexandertitus/Documents/github/DNAm_data_generation/results/anno450K_node58.csv')
table(node58$Relation_to_Island)
table(node58$Relation_to_Island, node58$Enhancer)

node86 = read.csv('/Users/alexandertitus/Documents/github/DNAm_data_generation/results/anno450K_node86.csv')
table(node86$Relation_to_Island)
table(node86$Relation_to_Island, node86$Enhancer)

node93 = read.csv('/Users/alexandertitus/Documents/github/DNAm_data_generation/results/anno450K_node93.csv')
table(node93$Relation_to_Island)
table(node93$Relation_to_Island, node93$Enhancer)

setwd('/Users/alexandertitus/Documents/github/DNAm_data_generation/results/')
png('nodeCpG_correllation_relationToisland.png', width = 800, height = 1200, res = 100)
par(mfrow=c(3,3))
barplot(table(node1$Relation_to_Island), las=2)
barplot(table(node13$Relation_to_Island), las=2)
barplot(table(node16$Relation_to_Island), las=2)
barplot(table(node40$Relation_to_Island), las=2)
barplot(table(node42$Relation_to_Island), las=2)
barplot(table(node56$Relation_to_Island), las=2)
barplot(table(node58$Relation_to_Island), las=2)
barplot(table(node86$Relation_to_Island), las=2)
barplot(table(node93$Relation_to_Island), las=2)
dev.off()


png('nodeCpG_correllation_chr.png', width = 1200, height = 1200, res = 100)
par(mfrow=c(3,3))
names = as.numeric(substr(names(table(node1$chr)), 4, 100))
barplot(table(node1$chr)[order(names)], las=2)

names = as.numeric(substr(names(table(node13$chr)), 4, 100))
barplot(table(node13$chr)[order(names)], las=2)

names = as.numeric(substr(names(table(node16$chr)), 4, 100))
barplot(table(node16$chr)[order(names)], las=2)

names = as.numeric(substr(names(table(node40$chr)), 4, 100)); names
barplot(table(node40$chr)[order(names)], las=2)

names = as.numeric(substr(names(table(node42$chr)), 4, 100)); names
barplot(table(node42$chr)[order(names)], las=2)

names = as.numeric(substr(names(table(node56$chr)), 4, 100)); names
barplot(table(node56$chr)[order(names)], las=2)

names = as.numeric(substr(names(table(node58$chr)), 4, 100)); names
barplot(table(node58$chr)[order(names)], las=2)

names = as.numeric(substr(names(table(node86$chr)), 4, 100)); names
barplot(table(node86$chr)[order(names)], las=2)

names = as.numeric(substr(names(table(node93$chr)), 4, 100)); names
barplot(table(node93$chr)[order(names)], las=2)
dev.off()