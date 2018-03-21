base = 'C:/Users/atitus/Documents/github/DNAm_data_generation/results/'

node1 = read.csv(paste(base, 'anno450K_node1.csv', sep = ''))
node1$node = 1
table(node1$Relation_to_Island)
table(node1$Relation_to_Island, node1$Enhancer)

node8 = read.csv(paste(base, 'anno450K_node8.csv', sep = ''))
node8$node = 8
table(node8$Relation_to_Island)
table(node8$Relation_to_Island, node8$Enhancer)

node13 = read.csv(paste(base, 'anno450K_node13.csv', sep = ''))
node13$node = 13
table(node13$Relation_to_Island)
table(node13$Relation_to_Island, node13$Enhancer)

node16 = read.csv(paste(base, 'anno450K_node16.csv', sep = ''))
node16$node = 16
table(node16$Relation_to_Island)
table(node16$Relation_to_Island, node16$Enhancer)

node33 = read.csv(paste(base, 'anno450K_node33.csv', sep = ''))
node33$node = 33
table(node33$Relation_to_Island)
table(node33$Relation_to_Island, node33$Enhancer)

node40 = read.csv(paste(base, 'anno450K_node40.csv', sep = ''))
node40$node = 40
table(node40$Relation_to_Island)
table(node40$Relation_to_Island, node40$Enhancer)

node48 = read.csv(paste(base, 'anno450K_node48.csv', sep = ''))
node48$node = 48
table(node48$Relation_to_Island)
table(node48$Relation_to_Island, node48$Enhancer)

node49 = read.csv(paste(base, 'anno450K_node49.csv', sep = ''))
node49$node = 49
table(node49$Relation_to_Island)
table(node49$Relation_to_Island, node49$Enhancer)

node55 = read.csv(paste(base, 'anno450K_node55.csv', sep = ''))
node55$node = 55
table(node55$Relation_to_Island)
table(node55$Relation_to_Island, node55$Enhancer)

node56 = read.csv(paste(base, 'anno450K_node56.csv', sep = ''))
node56$node = 56
table(node56$Relation_to_Island)
table(node56$Relation_to_Island, node56$Enhancer)

node58 = read.csv(paste(base, 'anno450K_node58.csv', sep = ''))
node58$node = 58
table(node58$Relation_to_Island)
table(node58$Relation_to_Island, node58$Enhancer)

node66 = read.csv(paste(base, 'anno450K_node66.csv', sep = ''))
node66$node = 66
table(node66$Relation_to_Island)
table(node66$Relation_to_Island, node66$Enhancer)

node68 = read.csv(paste(base, 'anno450K_node68.csv', sep = ''))
node68$node = 68
table(node68$Relation_to_Island)
table(node68$Relation_to_Island, node68$Enhancer)

node86 = read.csv(paste(base, 'anno450K_node86.csv', sep = ''))
node86$node = 86
table(node86$Relation_to_Island)
table(node86$Relation_to_Island, node86$Enhancer)

node93 = read.csv(paste(base, 'anno450K_node93.csv', sep = ''))
node93$node = 93
table(node93$Relation_to_Island)
table(node93$Relation_to_Island, node93$Enhancer)


full.data = rbind(node1, node8, node13, node16, node33,
                  node40, node48, node49, node55, node56,
                  node58, node66, node68, node86, node93)

setwd('/Users/alexandertitus/Documents/github/DNAm_data_generation/results/')
png('nodeCpG_correllation_relationToisland.png', width = 1600, height = 2400, res = 100)
par(mfrow=c(3,5))
barplot(table(node1$Relation_to_Island), las=2, main = 'Node 1')
barplot(table(node8$Relation_to_Island), las=2, main = 'Node 8')
barplot(table(node13$Relation_to_Island), las=2, main = 'Node 13')
barplot(table(node16$Relation_to_Island), las=2, main = 'Node 16')
barplot(table(node33$Relation_to_Island), las=2, main = 'Node 33')
barplot(table(node40$Relation_to_Island), las=2, main = 'Node 40')
barplot(table(node48$Relation_to_Island), las=2, main = 'Node 48')
barplot(table(node49$Relation_to_Island), las=2, main = 'Node 49')
barplot(table(node55$Relation_to_Island), las=2, main = 'Node 55')
barplot(table(node56$Relation_to_Island), las=2, main = 'Node 56')
barplot(table(node58$Relation_to_Island), las=2, main = 'Node 58')
barplot(table(node66$Relation_to_Island), las=2, main = 'Node 66')
barplot(table(node68$Relation_to_Island), las=2, main = 'Node 68')
barplot(table(node86$Relation_to_Island), las=2, main = 'Node 86')
barplot(table(node93$Relation_to_Island), las=2, main = 'Node 93')
dev.off()

## Plot CpGs by island context
n1 = ggplot(node1, aes(x = Relation_to_Island)) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 1") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per CpG island context") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

n8 = ggplot(node8, aes(x = Relation_to_Island)) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 8") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per CpG island context") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

n13 = ggplot(node13, aes(x = Relation_to_Island)) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 13") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per CpG island context") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

n16 = ggplot(node16, aes(x = Relation_to_Island)) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 16") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per CpG island context") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

n33 = ggplot(node33, aes(x = Relation_to_Island)) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 33") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per CpG island context") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

n40 = ggplot(node40, aes(x = Relation_to_Island)) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 40") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per CpG island context") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

n48 = ggplot(node48, aes(x = Relation_to_Island)) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 48") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per CpG island context") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

n49 = ggplot(node49, aes(x = Relation_to_Island)) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 49") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per CpG island context") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

n55 = ggplot(node55, aes(x = Relation_to_Island)) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 55") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per CpG island context") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

n56 = ggplot(node56, aes(x = Relation_to_Island)) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 56") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per CpG island context") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

n58 = ggplot(node58, aes(x = Relation_to_Island)) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 58") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per CpG island context") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

n66 = ggplot(node66, aes(x = Relation_to_Island)) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 66") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per CpG island context") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

n68 = ggplot(node68, aes(x = Relation_to_Island)) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 68") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per CpG island context") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

n86 = ggplot(node86, aes(x = Relation_to_Island)) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 86") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per CpG island context") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

n93 = ggplot(node93, aes(x = Relation_to_Island)) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 93") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per CpG island context") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

g = grid.arrange(n1, n8, n13, n16, n33, 
                 n40, n48, n49, n55, n56,
                 n58, n66, n68, n86, n93, ncol=4)

ggsave(g, file="nodeCpG_correllation_relationToisland.png", width=16, height=16)


## Plot CpGs by chromosome
library(scales)
library(gridExtra)

names1 = as.numeric(substr(node1$chr, 4, 100))
n1 = ggplot(node1, aes(x = reorder(chr, names1))) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 1") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per chromosome") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

names8 = as.numeric(substr(node8$chr, 4, 100))
n8 = ggplot(node8, aes(x = reorder(chr, names8))) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 8") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per chromosome") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

names13 = as.numeric(substr(node13$chr, 4, 100))
n13 = ggplot(node13, aes(x = reorder(chr, names13))) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 13") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per chromosome") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

names16 = as.numeric(substr(node16$chr, 4, 100))
n16 = ggplot(node16, aes(x = reorder(chr, names16))) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 16") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per chromosome") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

names33 = as.numeric(substr(node33$chr, 4, 100))
n33 = ggplot(node33, aes(x = reorder(chr, names33))) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 33") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per chromosome") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

names40 = as.numeric(substr(node40$chr, 4, 100))
n40 = ggplot(node40, aes(x = reorder(chr, names40))) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 40") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per chromosome") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

names48 = as.numeric(substr(node48$chr, 4, 100))
n48 = ggplot(node48, aes(x = reorder(chr, names48))) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 48") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per chromosome") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

names49 = as.numeric(substr(node49$chr, 4, 100))
n49 = ggplot(node49, aes(x = reorder(chr, names49))) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 49") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per chromosome") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

names55 = as.numeric(substr(node55$chr, 4, 100))
n55 = ggplot(node55, aes(x = reorder(chr, names55))) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 55") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per chromosome") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

names56 = as.numeric(substr(node56$chr, 4, 100))
n56 = ggplot(node56, aes(x = reorder(chr, names56))) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 56") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per chromosome") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

names58 = as.numeric(substr(node58$chr, 4, 100))
n58 = ggplot(node58, aes(x = reorder(chr, names58))) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 58") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per chromosome") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

names66 = as.numeric(substr(node66$chr, 4, 100))
n66 = ggplot(node66, aes(x = reorder(chr, names66))) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 66") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per chromosome") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

names68 = as.numeric(substr(node68$chr, 4, 100))
n68 = ggplot(node68, aes(x = reorder(chr, names68))) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 68") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per chromosome") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

names86 = as.numeric(substr(node86$chr, 4, 100))
n86 = ggplot(node86, aes(x = reorder(chr, names86))) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 86") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per chromosome") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

names93 = as.numeric(substr(node93$chr, 4, 100))
n93 = ggplot(node93, aes(x = reorder(chr, names93))) + geom_bar(aes(y=..count../sum(..count..))) +
     scale_y_continuous(labels=percent_format()) + ggtitle("Node 93") +
     xlab("Chromosome") + ylab("Proportion of associated CpGs per chromosome") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

g = grid.arrange(n1, n8, n13, n16, n33, 
             n40, n48, n49, n55, n56,
             n58, n66, n68, n86, n93, ncol=4)

ggsave(g, file="nodeCpG_correllation_chr.png", width=16, height=16)
