mainDir = '/Users/linwang/Work/dingcheng lineage tracing/amplican crispr analysis/result all cell 1'
subDir <- "network"
subDir1 <- 'd5'
dir.create(file.path(mainDir, subDir, subDir1), showWarnings = FALSE)
setwd(file.path(mainDir, subDir, subDir1))

load('network.direction.d5A21.RData')
load('network.direction.d5B25.RData')

mainDir = '/Users/linwang/Work/dingcheng lineage tracing/amplican crispr analysis/result all cell 1'
subDir <- "network"
subDir1 <- 'd14'
dir.create(file.path(mainDir, subDir, subDir1), showWarnings = FALSE)
setwd(file.path(mainDir, subDir, subDir1))

load('network.direction.d14A21.RData')
load('network.direction.d14B25.RData')


############################################################## combine network A21
node.d5A21
node.d14A21

node.cellline = full_join(node.d5A21, node.d14A21, by  =c('node','node1'))

node.cellline[is.na(node.cellline[,3]),3] = 0
node.cellline[is.na(node.cellline[,4]),4] = 0
colnames(node.cellline)[3:4] = c('cell.number.d5A21', 'cell.number.d4A21')

network.direction.d5A21$sample = 'd5A21'
network.direction.d14A21$sample = 'd14A21'

network.direction.cellline = full_join(network.direction.d5A21[,c(1,2,4)], network.direction.d14A21[,c(1,2,4)], by  =c('from','to'))
colnames(network.direction.cellline)[3:4] = c('sample.d5A21','sample.d14A21')

mainDir = '/Users/linwang/Work/dingcheng lineage tracing/amplican crispr analysis/result all cell 1'
subDir <- "network"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))

write.table(node.cellline, file = 'node.cellline.A21.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(network.direction.cellline, file = 'network.direction.cellline.A21.txt', quote = FALSE, sep = '\t', row.names = FALSE)

############################################################## combine network B25
node.d5B25
node.d14B25

node.cellline = full_join(node.d5B25, node.d14B25, by  =c('node','node1'))

node.cellline[is.na(node.cellline[,3]),3] = 0
node.cellline[is.na(node.cellline[,4]),4] = 0
colnames(node.cellline)[3:4] = c('cell.number.d5B25', 'cell.number.d4B25')

network.direction.d5B25$sample = 'd5B25'
network.direction.d14B25$sample = 'd14B25'

network.direction.cellline = full_join(network.direction.d5B25[,c(1,2,4)], network.direction.d14B25[,c(1,2,4)], by  =c('from','to'))
colnames(network.direction.cellline)[3:4] = c('sample.d5B25','sample.d14B25')

mainDir = '/Users/linwang/Work/dingcheng lineage tracing/amplican crispr analysis/result all cell 1'
subDir <- "network"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))

write.table(node.cellline, file = 'node.cellline.B25.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(network.direction.cellline, file = 'network.direction.cellline.B25.txt', quote = FALSE, sep = '\t', row.names = FALSE)
############################################################## combine network
node.d5
node.d14

node.cellline = full_join(node.d5, node.d14, by  =c('node','node1'))

node.cellline[is.na(node.cellline[,3]),3] = 0
node.cellline[is.na(node.cellline[,4]),4] = 0
colnames(node.cellline)[3:4] = c('cell.number.d5', 'cell.number.d4')

network.direction.d5$sample = 'd5'
network.direction.d14$sample = 'd14'

network.direction.cellline = full_join(network.direction.d5[,c(1,2,4)], network.direction.d14[,c(1,2,4)], by  =c('from','to'))
colnames(network.direction.cellline)[3:4] = c('sample.d5','sample.d14')

mainDir = '/Users/linwang/Work/dingcheng lineage tracing/amplican crispr analysis/result all cell 1'
subDir <- "network"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))

write.table(node.cellline, file = 'node.cellline.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(network.direction.cellline, file = 'network.direction.cellline.txt', quote = FALSE, sep = '\t', row.names = FALSE)

################################# hierarchical cluster
# unique.aligh = unique(data$align.short)
# 
# smilarity.matrix = matrix(,nrow = length(unique.aligh), ncol = length(unique.aligh))
# 
# for (i in 1 : length(unique.aligh)){
#   smilarity.matrix[i,] = levenshteinSim(unique.aligh[i],unique.aligh)
# }
# 
# as.dist(smilarity.matrix)
# tree = hclust(as.dist(smilarity.matrix))
# plot(tree)
# cutree(tree, k = 6)
# 
# unique.aligh[which(cutree(tree, k = 8) == 1)]
# 
# tmp =  data %>% count(align.short)
# 
# tmp[order(tmp[,2],decreasing = TRUE),]
# 
# sum(tmp[,2])
# 
# tmp1 = data[which(data$align.short == 'GTTCCCGTCCAGTAATCGTGGGG'),]
# 
# nrow(tmp1)
# 
# data[order(data$num.reads, decreasing = TRUE),]


#################################################hierarchical clustering
