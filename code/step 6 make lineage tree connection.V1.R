setwd('/Users/linwang/Work/dingcheng lineage tracing/amplican crispr analysis/result all cell 1')
setwd('./cell.read')

setwd('./d14')

library(dplyr)
library(stringr)
library(RecordLinkage)

load('cell.read.d14.RData')

data = cell.read.d14A21
data = cell.read.d14B25

###############################################################
data = data[data$number.of.notinread ==0,]
# data = data[which(data$number.of.different <= 23),]

data = data[!is.na(data$align.short),]

table(data$number.of.different)
table(data$number.of.mutation)

table(data$number.of.insertion)
table(data$number.of.deletion)

table(data$length.of.insertion)
table(data$length.of.deletion)

table(data$number.of.mismatch)

nrow(unique(data[,c('mutation.correction1','align.correction1','mutation.report')]))

data <-
  data %>% 
  group_by(mutation.correction1,align.correction1,mutation.report) %>% 
  dplyr::summarise(n = n())
data = as.data.frame(data)

nrow(data)
# 198 unique cell.read.d 14B25
# 450 unique cell.read.d 14A21
# 239 unique cell.read.d 5A21

##########################
network.direction.d14A21 = make.network.direction(data)

node.d14A21 = cbind(data[,c('mutation.correction1')], data[,c('mutation.report')], data[,'n'])

node.d14A21 = as.data.frame(node.d14A21)

colnames(node.d14A21) = c('node','node1','cell.number')

network.direction.d14B25 = make.network.direction.B25(data)

node.d14B25 = cbind(data[,c('mutation.correction1')], data[,c('mutation.report')], data[,'n'])

node.d14B25 = as.data.frame(node.d14B25)

colnames(node.d14B25) = c('node','node1','cell.number')

mainDir = '/Users/linwang/Work/dingcheng lineage tracing/amplican crispr analysis/result all cell 1'
subDir <- "network"
subDir1 <- 'd14'
dir.create(file.path(mainDir, subDir, subDir1), showWarnings = FALSE)
setwd(file.path(mainDir, subDir, subDir1))

save(network.direction.d14A21, node.d14A21, file = 'network.direction.d14A21.RData')
save(network.direction.d14B25, node.d14B25, file = 'network.direction.d14B25.RData')

write.table(network.direction.d14A21, file = 'network.direction.d14A21.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(node.d14A21, file = 'node.d14A21.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(network.direction.d14B25, file = 'network.direction.d14B25.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(node.d14B25, file = 'node.d14B25.txt', quote = FALSE, sep = '\t', row.names = FALSE)

###############################################################
data = d14A21
data1 = d14B25
data = data[!duplicated(data$ID),]
data1 = data1[!duplicated(data1$ID),]

data = full_join(data[,c('ID','mutation.correction1','align.correction1','mutation.report')], 
                 data1[,c('ID','mutation.correction1','align.correction1','mutation.report')], by  ='ID')


data <-
  data %>% 
  group_by(mutation.correction1.x,align.correction1.x,mutation.report.x,mutation.correction1.y,align.correction1.y,mutation.report.y) %>% 
  dplyr::summarise(n = n())
data = as.data.frame(data)

network.direction.d14 = make.network.direction.combine(data,2)

node.d14 = cbind(paste(data[,c('mutation.correction1.x')], data[,c('mutation.correction1.y')], sep = ','), 
                 paste(data[,c('mutation.report.x')], data[,c('mutation.report.y')], sep = ' && '), data[,'n'])

node.d14 = as.data.frame(node.d14)
colnames(node.d14) = c('node','node1','cell.number')

mainDir = '/Users/linwang/Work/dingcheng lineage tracing/amplican crispr analysis/result all cell 1'
subDir <- "network"
subDir1 <- 'd14'
dir.create(file.path(mainDir, subDir, subDir1), showWarnings = FALSE)
setwd(file.path(mainDir, subDir, subDir1))

save(network.direction.d14, node.d14, file = 'network.direction.d14.RData')

write.table(node.d14, file = 'node.d14.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(network.direction.d14, file = 'network.direction.d14.txt', quote = FALSE, sep = '\t', row.names = FALSE)

############################################################# function
################################### make.mutation.type.matrix
make.mutation.type.matrix <- function(mutation.report.matrix.list){
	mutation.type.matrix = matrix(,ncol = 47,nrow = length(mutation.report.matrix.list))
	for (i in 1 : nrow(mutation.type.matrix)){
		  mutation.type.matrix[i,] = mutation.report.matrix.list[[i]][,1]
	}
	tmp = mutation.type.matrix
	tmp[is.na(mutation.type.matrix)]  = 0
	tmp[mutation.type.matrix ==  'insert' ] = 1
	tmp[mutation.type.matrix ==  'mismatch'] =  0.5
	tmp[mutation.type.matrix ==  'miss'] = 1
	tmp = matrix(as.numeric(tmp), ncol = ncol(tmp))
	return(tmp)
	return(mutation.type.matrix)
}

make.mutation.type.matrix.B25 <- function(mutation.report.matrix.list){
  mutation.type.matrix = matrix(,ncol = 55,nrow = length(mutation.report.matrix.list))
  for (i in 1 : nrow(mutation.type.matrix)){
      mutation.type.matrix[i,] = mutation.report.matrix.list[[i]][,1]
  }
  tmp = mutation.type.matrix
  tmp[is.na(mutation.type.matrix)]  = 0
  tmp[mutation.type.matrix ==  'insert' ] = 1
  tmp[mutation.type.matrix ==  'mismatch'] =  0.5
  tmp[mutation.type.matrix ==  'miss'] = 1
  tmp = matrix(as.numeric(tmp), ncol = ncol(tmp))
  return(tmp)
}
################################### make.mutation.type.matrix.plus
make.mutation.type.matrix.plus <- function(mutation.report.matrix.list){
	mutation.type.matrix = matrix(,ncol = 47,nrow = length(mutation.report.matrix.list))
	for (i in 1 : nrow(mutation.type.matrix)){
		    mutation.type.matrix[i,] = mutation.report.matrix.list[[i]][,3]
	}
	mutation.type.matrix = matrix(as.numeric(mutation.type.matrix), ncol = ncol(mutation.type.matrix))
	return(mutation.type.matrix)
}

make.mutation.type.matrix.plus.B25 <- function(mutation.report.matrix.list){
  mutation.type.matrix = matrix(,ncol = 55,nrow = length(mutation.report.matrix.list))
  for (i in 1 : nrow(mutation.type.matrix)){
      mutation.type.matrix[i,] = mutation.report.matrix.list[[i]][,3]
  }
  mutation.type.matrix = matrix(as.numeric(mutation.type.matrix), ncol = ncol(mutation.type.matrix))
  return(mutation.type.matrix)
}




################################### make.network.direction

make.network.direction <- function(data){
  
  network.direction = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(network.direction) <- c("from", "to",'score')
  
  mutation.report.matrix = make.mutation.report.matrix(data[,'mutation.report'])
  mutation.type.matrix = make.mutation.type.matrix(mutation.report.matrix)
  mutation.type.matrix.plus = make.mutation.type.matrix.plus(mutation.report.matrix)
  
  for( i in 1 : nrow(data)){
    tmp.ind = c()
    for (j in setdiff(1:nrow(data),i)){
      if(min(mutation.type.matrix[i,] - mutation.type.matrix[j,]) >=0  ){
        if ( !identical(mutation.type.matrix[i,],  mutation.type.matrix[j,])){
          tmp.ind = c(tmp.ind,j)
        }
      }
    }
    
    if(is.null(tmp.ind)) {next}
    
    score1 = rowSums(mutation.type.matrix[tmp.ind,, drop = FALSE])
    # score2 = rowSums(mutation.type.matrix.plus[tmp.ind,, drop = FALSE])
    # 
    # tmp <-
    #   data.frame(tmp.ind, score1) %>% 
    #   arrange(desc(score1)) 
    tmp = tmp.ind[which(score1 == max(score1))]
    
    network.direction[nrow(network.direction) + 1,] = c(data[tmp[1],'mutation.report'],data[i,'mutation.report'],data[i,'n'])
    
  }
  
  return(network.direction)
}



make.network.direction.B25 <- function(data){
  
  network.direction = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(network.direction) <- c("from", "to",'score')
  
  mutation.report.matrix = make.mutation.report.matrix.B25(data[,'mutation.report'])
  mutation.type.matrix = make.mutation.type.matrix.B25(mutation.report.matrix)
  mutation.type.matrix.plus = make.mutation.type.matrix.plus.B25(mutation.report.matrix)
 
  for( i in 1 : nrow(data)){
    tmp.ind = c()
    for (j in setdiff(1:nrow(data),i)){
      if(min(mutation.type.matrix[i,] - mutation.type.matrix[j,]) >=0  ){
        if (!identical(mutation.type.matrix[i,],  mutation.type.matrix[j,]) ){
          tmp.ind = c(tmp.ind,j)
        }
      }
    }
    
    if(is.null(tmp.ind)) {next}
    
    score1 = rowSums(mutation.type.matrix[tmp.ind,, drop = FALSE])
    # score2 = rowSums(mutation.type.matrix.plus[tmp.ind,, drop = FALSE])
    # 
    #   tmp <-
    #   data.frame(tmp.ind, score1,score2) %>% 
    #   arrange(desc(score1),score2) 
    tmp = tmp.ind[which(score1 == max(score1))]
    
    network.direction[nrow(network.direction) + 1,] = c(data[tmp[1],'mutation.report'],data[i,'mutation.report'],data[i,'n'])

  }
  
  return(network.direction)
}
  
# make.network.direction.combine <- function(data,num){
#   
#   network.direction = data.frame(matrix(ncol = 3, nrow = 0))
#   colnames(network.direction) <- c("from", "to",'score')
#   
#   mutation.report.matrix.x = make.mutation.report.matrix(data[,'mutation.report.x'])
#   mutation.type.matrix.x = make.mutation.type.matrix(mutation.report.matrix.x)
#   mutation.type.matrix.plus.x = make.mutation.type.matrix.plus(mutation.report.matrix.x)
#   
#   mutation.report.matrix.y = make.mutation.report.matrix.B25(data[,'mutation.report.y'])
#   mutation.type.matrix.y = make.mutation.type.matrix.B25(mutation.report.matrix.y)
#   mutation.type.matrix.plus.y = make.mutation.type.matrix.plus.B25(mutation.report.matrix.y)
#   
#   for( i in 1 : nrow(data)){
#     tmp.ind = c()
#     for (j in setdiff(1:nrow(data),i)){
#       if(min(mutation.type.matrix.x[i,] - mutation.type.matrix.x[j,]) >=0  & min(mutation.type.matrix.plus.x[i,] - mutation.type.matrix.plus.x[j,]) >=0
#          & min(mutation.type.matrix.y[i,] - mutation.type.matrix.y[j,]) >=0  & min(mutation.type.matrix.plus.y[i,] - mutation.type.matrix.plus.y[j,]) >=0){
#         if (!identical(mutation.type.matrix.x[i,],  mutation.type.matrix.x[j,]) | !identical(mutation.type.matrix.plus.x[i,] , mutation.type.matrix.plus.x[j,])
#             | !identical(mutation.type.matrix.y[i,],  mutation.type.matrix.y[j,]) | !identical(mutation.type.matrix.plus.y[i,] , mutation.type.matrix.plus.y[j,])){
#           tmp.ind = c(tmp.ind,j)
#         }
#       }
#     }
#     if(is.null(tmp.ind)) {next}
#     
#     score1 = rowSums(mutation.type.matrix.x[tmp.ind,, drop = FALSE]) + rowSums(mutation.type.matrix.y[tmp.ind,, drop = FALSE])
#     score2 = rowSums(mutation.type.matrix.plus.x[tmp.ind,, drop = FALSE]) + rowSums(mutation.type.matrix.plus.y[tmp.ind,, drop = FALSE])
#     
#     tmp <-
#       data.frame(tmp.ind, score1,score2) %>% 
#       arrange(desc(score1),score2) 
#     
#     network.direction[nrow(network.direction) + 1,] = c(paste(data[tmp[1,1],'mutation.correction1.x'], data[tmp[1,1],'mutation.correction1.y'], sep = ','),
#                                                         paste(data[i,'mutation.correction1.x'], data[i,'mutation.correction1.y'], sep = ','),data[i,'n'])
#     
#   }
#   return(network.direction)
# }


make.network.direction.combine <- function(data,num){
  
  network.direction = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(network.direction) <- c("from", "to",'score')
  
  mutation.report.matrix.x = make.mutation.report.matrix(data[,'mutation.report.x'])
  mutation.type.matrix.x = make.mutation.type.matrix(mutation.report.matrix.x)
  mutation.type.matrix.plus.x = make.mutation.type.matrix.plus(mutation.report.matrix.x)
  
  mutation.report.matrix.y = make.mutation.report.matrix.B25(data[,'mutation.report.y'])
  mutation.type.matrix.y = make.mutation.type.matrix.B25(mutation.report.matrix.y)
  mutation.type.matrix.plus.y = make.mutation.type.matrix.plus.B25(mutation.report.matrix.y)
  
  for( i in 1 : nrow(data)){
    tmp.ind = c()
    for (j in setdiff(1:nrow(data),i)){
      if(min(mutation.type.matrix.x[i,] - mutation.type.matrix.x[j,]) >=0  
         & min(mutation.type.matrix.y[i,] - mutation.type.matrix.y[j,]) >=0  ){
        if (!identical(mutation.type.matrix.x[i,],  mutation.type.matrix.x[j,]) 
            | !identical(mutation.type.matrix.y[i,],  mutation.type.matrix.y[j,]) ){
          tmp.ind = c(tmp.ind,j)
        }
      }
    }
    if(is.null(tmp.ind)) {next}
    
    score1 = rowSums(mutation.type.matrix.x[tmp.ind,, drop = FALSE]) + rowSums(mutation.type.matrix.y[tmp.ind,, drop = FALSE])
    # score2 = rowSums(mutation.type.matrix.plus.x[tmp.ind,, drop = FALSE]) + rowSums(mutation.type.matrix.plus.y[tmp.ind,, drop = FALSE])
    
    # tmp <-
    #   data.frame(tmp.ind, score1,score2) %>% 
    #   arrange(desc(score1),score2) 
    
    tmp = tmp.ind[which(score1 == max(score1))]
    
    # network.direction[nrow(network.direction) + 1,] = c(paste(data[tmp[1],'mutation.correction1.x'], data[tmp[1],'mutation.correction1.y'], sep = ','),
    #                                                     paste(data[i,'mutation.correction1.x'], data[i,'mutation.correction1.y'], sep = ','),data[i,'n'])
    network.direction[nrow(network.direction) + 1,] = c(paste(data[tmp[1],'mutation.report.x'], data[tmp[1],'mutation.report.y'], sep = ' && '),
                                                        paste(data[i,'mutation.report.x'], data[i,'mutation.report.y'], sep = ' && '),data[i,'n'])
    
  }
  return(network.direction)
}






  	