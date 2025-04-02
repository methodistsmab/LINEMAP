setwd('/Users/linwang/Work/dingcheng lineage tracing/amplican crispr analysis/result all cell 1')
setwd('./cell.read.matrix')

library(dplyr)

setwd('./d14')

cell.read.matrix.d14A21.1 = cell.read.matrix.d14A21[which(ind1 >0),]
cell.read.matrix.d14B25.1 = cell.read.matrix.d14B25[which(ind2 >0),]

tmp <-
  cell.read.matrix.d14A21.1 %>% 
  group_by(ID) %>% 
  dplyr::summarise(n = n())
tmp = as.data.frame(tmp)

cell.read.matrix.d14A21.1 = make.mutation.merge(cell.read.matrix.d14A21.1)

for( i in 1 : nrow(cell.read.matrix.d14A21.1)){
  
  result = make.mutation.report(cell.read.matrix.d14A21.1[i,'mutation.correction1'], cell.read.matrix.d14A21.1[i,"align.correction1"])
  cell.read.matrix.d14A21.1$number.of.different[i] = result[[1]]
  cell.read.matrix.d14A21.1$mutation.report[i] = result[[2]]
}

cell.read.matrix.d14A21.1$number.of.mutation = str_count(cell.read.matrix.d14A21.1$mutation.report, ";")

cell.read.matrix.d14A21.1$complexity.score = cell.read.matrix.d14A21.1$len.scaffold - cell.read.matrix.d14A21.1$number.of.different.scaffold -
  cell.read.matrix.d14A21.1$number.of.different - cell.read.matrix.d14A21.1$number.of.mutation

cell.read.matrix.d14A21.1$number.of.mismatch = lengths(regmatches(cell.read.matrix.d14A21.1[,"mutation.report"], gregexpr("mismatch", cell.read.matrix.d14A21.1[,"mutation.report"])))

cell.read.matrix.d14A21.1$length.of.insertion = lengths(regmatches(cell.read.matrix.d14A21.1[,"align.correction1"], gregexpr("-", cell.read.matrix.d14A21.1[,"align.correction1"])))

cell.read.matrix.d14A21.1$number.of.insertion = lengths(regmatches(cell.read.matrix.d14A21.1[,"mutation.report"], gregexpr("insert", cell.read.matrix.d14A21.1[,"mutation.report"])))

cell.read.matrix.d14A21.1$length.of.deletion = lengths(regmatches(cell.read.matrix.d14A21.1[,"mutation.correction1"], gregexpr("-", cell.read.matrix.d14A21.1[,"mutation.correction1"])))

cell.read.matrix.d14A21.1$number.of.deletion = lengths(regmatches(cell.read.matrix.d14A21.1[,"mutation.report"], gregexpr("miss", cell.read.matrix.d14A21.1[,"mutation.report"])))


cell.read.d14A21 <- cell.read.matrix.d14A21.1 %>% 
group_by(ID) %>%
filter(complexity.score == max(complexity.score)) 
cell.read.d14A21 = as.data.frame(cell.read.d14A21)

cell.read.d14A21 = cell.read.d14A21[!duplicated(cell.read.d14A21[,c('ID','mutation.correction1','align.correction1')]),]

cell.read.d14A21 <- cell.read.d14A21 %>% 
  group_by(ID) %>%
  filter(number.of.mismatch == min(number.of.mismatch)) 
cell.read.d14A21 = as.data.frame(cell.read.d14A21)

cell.read.d14A21 <- cell.read.d14A21 %>% 
  group_by(ID) %>%
  filter(number.of.different == min(number.of.different)) 
cell.read.d14A21 = as.data.frame(cell.read.d14A21)


cell.read.d14A21 <- cell.read.d14A21 %>% 
  group_by(ID) %>%
  filter(count == max(count)) 
cell.read.d14A21 = as.data.frame(cell.read.d14A21)


# cell.read.d14A21[which(cell.read.d14A21$ID == 'GTGCAGCAGATGAGGG'),]
# cell.read.d14B25[which(cell.read.d14B25$ID == 'GTGCAGCAGATGAGGG'),]

tmp <-
  cell.read.matrix.d14B25.1 %>% 
  group_by(ID) %>% 
  dplyr::summarise(n = n())
tmp = as.data.frame(tmp)

cell.read.matrix.d14B25.1 = make.mutation.merge.B25(cell.read.matrix.d14B25.1)

for( i in 1 : nrow(cell.read.matrix.d14B25.1)){
  
  result = make.mutation.report.B25(cell.read.matrix.d14B25.1[i,'mutation.correction1'], cell.read.matrix.d14B25.1[i,"align.correction1"])
  cell.read.matrix.d14B25.1$number.of.different[i] = result[[1]]
  cell.read.matrix.d14B25.1$mutation.report[i] = result[[2]]
}

cell.read.matrix.d14B25.1$number.of.mutation = str_count(cell.read.matrix.d14B25.1$mutation.report, ";")

cell.read.matrix.d14B25.1$complexity.score = cell.read.matrix.d14B25.1$len.scaffold - cell.read.matrix.d14B25.1$number.of.different.scaffold -
  cell.read.matrix.d14B25.1$number.of.different - cell.read.matrix.d14B25.1$number.of.mutation

cell.read.matrix.d14B25.1$number.of.mismatch = lengths(regmatches(cell.read.matrix.d14B25.1[,"mutation.report"], gregexpr("mismatch", cell.read.matrix.d14B25.1[,"mutation.report"])))

cell.read.matrix.d14B25.1$length.of.insertion = lengths(regmatches(cell.read.matrix.d14B25.1[,"align.correction1"], gregexpr("-", cell.read.matrix.d14B25.1[,"align.correction1"])))

cell.read.matrix.d14B25.1$number.of.insertion = lengths(regmatches(cell.read.matrix.d14B25.1[,"mutation.report"], gregexpr("insert", cell.read.matrix.d14B25.1[,"mutation.report"])))

cell.read.matrix.d14B25.1$length.of.deletion = lengths(regmatches(cell.read.matrix.d14B25.1[,"mutation.correction1"], gregexpr("-", cell.read.matrix.d14B25.1[,"mutation.correction1"])))

cell.read.matrix.d14B25.1$number.of.deletion = lengths(regmatches(cell.read.matrix.d14B25.1[,"mutation.report"], gregexpr("miss", cell.read.matrix.d14B25.1[,"mutation.report"])))



cell.read.d14B25 <- cell.read.matrix.d14B25.1 %>% 
  group_by(ID) %>%
  filter(complexity.score == max(complexity.score)) 
cell.read.d14B25 = as.data.frame(cell.read.d14B25)

cell.read.d14B25 = cell.read.d14B25[!duplicated(cell.read.d14B25[,c('ID','mutation.correction1','align.correction1')]),]

cell.read.d14B25 <- cell.read.d14B25 %>% 
  group_by(ID) %>%
  filter(number.of.mismatch == min(number.of.mismatch)) 
cell.read.d14B25 = as.data.frame(cell.read.d14B25)

cell.read.d14B25 <- cell.read.d14B25 %>% 
  group_by(ID) %>%
  filter(number.of.different == min(number.of.different)) 
cell.read.d14B25 = as.data.frame(cell.read.d14B25)

cell.read.d14B25 <- cell.read.d14B25 %>% 
  group_by(ID) %>%
  filter(count == max(count)) 
cell.read.d14B25 = as.data.frame(cell.read.d14B25)


length(unique(cell.read.d14A21$ID))
length(unique(cell.read.d14B25$ID))

nrow(cell.read.d14A21)
nrow(cell.read.d14B25)


cell.read.d14A21[duplicated(cell.read.d14A21$ID) | duplicated(cell.read.d14A21$ID, fromLast=TRUE),]
cell.read.d14B25[duplicated(cell.read.d14B25$ID) | duplicated(cell.read.d14B25$ID, fromLast=TRUE),]

common.cell.ID = intersect(unique(cell.read.d14A21$ID), unique(cell.read.d14B25$ID))
d14A21 = cell.read.d14A21[which(cell.read.d14A21$ID %in% common.cell.ID),]
d14B25 = cell.read.d14B25[which(cell.read.d14B25$ID %in% common.cell.ID),]

mainDir = '/Users/linwang/Work/dingcheng lineage tracing/amplican crispr analysis/result all cell 1'
subDir <- "cell.read"
subDir1 <- 'd14'
dir.create(file.path(mainDir, subDir, subDir1), showWarnings = FALSE)
setwd(file.path(mainDir, subDir, subDir1))

save(cell.read.d14A21,cell.read.d14B25,d14A21, d14B25, file = 
       'cell.read.d14.RData')

#########################################################################

make.mutation.merge <- function(read.matrix){
  tmp = unique(read.matrix[,c('mutation.correction','align.correction','mutation.report')])
  tmp1 = make.mutation.report.matrix(tmp[,3])
  
  tmp$mutation.correction1 = tmp$mutation.correction
  tmp$align.correction1 = tmp$align.correction
  
  tmp.ind = c()
  for (i in 1 : nrow(tmp)){
    if(!is.na(tmp1[[i]][24,1])){
      tmp.ind = c(tmp.ind,i)
      substr(tmp[i,'mutation.correction1'], 1 , as.numeric(tmp1[[i]][24,3])) = paste(rep('*',as.numeric(tmp1[[i]][24,3])), collapse = "")
      substr(tmp[i,'align.correction1'], 1 , as.numeric(tmp1[[i]][24,3])) = paste(rep('*',as.numeric(tmp1[[i]][24,3])), collapse = "")
      
      tmp[i,'mutation.correction1'] = gsub('\\*','', tmp[i,'mutation.correction1'] )
      tmp[i,'align.correction1'] = gsub('\\*','', tmp[i,'align.correction1'] )
    }
  }
  
  tmp.ind = which(tmp[,1] == 'GGGTTCCCAGT------------GGG')
  if(length(tmp.ind) >0){
    tmp[tmp.ind,'mutation.correction1'] = 'GGGTTCCC----AGT-------GGG'
    tmp[tmp.ind,'align.correction1'] =    'G--TTCCCGTCCAGTAATCGTGGGG'
    tmp[tmp.ind,'mutation.report'] = make.mutation.report('GGGTTCCC----AGT-------GGG','G--TTCCCGTCCAGTAATCGTGGGG')[[2]]
    tmp1[[tmp.ind]] = make.mutation.report.matrix(tmp[tmp.ind,'mutation.report'])[[1]]
  }

  tmp.ind = c()
  for (i in 1 : nrow(tmp)){
    if(!is.na(tmp1[[i]][25,1])){
      tmp.ind = c(tmp.ind,i)
      substr(tmp[i,'mutation.correction1'], 2 , (as.numeric(tmp1[[i]][25,3]) +1)) = paste(rep('*',as.numeric(tmp1[[i]][25,3])), collapse = "")
      substr(tmp[i,'align.correction1'], 2 , (as.numeric(tmp1[[i]][25,3]) +1)) = paste(rep('*',as.numeric(tmp1[[i]][25,3])), collapse = "")
      
      tmp[i,'mutation.correction1'] = gsub('\\*','', tmp[i,'mutation.correction1'] )
      tmp[i,'align.correction1'] = gsub('\\*','', tmp[i,'align.correction1'] )
    }
  }
  
  tmp.ind = c()
  for (i in 1 : nrow(tmp)){
    if(!is.na(tmp1[[i]][26,1])){
      tmp.ind = c(tmp.ind,i)
      substr(tmp[i,'mutation.correction1'], 3 , (as.numeric(tmp1[[i]][26,3]) +2)) = paste(rep('*',as.numeric(tmp1[[i]][26,3])), collapse = "")
      substr(tmp[i,'align.correction1'], 3 , (as.numeric(tmp1[[i]][26,3]) +2)) = paste(rep('*',as.numeric(tmp1[[i]][26,3])), collapse = "")

      tmp[i,'mutation.correction1'] = gsub('\\*','', tmp[i,'mutation.correction1'] )
      tmp[i,'align.correction1'] = gsub('\\*','', tmp[i,'align.correction1'] )
    }
  }
  
  read.matrix = left_join(read.matrix, tmp[c(1,2,4,5)], by=c('mutation.correction','align.correction'))
  
  return(read.matrix)

}
  

make.mutation.merge.B25 <- function(read.matrix){
  tmp = unique(read.matrix[,c('mutation.correction','align.correction','mutation.report')])
  tmp1 = make.mutation.report.matrix.B25(tmp[,3])
  
  tmp$mutation.correction1 = tmp$mutation.correction
  tmp$align.correction1 = tmp$align.correction
  
  tmp.ind = c()
  for (i in 1 : nrow(tmp)){
    if(!is.na(tmp1[[i]][28,1])){
      tmp.ind = c(tmp.ind,i)
      substr(tmp[i,'mutation.correction1'], 1 , as.numeric(tmp1[[i]][28,3])) = paste(rep('*',as.numeric(tmp1[[i]][28,3])), collapse = "")
      substr(tmp[i,'align.correction1'], 1 , as.numeric(tmp1[[i]][28,3])) = paste(rep('*',as.numeric(tmp1[[i]][28,3])), collapse = "")
      
      tmp[i,'mutation.correction1'] = gsub('\\*','', tmp[i,'mutation.correction1'] )
      tmp[i,'align.correction1'] = gsub('\\*','', tmp[i,'align.correction1'] )
    }
  }
  
  
  tmp.ind = c()
  for (i in 1 : nrow(tmp)){
    if(!is.na(tmp1[[i]][29,1])){
      tmp.ind = c(tmp.ind,i)
      substr(tmp[i,'mutation.correction1'], 2 , (as.numeric(tmp1[[i]][29,3]) +1)) = paste(rep('*',as.numeric(tmp1[[i]][29,3])), collapse = "")
      substr(tmp[i,'align.correction1'], 2 , (as.numeric(tmp1[[i]][29,3]) +1)) = paste(rep('*',as.numeric(tmp1[[i]][29,3])), collapse = "")
      
      tmp[i,'mutation.correction1'] = gsub('\\*','', tmp[i,'mutation.correction1'] )
      tmp[i,'align.correction1'] = gsub('\\*','', tmp[i,'align.correction1'] )
    }
  }
  
  tmp.ind = c()
  for (i in 1 : nrow(tmp)){
    if(!is.na(tmp1[[i]][30,1])){
      tmp.ind = c(tmp.ind,i)
      substr(tmp[i,'mutation.correction1'], 3 , (as.numeric(tmp1[[i]][30,3]) +2)) = paste(rep('*',as.numeric(tmp1[[i]][30,3])), collapse = "")
      substr(tmp[i,'align.correction1'], 3 , (as.numeric(tmp1[[i]][30,3]) +2)) = paste(rep('*',as.numeric(tmp1[[i]][30,3])), collapse = "")
      
      tmp[i,'mutation.correction1'] = gsub('\\*','', tmp[i,'mutation.correction1'] )
      tmp[i,'align.correction1'] = gsub('\\*','', tmp[i,'align.correction1'] )
    }
  }
  
  read.matrix = left_join(read.matrix, tmp[c(1,2,4,5)], by=c('mutation.correction','align.correction'))
  
  return(read.matrix)
  
}

  
  