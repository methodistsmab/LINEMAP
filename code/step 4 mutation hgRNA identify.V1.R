setwd('/Users/linwang/Work/dingcheng lineage tracing/amplican crispr analysis/result all cell 1')
setwd('./cell.read.matrix')

library(dplyr)

setwd('./d14')

load('cell.read.matrix.d14.RData')
cell.read.matrix.d14A21 = cell.read.matrix
load('mutation.report.matrix.d14.Rdata')
mutation.report.matrix.d14A21 = mutation.report.matrix

load('cell.read.matrix.B25.d14.RData')
cell.read.matrix.d14B25 = cell.read.matrix
load('mutation.report.matrix.B25.d14.Rdata')
mutation.report.matrix.d14B25 = mutation.report.matrix

##########################################################################


ind1 = rep(0,nrow(cell.read.matrix.d14A21))
ind2 = rep(0,nrow(cell.read.matrix.d14B25))

for (i in 1 : length(ind1)){
  if(cell.read.matrix.d14A21$number.of.different[i] < cell.read.matrix.d14B25$number.of.different[i] &
     cell.read.matrix.d14A21$complexity.score[i] > cell.read.matrix.d14B25$complexity.score[i] &
     cell.read.matrix.d14A21$len.substring.A21[i] >= cell.read.matrix.d14B25$len.substring.B25[i]){
    ind1[i] =1
  }
  
  if(cell.read.matrix.d14A21$number.of.different[i] > cell.read.matrix.d14B25$number.of.different[i] &
     cell.read.matrix.d14A21$complexity.score[i] < cell.read.matrix.d14B25$complexity.score[i] &
     cell.read.matrix.d14A21$len.substring.A21[i] <= cell.read.matrix.d14B25$len.substring.B25[i]){
    ind2[i] =1
  }
  
  if(ind1[i] == 0 & ind2[i] == 0){
    if(cell.read.matrix.d14A21$len.substring.A21[i] - cell.read.matrix.d14B25$len.substring.B25[i] >=4){
      ind1[i] = 0.5
    }
    if(cell.read.matrix.d14A21$len.substring.A21[i] - cell.read.matrix.d14B25$len.substring.B25[i] <= -4){
      ind2[i] = 0.5
    }
  }
  
}


ind1[cell.read.matrix.d14A21$number.of.different.scaffold/cell.read.matrix.d14A21$len.scaffold > 0.1] = -1
ind1[grepl('notinread',cell.read.matrix.d14A21$mutation.report)] = -1
ind1[cell.read.matrix.d14A21$length.of.deletion >17] = -1
ind1[cell.read.matrix.d14A21$number.of.mismatch/(23-cell.read.matrix.d14A21$length.of.deletion) >= 0.2] = -1


ind2[cell.read.matrix.d14B25$number.of.different.scaffold/cell.read.matrix.d14B25$len.scaffold > 0.1] = -1
ind2[grepl('notinread',cell.read.matrix.d14B25$mutation.report)] = -1
ind2[cell.read.matrix.d14B25$length.of.deletion >17] = -1
ind2[cell.read.matrix.d14B25$number.of.mismatch/(27-cell.read.matrix.d14B25$length.of.deletion) >= 0.2] = -1


table(ind1,ind2)


ind = rep(0,nrow(cell.read.matrix.d14A21))
ind[which(ind1 >0)] =1
ind[which(ind2 >0)] =2


########################################################### checking the identification
plot(cell.read.matrix.d14A21$len.substring.A21,cell.read.matrix.d14B25$len.substring.B25,col = ind, pch = ind)
plot(cell.read.matrix.d14A21$complexity.score,cell.read.matrix.d14B25$complexity.score, col = ind, pch = ind)
plot(cell.read.matrix.d14A21$number.of.different,cell.read.matrix.d14B25$number.of.different, col = ind, pch = ind)


table(cell.read.matrix.d14A21$number.of.different[ind1==0])
table(cell.read.matrix.d14B25$number.of.different[ind2==0])

table(cell.read.matrix.d14A21$complexity.score[ind1==0])
table(cell.read.matrix.d14B25$complexity.score[ind2==0])

table(cell.read.matrix.d14A21$len.substring.A21[ind1<=0])
table(cell.read.matrix.d14B25$len.substring.B25[ind2==0])


length(unique(cell.read.matrix.d14A21$ID))
length(unique(cell.read.matrix.d14A21$ID[ind !=0]))
length(unique(cell.read.matrix.d14A21$ID[ind ==0]))


plot(cell.read.matrix.d14A21$len.substring.A21[ind==0],cell.read.matrix.d14B25$len.substring.B25[ind==0])

plot(cell.read.matrix.d14A21$complexity.score[ind==0],cell.read.matrix.d14B25$complexity.score[ind==0])
plot(cell.read.matrix.d14A21$number.of.different[ind==0],cell.read.matrix.d14B25$number.of.different[ind==0])


