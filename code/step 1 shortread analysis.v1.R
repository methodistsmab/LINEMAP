setwd('/Users/linwang/Work/dingcheng lineage tracing')
setwd('amplican crispr analysis')
setwd('./fastq/')

library(ggplot2)
library(ShortRead)
library("dplyr")

methods(readFastq)
methods(writeFastq)
methods(countFastq)

# rm(list=ls())
# rm()


# fqr1 = readFastq('Lung_G1_CRISPR_S19_L002_R1_001.fastq.gz')
# fqr2 = readFastq('Lung_G1_CRISPR_S19_L002_R2_001.fastq.gz')
# 88381176
# 23785595


# fqr1 = readFastq('G1_CRISPR_L004_R1_001.fastq.gz')
# fqr2 = readFastq('G1_CRISPR_L004_R2_001.fastq.gz')
# 76225233
# 15632405

fqr1 = readFastq('G2_CRISPR_L004_R1_001.fastq.gz')
fqr2 = readFastq('G2_CRISPR_L004_R2_001.fastq.gz')

# 77120444
# 13541138

fqr1 = readFastq('Lung_G2_CRISPR_S20_L002_R1_001.fastq.gz')
fqr2 = readFastq('Lung_G2_CRISPR_S20_L002_R2_001.fastq.gz')

# 136522750
# 73880836

fqr1 = readFastq('./cell_line/D5_CRISPR_S5_L003_R1_001.fastq.gz')
fqr2 = readFastq('./cell_line/D5_CRISPR_S5_L003_R2_001.fastq.gz')
# 96126978
# 

fqr1 = readFastq('./cell_line/D14_CRISPR_S6_L003_R1_001.fastq.gz')
fqr2 = readFastq('./cell_line/D14_CRISPR_S6_L003_R2_001.fastq.gz')

############################################################## fastq file quality check
# quals = quality(fqr1)
# head(quals)
# numscores = as(quals,'matrix')
# which(is.na(numscores))
# 
# avgscore = rowMeans(numscores, na.rm = T)
# avgscore = as.data.frame(avgscore)
# 
# ggplot(avgscore) + geom_histogram(aes(x = avgscore),binwidth = 2)
# 
# which(avgscore)
# 
# 
# quals = quality(fqr2)
# head(quals)
# numscores = as(quals,'matrix')
# which(is.na(numscores))
# 
# avgscore = rowMeans(numscores, na.rm = T)
# avgscore = as.data.frame(avgscore)
# 
# ggplot(avgscore) + geom_histogram(aes(x = avgscore),binwidth = 2)
# 
# which(avgscore <20) 
#####################################################################

barcodelist = substring(as.character(sread(fqr1)),1,16)

umilist = substring(as.character(sread(fqr1)),17,28)
rm(fqr1)

read.string = as.character(sread(fqr2))

rm(fqr2)

mainDir = '/Users/linwang/Work/dingcheng lineage tracing/amplican crispr analysis/result all cell 1/fastqinformation'
subDir <- "d14"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
save(barcodelist,umilist,read.string, file = 'fastq.abstract.RData')
# load('fastq.abstract.RData')

#####################################################################

N.index = which(grepl('N',read.string,fixed = TRUE))

unique.read.string = unique(read.string)

len = length(unique.read.string)
num = 0
longest.string = c()

for (i in 1 : len){
  if(i %% 1000 ==0){
    print(i)
  }
  longest.string.tmp = find.longest.substring('TTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTC',unique.read.string[i])
  longest.string = c(longest.string,longest.string.tmp)
  if (nchar(longest.string.tmp) >10){
    num = num + 1
    print(paste0('number of read >10: ',num))
  }
}
options(scipen = 100)
write.table(longest.string, file = 'longest.string.0_5000000_lung1.txt', quote = FALSE, sep = '\t', col.names = FALSE, row.names = FALSE)


######################################################################################################## function 
# find all the matched substrings
##########################
find.substring.lacation <- function(substring, string){
	len.string = nchar(string)
	len.substring = nchar(substring)
	first = as.numeric(NA)
	last = as.numeric(NA)
	for( i in 1 : (len.string - len.substring +1)){
		if (identical(substring, substring(string, first = i, last = (i + len.substring -1)))){
			first = i
			last = i + len.substring -1
			break
		}
	}
	return(c(first,last))
}

find.longest.substring <- function(string1, string2.list){
  stringreturn = c()
  for (i in 1 : length(string2.list)){
    string2 = string2.list[i]
    if (nchar(string1) < nchar(string2)){
      a = string1
      b = string2
    }else{
      a = string2
      b = string1
    }
    
    len1 = nchar(a)
    len2 = nchar(b)
    sub.a = c()
    for (i in 1 : len1){
      sub.sub = substring(a, first = i, last = i:len1)
      sub.a = c(sub.a,sub.sub)
    }
    sub.a = sub.a[order(nchar(sub.a),sub.a, decreasing = TRUE)]
    
    length.subset.a = length(sub.a)
    
    for (i in 1 : length.subset.a){
      if (grepl(sub.a[i],b, fixed = TRUE)){
        break
      }
    }
    stringreturn = c(stringreturn,sub.a[i])
  }
    	return(stringreturn)
}

########################################################################################################




















