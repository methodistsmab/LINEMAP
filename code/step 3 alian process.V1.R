setwd('/Users/linwang/Work/dingcheng lineage tracing')
setwd('./amplican crispr analysis')
setwd('./result all cell 1')

library("stringr")
library('RecordLinkage')

mainDir = '/Users/linwang/Work/dingcheng lineage tracing/amplican crispr analysis/result all cell 1'
subDir <- "fastqinformation"
subDir1 <- 'd14'
path = file.path(mainDir,subDir,subDir1)
setwd(path)
tmp = load('fastqinformation.RData')

subDir <- "aligment"
subDir1 <- 'd14/V3'
path = file.path(mainDir,subDir,subDir1)
setwd(path)

d14 = read.delim('alignments.txt', header = FALSE)

cell.read.align = d14

cell.read.matrix = make.align.matrix(cell.read.align)

cell.read.matrix = add.umi(cell.read.matrix, count.table)

cell.read.matrix = cell.read.matrix[!is.na(cell.read.matrix$no.umi),]

# cell.read.matrix.choose = choose.cell.read.matrix(cell.read.matrix)

length(unique(cell.read.matrix[,4]))
# 6562 d14
# 5011 d14
# 3350 d14
length(unique(cell.read.matrix[,5]))
# 1028 d14
# 801 d14
# 257 d14


cell.read.matrix$read.short = NA
cell.read.matrix$align.short = NA
cell.read.matrix$scaffold.read = NA
cell.read.matrix$scaffold.align = NA
cell.read.matrix$number.of.different = NA
cell.read.matrix$mutation.report = NA

for( i in 1 : nrow(cell.read.matrix)){
  result = read.separate(cell.read.matrix[i,'read'], cell.read.matrix[i,'align'])
  if(length(result) ==4){
    cell.read.matrix$read.short[i] = result[1]
    cell.read.matrix$align.short[i] = result[2]
    cell.read.matrix$scaffold.read[i] = result[3]
    cell.read.matrix$scaffold.align[i] = result[4]
    print(i)
  }
  
  result = make.mutation.report(cell.read.matrix[i,'read.short'], cell.read.matrix[i,"align.short"])
  cell.read.matrix$number.of.different[i] = result[[1]]
  cell.read.matrix$mutation.report[i] = result[[2]]
}

cell.read.matrix = mutation.correction(cell.read.matrix)

for( i in 1 : nrow(cell.read.matrix)){
  
  result = make.mutation.report(cell.read.matrix[i,'mutation.correction'], cell.read.matrix[i,"align.correction"])
  cell.read.matrix$number.of.different[i] = result[[1]]
  cell.read.matrix$mutation.report[i] = result[[2]]
}

cell.read.matrix$len.scaffold = nchar(cell.read.matrix$scaffold.read)

cell.read.matrix$number.of.different.scaffold = NA
for (i in 1 : nrow(cell.read.matrix)){
  cell.read.matrix$number.of.different.scaffold[i]= 
    check.string.different(cell.read.matrix$scaffold.read[i],cell.read.matrix$scaffold.align[i])
}

cell.read.matrix$number.of.mutation = str_count(cell.read.matrix$mutation.report, ";")

cell.read.matrix$complexity.score = cell.read.matrix$len.scaffold - cell.read.matrix$number.of.different.scaffold -
  cell.read.matrix$number.of.different - cell.read.matrix$number.of.mutation

cell.read.matrix$number.of.mismatch = lengths(regmatches(cell.read.matrix[,"mutation.report"], gregexpr("mismatch", cell.read.matrix[,"mutation.report"])))

cell.read.matrix$length.of.insertion = lengths(regmatches(cell.read.matrix[,"align.correction"], gregexpr("-", cell.read.matrix[,"align.correction"])))

cell.read.matrix$number.of.insertion = lengths(regmatches(cell.read.matrix[,"mutation.report"], gregexpr("insert", cell.read.matrix[,"mutation.report"])))

cell.read.matrix$length.of.deletion = lengths(regmatches(cell.read.matrix[,"mutation.correction"], gregexpr("-", cell.read.matrix[,"mutation.correction"])))

cell.read.matrix$number.of.deletion = lengths(regmatches(cell.read.matrix[,"mutation.report"], gregexpr("miss", cell.read.matrix[,"mutation.report"])))

cell.read.matrix$number.of.notinread = lengths(regmatches(cell.read.matrix[,"mutation.report"], gregexpr("\\*", cell.read.matrix[,"mutation.report"])))



A21 = 'GTTCCCGTCCAGTAATCGTG'
B25 = 'GTCGTTGTAGCAACCTATCGGGTG'
for (i in 1 : nrow(cell.read.matrix)){
  if(i %% 1000 ==0){
    print(i)
  }
  short.read = gsub('-','',cell.read.matrix[i,'mutation.correction'])
  cell.read.matrix$substring.A21[i] = find.longest.substring(A21,short.read)
  cell.read.matrix$len.substring.A21[i] = nchar(cell.read.matrix$substring.A21[i])
}


subDir <- "cell.read.matrix"
subDir1 <- 'd14'
dir.create(file.path(mainDir, subDir, subDir1), showWarnings = FALSE)
setwd(file.path(mainDir, subDir, subDir1))

write.table(cell.read.matrix, file= 'cell.read.matrix.d14.txt',quote = FALSE, sep = '\t',row.names = FALSE)
save(cell.read.matrix, file = 'cell.read.matrix.d14.RData')


mutation.report.matrix = make.mutation.report.matrix(cell.read.matrix[,'mutation.report'])

save(mutation.report.matrix, file = 'mutation.report.matrix.d14.RData')

############################################

subDir <- "aligment"
subDir1 <- 'd14/v3B25'
path = file.path(mainDir,subDir,subDir1)
setwd(path)

d14 = read.delim('alignments.txt', header = FALSE)

cell.read.align = d14

cell.read.matrix = make.align.matrix(cell.read.align)

cell.read.matrix = add.umi(cell.read.matrix, count.table)

cell.read.matrix = cell.read.matrix[!is.na(cell.read.matrix$no.umi),]

length(unique(cell.read.matrix[,4]))
# 6407 d14
# 4930 d14
# 3297 d14
length(unique(cell.read.matrix[,5]))
# 923 d14
# 677 d14
# 136 d14


cell.read.matrix$read.short = NA
cell.read.matrix$align.short = NA
cell.read.matrix$scaffold.read = NA
cell.read.matrix$scaffold.align = NA
cell.read.matrix$number.of.different = NA
cell.read.matrix$mutation.report = NA

for( i in 1 : nrow(cell.read.matrix)){
  result = read.separate.B25(cell.read.matrix[i,'read'], cell.read.matrix[i,'align'])
  if(length(result) ==4){
    cell.read.matrix$read.short[i] = result[1]
    cell.read.matrix$align.short[i] = result[2]
    cell.read.matrix$scaffold.read[i] = result[3]
    cell.read.matrix$scaffold.align[i] = result[4]
    print(i)
  }
  
  result = make.mutation.report.B25(cell.read.matrix[i,'read.short'], cell.read.matrix[i,"align.short"])
  cell.read.matrix$number.of.different[i] = result[[1]]
  cell.read.matrix$mutation.report[i] = result[[2]]
}

cell.read.matrix = mutation.correction.B25(cell.read.matrix)

for( i in 1 : nrow(cell.read.matrix)){
  
  result = make.mutation.report.B25(cell.read.matrix[i,'mutation.correction'], cell.read.matrix[i,"align.correction"])
  cell.read.matrix$number.of.different[i] = result[[1]]
  cell.read.matrix$mutation.report[i] = result[[2]]
}

cell.read.matrix$len.scaffold = nchar(cell.read.matrix$scaffold.read)

cell.read.matrix$number.of.different.scaffold = NA
for (i in 1 : nrow(cell.read.matrix)){
  cell.read.matrix$number.of.different.scaffold[i]= 
    check.string.different(cell.read.matrix$scaffold.read[i],cell.read.matrix$scaffold.align[i])
}


cell.read.matrix$number.of.mutation = str_count(cell.read.matrix$mutation.report, ";")

cell.read.matrix$complexity.score = cell.read.matrix$len.scaffold - cell.read.matrix$number.of.different.scaffold - 
  cell.read.matrix$number.of.different - cell.read.matrix$number.of.mutation

cell.read.matrix$number.of.mismatch = lengths(regmatches(cell.read.matrix[,"mutation.report"], gregexpr("mismatch", cell.read.matrix[,"mutation.report"])))

cell.read.matrix$length.of.insertion = lengths(regmatches(cell.read.matrix[,"align.correction"], gregexpr("-", cell.read.matrix[,"align.correction"])))

cell.read.matrix$number.of.insertion = lengths(regmatches(cell.read.matrix[,"mutation.report"], gregexpr("insert", cell.read.matrix[,"mutation.report"])))

cell.read.matrix$length.of.deletion = lengths(regmatches(cell.read.matrix[,"mutation.correction"], gregexpr("-", cell.read.matrix[,"mutation.correction"])))

cell.read.matrix$number.of.deletion = lengths(regmatches(cell.read.matrix[,"mutation.report"], gregexpr("miss", cell.read.matrix[,"mutation.report"])))

cell.read.matrix$number.of.notinread = lengths(regmatches(cell.read.matrix[,"mutation.report"], gregexpr("\\*", cell.read.matrix[,"mutation.report"])))

A21 = 'GTTCCCGTCCAGTAATCGTG'
B25 = 'GTCGTTGTAGCAACCTATCGGGTG'
for (i in 1 : nrow(cell.read.matrix)){
  if(i %% 1000 ==0){
    print(i)
  }
  short.read = gsub('-','',cell.read.matrix[i,'mutation.correction'])
  cell.read.matrix$substring.B25[i] = find.longest.substring(B25,short.read)
  cell.read.matrix$len.substring.B25[i] = nchar(cell.read.matrix$substring.B25[i])
  
}

subDir <- "cell.read.matrix"
subDir1 <- 'd14'
dir.create(file.path(mainDir, subDir, subDir1), showWarnings = FALSE)
setwd(file.path(mainDir, subDir, subDir1))

write.table(cell.read.matrix, file= 'cell.read.matrix.B25.d14.txt',quote = FALSE, sep = '\t',row.names = FALSE)
save(cell.read.matrix, file = 'cell.read.matrix.B25.d14.RData')

mutation.report.matrix = make.mutation.report.matrix.B25(cell.read.matrix[,'mutation.report'])

save(mutation.report.matrix, file = 'mutation.report.matrix.B25.d14.RData')


################################################# functions
make.align.matrix <- function(cell.read.align){
	cell.read.info = strsplit(cell.read.align[seq(1,nrow(cell.read.align),3),1], "\\s+")
    cell.read.info = matrix(unlist(cell.read.info), ncol = 6, nrow = nrow(cell.read.align)/3 , byrow = TRUE)
    cell.read.info = data.frame(cell.read.info)
    cell.read.info[,6] = as.numeric(cell.read.info[,6])
    
    cell.read.matrix = data.frame(matrix(ncol = 5, nrow = 0))
    colnames(cell.read.matrix) <- x <- c("ID", "order", "count", "read", "align")
    
    for(i in 1: nrow(cell.read.info)){
    	cell.read.matrix[i,"ID"] = cell.read.info[i,2]
      	cell.read.matrix[i,"order"] = cell.read.info[i,4]
      	cell.read.matrix[i,"count"] = cell.read.info[i,6]
      	cell.read.matrix[i,"read"] = cell.read.align[(i*3-1),1]
      	cell.read.matrix[i,"align"] = cell.read.align[(i*3),1]
      	print(i)
    }
    return(cell.read.matrix)
}


add.umi <- function(cell.read.matrix, count.table){
  cell.read.matrix$no.umi = 0
  no.row = nrow(cell.read.matrix)
  for(i in 1 : no.row){
    ind1 = which(count.table[,'barcodelist1'] %in% cell.read.matrix$ID[i])
    read.original = gsub("-", "", cell.read.matrix$read[i]) 
    ind2 = which(grepl(read.original, count.table[,'read.string1'], fixed = TRUE))
    ind = intersect(ind1,ind2)
    
    if (length(ind) == 1){
      cell.read.matrix$no.umi[i] = count.table[ind,'no.umi']
    }else{
      cell.read.matrix$no.umi[i] = count.table[ind,'no.umi'][match(cell.read.matrix[i,'count'], count.table[ind,"count"])]
    }
  }
  return(cell.read.matrix)
}

check.string.different <- function(string1, string2){
  tmp1 = strsplit(string1,'')[[1]]
  tmp2 = strsplit(string2,'')[[1]]
  return(length(which(tmp1 != tmp2)))
}

strat.end.of.underscore <- function(sequence){
	start.end.vector = c()
	sequence = strsplit(sequence,'')[[1]]
	
	ind = which(sequence == '-')
	diff.ind = which(diff(ind) >1)
	
	if (length(ind) ==0){
		return()
	}
	
	start.vector =c(ind[1],ind[(diff.ind+1)])
	end.vector = c(ind[diff.ind],ind[length(ind)])
	start.end.vector = 	cbind(start.vector,end.vector)
	return(start.end.vector)		
}

# read =  'AAGCAGTGGTATCAACGCAGAGTACATGGG------------TGAGTGTGTGTAAAAAGAAATAATTA------------------AATC'
# align = 'AAGCAGTGGTATCAACGCAGAGTACATGGGGTTCCCGTCCAGTAATCGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTC'


read.separate <- function(read, align){
  len1 = nchar(read)
  len2 = nchar(align)
  if (len1 != len2){
    print('different length')
    return()
  }
  
  if(substr(align, 1, 30) != 'AAGCAGTGGTATCAACGCAGAGTACATGGG'){
    print('primer wrong')
    return()
  }
  
  for ( j in 53: len2){
    target.align = substr(align, 31,j)
    if ((nchar(target.align)-str_count(target.align,'-')) ==23 ){
      break
    }
  }
  
  if (nchar(target.align) == 0){
    print('all hgRNA missed')
    return()
  }
  
  
  len2 = nchar(target.align)
  target.read = substr(read,31,(30+len2))
  
  scaffold.align =  substr(align, (30+ len2 + 1),nchar(align))
  scaffold.read =  substr(read, (30+ len2 + 1),nchar(align))
  
  result = c(target.read, target.align,scaffold.read,scaffold.align)
  return(result)
}


read.separate.B25 <- function(read, align){
  len1 = nchar(read)
  len2 = nchar(align)
  if (len1 != len2){
    print('different length')
    return()
  }
  
  if( substr(align, 1, 30) != 'AAGCAGTGGTATCAACGCAGAGTACATGGG'){
    print('primer wrong')
    return()
  }
  
  for ( j in 57: len2){
    target.align = substr(align, 31,j)
    if ((nchar(target.align)-str_count(target.align,'-')) ==27 ){
      break
    }
  }
  
  if (nchar(target.align) == 0){
    print('all hgRNA missed')
    return()
  }
  
  
  len2 = nchar(target.align)
  target.read = substr(read,31,(30+len2))
  
  scaffold.align =  substr(align, (30+ len2 + 1),nchar(align))
  scaffold.read =  substr(read, (30+ len2 + 1),nchar(align))
  
  result = c(target.read, target.align,scaffold.read,scaffold.align)
  return(result)
}



make.mutation.report <- function(target.read, target.align){
    
    align.start.end.underscore = strat.end.of.underscore(target.align)
    read.start.end.underscore = strat.end.of.underscore(target.read)
    
    target.align.list = strsplit(target.align,'')[[1]]
    target.read.list = strsplit(target.read,'')[[1]]   
     
    align.underscore_index = which(target.align.list =="-")
    read.underscore_index = which(target.read.list =="-")

    base.number = rep(0,nchar(target.align))
    base.number[which(target.align.list !='-')] = 1:length(which(target.align.list !='-'))
    
    
    mutation.report = ''
    i = 1

    while(i <= nchar(target.align) ){
    	
    	if( !i %in% align.underscore_index  & !i %in% read.underscore_index){
    		if (target.align.list[i] != target.read.list[i]){
    			mutation.report = paste(mutation.report, base.number[i],'mismatch',target.align.list[i], target.read.list[i])
    			mutation.report = paste(mutation.report,';',sep = '')
    		}
    		i = i + 1

    	}else if(i %in% align.underscore_index){
    		row.index = which(align.start.end.underscore[,1] %in% i)
    		align.underscore = substr(target.align,align.start.end.underscore[row.index,1],align.start.end.underscore[row.index,2])
    		insert.in.read = substr(target.read,align.start.end.underscore[row.index,1],align.start.end.underscore[row.index,2])
    		mutation.report = paste(mutation.report, paste((base.number[i-1]),'<',(align.start.end.underscore[row.index,2] - align.start.end.underscore[row.index,1]+1),sep = ''),'insert',align.underscore,insert.in.read)
    		mutation.report = paste(mutation.report,';',sep = '')
    		i = align.start.end.underscore[row.index,2] + 1 

    	}else if (i %in% read.underscore_index){
    		row.index = which( read.start.end.underscore[,1] %in% i)
    		should.in.align = substr(target.align,read.start.end.underscore[row.index,1],read.start.end.underscore[row.index,2])
    		missed.in.read = substr(target.read,read.start.end.underscore[row.index,1],read.start.end.underscore[row.index,2])
    		mutation.report = paste(mutation.report, paste(base.number[i],'-',base.number[read.start.end.underscore[row.index,2]],sep = ''),'miss',should.in.align,missed.in.read)
    		mutation.report = paste(mutation.report,';',sep = '')
            i = read.start.end.underscore[row.index,2] + 1
    	}
    	
    }


    number.of.different = length(which(target.align.list !=target.read.list))
    
    if ((nchar(target.align)-str_count(target.align,'-')) <23){
    	not.in.align = substr('GTTCCCGTCCAGTAATCGTGGGG',((nchar(target.align)-str_count(target.align,'-'))+1), 23)
    	mutation.report = paste(mutation.report, paste((nchar(target.align)-str_count(target.align,'-') +1), '-23', sep = ''), 'notinread', not.in.align,strrep('*',(23+str_count(target.align,'-')-nchar(target.align))))
        mutation.report = paste(mutation.report,';',sep = '')
    }
    number.of.different = number.of.different + (23+str_count(target.align,'-')-nchar(target.align))
    
    return(list(number.of.different, mutation.report))

}

make.mutation.report.B25 <- function(target.read, target.align){
  
  align.start.end.underscore = strat.end.of.underscore(target.align)
  read.start.end.underscore = strat.end.of.underscore(target.read)
  
  target.align.list = strsplit(target.align,'')[[1]]
  target.read.list = strsplit(target.read,'')[[1]]   
  
  align.underscore_index = which(target.align.list =="-")
  read.underscore_index = which(target.read.list =="-")
  
  base.number = rep(0,nchar(target.align))
  base.number[which(target.align.list !='-')] = 1:length(which(target.align.list !='-'))
  
  
  mutation.report = ''
  i = 1
  
  while(i <= nchar(target.align) ){
    
    if( !i %in% align.underscore_index  & !i %in% read.underscore_index){
      if (target.align.list[i] != target.read.list[i]){
        mutation.report = paste(mutation.report, base.number[i],'mismatch',target.align.list[i], target.read.list[i])
        mutation.report = paste(mutation.report,';',sep = '')
      }
      i = i + 1
      
    }else if(i %in% align.underscore_index){
      row.index = which(align.start.end.underscore[,1] %in% i)
      align.underscore = substr(target.align,align.start.end.underscore[row.index,1],align.start.end.underscore[row.index,2])
      insert.in.read = substr(target.read,align.start.end.underscore[row.index,1],align.start.end.underscore[row.index,2])
      mutation.report = paste(mutation.report, paste((base.number[i-1]),'<',(align.start.end.underscore[row.index,2] - align.start.end.underscore[row.index,1]+1),sep = ''),'insert',align.underscore,insert.in.read)
      mutation.report = paste(mutation.report,';',sep = '')
      i = align.start.end.underscore[row.index,2] + 1 
      
    }else if (i %in% read.underscore_index){
      row.index = which( read.start.end.underscore[,1] %in% i)
      should.in.align = substr(target.align,read.start.end.underscore[row.index,1],read.start.end.underscore[row.index,2])
      missed.in.read = substr(target.read,read.start.end.underscore[row.index,1],read.start.end.underscore[row.index,2])
      mutation.report = paste(mutation.report, paste(base.number[i],'-',base.number[read.start.end.underscore[row.index,2]],sep = ''),'miss',should.in.align,missed.in.read)
      mutation.report = paste(mutation.report,';',sep = '')
      i = read.start.end.underscore[row.index,2] + 1
    }
    
  }
  
  
  number.of.different = length(which(target.align.list !=target.read.list))
  
  if ((nchar(target.align)-str_count(target.align,'-')) <27){
    not.in.align = substr('GTTCCCGTCCAGTAATCGTGGGG',((nchar(target.align)-str_count(target.align,'-'))+1), 27)
    mutation.report = paste(mutation.report, paste((nchar(target.align)-str_count(target.align,'-') +1), '-27', sep = ''), 'notinread', not.in.align,strrep('*',(27+str_count(target.align,'-')-nchar(target.align))))
    mutation.report = paste(mutation.report,';',sep = '')
  }
  number.of.different = number.of.different + (27+str_count(target.align,'-')-nchar(target.align))
  
  return(list(number.of.different, mutation.report))
  
}


make.mutation.report.matrix <- function(mutation.report){
	
	number.of.report = length(mutation.report)
	list.of.matrix = list()
	for (i in 1: number.of.report){
	  mutation.report.matrix = matrix(, nrow = 47, ncol = 5)
	  colnames(mutation.report.matrix) = c('mutation.type', 'location', 'number.of.base', 'original','actual')
	  mutation.report.matrix[,2] = c(1:23,0:23)
	  mutation.report.matrix[,3] = 0
	  
		if (is.na(mutation.report[i])){
			list.of.matrix[[i]] = mutation.report.matrix
			next
		}
	  if (mutation.report[i] == ''){
	    list.of.matrix[[i]] = mutation.report.matrix
	    next
	  }
	  
		report1 = mutation.report[i]
		report1 = strsplit(report1, split = ';')[[1]]
		report1 = str_trim(report1,'both')
			
		for(j in 1 : length(report1)){
			report2 = strsplit(report1[j],split = ' ')[[1]]
			if (report2[2] == 'insert'){
				tmp = as.numeric(strsplit(report2[1], split = '<')[[1]])
				if(is.na(tmp[1])){
					tmp[1] = 0
				}
				mutation.report.matrix[(23+tmp[1] +1),'mutation.type'] = 'insert'
				mutation.report.matrix[(23+tmp[1] +1),'number.of.base'] = tmp[2]
				mutation.report.matrix[(23+tmp[1] +1),'original'] = report2[3]
				mutation.report.matrix[(23+tmp[1] +1),'actual'] = report2[4]
			}else if (report2[2] == 'miss'){
				tmp = as.numeric(strsplit(report2[1], split = '-')[[1]])
				for(k in tmp[1] : tmp[2]){
					mutation.report.matrix[k,'mutation.type'] = 'miss'
				    mutation.report.matrix[k,'number.of.base'] = 1
				    mutation.report.matrix[k,'original'] = substring(report2[3],(k-tmp[1] + 1), (k-tmp[1] + 1))
				    mutation.report.matrix[k,'actual'] = '-'
				}
			} else if (report2[2] == 'notinread'){
				tmp = as.numeric(strsplit(report2[1], split = '-')[[1]])
				for (k in tmp[1] : tmp[2]){
					mutation.report.matrix[k,'mutation.type'] = 'miss*'
				    mutation.report.matrix[k,'number.of.base'] = 1
				    mutation.report.matrix[k,'original'] = substring(report2[3],(k-tmp[1] + 1), (k-tmp[1] + 1))
				    mutation.report.matrix[k,'actual'] = '*'
				} 
			}else {
				    mutation.report.matrix[as.numeric(report2[1]),'mutation.type'] = report2[2]
				    mutation.report.matrix[as.numeric(report2[1]),'number.of.base'] = 1
				    mutation.report.matrix[as.numeric(report2[1]),'original'] = report2[3]
				    mutation.report.matrix[as.numeric(report2[1]),'actual'] = report2[4]
			}
		}
		list.of.matrix[[i]] = mutation.report.matrix
	}
	return(list.of.matrix)
}


make.mutation.report.matrix.B25 <- function(mutation.report){
  
  number.of.report = length(mutation.report)
  list.of.matrix = list()
  for (i in 1: number.of.report){
    mutation.report.matrix = matrix(, nrow = 55, ncol = 5)
    colnames(mutation.report.matrix) = c('mutation.type', 'location', 'number.of.base', 'original','actual')
    mutation.report.matrix[,2] = c(1:27,0:27)
    mutation.report.matrix[,3] = 0
    
    if (is.na(mutation.report[i])){
      list.of.matrix[[i]] = mutation.report.matrix
      next
    }
    if (mutation.report[i] == ''){
      list.of.matrix[[i]] = mutation.report.matrix
      next
    }
    
    report1 = mutation.report[i]
    report1 = strsplit(report1, split = ';')[[1]]
    report1 = str_trim(report1,'both')
    
    
    for(j in 1 : length(report1)){
      report2 = strsplit(report1[j],split = ' ')[[1]]
      if (report2[2] == 'insert'){
        tmp = as.numeric(strsplit(report2[1], split = '<')[[1]])
        if(is.na(tmp[1])){
          tmp[1] = 0
        }
        mutation.report.matrix[(27+tmp[1] +1),'mutation.type'] = 'insert'
        mutation.report.matrix[(27+tmp[1] +1),'number.of.base'] = tmp[2]
        mutation.report.matrix[(27+tmp[1] +1),'original'] = report2[3]
        mutation.report.matrix[(27+tmp[1] +1),'actual'] = report2[4]
      }else if (report2[2] == 'miss'){
        tmp = as.numeric(strsplit(report2[1], split = '-')[[1]])
        for(k in tmp[1] : tmp[2]){
          mutation.report.matrix[k,'mutation.type'] = 'miss'
          mutation.report.matrix[k,'number.of.base'] = 1
          mutation.report.matrix[k,'original'] = substring(report2[3],(k-tmp[1] + 1), (k-tmp[1] + 1))
          mutation.report.matrix[k,'actual'] = '-'
        }
      } else if (report2[2] == 'notinread'){
        tmp = as.numeric(strsplit(report2[1], split = '-')[[1]])
        for (k in tmp[1] : tmp[2]){
          mutation.report.matrix[k,'mutation.type'] = 'miss*'
          mutation.report.matrix[k,'number.of.base'] = 1
          mutation.report.matrix[k,'original'] = substring(report2[3],(k-tmp[1] + 1), (k-tmp[1] + 1))
          mutation.report.matrix[k,'actual'] = '*'
        } 
      }else {
        mutation.report.matrix[as.numeric(report2[1]),'mutation.type'] = report2[2]
        mutation.report.matrix[as.numeric(report2[1]),'number.of.base'] = 1
        mutation.report.matrix[as.numeric(report2[1]),'original'] = report2[3]
        mutation.report.matrix[as.numeric(report2[1]),'actual'] = report2[4]
      }
    }
    list.of.matrix[[i]] = mutation.report.matrix
  }
  return(list.of.matrix)
}


mutation.correction <- function(cell.matrix){
  # ind = rep(0, nrow(cell.matrix))
  # ind[which(cell.matrix$number.of.different.scaffold/cell.matrix$len.scaffold > 0.1)] = -1
  # ind[which(grepl('notinread',cell.matrix$mutation.report))] = -1
  # ind[which(cell.matrix$number.of.mismatch >3)] = -1
  
  tmp = unique(cell.matrix[,c("read.short","align.short",'mutation.report')])
  tmp1 = make.mutation.report.matrix(tmp[,3])
  
  tmp$mutation.correction = tmp$read.short
  tmp$align.correction = tmp$align.short
  
  tmp.ind = c()
  for (i in 1 : nrow(tmp)){
    if ( is.na(tmp1[[i]][1,1]) &  is.na(tmp1[[i]][2,1]) & is.na(tmp1[[i]][3,1])){
      if( !is.na(tmp1[[i]][24,1])){
        tmp.ind = c( tmp.ind,i)
        substr(tmp[i,'mutation.correction'], 1 , as.numeric(tmp1[[i]][24,3])) = paste(rep('*',as.numeric(tmp1[[i]][24,3])), collapse = "")
        substr(tmp[i,'align.correction'], 1 , as.numeric(tmp1[[i]][24,3])) = paste(rep('*',as.numeric(tmp1[[i]][24,3])), collapse = "")
        
        tmp[i,'mutation.correction'] = gsub('\\*','', tmp[i,'mutation.correction'] )
        tmp[i,'align.correction'] = gsub('\\*','', tmp[i,'align.correction'] )
      }
    }
  }
  
  tmp.ind = c()
  for (i in 1 : nrow(tmp)){
    if ( is.na(tmp1[[i]][1,1]) &  is.na(tmp1[[i]][2,1]) & is.na(tmp1[[i]][3,1]) & is.na(tmp1[[i]][4,1])){
      if( !is.na(tmp1[[i]][25,1]) & substr(tmp1[[i]][25,5], nchar(tmp1[[i]][25,5]), nchar(tmp1[[i]][25,5])) == 'G' ){
        tmp.ind = c( tmp.ind,i)
        substr(tmp[i,'mutation.correction'], 2 , (as.numeric(tmp1[[i]][25,3]) +1)) = paste(rep('*',as.numeric(tmp1[[i]][25,3])), collapse = "")
        substr(tmp[i,'align.correction'], 2 , (as.numeric(tmp1[[i]][25,3]) +1)) = paste(rep('*',as.numeric(tmp1[[i]][25,3])), collapse = "")
        
        tmp[i,'mutation.correction'] = gsub('\\*','', tmp[i,'mutation.correction'] )
        tmp[i,'align.correction'] = gsub('\\*','', tmp[i,'align.correction'] )
      }
    }
  }
  
  
  
  tmp.ind = c()
  for (i in 1 : nrow(tmp)){
    if ( is.na(tmp1[[i]][1,1]) &  is.na(tmp1[[i]][2,1]) & is.na(tmp1[[i]][3,1]) & is.na(tmp1[[i]][4,1]) & is.na(tmp1[[i]][5,1])){
      if( !is.na(tmp1[[i]][26,1]) & substr(tmp1[[i]][26,5], nchar(tmp1[[i]][26,5])-1, nchar(tmp1[[i]][26,5])) == 'GT' ){
        tmp.ind = c( tmp.ind,i)
        substr(tmp[i,'mutation.correction'], 3 , (as.numeric(tmp1[[i]][26,3]) +2)) = paste(rep('*',as.numeric(tmp1[[i]][26,3])), collapse = "")
        substr(tmp[i,'align.correction'], 3 , (as.numeric(tmp1[[i]][26,3]) +2)) = paste(rep('*',as.numeric(tmp1[[i]][26,3])), collapse = "")
        
        tmp[i,'mutation.correction'] = gsub('\\*','', tmp[i,'mutation.correction'] )
        tmp[i,'align.correction'] = gsub('\\*','', tmp[i,'align.correction'] )
      }
    }
  }
  
  cell.matrix = left_join(cell.matrix, tmp[c(1,2,4,5)], by=c('read.short','align.short'))
  
  return(cell.matrix)
  
}


mutation.correction.B25 <- function(cell.matrix){
  
  tmp = unique(cell.matrix[,c("read.short","align.short",'mutation.report')])
  tmp1 = make.mutation.report.matrix.B25(tmp[,3])
  
  tmp$mutation.correction = tmp$read.short
  tmp$align.correction = tmp$align.short
  
  tmp.ind = c()
  for (i in 1 : nrow(tmp)){
    if ( is.na(tmp1[[i]][1,1]) &  is.na(tmp1[[i]][2,1]) & is.na(tmp1[[i]][3,1])){
      if( !is.na(tmp1[[i]][28,1])){
        tmp.ind = c( tmp.ind,i)
        substr(tmp[i,'mutation.correction'], 1 , as.numeric(tmp1[[i]][28,3])) = paste(rep('*',as.numeric(tmp1[[i]][28,3])), collapse = "")
        substr(tmp[i,'align.correction'], 1 , as.numeric(tmp1[[i]][28,3])) = paste(rep('*',as.numeric(tmp1[[i]][28,3])), collapse = "")
        
        tmp[i,'mutation.correction'] = gsub('\\*','', tmp[i,'mutation.correction'] )
        tmp[i,'align.correction'] = gsub('\\*','', tmp[i,'align.correction'] )
      }
    }
  }
  
  tmp.ind = c()
  for (i in 1 : nrow(tmp)){
    if ( is.na(tmp1[[i]][1,1]) &  is.na(tmp1[[i]][2,1]) & is.na(tmp1[[i]][3,1]) & is.na(tmp1[[i]][4,1])){
      if( !is.na(tmp1[[i]][29,1]) & substr(tmp1[[i]][29,5], nchar(tmp1[[i]][29,5]), nchar(tmp1[[i]][29,5])) == 'G' ){
        tmp.ind = c( tmp.ind,i)
        substr(tmp[i,'mutation.correction'], 2 , (as.numeric(tmp1[[i]][29,3]) +1)) = paste(rep('*',as.numeric(tmp1[[i]][29,3])), collapse = "")
        substr(tmp[i,'align.correction'], 2 , (as.numeric(tmp1[[i]][29,3]) +1)) = paste(rep('*',as.numeric(tmp1[[i]][29,3])), collapse = "")
        
        tmp[i,'mutation.correction'] = gsub('\\*','', tmp[i,'mutation.correction'] )
        tmp[i,'align.correction'] = gsub('\\*','', tmp[i,'align.correction'] )
      }
    }
  }
  
  
  
  tmp.ind = c()
  for (i in 1 : nrow(tmp)){
    if ( is.na(tmp1[[i]][1,1]) &  is.na(tmp1[[i]][2,1]) & is.na(tmp1[[i]][3,1]) & is.na(tmp1[[i]][4,1]) & is.na(tmp1[[i]][5,1])){
      if( !is.na(tmp1[[i]][30,1]) & substr(tmp1[[i]][30,5], nchar(tmp1[[i]][30,5])-1, nchar(tmp1[[i]][30,5])) == 'GT' ){
        tmp.ind = c( tmp.ind,i)
        substr(tmp[i,'mutation.correction'], 3 , (as.numeric(tmp1[[i]][30,3]) +2)) = paste(rep('*',as.numeric(tmp1[[i]][30,3])), collapse = "")
        substr(tmp[i,'align.correction'], 3 , (as.numeric(tmp1[[i]][30,3]) +2)) = paste(rep('*',as.numeric(tmp1[[i]][30,3])), collapse = "")
        
        tmp[i,'mutation.correction'] = gsub('\\*','', tmp[i,'mutation.correction'] )
        tmp[i,'align.correction'] = gsub('\\*','', tmp[i,'align.correction'] )
      }
    }
  }
  
  cell.matrix = left_join(cell.matrix, tmp[c(1,2,4,5)], by=c('read.short','align.short'))
  
  return(cell.matrix)
  
}





