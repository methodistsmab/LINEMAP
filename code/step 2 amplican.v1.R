setwd('/Users/linwang/Work/dingcheng lineage tracing')
setwd('./amplican crispr analysis')

library(amplican)

#########################################

quality.sequence = 'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF'

read.construct = function(id, sequence){
		fastqread = matrix('',4*length(id),1)
	for (i in 1 : length(id)){
		fastqread[((i-1)*4+1),1] = id[i]
	    fastqread[((i-1)*4+2),1] = sequence[i]
	    fastqread[((i-1)*4+3),1] = '+'
	    fastqread[((i-1)*4+4),1] = quality.sequence
	}
	return(fastqread)
}



make.config = function(ID,Forward_Reads, Group,GuideRNA, Forward_Primer, Amplicon){
	config = data.frame(ID= '', Barcode = 'Barcode_1', Forward_Reads = '',
	Reverse_Reads = '', Group = '', Control = 0, guideRNA = '',
	Forward_Primer = '', Reverse_Primer = '', Direction = 0, Amplicon = '',Donor = '')
	config['ID'] = ID
	config['Forward_Reads'] = Forward_Reads
	config['Group'] = Group
	config['guideRNA'] = GuideRNA
	config['Forward_Primer'] = Forward_Primer
	config['Amplicon'] = Amplicon
	return(config)
	
}

ID = 'ss'
Forward_Reads = 'py1read.fastq.gz'
Group = c('PY','LU')
GuideRNA = c('GTTCCCGTCCAGTAATCGTG', 'GTCGTTGTAGCAACCTATCGGGTG')
Forward_Primer = c('AAGCAGTGGTATCAACGCAGAGTACATGGG','AAGCAGTGGTATCAACGCAGAGTACATGGG')
Amplicon = c('aagcagtggtatcaacgcagagtacatgggGTTCCCGTCCAGTAATCGTGgggttagagctagaaatagcaagttaacctaaggctagtc',
'aagcagtggtatcaacgcagagtacatgggGTCGTTGTAGCAACCTATCGGGTGgggttagagctagaaatagcaagttaacctaaggct')

make.config(ID,Forward_Reads, Group[1],GuideRNA[1], Forward_Primer[1], Amplicon[1])

mainDir = '/Users/linwang/Work/dingcheng lineage tracing/amplican crispr analysis/result all cell 1/aligment'
subDir <- "lu2"
dir.create(file.path(mainDir, subDir))
subDir1 <- 'B25'
dir.create(file.path(mainDir, subDir,subDir1))

setwd(file.path(mainDir, subDir,subDir1))

unlink("./alignments.txt",recursive = TRUE)
unlink("./events_filtered_shifted_normalized.txt",recursive = TRUE)
unlink("./events_filtered_shifted.txt",recursive = TRUE)
unlink("./raw_events.txt",recursive = TRUE)
unlink("./unassigned_reads.txt",recursive = TRUE)

cell.number = length(names.cell.barcode)

for(i in 1: cell.number){
	index.tmp = which(indexlist == i)
	all.sequence = read.string1[index.tmp]
	if (length(all.sequence) == 0){
		next
	}	
	if(!any(grepl(Forward_Primer[2],all.sequence,fixed = TRUE))){
	  next
	}
	
	unique.all.sequence = unique(all.sequence[which(grepl(Forward_Primer[2],all.sequence,fixed = TRUE))])
	if(length(unique.all.sequence) == 1){
	  if(unique.all.sequence == 'AAGCAGTGGTATCAACGCAGAGTACATGGGGTCGTTGTAGCAACCTATCGGGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCT'){
	    all.sequence = c(all.sequence,'AAGCAGTGGTATCAACGCAGAGTACATGGGGGTCGTTGTAGCAACCTATCGGGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGC')
	  }
	  if(unique.all.sequence == 'AAGCAGTGGTATCAACGCAGAGTACATGGGGTTCCCGTCCAGTAATCGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTC'){
	    all.sequence = c(all.sequence,'AAGCAGTGGTATCAACGCAGAGTACATGGGGGTTCCCGTCCAGTAATCGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGT')
	  }
	}
	
	d = length(all.sequence)
	id = rep('',d)
	for (j in 1 :d){
		id[j] = paste('@',names.cell.barcode[i],'_',j,sep = '')
	}
	
	tmp.read = read.construct(id,all.sequence)
	
	unlink("./tmp.read.fastq.gz",recursive = TRUE)
	write.table(tmp.read,file = './tmp.read.fastq.gz',quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

    ID = names.cell.barcode[i]
    Forward_Reads = 'tmp.read.fastq.gz'
    tmp.config <- make.config(ID,Forward_Reads, Group[1],GuideRNA[2], Forward_Primer[2], Amplicon[2])
    unlink("./tmp.config.csv",recursive = TRUE)
    write.csv(tmp.config,file = 'tmp.config.csv',quote = FALSE,row.names = FALSE)

	unlink("./results_folder_tmp",recursive = TRUE)
    dir.create("./results_folder_tmp") 
    config <- './tmp.config.csv'
    fastq_folder <- './'
    results_folder <- './results_folder_tmp'
    amplicanPipeline(config, fastq_folder, results_folder, fastqfiles=1,knit_reports = FALSE)
    
    alignments = read.csv('./results_folder_tmp/alignments/alignments.txt', header = FALSE)
    write.table(alignments,file = 'alignments.txt', quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
    
    events_filtered_shifted_normalized = read.csv('./results_folder_tmp/alignments/events_filtered_shifted_normalized.csv')
    write.table(events_filtered_shifted_normalized,file = 'events_filtered_shifted_normalized.txt', quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
    
     events_filtered_shifted = read.csv('./results_folder_tmp/alignments/events_filtered_shifted.csv')
    write.table(events_filtered_shifted,file = 'events_filtered_shifted.txt', quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
    
     raw_events = read.csv('./results_folder_tmp/alignments/raw_events.csv')
     write.table(raw_events,file = 'raw_events.txt', quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
     
     if(file.exists("./results_folder_tmp/alignments/unassigned_reads.csv")){
       unassigned_reads = read.csv('./results_folder_tmp/alignments/unassigned_reads.csv')
       write.table(unassigned_reads,file = 'unassigned_reads.txt', quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
     }

       print(i)

}


