
args = commandArgs(trailingOnly = TRUE)
antisense.transcripts = args[1]
annotation = args[2]

antisense.transcripts <- read.table(antisense.transcripts, header = FALSE, quote = "#", sep = "\t")
annotation <- read.table(annotation, header = FALSE, quote = "#", sep = "\t")

print(head(antisense.transcripts))
print(head(annotation))
nrow(annotation)
nrow(antisense.transcripts)
i <- 2

retain.transcript <- vector(mode = "logical",length = nrow(antisense.transcripts))

for(i in 1:nrow(antisense.transcripts))
{
  print(i)
  current.chr <- antisense.transcripts[i,1]
  current.start <- antisense.transcripts[i,4]
  current.end <- antisense.transcripts[i,5]
  
  overlap.genes <- subset(annotation, 
                          V1 == current.chr & ( 
                            V4 >= current.start & V4 <= current.end | 
                              V5 >= current.start & V5 <= current.end))[,c(1,4,5)]
  
  overlap.genes <- overlap.genes[!duplicated(overlap.genes),]
  retain.transcript[i] <- nrow(overlap.genes) == 1
}

head(retain.transcript)
sum(retain.transcript)

sum(retain.transcript)/7668

filtered.antisense.transcripts <- antisense.transcripts[retain.transcript,]
nrow(filtered.antisense.transcripts)

write.table(x = filtered.antisense.transcripts,file = "annot_merged_filtered.gtf",quote = F,sep = "\t",row.names = F,col.names = F)
