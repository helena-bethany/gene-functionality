working_d <- getwd()
setwd(working_d)
library(PopGenome)

genome.class <- readData("FASTA", format="VCF",gffpath=FALSE)
genome.class <- neutrality.stats(genome.class)
dataset <- get.neutrality(genome.class)[[1]]
write.csv(dataset, "name.csv")     
