##Double cut table maker##

#Make the post-annotation DC master table


library(stringr)
library(plyr)
A_master <- NULL
files = list.files(pattern="_DC_annotated_table.txt") 
for (j in 1:length(files)) {
  data <- read.table(files[j], sep = "\t", header=TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, fill = TRUE)
  #data <-data[,-64]
  A_master <- rbind.fill(A_master, data)
}

write.table(A_master, file="DC_master_table.txt", col.names =TRUE, row.names = FALSE, quote = FALSE, sep="\t")
