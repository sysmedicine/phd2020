# Load Library
library('readr')

#load File
temp=read_delim('../data/deseq.txt', delim="\t",col_names = T)

# head the dataframe
head(temp)

# draw histogram of pvalue
hist(temp$pvalue)
