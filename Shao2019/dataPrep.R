# prepare data

library(tidyverse)

df = read.csv("./Data/count-table.csv") %>% as_tibble()
taxa = read.csv("./Data/speciesMetadata.csv") %>% as_tibble()
sampleMeta = read.csv("./Data/sampleMetadata.csv", skip=2) %>% as_tibble()
mapping = read.csv("./Data/mapping.tsv", sep="\t") %>% as_tibble()

temp= sampleMeta %>% left_join(mapping, by=c("Accession"="secondary_sample_accession"))
colnames(df) = str_split_fixed(colnames(df),"_",3)[,1]
order = temp %>% filter(run_accession %in% colnames(df)) %>% select(run_accession)

df = df %>% select(-X)
df = df %>% select(order %>% pull())

write.table(df %>% t() %>% as_tibble(), "./Data/count-table_clean.csv", row.names=FALSE, col.names=FALSE)
write.table(temp %>% filter(run_accession %in% colnames(df)), "./Data/sampleMetadata_clean.csv", row.names=FALSE, col.names=FALSE)
write.table(taxa, "./Data/speciesMetadata_clean.csv", row.names=FALSE, col.names=FALSE)
