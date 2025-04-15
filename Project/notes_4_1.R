library(rentrez)
library(tidyverse)
library(ape)
library(phangorn)
library(ggtree)

vignette("Trees", package="phangorn")

#install ape
#install phytools
#phangorn
#install DECIPHER
#install ggtree


entrez_db_summary("protein")

dat <- read_csv("phi-base_current.csv")

View(dat)

#removing the first row

#need to pull out the rip ones or whatever can make a tree
dat1 <- 
  dat[-1,] %>% 
  filter(Pathogen_species=="Ralstonia solanacearum") %>%
  filter(!is.na(GeneLocusID)) %>% filter(!is.na(Disease_name)) %>% filter(!is.na(Experimental_strain)) %>%
  filter(grepl("^rip", Gene_name, ignore.case = TRUE)) %>%
  select("ProteinID","GeneLocusID","Gene_name","Pathogen_species","Disease_name","Phenotype_of_mutant","Experimental_host_species","Host_descripton", "Experimental_strain")
  
View(dat1)
#select ProtienID, GeneLocusID, Genomic Sequence Providing Strain, Gene_name, Disease_name, 
#raldat <- dat1 %>%
 # filter(Pathogen_species=="Ralstonia solanacearum") %>%
#  pivot_longer()
#View(raldat)
