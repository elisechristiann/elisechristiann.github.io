library(tidyverse)
library(ape)
library(phangorn)
library(ggtree)
library(seqinr)
library(msa)



dat <- read_csv("phi-base_current.csv")

View(dat)


##Cleaning the data
#removing the first row

#need to pull out the rip ones or whatever can make a tree
dat1 <- 
  dat[-1,] %>% 
  filter(Pathogen_species=="Ralstonia solanacearum") %>%
  filter(!is.na(GeneLocusID)) %>% filter(!is.na(Disease_name)) %>% filter(!is.na(Experimental_strain)) %>%
  filter(grepl("^rip", Gene_name, ignore.case = TRUE)) %>%
  select("ProteinID","GeneLocusID","Gene_name","Pathogen_species","Disease_name","Phenotype_of_mutant","Experimental_host_species","Host_descripton", "Experimental_strain") %>%
  mutate(Gene_name = as.character(gsub("Rip", "rip", Gene_name)))

write.csv(dat1, "C:/Users/elise/OneDrive/Desktop/elisechristiann.github.io/Project/raldat.csv", row.names = FALSE)

#making a list of acession numbers

#need to remove all duplicate rows

df <- read_csv('raldat.csv')

raldat_norepeats <- unique(df["GeneLocusID"])

#getting accession numbers into a single vector

ral_acession_numbers <- raldat_norepeats$GeneLocusID

#downloading DNAbin seq

ral_sequences <- read.GenBank(ral_acession_numbers,type="AA")

attributes(ral_sequences)

#giving them better names

name_data <- attr(ral_sequences, "description")

# applying the names
new_raldatseq <- paste(name_data, names(ral_sequences), sep = "_")


#aligning it 

#  temporary FASTA file
fasta_file <- tempfile(fileext = ".fasta")
write.dna(new_raldatseq, file = fasta_file, format = "fasta")

# multiple sequence alignment 
alignment <- msaMuscle(fasta_file, type="protein")  # or msaClustalW(), msaClustalOmega()

# Convert the alignment to ape format
alignment_ape <- msaConvert(alignment, type = "ape::AAbin")  

phydat_alignment <- phyDat(alignment_ape, type = "AA")

# Distance matrix
dist_matrix <- dist.ml(phydat_alignment)

# Try NJ tree first
start_tree <- tryCatch(NJ(dist_matrix), error = function(e) NULL)

# If NJ fails, use UPGMA
if (is.null(start_tree)) {
  start_tree <- upgma(dist_matrix)
}

# Fix tip labels
start_tree$tip.label <- names(phydat_alignment)

# Build model
fit <- pml(start_tree, data = phydat_alignment)

# Optimize
fit_opt <- optim.pml(fit, model = "JTT", optInv = TRUE, optGamma = TRUE, rearrangement = "stochastic")

# Plot final tree
fit %>%
plot(tree, main = "Maximum Likelihood Phylogenetic Tree")





