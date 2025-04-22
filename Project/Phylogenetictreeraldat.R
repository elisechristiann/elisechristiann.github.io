library(tidyverse)
library(ape)
library(msa)
library(phangorn)
library(rentrez)
library(seqinr)
library(ggmsa)
library(Biostrings)
library(janitor)

#getting the sequences from NCBI

effect <- read_csv("effectors_entire_spreadsheet.csv")

#cleaning the data again
#now there are only proteins
clean_effect <- effect %>%
  select("#Species","accession","type","name","description")%>%
  filter(effect$accession !="NEW",) %>%
  filter(type== "protein")

#making a new csv
write.csv(clean_effect,"C:/Users/elise/OneDrive/Desktop/elisechristiann.github.io/Project/clean_effect.csv", row.names = FALSE)

#now to get the sequences again pls ape

ral_df <- read_csv("clean_effect.csv")

ral_dfaccession <- ral_df$accession

#getting sequences
# hello chat
# Your list of accession numbers

# List of IDs
ids <- ral_dfaccession

# Function to split into batches of 20
batch_size <- 20
id_batches <- split(ids, ceiling(seq_along(ids) / batch_size))

all_seqs <- list()

for (i in seq_along(id_batches)) {
  cat("Downloading batch", i, "\n")
  
  # Fetch the FASTA sequences
  fasta <- entrez_fetch(db="protein", id=id_batches[[i]], rettype="fasta", retmode="text")
  
  # Parse the FASTA into individual sequences
  temp <- tempfile()
  writeLines(fasta, temp)
  seqs <- read.fasta(temp, seqtype="AA", as.string=TRUE)
  
  all_seqs <- c(all_seqs, seqs)
}

# Save all sequences into a file
write.fasta(sequences = all_seqs,
            names = names(all_seqs),
            file.out = "ralproteins.fasta")

#now downloaded, can continue
#renaming so not accession numbers
sequences <- read.FASTA("ralproteins.fasta", type="AA")
#accession numbers as sequence names
names_df <- data.frame(accession = names(sequences))
#merge tables to find where differences and get rid of different ones
merged_df <- merge(names_df, ral_df, by = "accession", all.x = TRUE)
# new names
merged_df$new_name <- paste0(merged_df$description, "_", seq_len(nrow(merged_df)))
# renamed sequences
names(sequences) <- merged_df$new_name
#saved the renamed sequences
write.fasta(sequences, names = names(sequences), file.out = "ral_renamed.fasta")

#wow that's a big file, let's select a few Rip proteins to look into

#some instructions from google on how the heck to subset a fasta file
# 1. Read in your FASTA file
sequences <- readAAStringSet("ral_renamed.fasta")

# 2. Look at the names
names(sequences)

# 3. Subset sequences where the name contains a keyword
# Example: keep sequences where the name contains "RipU"
subset_sequences <- sequences[grep("RipU", names(sequences))]

# 4. Check
subset_sequences

# 5. (Optional) Save them to a new FASTA file
writeXStringSet(subset_sequences, filepath = "subset_RipU_sequences.fasta")
#visualizing
ggmsa(subset_sequences, start=50, end= 100, color = "Clustal", 
      font = "helvetical", char_width = 0.5 )


#alignment of RipU
alignment_ripU <- msa(subset_sequences, method="Muscle", type="protein" )

# Convert to AAStringSet
alignment_aa_ripU <- msaConvert(alignment_ripU, type = "ape::AAbin")
alignment_char_ripU <- as.character(alignment_aa)

# Now collapse each character vector into a single string
alignment_strings_ripU <- sapply(alignment_char_ripU, paste0, collapse = "")

# Now make an AAStringSet!
alignment_AAStringSet_ripU <- AAStringSet(alignment_strings_ripU)

# Now, convert manually to phyDat (for amino acids)
phang_align_ripU <- phyDat(alignment_aa, type = "AA")

dm <- dist.ml(phang_align_ripU)

#making a tree???
treeNJ <- NJ(dm)
fit <- pml(treeNJ, phang_align_ripU)
fit_optimized <- optim.pml(fit, model = "JTT")

plot(fit_optimized$tree, main = "Phylogenetic Tree")







