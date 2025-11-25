library(ape)   
library(seqinr)
library(Biostrings)
library(phytools)
library(phangorn)

# To contruct a phylogenetic tree from aligned sequences 
make_tree <- function(input_file, output_file) {

  #input file: aligned sequences in fasta format
  #output file: .phy phylogenetic tree as newick string 
  
  aligned_protein_sequences <- read.alignment(file = input_file, format = "fasta")
  
  aligned_protein_phyDat <- as.phyDat(aligned_protein_sequences)
  
  dist_matrix_protein <- dist.alignment(aligned_protein_sequences, "identity")
  
  protein_tree <- nj(dist_matrix_protein)
  
  plot(protein_tree)
  
  newick_string <- write.tree(protein_tree, file = "")
  
  write(newick_string, file = output_file)
}
