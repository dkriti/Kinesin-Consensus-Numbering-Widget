library(bio3d)
require(dplyr)
library(splitstackshape)
library(DT)
library(shiny)

load("D:/raw/kinesin_subspace_3Apr12.RData")

#file_fastaa = "D:/aln_combine_Trimed_ID10.fa"
file_fastaa = "D:/final_alignment-h4-5_tempus2_barry_version.fa"
seqq <- read.fasta(file_fastaa, rm.dup = TRUE, to.upper = FALSE, to.dash=TRUE)

common_ids <- pdbs$id

for(i in 1:8){common_ids <- common_ids[-115]}
#remove this later
for(i in 1:3){common_ids <- common_ids[-15]}
common_ids <- basename.pdb(common_ids)
common_ids <- as.character(common_ids)
k <- length(common_ids)

aln <- seqq$ali
seqq_ids <- as.character(basename.pdb(seqq$id))
common_ids %in% seqq_ids
aln_new <- cbind(seqq_ids, seqq$ali)
pdb <- pdbs$ali
pdb_id <- basename.pdb(pdbs$id)
pdb_id <- as.character(pdb_id)
pdb_id <- as.matrix(pdb_id)
pdb_new <- cbind(pdb_id, pdbs$ali)

sse_tab_all1 <- NULL
sse_tab_all <- NULL

for(jh in 1:length(common_ids)){
  for(v in 1:nrow(aln_new)){
    if(common_ids[jh] == aln_new[v]){
      aln_seq <- aln[v,]
    }
  }
  for(v in 1:nrow(pdb_new)){
    if(common_ids[jh] == pdb_new[v]){
      pdb_seq1 <- pdbs$ali[v,]
      pdb_seq <- pdb_id[v]
      
    }
  }
  
  pdb_seq <- strsplit(pdb_seq, "_")
  pdb_seq <- lapply(pdb_seq,function(x){x[1]})
  pdb_seq <- unlist(pdb_seq)
  
  
  pdb <- read.pdb(pdb_seq)
  seq <- aln_seq
  
  sse <- dssp(pdb, exefile="D:/dssp.exe", resno=TRUE)
 
  res_no <- unbound(sse$helix$start, sse$helix$end)
  helix_table <- as.data.frame(res_no)
  helix_table <- cbind(RegionName = "H", helix_table)
  
  res_no <- unbound(sse$sheet$start, sse$sheet$end)
  sheet_table <- as.data.frame(res_no)
  sheet_table <- cbind(RegionName = "S", sheet_table)
  
  sse_table <- rbind(helix_table, sheet_table)
  
  helix_residues <- unbound(sse$helix$start, sse$helix$end)
  sheet_residues <- unbound(sse$sheet$start, sse$sheet$end)
  
  ca.inds <- atom.select(pdb, string="calpha")
  last.residue <- max(pdb$atom[ca.inds$atom,"resno"])
  first.residue <- min(pdb$atom[ca.inds$atom,"resno"]) 
  
  everything <- sort(unique(pdb$atom[,7]))
  everything_new <- subset(everything, everything <= last.residue)
  
  everything_table <- as.data.frame(everything_new)
  helix_res_rm_table <- as.data.frame(helix_residues)
  sheet_res_rm_table <- as.data.frame(sheet_residues)
  colnames(everything_table) <- c("v1")
  colnames(helix_res_rm_table) <- c("v1")
  colnames(sheet_res_rm_table) <- c("v1")
  everything_less_helix <- anti_join(everything_table,helix_res_rm_table)
  loop_residues <- anti_join(everything_less_helix, sheet_res_rm_table)
  res_no <- sort(loop_residues[,1])
  loop_table <- as.data.frame(res_no)
  loop_table <- cbind(RegionName = "L", loop_table)

  sse_table <- rbind(sse_table, loop_table)
  
  seq_tab <- as.data.frame(seq)
  seq_tab["Res_No"] <- NA
  Residues_numbers <- as.data.frame(pdb$atom[ca.inds$atom,"resno"])
  m <- nrow(seq_tab)
  n <- nrow(Residues_numbers)
  s <- 1
  for(t in 1:m){
    if(seq_tab[t,1]!="-"){
        seq_tab[t,2] <- Residues_numbers[s,1]
        s = s+1
    }
  }
  seq_tab["Aln_No"] <- rownames(seq_tab)
  
  sse_table$aln_no = seq_tab[match(sse_table$res_no, seq_tab$Res_No),"Aln_No"]
  sse_table$residue = seq_tab[match(sse_table$res_no, seq_tab$Res_No),"seq"]
  
  gaps_info <- gap.inspect(seq)
  gaps_index <- gaps_info$t.inds
  gaps_tabb <- as.data.frame(gaps_index)
  gaps_tabb$aln_prev <- gaps_tabb$gaps_index - 1
  gaps_tabb$aln_after <- gaps_tabb$gaps_index + 1

  sse_tab_new <- data.frame(index = seq(m))
  sse_tab_new1 <- sse_table[match(sse_tab_new$index, sse_table$aln_no),"res_no"]
  sse_tab_new <- sse_table[match(sse_tab_new$index, sse_table$aln_no),"RegionName"] 
  sse_tab_new1 <- t(as.matrix(sse_tab_new1))
  sse_tab_new <- t(as.matrix(sse_tab_new))
  sse_tab_all1 <- rbind(sse_tab_all1,sse_tab_new1)
  sse_tab_all <- rbind(sse_tab_all,sse_tab_new)
  
  
}

write.table(sse_tab_all, "D:/sse_percentages_again_guido.xls", row.names=T, col.names=T, quote=F, sep="\t")
write.table(sse_tab_all1, "D:/sse_residue_nos_guido.xls", row.names=T, col.names=T, quote=F, sep="\t")
final_tab_reg <- read.table("D:/sse_percentages_values.txt", header = FALSE, sep = "\t")

#file_fastaa = "D:/aln_combine_Trimed_ID10.fa"
file_fastaa = "D:/final_alignment-h4-5_tempus2_barry_version.fa"
full_alignment <- read.fasta(file_fastaa, rm.dup = TRUE, to.upper = FALSE, to.dash=TRUE)
full_aln <- as.data.frame(full_alignment$ali)
new.row <- final_tab_reg[1,]
last_tab <- rbind(new.row, full_aln)

write.table(last_tab, "D:/consensus_table.xls", row.names=T, col.names=T, quote=F, sep="\t")
consensus_table <- read.table("D:/consensus_table.xls", header = FALSE, sep = "\t")


#perform seq2aln in a for loop

fasta_to_add = "D:/aln_combine_Trimed_ID10_new.fa"
file_fasta_to_add <- read.fasta(fasta_to_add, rm.dup = TRUE, to.upper = FALSE, to.dash=TRUE)
fasta_file_new = "D:/final_alignment-h4-5_tempus2_barry_version.fa"
file_fasta_file_new <- read.fasta(fasta_file_new1, rm.dup = TRUE, to.upper = FALSE, to.dash=TRUE)

for(i in 5878:6001){
  tmp <- file_fasta_to_add$ali[i,]
  id_added <- file_fasta_to_add$id[i]
  remove <- which(tmp=="-")
  tmp <- tmp[-remove]
  fasta_file_new = paste("D:/new_alignment","_",i-1,".fa",sep="")
  file_fasta_file_new <- read.fasta(fasta_file_new, rm.dup = TRUE, to.upper = FALSE, to.dash=TRUE)
  seq2aln(tmp, file_fasta_file_new, id = id_added, file = paste("D:/new_alignment","_",i,".fa",sep=""))
}

final_aln <- "D:/new_alignment_with all kinesins.fa"
file_final_aln <- read.fasta(final_aln, rm.dup = TRUE, to.upper = FALSE, to.dash=TRUE)

#conservation scores
conservation <- conserv(x=file_fasta_file_new$ali, method="similarity", sub.matrix="bio3d")


