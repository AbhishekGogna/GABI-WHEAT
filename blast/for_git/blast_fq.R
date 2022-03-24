# Packages ----------------------------------------------------------------
from_cran <- c("tidyverse",  # 1.3.1
              "readxl", # 1.3.1
              "hablar", # 0.3.0
              "stringr", # 1.4.0
              "data.table", # 1.14.0
              "future", # 1.22.1
              "future.apply", # 1.8.1
              "qs") # 0.25.1 

for (i in from_cran){
  if (!requireNamespace(i, quietly = TRUE)) 
    install.packages(i)}

if(sum(unlist(lapply(from_cran, require, character.only = TRUE))) == length(from_cran)) {
  print("All required packages are present and are loaded. Version check was not done.")
} else {
  print("Some packages were not loaded or are not installed. Please install and load packages manually by using devtools::install_version()")
}


# Folder structure --------------------------------------------------------


# functions ---------------------------------------------------------------
`%!in%` = Negate(`%in%`) # negates %in%

query <- function(x) { # to process raw xlsx data sequences
  a <- gsub("(\\S+)?\\[(\\S)\\/(\\S)\\](\\S+)?", "\\1", x, perl = T) # sequence before the SNP
  b <- gsub("(\\S+)?\\[(\\S)\\/(\\S)\\](\\S+)?", paste0("\\", sample(2:3, 1, prob = c(0.5, 0.5))), x, perl = T) # sampling SNPs since only one can go in at this position in th fasta file
  c <- gsub("(\\S+)?\\[(\\S)\\/(\\S)\\](\\S+)?", "\\4", x, perl = T) # sequence after the SNP
  seq <- str_c(a,b,c)
  return(seq)
}

blast_search <- function(raw_data, 
                         index, 
                         db = "DNA", 
                         markers = NULL,
                         path = "/proj/GABI/blast/for_git/")
{
  #Filter data
  
  if (is.null(markers) == T){
    data <- raw_data %>% filter(chip == index)
  } else if(is.null(markers) == F){
    data <- raw_data %>% filter(chip == index, Locus_Name %in% markers)
  }
  
  ##Write fasta
  system.time(for (i in 1:nrow(data)){
    
    if(i == 1) {
      cat(paste0(">", data[i, "Locus_Name"], "\n"), file = paste0(path, index, "_blast" ,"_", db,".fa"))
    } else {
      cat(paste0(">", data[i, "Locus_Name"], "\n"), file = paste0(path, index, "_blast" ,"_", db,".fa"), append = T)
    }
    
    if(i < nrow(data)){
      cat(paste0(data[i, "Sequence"], "\n"), file = paste0(path, index, "_blast" ,"_", db,".fa"), append = T)
    } else if (i == nrow(data)){
      cat(paste0(data[i, "Sequence"], "\n"), file = paste0(path, index, "_blast" ,"_", db,".fa"), append = T)
    }
  })
  message("Fasta file written!")
  
  # Run blast
  
  # I put the code out since in my automated pipeline blast happens on its own. But to make things simple i modified the code so that this function only writes the fasta files.
  # Please set up blast in bash 
  # guide at - https://www.ncbi.nlm.nih.gov/books/NBK52640/ # i used https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.12.0/
  # download databases. This is chinese spring v1. But i guess you could also download latest ones.
  #makeblastdb -in Triticum_aestivum.IWGSC.cdna.all.fa -dbtype nucl -out ../databases/wheat_e51_cDNA/wheat_e51_cDNA # http://ftp.ensemblgenomes.org/pub/plants/release-51/fasta/triticum_aestivum/cdna/Triticum_aestivum.IWGSC.cdna.all.fa.gz
  #makeblastdb -in Triticum_aestivum.IWGSC.dna.toplevel.fa -dbtype nucl -out ../databases/wheat_e51_DNA/wheat_e51.blastdb # http://ftp.ensemblgenomes.org/pub/plants/release-51/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz
  
  message("Please use blast_custom.sh to run a blast search on the fasta file.")
}

read_res <- function(index, db, path = "/proj/GABI/blast/results/"){
  input_res <- fread(paste0(path, index,"_", db,"_blast_out"),
                     col.names = c("query_id", "query_length", "q_start", "q_end", 
                                   "subject_id", "subject_length", "s_start", "s_end", 
                                   "alignment_length", "pcent_identity", "mismatches", 
                                   "gaps", "evalue", "bit_score", "subject_strand"))
  input_res <- setDF(input_res)
  return(input_res)
}

filtering_data_dna <- function(data){
  
  out <- list()
  
  out[["alignment_test_yes"]] <- data %>% filter(alignment_length >= 0.9*query_length)
  
  out[["alignment_test_no"]] <- data %>% filter(alignment_length < 0.9*query_length)
  
  data_s1 <- data %>% filter(alignment_length >= 0.9*query_length)
  evals <- data_s1 %>% pull(evalue)
  mismatch <- data_s1 %>% pull(mismatches)
  pcent_identity <- data_s1 %>% pull(pcent_identity)
  
  if (nrow(data_s1)==1){
    out[["uniq_map"]] <- data_s1
  } else if (length(unique(evals)) >1){
    out[["uniq_map"]] <- data_s1 %>%
      slice(which.min(evalue))
  } else if (length(unique(mismatch)) >1){
    mis_match <- data_s1 %>% filter(mismatch < 2) %>%
      slice(which.min(mismatch))
    if (nrow(mis_match) != 0){
      out[["uniq_map"]] <- mis_match} 
    else {
      out[["uniq_map"]] <- NA
    }
  } else if (length(unique(pcent_identity)) >1){
    out[["uniq_map"]] <- data_s1 %>%
      slice(which.max(pcent_identity))
  } else {
    out[["uniq_map"]] <- NA
  }
  
  return(out)
}

filtering_data_cdna <- function(data){
  
  out <- list()
  
  out[["alignment_test_yes"]] <- data %>% filter(alignment_length >= 0.9*query_length)
  
  out[["alignment_test_no"]] <- data %>% filter(alignment_length < 0.9*query_length)
  
  data_s1 <- data %>% filter(alignment_length >= 0.9*query_length)
  evals <- data_s1 %>% pull(evalue)
  mismatch <- data_s1 %>% pull(mismatches)
  pcent_identity <- data_s1 %>% pull(pcent_identity)
  gene <- data_s1 %>% 
    mutate(subject_id = gsub("(\\S+)\\.\\d", "\\1", subject_id, perl = T)) %>%
    pull(subject_id)
  gene_variants <- data_s1 %>% pull(subject_id)# This helps ignore alternative splicing situation
  
  if (nrow(data_s1)==1){
    out[["uniq_map"]] <- data_s1 %>% mutate(variants = NA)
  } else if (length(unique(evals)) >1){
    out[["uniq_map"]] <- data_s1 %>%
      slice(which.min(evalue)) %>% mutate(variants = NA)
  } else if (length(unique(mismatch)) >1){
    mis_match <- data_s1 %>% filter(mismatch < 2) %>%
      slice(which.min(mismatch)) %>% mutate(variants = NA)
    if (nrow(mis_match) != 0){
      out[["uniq_map"]] <- mis_match} 
    else {
      out[["uniq_map"]] <- NA 
    }
  } else if (length(unique(pcent_identity)) >1){
    out[["uniq_map"]] <- data_s1 %>%
      slice(which.max(pcent_identity)) %>% mutate(variants = NA)
  } else if (length(unique(gene)) == 1){
    out[["uniq_map"]] <- data_s1[1,] %>% 
      mutate(subject_id = gsub("(\\S+)\\.\\d", "\\1", subject_id, perl = T)) %>%
      mutate(variants = length(unique(gene_variants)))
  } else {
    out[["uniq_map"]] <- NA
  }
  
  return(out)
}

process_data <- function(index, db){
  data <- read_res(index, db)
  
  message("File read")
  
  hits_per_query <- as_tibble(data) %>% count(query_id) %>% arrange(desc(n)) %>% pull(n)
  
  plan(cluster, workers = 20)
  
  split_data <- data %>%
    mutate(marker = query_id) %>%
    group_by(query_id) %>%
    group_split()
  names(split_data) <- data %>% count(query_id) %>% pull(query_id)
  
  message("Filtering started")
  
  if (db == "DNA"){
    system.time(check <- future_lapply(split_data, filtering_data_dna)) # around 3 minutes
    check_2 <- lapply(check, function(x) {
      test <- is.na(x[["uniq_map"]])
      if (test[1] == T){
        out <- T
      } else {
        out <- F
      }
      return(out)
    })
    
    no_uniq_pos <- check[which(check_2 == T)]
    
    uniq_pos <- do.call(rbind, lapply(check[which(check_2 == F)], function(x) x[["uniq_map"]])) %>%
      select(-marker) %>% 
      left_join(data_probe_raw %>% filter(chip == index), by = c("query_id" = "Locus_Name")) %>%
      filter(length_pre > 10 & length_post > 10) %>%
      mutate(pos = ifelse(subject_strand == "plus", (s_start+length_pre), (s_end+length_post))) %>%
      select(query_id, subject_id, pos, subject_strand)
  } else if (db == "cDNA"){
    system.time(check <- future_lapply(split_data, filtering_data_cdna)) 
    
    check_2 <- lapply(check, function(x) {
      test <- is.na(x[["uniq_map"]])
      if (test[1] == T){
        out <- T
      } else {
        out <- F
      }
      return(out)
    })
    
    no_uniq_pos <- check[which(check_2 == T)]
    uniq_pos <- do.call(rbind, lapply(check[which(check_2 == F)], function(x) x[["uniq_map"]])) %>%
      select(query_id, subject_id, subject_length, s_start, s_end, subject_strand, variants)
  }
  
  message("Filtering done")
  
  #Generate output
  output <- list()
  output[["no_uniq_pos"]] <- no_uniq_pos
  output[["uniq_pos"]] <- uniq_pos
  return(output)
}

# Input data --------------------------------------------------------------

data_probe <- data_probe_raw <- qread("/proj/GABI/blast/for_git/chip_90k_raw_data.qs")

data_probe$Sequence <- query(data_probe$Sequence)

# Process data ------------------------------------------------------------

system.time(res_90k <- blast_search(data_probe, 
                                    index = "90k"))

data_90k <- process_data("90k", "DNA")

# Process data - 2; Check markers with cDNA data --------------------------

system.time(res_90k_c <- blast_search(data_probe, 
                                      index = "90k", 
                                      db = "cDNA", 
                                      markers = names(data_90k$no_uniq_pos)))

data_90k_cDNA <- process_data("90k", "cDNA")