# R version 4.0.2

# Packages
from_cran <- c("tidyverse",  # v1.3.1 
               "asreml", # v4.1.0.110
               "hablar") # v0.3.0

for (i in from_cran){
  if (!requireNamespace(i, quietly = TRUE)) 
    install.packages(i)}

if(sum(unlist(lapply(from_cran, require, character.only = TRUE))) == length(from_cran)) {
  print("All required packages are present and are loaded. Version check was not done.")
} else {
    print("Some packages were not loaded or are not installed. Please install and load packages manually")
  }

# Functions
`%!in%` = Negate(`%in%`)

get_check_list <- function(data){
  split_by_rep <- data %>% 
    group_by(id, rep) %>% 
    group_split()
  
  miss_indices <- unlist(lapply(split_by_rep, function(x) sum(is.na(x$value))))
  
  check_list_miss <- cbind(data %>% 
                             group_by(id, rep) %>% 
                             summarise(n = n(), .groups = "drop"), 
                           missing = miss_indices)
  
  to_omit <- check_list_miss %>% 
    filter(missing >= 0.95*n) # where missing data is > 95 % of available data
  
  data <- data %>% 
    anti_join(to_omit[, c("id", "rep")], by = c("id", "rep"))
  
  split_by_id <- data %>% 
    group_by(id) %>% 
    group_split() 
  
  rep_indices <- unlist(lapply(lapply(split_by_id, 
                                      function(x) unique(x[,"rep"])), 
                               function(x) length(unlist(x))))
  
  check_list <- cbind(data %>% 
                        group_by(id) %>% 
                        summarise(n = n()), 
                      rep_indices = rep_indices)
  
  return(check_list)
}

BLUEs <- function(data, 
                  pheno, 
                  heritability = T, 
                  heritability_plot = T, 
                  repeatab = T, 
                  BLUE = T){
  
  if ((pheno %in% traits_trial_1)[1]){
    n_reps <- 2 # i put n_reps as fixed for calculation of repeatabilites / hritabilities. You could choose to put mean(check_list$rep_indices) in its place. 
    
    trait <- data %>% 
      convert(fct(Coding, block, rep, year, loc)) %>%
      mutate(id = paste0(loc, ".", year)) # Convert data columns to required data classes
    
    #trait %>% group_by(id) %>% summarise(n = n()) # shows values per environment
    
    #Generate information for repetitions
    
    check_list <- get_check_list(trait)
    
    env <- check_list$id[which(check_list$rep_indices==2)]
    
    env.un <- check_list$id[which(check_list$rep_indices==1)]
    
    # set asreml options
    asreml.options(maxit = 50,             
                   workspace = "128mb",       
                   pworkspace = "128mb",
                   trace=F,
                   extra = 10)
    
    #For replicated environment follows thoery from the original publication
    
    if (length(env) > 0){
      repeatability <- rep(NA, length(env))
      residual <- rep(NA,length(env))
      
      for (i in 1:length(env))
      {
        field.red <- trait[which(trait$id == env[i]),] %>% droplevels()
        
        env.1 <- asreml(fixed = value ~  1, 
                        random = ~ Coding + rep + block:rep,
                        data = field.red # from methods formula (1)
        )
        residual[i] <- summary(env.1)$varcomp["units!R", "component"]
        repeatability[i] <- summary(env.1)$varcomp["Coding","component"]/(summary(env.1)$varcomp["Coding","component"] + summary(env.1)$varcomp["units!R","component"]/n_reps)  # from methods formula (2)
        
        repeatability_df <- data.frame("Environment" = env, "Repeatibility" = repeatability)
        
        env.2 <- asreml(fixed = value ~  Coding, 
                        random = ~ rep + block:rep,
                        data = field.red  # from methods formula (1)
        )
        Geno <- predict.asreml(env.2, classify = "Coding",)$pvals[,1:2]
        colnames(Geno)[2] <- env[i]
        if (i == 1) {BLUES <- Geno} else {BLUES <- merge(BLUES,Geno)} 
      }
      
      BLUES[,1] <- as.character(BLUES[,1])
      
      BLUEs_rep_wide <- BLUES
      
      BLUEs_rep <-  pivot_longer(BLUEs_rep_wide, !Coding, names_to = "env", values_to = "trait") %>% convert(fct(env, Coding))
    }
    
    # for unreplicated environments
    if (pheno %!in% c("ZEL", "SDS")){
      if (length(env.un) > 0){
        for (i in 1:length(env.un))
        {
          field.red <- trait[which(trait$id == env.un[i]),] %>% droplevels()
          
          #env.1 <- asreml(fixed = value ~  1, 
          #                random = ~ Coding + block,
          #                data = field.red
          #) # not needed hence put out
          
          env.2 <- asreml(fixed = value ~  Coding, 
                          random = ~ block,
                          data = field.red
          )
          Geno <- predict.asreml(env.2, classify = "Coding",)$pvals[,1:2]
          colnames(Geno)[2] <- env.un[i]
          if (i == 1) {BLUES <- Geno} else {BLUES <- merge(BLUES,Geno)} 
        }
        
        BLUES[,1] <- as.character(BLUES[,1])
        
        BLUEs_unrep_wide <- BLUES
        
        BLUEs_unrep <-  pivot_longer(BLUES, !Coding, names_to = "env", values_to = "trait") %>% 
          convert(fct(env, Coding))
      }
    } else {
      if (length(env.un) > 0){
        BLUEs_unrep <- NULL
        for (i in env.un){
          unrep <- (tapply(unlist(trait[which(trait$id == i)[1:400], "value"]), unlist(trait[which(trait$id == i)[1:400], "Coding"]), mean)) # 1:400 are entries in a given replication
          unrep <- data.frame(names(unrep),as.vector(unrep))
          unrep$env <- i
          colnames(unrep)[1:3]<-c("Coding", "trait" ,"env")
          BLUEs_unrep <- rbind(BLUEs_unrep, unrep[,c("Coding", "env", "trait")])
        }
        BLUEs_unrep_wide <- pivot_wider(BLUEs_unrep, names_from= env, values_from = trait)
        #if (length(env.un) == 1){colnames(BLUEs_unrep_wide <- c("Coding", env.un))}
      }
    }
    # merge the two datasets
    
    if (length(env.un)>0){
      if (length(env)>0){
        com_dta <- rbind(BLUEs_rep, BLUEs_unrep)}else{
          com_dta <- BLUEs_unrep %>% convert(fct(Coding, env))
        }
    }else(com_dta <- BLUEs_rep)
    
    # across env. BLUEs
    
    com_dta <- com_dta %>% arrange(env, Coding)
    
    env.1 <- asreml(fixed = trait ~ 1 , 
                    random = ~ Coding + env,
                    data = com_dta
    )
    
    env.2 <- asreml(fixed = trait ~ Coding , 
                    random = ~ env,
                    data = com_dta
    )
    
    B.acr <- predict.asreml(env.2, classify = "Coding",)$pvals[, c(1, 2)]
    colnames(B.acr) <- c("Coding", "Blues.acr")
    
    # wide form blues
    
    if (length(env.un)>0){
      if (length(env)>0){
        message(paste0("There are unreplicated entries in ", pheno, " data!"))
        com <- merge(BLUEs_rep_wide, BLUEs_unrep_wide)
        com_BLUEs <- merge(com, B.acr)
      }else{
        com_BLUEs <- merge(BLUEs_unrep_wide, B.acr)
      }
    }else{
      com_BLUEs <- merge(BLUEs_rep_wide, B.acr)
    }
    
  } else if ((pheno %in% traits_trial_2)[1]){
    # put different from traits_trial_1 since these do not have block structure 
    
    n_reps <- 3 # i put n_reps as fixed for calculation of repeatabilites / hritabilities. You could choose to put mean(check_list$rep_indices) in its place. 
    
    trait <- data %>% 
      convert(fct(Coding,rep, year, loc)) %>%
      mutate(id = paste0(loc, ".", year))
    
    #trait %>% group_by(id) %>% summarise(n = n()) # shows values per environment
    
    #Generate information for repetitions
    
    check_list <- get_check_list(trait)
    
    env <- check_list$id[which(check_list$rep_indices > 2)]
    env.un <- check_list$id[which(check_list$rep_indices==1)]
    
    repeatability <- rep(NA, length(env))
    residual <- rep(NA,length(env))
    
    #BLUEs within env
    
    for (i in 1:length(env)){
      env_n <- env[i] %>% as.character()
      field.red <- trait %>% filter(id == env_n) %>% droplevels()
      
      env.1 <- asreml(fixed = value ~  1, 
                      random = ~ Coding + rep,
                      data = field.red  # from methods formula (6)
      )
      
      residual[i] <- summary(env.1)$varcomp["units!R", "component"]
      repeatability[i] <- summary(env.1)$varcomp["Coding","component"]/(summary(env.1)$varcomp["Coding","component"] + summary(env.1)$varcomp["units!R","component"]/n_reps)
      
      repeatability_df <- data.frame("Environment" = env, "Repeatibility" = repeatability)
      
      env.2 <- asreml(fixed = value ~  Coding, 
                      random = ~ rep,
                      data = field.red
      )  # from methods formula (7)
      
      Geno <- predict.asreml(env.2, classify = "Coding",)$pvals[,1:2]
      
      colnames(Geno)[2] <- paste0(as.character(env[i]), ".", pheno)
      
      if (i == 1) {BLUES <- Geno} else {BLUES <- merge(BLUES,Geno)} 
      
    }
    
    BLUES[,1] <- as.character(BLUES[,1])
    
    dta <-  pivot_longer(BLUES, !Coding, names_to = "env", values_to = "trait") %>% 
      convert(fct(env, Coding)) %>%
      arrange(env, Coding)
    
    #BLUEs across env
    
    env.1 <- asreml(fixed = trait ~ 1 , 
                    random = ~ Coding + env,
                    data = dta
    )  # from methods formula (3)
    
    env.2 <- asreml(fixed = trait ~  Coding, 
                    random = ~ env,
                    data = dta)   # from methods formula (3)
    Blues.acr <- predict.asreml(env.2, classify = "Coding",)$pvals[,2]
    
    # wide form blues
    
    com_BLUEs <- cbind(BLUES, Blues.acr)
  }
  
  # Calculating heritabilities
  
  sigma.g <- summary(env.1)$varcomp["Coding","component"]
  
  if (length(env.un)>0){
    if (length(env)>0){
      sigma.g.e <- summary(env.1)$varcomp["units!R","component"]-(mean(residual)/n_reps)
      sigma.e <- mean(residual)
    }else{
      sigma.g.e <- NA
      sigma.e <- summary(env.1)$varcomp["units!R","component"]
    }
  }else{sigma.g.e <- summary(env.1)$varcomp["units!R","component"]-(mean(residual)/n_reps)
  sigma.e <- mean(residual)
  }  
  
  if(is.na(sigma.g.e)){
    
    herit <-  sigma.g/(sigma.g + sigma.e/length(env.un))  
    message(paste0("Since there are no replicates the term for sigma.g.e was omitted"))
    message( paste0("heritability_entry_mean_based = ", herit))
    
    herit_plot<- sigma.g/(sigma.g + sigma.e) 
    message(paste0("heritability_plot_based = ", herit_plot))
    
  }else{

    herit <- sigma.g/(sigma.g + sigma.g.e/length(env) + sigma.e/(n_reps*length(env)))    # from methods formula (5)
    message(paste0("heritability_entry_mean_based = ", herit))
    
    herit_plot<- sigma.g/(sigma.g + sigma.g.e + sigma.e)   # from methods formula (4)
    message(paste0("heritability_plot_based = ", herit_plot))
  }
  
  # Producing outputs
  
  output <- list()
  
  if (repeatab == T) {
    if (length(env) > 0){
      output["Repeatability"] <- list(repeatability_df)
    }else {
      repeatability <- NA
      output["Repeatability"] <- repeatability
    }
  } 
  
  if (heritability == T) {
    output["Heritabiliy_entry"] <- herit
  }
  
  if (heritability_plot == T) {
    output["Heritabiliy_plot"] <- herit_plot
  }
  
  if (BLUE == T) {
    output["BLUEs"] <- list(com_BLUEs)
  }
  
  return(output) 
}

process_isa_files <- function(study_file_path, assay_file_path){
  # Read files
  study <- read.delim(study_file_path) 
  assay <- read.delim(assay_file_path)
  
  # Process files
  assay_ana <- assay %>% select("Sample.Name", grep("Parameter.Value", colnames(assay), value = T))
  
  if(length(grep("block", colnames(study), ignore.case = T, value = T)) !=0){
    study_cols <- c("Source.Name", "Characteristics.Original.Coding.", 
                    "Sample.Name", "Factor.Value.year.", "Factor.Value.loc.", 
                    "Factor.Value.rep.", "Factor.Value.block.")
  } else if (length(grep("block", colnames(study), ignore.case = T, value = T)) ==0){
    study_cols <- c("Source.Name", "Characteristics.Original.Coding.", 
                    "Sample.Name", "Factor.Value.year.", "Factor.Value.loc.", 
                    "Factor.Value.rep.")
  }
  
  phenodata <- study %>% select(all_of(study_cols)) %>%
    left_join(assay_ana, by = "Sample.Name") %>%
    mutate(original_index = as.numeric(gsub("\\S+\\_(\\d+)", "\\1", Sample.Name, perl = T))) %>%
    arrange(original_index) %>% select(-Sample.Name, -original_index, -Source.Name)
  
  colnames(phenodata) <- gsub("\\S+\\.\\S+\\.(\\S+)\\.", "\\1", colnames(phenodata), perl = T)
  
  # generate output
  
  return(phenodata)
}

plot_herit <- function(list) {
  output <- list()
  for (i in names(list)){
    output[[i]] <- list[[i]]["Heritabiliy_plot"]
  }
  output_df <- data.frame("Trait" = gsub("(\\S+)\\.\\S+\\_.*", "\\1"
                                         , names(unlist(output)), perl = T)
                          , "Hertitability_plot" = unname(unlist(output))) 
  return(output_df)
}

entry_herit <- function(list) {
  output <- list()
  for (i in names(list)){
    output[[i]] <- list[[i]]["Heritabiliy_entry"]
  }
  output_df <- data.frame("Trait" = gsub("(\\S+)\\.\\S+\\_.*", "\\1"
                                         , names(unlist(output)), perl = T)
                          , "Hertitability_entry" = unname(unlist(output))) %>%
    convert(fct(Trait))
  return(output_df)
}

# Load data for two trials

#the files may be downloaded at https://doi.ipk-gatersleben.de/DOI/6c9a091d-2de4-4db5-a862-393a98183691/52f2f989-a332-4aa9-a461-168833766757/2/1847940088


my_data <- list() #  creates empty list to store loaded data

my_data[["trial_1"]] <- process_isa_files(study_file_path = "/path to files/isa_files/trial_1/s_alpha_lattice.txt", 
                                          assay_file_path = "/path to files/isa_files/trial_1/a_alpha_lattice.txt") # loads data from trial 1

my_data[["trial_2"]] <- process_isa_files(study_file_path = "/path to files/isa_files/trial_2/s_RCBD.txt", 
                                          assay_file_path = "/path to files/isa_files/trial_2/a_RCBD.txt") # loads data from trial 2

traits_trial_1 <- c("HD", "PH", "TKW" , "EW", "GPE", "GY", "SW", "GH", "STC", "PC", "SDS", "HAG", "ZEL") 
traits_trial_2 <- c( "FHB", "DTR", "SEP") 

# Calculate BLUES

store_all_output <- list()

for (i in traits_trial_1){
  
  trait <- i

  print(i)

  raw_data <- my_data$trial_1[, c("year", "loc", "rep", "block", "Coding", trait)]

  colnames(raw_data)[length(raw_data)] <- "value"

  store_all_output[[i]] <- BLUEs(raw_data, pheno = trait)
}

for (i in traits_trial_2){
  trait <- i
  
  print(i)
  
  raw_data <- my_data$trial_2[, c("year", "loc", "rep", "Coding", trait)]
  
  colnames(raw_data)[length(raw_data)] <- "value"
  
  store_all_output[[i]] <- BLUEs(raw_data, pheno = trait)
}

# check heritablilities

herit <- plot_herit(store_all_output) %>% 
  left_join(entry_herit(store_all_output), by = "Trait")

# values i produced

#Trait Hertitability_plot Hertitability_entry
#HD          0.8495462           0.9832407
#PH          0.8663641           0.9858973
#TKW          0.7022659           0.9445112
#EW          0.3509120           0.6540255
#GPE          0.3005748           0.4958687
#GY          0.4367285           0.8850901
#SW          0.6915210           0.9458670
#GH          0.7488045           0.9226237
#STC          0.6758745           0.8839153
#PC          0.5156306           0.8325704
#SDS          0.7358135           0.9176333
#HAG          0.3475953           0.7478804
#ZEL          0.7735023           0.9483337
#FHB          0.5584156           0.9007160
#DTR          0.1106081           0.2924300
#SEP          0.4366014           0.7284495


