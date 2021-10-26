#Dummy R script
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

model <- function(fix_term, rand_term, data){
  asreml.options(maxit = 50,             
                 workspace = "128mb",       
                 pworkspace = "128mb",
                 trace=F)
  
  fit_model <- asreml(fixed = as.formula(fix_term),
                      random = as.formula(paste0("~ ", 
                                                 paste(rand_term, 
                                                       collapse = " + ")
                                                 )
                                          ),
                      data = data)
  
  return(fit_model)
}

get_variance_components <- function(asreml_obj, rand_term){
  out <- list()
  for (i in rand_term){
    out[[i]] <- summary(asreml_obj)$varcomp[i , "component"]
  }
  out[["units!R"]] <- summary(asreml_obj)$varcomp["units!R" , "component"]
  return(out)
}

BLUEs <- function(data, 
                  pheno, 
                  heritability = T, 
                  heritability_plot = T, 
                  repeatab = T, 
                  BLUE = T){
  library(tidyverse)
  library(asreml)
  library(hablar)
  `%!in%` = Negate(`%in%`)
  
  if ((pheno %in% traits_trial_1)[1])
    {
  trait <- data %>% 
    convert(fct(Coding, block, rep, year, loc)) %>%
    mutate(id = paste0(loc, ".", year))
  
  #trait %>% group_by(id) %>% summarise(n = n()) # shows values per environment
  
  #Generate information for repetitions
  
  check_list <- get_check_list(trait)
  
  env <- check_list$id[which(check_list$rep_indices==2)]
  
  env.un <- check_list$id[which(check_list$rep_indices==1)]
  
  #For replicated environment
  
  if (length(env) > 0){
    repeatability <- rep(NA, length(env))
    residual <- rep(NA,length(env))
    
    for (i in 1:length(env))
    {
      field.red <- trait[which(trait$id == env[i]),] %>% droplevels()
      asreml.options(maxit = 50,             
                     workspace = "128mb",       
                     pworkspace = "128mb",
                     trace=F)
      env.1 <- asreml(fixed = value ~  1, 
                      random = ~ Coding + rep + block:rep,
                      data = field.red
      )
      residual[i] <- summary(env.1)$varcomp["units!R", "component"]
      repeatability[i] <- summary(env.1)$varcomp["Coding","component"]/(summary(env.1)$varcomp["Coding","component"] + summary(env.1)$varcomp["block:rep","component"]/2)
      
      repeatability_df <- data.frame("Environment" = env, "Repeatibility" = repeatability)
      
      env.2 <- asreml(fixed = value ~  Coding, 
                      random = ~ rep + block:rep,
                      data = field.red
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
        asreml.options(maxit = 50,             
                       workspace = "128mb",       
                       pworkspace = "128mb",
                       trace=F)
        env.1 <- asreml(fixed = value ~  1, 
                        random = ~ Coding + block,
                        data = field.red
        )
        
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
        unrep <- (tapply(unlist(trait[which(trait$id == i)[1:400], "value"]), unlist(trait[which(trait$id == i)[1:400], "Coding"]), mean))
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
  
  # Calculating heritabilities
  
  sigma.g <- summary(env.1)$varcomp["Coding","component"]
  
  if (length(env.un)>0){
    if (length(env)>0){
      sigma.g.e <- summary(env.1)$varcomp["units!R","component"]-mean(residual)/mean(check_list$rep_indices)
      print(cbind(pheno, sigma.g.e))
      sigma.e <- mean(residual)
    }else{
      sigma.g.e <- NA
      sigma.e <- summary(env.1)$varcomp["units!R","component"]
    }
  }else{sigma.g.e <- summary(env.1)$varcomp["units!R","component"]-(mean(residual)/2)
  print(cbind(pheno, sigma.g.e))
  sigma.e <- mean(residual)
  }  
  
  if(is.na(sigma.g.e)){
    
    herit <-  sigma.g/(sigma.g + sigma.e/length(env.un)) 
    message(paste0("Since there are no replicates the term for sigma g_e was omitted"))
    message( paste0("heritability_entry_mean_based = ", herit))
    
    herit_plot<- sigma.g/(sigma.g + sigma.e) 
    message(paste0("heritability_plot_based = ", herit_plot))
    
  }else{
    
    herit <- sigma.g/(sigma.g + sigma.g.e/length(env) + sigma.e/(2*length(env))) 
    message(paste0("heritability_entry_mean_based = ", herit))
    
    herit_plot<- sigma.g/(sigma.g + sigma.g.e + sigma.e) 
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
  } else if ((pheno %in% traits_trial_2)[1])
    {
    
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
    
    for (i in 1:length(env))
    {
      env_n <- env[i] %>% as.character()
      field.red <- trait %>% filter(id == env_n) %>% droplevels()
      
      asreml.options(maxit = 50,             
                     workspace = "128mb",       
                     pworkspace = "128mb",
                     trace=F)
      env.1 <- asreml(fixed = value ~  1, 
                      random = ~ Coding + rep,
                      data = field.red
      )
      
      residual[i] <- summary(env.1)$varcomp["units!R", "component"]
      repeatability[i] <- summary(env.1)$varcomp["Coding","component"]/(summary(env.1)$varcomp["Coding","component"] + summary(env.1)$varcomp["units!R","component"]/3)
      
      repeatability_df <- data.frame("Environment" = env, "Repeatibility" = repeatability)
      
      env.1 <- asreml(fixed = value ~  Coding, 
                      random = ~ rep,
                      data = field.red
      )
      
      Geno <- predict.asreml(env.1, classify = "Coding",)$pvals[,1:2]
      
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
    )
    
    env.2 <- asreml(fixed = trait ~  Coding, 
                    random = ~ env,
                    data = dta)
    
    sigma.g <- summary(env.1)$varcomp["Coding","component"]
    sigma.g.e <- summary(env.1)$varcomp["units!R","component"] - mean(residual)/mean(check_list$rep_indices)
    print(cbind(pheno, sigma.g.e))
    sigma.e <- mean(residual)
    
    herit_entry_mean <- sigma.g/(sigma.g + sigma.g.e/length(env) + sigma.e/(2*length(env)))
    message(paste0("heritability_entry_based = ", herit_entry_mean))
    
    herit_plot <- sigma.g/(sigma.g + sigma.g.e + sigma.e)
    message(paste0("heritability_plot_based = ", herit_plot))
    
    Blues.acr <- predict.asreml(env.2, classify = "Coding",)$pvals[,2]
    
    com_BLUEs <- cbind(BLUES, Blues.acr)
    
    # Producing outputs
    
    output <- list()
    
    output["Repeatability"] <- list(repeatability_df)
    
    output["Heritabiliy_entry"] <- herit_entry_mean
    
    output["Heritabiliy_plot"] <- herit_plot
    
    output["BLUEs"] <- list(com_BLUEs)
  
    return(output) 
  }
}

my_data <- readRDS("~/GABI/output_data/isa_data.Rdata")

traits_trial_1 <- c("HD", "PH", "TKW" , "EW", "GPE", "GY", "SW", "GH", "STC", "PC", "SDS", "HAG", "ZEL") 
traits_trial_2 <- c( "FHB", "DTR", "SEP") 

for (i in traits_trial_1){
  
trait <- i

print(i)

raw_data <- my_data$trial_1[, c("year", "loc", "rep", "block", "Coding", trait)]

colnames(raw_data)[length(raw_data)] <- "value"

check <- BLUEs(raw_data, pheno = trait)
}

for (i in traits_trial_2){
  trait <- i
  
  print(i)
  
  raw_data <- my_data$trial_2[, c("year", "loc", "rep", "Coding", trait)]
  
  colnames(raw_data)[length(raw_data)] <- "value"
  
  check <- BLUEs(raw_data, pheno = trait)
}

#herit_p
#Trait Hertitability_plot
#HD          0.8495466
#PH          0.8663639
#YIE          0.4367276
#SW          0.6915210
#HAG          0.3475945
#STC          0.6758751
#EW          0.3509110
#TKW          0.7174052
#PC          0.5278772
#ZEL          0.7775175
#GPE          0.3059756
#GH          0.7488266
#SDS          0.6870338
#FHB          0.5584159
#DTR          0.1106060
#SEP          0.4366028

#herit_e
#Trait Hertitability_entry
#HD                  0.983
#PH                  0.986
#YIE                 0.885
#SW                  0.946
#HAG                 0.748
#STC                 0.884
#EW                  0.654
#TKW                 0.950
#PC                  0.840
#ZEL                 0.950
#GPE                 0.503
#GH                 NA    
#SDS                NA    
#FHB                 0.883
#DTR                 0.262
#SEP                 0.694

B.acr <- check$BLUEs[, c("Coding", "Blues.acr")] %>% left_join(my_data$BLUEs[, c("Coding", trait)], by = "Coding")
cor(B.acr$Blues.acr, B.acr[, trait], use = "complete.obs")
