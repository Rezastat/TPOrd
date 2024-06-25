##################################
## Adjacent Category Simulation ##
##################################

library(parallel)
library(TPOrd)


file_path <- "/scratch/a/amitani/fahimi/MSC_Thesis/Experiment_1/log/"



array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset=1))



m <- as.numeric(Sys.getenv("ITERATIONS", "10")) # Number of iterations
mit <- as.numeric(Sys.getenv("MIT", "3"))
num_cores = as.numeric(Sys.getenv("NUM_CORES", "10"))

if (m <= mit) {
  iterations_to_log <- 1:m
} else {
  iterations_to_log <- sample(1:m, mit)
}

iterations_to_log <- sort(iterations_to_log)

target_cor <- as.numeric(Sys.getenv("TARGET_COR", "0.75"))
fractions_cat1 <- as.numeric(strsplit(Sys.getenv("FRACTIONS_CAT1", "0.35,0.25,0.2"), ",")[[1]])
fractions_cat2 <- as.numeric(strsplit(Sys.getenv("FRACTIONS_CAT2", "0.30,0.25,0.25,0.1"), ",")[[1]])
N <- as.numeric(Sys.getenv("N", "5000"))
cor_YZ_break <- as.numeric(Sys.getenv("COR_YZ_BREAK", "0.85"))
num_categories <- as.numeric(Sys.getenv("NUM_CATEGORIES", "5"))
n2_vector <- as.numeric(strsplit(Sys.getenv("N2_VECTOR", "3500,4000"), ",")[[1]])
model_type <- "Adjacent_Category"
Beta0 <- as.numeric(strsplit(Sys.getenv("BETA0", "-1.007578022,-0.510263682,2.252233149,2.27519740"), ",")[[1]])
Beta1 <- as.numeric(Sys.getenv("BETA1", "-0.2068461819"))
fam <- acat(reverse =T, parallel = T)



run_simulation <- function(params) {
  iteration <- params$iter
  n2 <- params$n2
  
  Sys.setenv(OMP_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1",
             OPENBLAS_NUM_THREADS = "1",
             BLAS_NUM_THREADS = "1",
             LAPACK_NUM_THREADS = "1",
             VECLIB_MAXIMUM_THREADS = "1",
             NUMEXPR_NUM_THREADS = "1")
  
  my_seed <- iteration + 1000 * array_id
  set.seed(my_seed)
  
  
  
  dat_sim <- sim_categorical_data(Beta0 = Beta0 ,# intercept
                                  Beta1 = Beta1, # effect size of the sequencing variant
                                  fractions_cat1 =fractions_cat1,
                                  fractions_cat2 =fractions_cat2,
                                  target_cor = target_cor,
                                  N = N, # phase 1  sample size
                                  n2 = n2, # phase 2 sample size
                                  cor_YZ_break = cor_YZ_break,
                                  num_categories = num_categories,
                                  model_type = model_type)
  dat_sim <- sample_dataframe(dat_sim, n2)
  dat_sim$fZ <- as.factor(dat_sim$Z)
  
  com_table <- table(dat_sim$Y, dat_sim$G1)
  R3_table <- table(dat_sim[dat_sim$R3==1,]$Y,dat_sim[dat_sim$R3==1,]$G1)
  R4_table <- table(dat_sim[dat_sim$R4==1,]$Y,dat_sim[dat_sim$R4==1,]$G1)
  
  vglmfit.com.G0 <- vglm(factor(Y, ordered = T) ~ G0+fZ, family = fam, data = dat_sim)
  vglmfit.com.Ho <- vglm(factor(Y, ordered = T) ~ fZ, family = fam,  data = dat_sim)
  out_com0.G0 <- c(coef(vglmfit.com.G0)[num_categories], diag(vcov(vglmfit.com.G0))[num_categories])
  Wcom.G0 <- out_com0.G0[1]^2/out_com0.G0[2]
  LRcom.G0 <- as.numeric(2*(logLik(vglmfit.com.G0) - logLik(vglmfit.com.Ho)))
  Scom.G0 <- Scom <- (score.stat(vglmfit.com.G0,orig.SE=F)^2)[1]
  out_com.G0 <- c( out_com0.G0, Wcom.G0, Scom.G0, LRcom.G0)
  names(out_com.G0)<-c("beta1","var_beta1","W","S","LR")
  
  vglmfit.com.G1 <- vglm(factor(Y, ordered = T) ~ G1+fZ, family = fam, data = dat_sim)
  
  out_com0.G1 <- c(coef(vglmfit.com.G1)[num_categories], diag(vcov(vglmfit.com.G1))[num_categories])
  Wcom.G1 <- out_com0.G1[1]^2/out_com0.G1[2]
  LRcom.G1 <- as.numeric(2*(logLik(vglmfit.com.G1) - logLik(vglmfit.com.Ho)))
  Scom.G1 <- Scom <- (score.stat(vglmfit.com.G1,orig.SE=F)^2)[1]
  out_com.G1 <- c( out_com0.G1, Wcom.G1, Scom.G1, LRcom.G1)
  names(out_com.G1)<-c("beta1","var_beta1","W","S","LR")
  #####################
  
  if (iteration %in% iterations_to_log) {
    # Log file naming includes the iteration number
    log_file_name <- file.path(file_path, sprintf("adcat_log_arrayid_%d_n2_%d_cores_%d_iteration_%d.txt", array_id, n2, num_cores, iteration))
    log_file <- file(log_file_name, open = "wt")
    
    data1.R3.G0 <- dat_sim[dat_sim$R3==1,c("Y","G0","Z","fZ")] # phase 2 data
    names(data1.R3.G0)[names(data1.R3.G0)=="G0"] <- "G"
    data0.R3.G0 <- dat_sim[dat_sim$R3==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R3.G0_time <- system.time({resHo.R3.G0 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R3.G0,data1.R3.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )})[3]
    cat(sprintf("%d %d resHo.R3.G0 %f\n", array_id, iteration, resHo.R3.G0_time), file = log_file, append = TRUE)
    resHa.R3.G0_time <- system.time({resHa.R3.G0 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R3.G0,data1.R3.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)})[3]
    cat(sprintf("%d %d resHa.R3.G0 %f\n", array_id, iteration, resHa.R3.G0_time), file = log_file, append = TRUE)
    index_of_G0 <- which(names(resHa.R3.G0$theta) == "G")
    res.R3.G0<-  c( resHa.R3.G0$theta[index_of_G0], resHa.R3.G0$var_theta[index_of_G0], resHa.R3.G0$Wobs, resHo.R3.G0$Sobs, 2*(resHa.R3.G0$ll-resHo.R3.G0$ll) )
    
    
    data1.R3.G1 <- dat_sim[dat_sim$R3==1,c("Y","G1","Z","fZ")] # phase 2 data
    names(data1.R3.G1)[names(data1.R3.G1)=="G1"] <- "G"
    data0.R3.G1 <- dat_sim[dat_sim$R3==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R3.G1_time <- system.time({resHo.R3.G1 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R3.G1,data1.R3.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )})[3]
    cat(sprintf("%d %d resHo.R3.G1 %f\n", array_id, iteration, resHo.R3.G1_time), file = log_file, append = TRUE)
    resHa.R3.G1_time <- system.time({resHa.R3.G1 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R3.G1,data1.R3.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)})[3]
    cat(sprintf("%d %d resHa.R3.G1 %f\n", array_id, iteration, resHa.R3.G1_time), file = log_file, append = TRUE)
    index_of_G1 <- which(names(resHa.R3.G1$theta) == "G")
    res.R3.G1 <-c( resHa.R3.G1$theta[index_of_G1], resHa.R3.G1$var_theta[index_of_G1], resHa.R3.G1$Wobs, resHo.R3.G1$Sobs, 2*(resHa.R3.G1$ll-resHo.R3.G1$ll) )
    
    data1.R4.G0 <- dat_sim[dat_sim$R4==1,c("Y","G0","Z","fZ")] # phase 2 data
    names(data1.R4.G0)[names(data1.R4.G0)=="G0"] <- "G"
    data0.R4.G0 <- dat_sim[dat_sim$R4==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R4.G0_time <- system.time({resHo.R4.G0 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R4.G0,data1.R4.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )})[3]
    cat(sprintf("%d %d resHo.R4.G0 %f\n", array_id, iteration, resHo.R4.G0_time), file = log_file, append = TRUE)
    resHa.R4.G0_time <- system.time({resHa.R4.G0 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R4.G0,data1.R4.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)})[3]
    cat(sprintf("%d %d resHa.R4.G0 %f\n", array_id, iteration, resHa.R4.G0_time), file = log_file, append = TRUE)
    index_of_G0 <- which(names(resHa.R4.G0$theta) == "G")
    res.R4.G0<-  c( resHa.R4.G0$theta[index_of_G0], resHa.R4.G0$var_theta[index_of_G0], resHa.R4.G0$Wobs, resHo.R4.G0$Sobs, 2*(resHa.R4.G0$ll-resHo.R4.G0$ll) )
    
    data1.R4.G1 <- dat_sim[dat_sim$R4==1,c("Y","G1","Z","fZ")] # phase 2 data
    names(data1.R4.G1)[names(data1.R4.G1)=="G1"] <- "G"
    data0.R4.G1 <- dat_sim[dat_sim$R4==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R4.G1_time <- system.time({resHo.R4.G1 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R4.G1,data1.R4.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type , num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)})[3]
    cat(sprintf("%d %d resHo.R4.G1 %f\n", array_id, iteration, resHo.R4.G1_time), file = log_file, append = TRUE)
    resHa.R4.G1_time <- system.time({resHa.R4.G1 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R4.G1,data1.R4.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)})[3]
    cat(sprintf("%d %d resHa.R4.G1 %f\n", array_id, iteration, resHa.R4.G1_time), file = log_file, append = TRUE)
    index_of_G1 <- which(names(resHa.R4.G1$theta) == "G")
    res.R4.G1 <-c( resHa.R4.G1$theta[index_of_G1], resHa.R4.G1$var_theta[index_of_G1], resHa.R4.G1$Wobs, resHo.R4.G1$Sobs, 2*(resHa.R4.G1$ll-resHo.R4.G1$ll) )
    
    
    # Close the log file
    close(log_file)
  } else {
    data1.R3.G0 <- dat_sim[dat_sim$R3==1,c("Y","G0","Z","fZ")] # phase 2 data
    names(data1.R3.G0)[names(data1.R3.G0)=="G0"] <- "G"
    data0.R3.G0 <- dat_sim[dat_sim$R3==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R3.G0 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R3.G0,data1.R3.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )
    resHa.R3.G0 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R3.G0,data1.R3.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
    index_of_G0 <- which(names(resHa.R3.G0$theta) == "G")
    res.R3.G0<-  c( resHa.R3.G0$theta[index_of_G0], resHa.R3.G0$var_theta[index_of_G0], resHa.R3.G0$Wobs, resHo.R3.G0$Sobs, 2*(resHa.R3.G0$ll-resHo.R3.G0$ll) )
    
    
    data1.R3.G1 <- dat_sim[dat_sim$R3==1,c("Y","G1","Z","fZ")] # phase 2 data
    names(data1.R3.G1)[names(data1.R3.G1)=="G1"] <- "G"
    data0.R3.G1 <- dat_sim[dat_sim$R3==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R3.G1 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R3.G1,data1.R3.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )
    resHa.R3.G1 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R3.G1,data1.R3.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
    index_of_G1 <- which(names(resHa.R3.G1$theta) == "G")
    res.R3.G1 <-c( resHa.R3.G1$theta[index_of_G1], resHa.R3.G1$var_theta[index_of_G1], resHa.R3.G1$Wobs, resHo.R3.G1$Sobs, 2*(resHa.R3.G1$ll-resHo.R3.G1$ll) )
    
    data1.R4.G0 <- dat_sim[dat_sim$R4==1,c("Y","G0","Z","fZ")] # phase 2 data
    names(data1.R4.G0)[names(data1.R4.G0)=="G0"] <- "G"
    data0.R4.G0 <- dat_sim[dat_sim$R4==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R4.G0 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R4.G0,data1.R4.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )
    resHa.R4.G0 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R4.G0,data1.R4.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
    index_of_G0 <- which(names(resHa.R4.G0$theta) == "G")
    res.R4.G0<-  c( resHa.R4.G0$theta[index_of_G0], resHa.R4.G0$var_theta[index_of_G0], resHa.R4.G0$Wobs, resHo.R4.G0$Sobs, 2*(resHa.R4.G0$ll-resHo.R4.G0$ll) )
    
    data1.R4.G1 <- dat_sim[dat_sim$R4==1,c("Y","G1","Z","fZ")] # phase 2 data
    names(data1.R4.G1)[names(data1.R4.G1)=="G1"] <- "G"
    data0.R4.G1 <- dat_sim[dat_sim$R4==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R4.G1 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R4.G1,data1.R4.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type , num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
    resHa.R4.G1 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R4.G1,data1.R4.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
    index_of_G1 <- which(names(resHa.R4.G1$theta) == "G")
    res.R4.G1 <-c( resHa.R4.G1$theta[index_of_G1], resHa.R4.G1$var_theta[index_of_G1], resHa.R4.G1$Wobs, resHo.R4.G1$Sobs, 2*(resHa.R4.G1$ll-resHo.R4.G1$ll) )
    
  }
  
  #############
  dat_nai.R3 <- dat_sim[dat_sim$R3==1,c("Y","G0", "G1", "Z","fZ")] # phase 2 data alone
  ordinalfit.nai.R3.G0 <- vglm(factor(Y,ordered = T) ~ G0+fZ, family = fam,  data = dat_nai.R3)
  ordinalfit.nai.Ho.R3 <- vglm(factor(Y,ordered = T) ~ fZ, family = fam,  data = dat_nai.R3)
  out_nai0.R3.G0 <- c(coef(ordinalfit.nai.R3.G0)[num_categories], diag(vcov(ordinalfit.nai.R3.G0))[num_categories])
  Wnai.R3.G0 <- out_nai0.R3.G0[1]^2/out_nai0.R3.G0[2]
  LRnai.R3.G0 <- as.numeric(2*(logLik(ordinalfit.nai.R3.G0) - logLik(ordinalfit.nai.Ho.R3)))
  Snai.R3.G0 <- (score.stat(ordinalfit.nai.R3.G0,orig.SE=T)^2)[1]
  out_nai.R3.G0 <- c( out_nai0.R3.G0, Wnai.R3.G0, Snai.R3.G0, LRnai.R3.G0)
  
  
  ordinalfit.nai.R3.G1 <- vglm(factor(Y,ordered = T) ~ G1+fZ, family = fam,  data = dat_nai.R3)
  out_nai0.R3.G1 <- c(coef(ordinalfit.nai.R3.G1)[num_categories], diag(vcov(ordinalfit.nai.R3.G1))[num_categories])
  Wnai.R3.G1 <- out_nai0.R3.G1[1]^2/out_nai0.R3.G1[2]
  LRnai.R3.G1 <- as.numeric(2*(logLik(ordinalfit.nai.R3.G1) - logLik(ordinalfit.nai.Ho.R3)))
  Snai.R3.G1 <- (score.stat(ordinalfit.nai.R3.G1,orig.SE=T)^2)[1]
  out_nai.R3.G1 <- c( out_nai0.R3.G1, Wnai.R3.G1, Snai.R3.G1, LRnai.R3.G1)
  
  
  dat_nai.R4 <- dat_sim[dat_sim$R4==1,c("Y","G0", "G1", "Z","fZ")] # phase 2 data alone
  ordinalfit.nai.R4.G0 <- vglm(factor(Y,ordered = T) ~ G0+fZ, family = fam,  data = dat_nai.R4)
  ordinalfit.nai.Ho.R4 <- vglm(factor(Y,ordered = T) ~ fZ, family = fam,  data = dat_nai.R4)
  out_nai0.R4.G0 <- c(coef(ordinalfit.nai.R4.G0)[num_categories], diag(vcov(ordinalfit.nai.R4.G0))[num_categories])
  Wnai.R4.G0 <- out_nai0.R4.G0[1]^2/out_nai0.R4.G0[2]
  LRnai.R4.G0 <- as.numeric(2*(logLik(ordinalfit.nai.R4.G0) - logLik(ordinalfit.nai.Ho.R4)))
  Snai.R4.G0 <- (score.stat(ordinalfit.nai.R4.G0,orig.SE=T)^2)[1]
  out_nai.R4.G0 <- c( out_nai0.R4.G0, Wnai.R4.G0, Snai.R4.G0, LRnai.R4.G0)
  
  
  ordinalfit.nai.R4.G1 <- vglm(factor(Y,ordered = T) ~ G1+fZ, family = fam,  data = dat_nai.R4)
  out_nai0.R4.G1 <- c(coef(ordinalfit.nai.R4.G1)[num_categories], diag(vcov(ordinalfit.nai.R4.G1))[num_categories])
  Wnai.R4.G1 <- out_nai0.R4.G1[1]^2/out_nai0.R4.G1[2]
  LRnai.R4.G1 <- as.numeric(2*(logLik(ordinalfit.nai.R4.G1) - logLik(ordinalfit.nai.Ho.R4)))
  Snai.R4.G1 <- (score.stat(ordinalfit.nai.R4.G1,orig.SE=T)^2)[1]
  out_nai.R4.G1 <- c( out_nai0.R4.G1, Wnai.R4.G1, Snai.R4.G1, LRnai.R4.G1)
  
  beta1_res <- data.frame()
  beta1_res <- data.frame(rbind(out_com.G0, res.R3.G0, out_nai.R3.G0, res.R4.G0, out_nai.R4.G0,
                                out_com.G1, res.R3.G1, out_nai.R3.G1, res.R4.G1, out_nai.R4.G1))
  G_vec <- c(rep("G0", 5), rep("G1", 5))
  Sample_vec <- rep(c("Complete", "Z-Balanced-EM", "Z-Balanced-Naive", "YZ-Balanced-EM", "YZ-Balanced-Naive" ), 2)
  True_Beta1 <- rep(c(0,Beta1), each= 5)
  
  beta1_adj <- cbind(Sample_vec, G_vec, True_Beta1, beta1_res)
  row.names(beta1_adj) <- NULL
  
  ####### Under construction!!!!
  sum.com.Ho <- summary(vglmfit.com.Ho)@coef3[,1:2]
  sum.com.Ho[,2] <- sum.com.Ho[,2]^2
  sum.com.G0 <- summary(vglmfit.com.G0)@coef3[,1:2]
  sum.com.G0[,2] <- sum.com.G0[,2]^2
  sum.com.G1 <- summary(vglmfit.com.G1)@coef3[,1:2]
  sum.com.G1[,2] <- sum.com.G1[,2]^2
  sum.nai.Ho.R3 <- summary(ordinalfit.nai.Ho.R3)@coef3[,1:2]
  sum.nai.Ho.R3[,2] <- sum.nai.Ho.R3[,2]^2
  sum.nai.R3.G0 <- summary(ordinalfit.nai.R3.G0)@coef3[,1:2]
  sum.nai.R3.G0[,2] <- sum.nai.R3.G0[,2]^2
  sum.nai.R3.G1 <- summary(ordinalfit.nai.R3.G1)@coef3[,1:2]
  sum.nai.R3.G1[,2] <- sum.nai.R3.G1[,2]^2
  sum.nai.Ho.R4 <- summary(ordinalfit.nai.Ho.R4)@coef3[,1:2]
  sum.nai.Ho.R4[,2] <- sum.nai.Ho.R4[,2]^2
  sum.nai.R4.G0 <- summary(ordinalfit.nai.R4.G0)@coef3[,1:2]
  sum.nai.R4.G0[,2] <- sum.nai.R4.G0[,2]^2
  sum.nai.R4.G1 <- summary(ordinalfit.nai.R4.G1)@coef3[,1:2]
  sum.nai.R4.G1[,2] <- sum.nai.R4.G1[,2]^2
  
  est.resHa.R3.G0 <- resHa.R3.G0[["theta"]]
  var.resHa.R3.G0 <- diag(solve(resHa.R3.G0[["FIM"]]))[1:length(est.resHa.R3.G0)]
  sum.resHa.R3.G0 <- cbind(est.resHa.R3.G0, var.resHa.R3.G0)
  est.resHo.R3.G1 <- resHo.R3.G1[["theta"]]
  var.resHo.R3.G1 <- diag(solve(resHo.R3.G1[["FIM"]]))[1:length(est.resHo.R3.G1)]
  sum.resHo.R3.G1 <- cbind(est.resHo.R3.G1, var.resHo.R3.G1)
  est.resHo.R4.G0 <- resHo.R4.G0[["theta"]]
  var.resHo.R4.G0 <- diag(solve(resHo.R4.G0[["FIM"]]))[1:length(est.resHo.R4.G0)]
  sum.resHo.R4.G0 <- cbind(est.resHo.R4.G0, var.resHo.R4.G0)
  est.resHo.R4.G1 <- resHo.R4.G1[["theta"]]
  var.resHo.R4.G1 <- diag(solve(resHo.R4.G1[["FIM"]]))[1:length(est.resHo.R4.G1)]
  sum.resHo.R4.G1 <- cbind(est.resHo.R4.G1, var.resHo.R4.G1)
  
  sum.mat<- rbind(sum.com.Ho ,
                  sum.com.G0 ,
                  sum.com.G1 ,
                  sum.nai.Ho.R3 ,
                  sum.nai.R3.G0 ,
                  sum.nai.R3.G1 ,
                  sum.nai.Ho.R4 ,
                  sum.nai.R4.G0 ,
                  sum.nai.R4.G1 ,
                  sum.resHa.R3.G0 ,
                  sum.resHo.R3.G1 ,
                  sum.resHo.R4.G0 ,
                  sum.resHo.R4.G1 )
  
  model_names <- c(rep("com.Ho", dim(sum.com.Ho)[1]) ,
                   rep("com.G0", dim(sum.com.G0)[1]) ,
                   rep("com.G1", dim(sum.com.G1)[1]) ,
                   rep("nai.Ho.R3", dim(sum.nai.Ho.R3)[1]) ,
                   rep("nai.R3.G0", dim(sum.nai.R3.G0)[1]) ,
                   rep("nai.R3.G1", dim(sum.nai.R3.G1)[1]) ,
                   rep("nai.Ho.R4", dim(sum.nai.Ho.R4)[1]) ,
                   rep ("nai.R4.G0", dim(sum.nai.R4.G0)[1]) ,
                   rep("nai.R4.G1", dim(sum.nai.R4.G1)[1]) ,
                   rep ("resHa.R3.G0", dim(sum.resHa.R3.G0)[1]) ,
                   rep ("resHo.R3.G1", dim(sum.resHo.R3.G1)[1]) ,
                   rep ("resHo.R4.G0", dim(sum.resHo.R4.G0)[1]) ,
                   rep ("resHo.R4.G1", dim(sum.resHo.R4.G1)[1]) )
  
  sum.dat <- data.frame(model_names= model_names,
                        coef_name = row.names(sum.mat),
                        est = sum.mat[,1],
                        var = sum.mat[,2])
  
  
  
  
  
  adj_dat <- list(iteration= iteration,
                  seed = my_seed,
                  sample_freq=list(com_table=com_table, R3_table=R3_table, R4_table=R4_table),
                  beta1_table = beta1_adj,
                  est_var = sum.dat)
  
  ###############
  
  return(adj_dat)
}


closeAllConnections()


cl <- makeCluster(num_cores)
clusterExport(cl, varlist = c(
  "array_id", "m", "num_cores",
  "target_cor", "fractions_cat1", "fractions_cat2",
  "N", "cor_YZ_break", "num_categories",
  "model_type", "Beta0", "Beta1",
  "fam", "iterations_to_log", "file_path", "run_simulation"
))
clusterEvalQ(cl, {
  library(TPOrd)
})

# Results initialization
all_results <- list()

# Run simulations for each n2 value
for (n2 in n2_vector) {
  tasks <- lapply(1:m, function(iter) list(iter=iter, n2=n2))
  
  total.time <- system.time({
    results <- parLapply(cl, tasks, run_simulation)
  })[3]
  
  # Save results
  save(total.time, results, file = sprintf("adcat_R3R4_arrayid_%d_n2_%d_cores_%d_iterations_%d.RData", array_id, n2, num_cores, m))
}

stopCluster(cl)
closeAllConnections()
