##################################
## Adjacent Category Simulation ##
##################################

library(parallel)
library(TPOrd)


file_path <- "/scratch/a/amitani/fahimi/MSC_Thesis/Experiment_2/log/"



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
  R1_table <- table(dat_sim[dat_sim$R1==1,]$Y,dat_sim[dat_sim$R1==1,]$G1)
  R2_table <- table(dat_sim[dat_sim$R2==1,]$Y,dat_sim[dat_sim$R2==1,]$G1)

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

    data1.R1.G0 <- dat_sim[dat_sim$R1==1,c("Y","G0","Z","fZ")] # phase 2 data
    names(data1.R1.G0)[names(data1.R1.G0)=="G0"] <- "G"
    data0.R1.G0 <- dat_sim[dat_sim$R1==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R1.G0_time <- system.time({resHo.R1.G0 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R1.G0,data1.R1.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )})[3]
    cat(sprintf("%d %d resHo.R1.G0 %f\n", array_id, iteration, resHo.R1.G0_time), file = log_file, append = TRUE)
    resHa.R1.G0_time <- system.time({resHa.R1.G0 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R1.G0,data1.R1.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)})[3]
    cat(sprintf("%d %d resHa.R1.G0 %f\n", array_id, iteration, resHa.R1.G0_time), file = log_file, append = TRUE)
    index_of_G0 <- which(names(resHa.R1.G0$theta) == "G")
    res.R1.G0<-  c( resHa.R1.G0$theta[index_of_G0], resHa.R1.G0$var_theta[index_of_G0], resHa.R1.G0$Wobs, resHo.R1.G0$Sobs, 2*(resHa.R1.G0$ll-resHo.R1.G0$ll) )


    data1.R1.G1 <- dat_sim[dat_sim$R1==1,c("Y","G1","Z","fZ")] # phase 2 data
    names(data1.R1.G1)[names(data1.R1.G1)=="G1"] <- "G"
    data0.R1.G1 <- dat_sim[dat_sim$R1==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R1.G1_time <- system.time({resHo.R1.G1 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R1.G1,data1.R1.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )})[3]
    cat(sprintf("%d %d resHo.R1.G1 %f\n", array_id, iteration, resHo.R1.G1_time), file = log_file, append = TRUE)
    resHa.R1.G1_time <- system.time({resHa.R1.G1 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R1.G1,data1.R1.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)})[3]
    cat(sprintf("%d %d resHa.R1.G1 %f\n", array_id, iteration, resHa.R1.G1_time), file = log_file, append = TRUE)
    index_of_G1 <- which(names(resHa.R1.G1$theta) == "G")
    res.R1.G1 <-c( resHa.R1.G1$theta[index_of_G1], resHa.R1.G1$var_theta[index_of_G1], resHa.R1.G1$Wobs, resHo.R1.G1$Sobs, 2*(resHa.R1.G1$ll-resHo.R1.G1$ll) )

    data1.R2.G0 <- dat_sim[dat_sim$R2==1,c("Y","G0","Z","fZ")] # phase 2 data
    names(data1.R2.G0)[names(data1.R2.G0)=="G0"] <- "G"
    data0.R2.G0 <- dat_sim[dat_sim$R2==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R2.G0_time <- system.time({resHo.R2.G0 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R2.G0,data1.R2.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )})[3]
    cat(sprintf("%d %d resHo.R2.G0 %f\n", array_id, iteration, resHo.R2.G0_time), file = log_file, append = TRUE)
    resHa.R2.G0_time <- system.time({resHa.R2.G0 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R2.G0,data1.R2.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)})[3]
    cat(sprintf("%d %d resHa.R2.G0 %f\n", array_id, iteration, resHa.R2.G0_time), file = log_file, append = TRUE)
    index_of_G0 <- which(names(resHa.R2.G0$theta) == "G")
    res.R2.G0<-  c( resHa.R2.G0$theta[index_of_G0], resHa.R2.G0$var_theta[index_of_G0], resHa.R2.G0$Wobs, resHo.R2.G0$Sobs, 2*(resHa.R2.G0$ll-resHo.R2.G0$ll) )

    data1.R2.G1 <- dat_sim[dat_sim$R2==1,c("Y","G1","Z","fZ")] # phase 2 data
    names(data1.R2.G1)[names(data1.R2.G1)=="G1"] <- "G"
    data0.R2.G1 <- dat_sim[dat_sim$R2==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R2.G1_time <- system.time({resHo.R2.G1 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R2.G1,data1.R2.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type , num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)})[3]
    cat(sprintf("%d %d resHo.R2.G1 %f\n", array_id, iteration, resHo.R2.G1_time), file = log_file, append = TRUE)
    resHa.R2.G1_time <- system.time({resHa.R2.G1 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R2.G1,data1.R2.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)})[3]
    cat(sprintf("%d %d resHa.R2.G1 %f\n", array_id, iteration, resHa.R2.G1_time), file = log_file, append = TRUE)
    index_of_G1 <- which(names(resHa.R2.G1$theta) == "G")
    res.R2.G1 <-c( resHa.R2.G1$theta[index_of_G1], resHa.R2.G1$var_theta[index_of_G1], resHa.R2.G1$Wobs, resHo.R2.G1$Sobs, 2*(resHa.R2.G1$ll-resHo.R2.G1$ll) )


    # Close the log file
    close(log_file)
  } else {
    data1.R1.G0 <- dat_sim[dat_sim$R1==1,c("Y","G0","Z","fZ")] # phase 2 data
    names(data1.R1.G0)[names(data1.R1.G0)=="G0"] <- "G"
    data0.R1.G0 <- dat_sim[dat_sim$R1==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R1.G0 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R1.G0,data1.R1.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )
    resHa.R1.G0 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R1.G0,data1.R1.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
    index_of_G0 <- which(names(resHa.R1.G0$theta) == "G")
    res.R1.G0<-  c( resHa.R1.G0$theta[index_of_G0], resHa.R1.G0$var_theta[index_of_G0], resHa.R1.G0$Wobs, resHo.R1.G0$Sobs, 2*(resHa.R1.G0$ll-resHo.R1.G0$ll) )


    data1.R1.G1 <- dat_sim[dat_sim$R1==1,c("Y","G1","Z","fZ")] # phase 2 data
    names(data1.R1.G1)[names(data1.R1.G1)=="G1"] <- "G"
    data0.R1.G1 <- dat_sim[dat_sim$R1==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R1.G1 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R1.G1,data1.R1.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )
    resHa.R1.G1 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R1.G1,data1.R1.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
    index_of_G1 <- which(names(resHa.R1.G1$theta) == "G")
    res.R1.G1 <-c( resHa.R1.G1$theta[index_of_G1], resHa.R1.G1$var_theta[index_of_G1], resHa.R1.G1$Wobs, resHo.R1.G1$Sobs, 2*(resHa.R1.G1$ll-resHo.R1.G1$ll) )

    data1.R2.G0 <- dat_sim[dat_sim$R2==1,c("Y","G0","Z","fZ")] # phase 2 data
    names(data1.R2.G0)[names(data1.R2.G0)=="G0"] <- "G"
    data0.R2.G0 <- dat_sim[dat_sim$R2==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R2.G0 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R2.G0,data1.R2.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )
    resHa.R2.G0 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R2.G0,data1.R2.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
    index_of_G0 <- which(names(resHa.R2.G0$theta) == "G")
    res.R2.G0<-  c( resHa.R2.G0$theta[index_of_G0], resHa.R2.G0$var_theta[index_of_G0], resHa.R2.G0$Wobs, resHo.R2.G0$Sobs, 2*(resHa.R2.G0$ll-resHo.R2.G0$ll) )

    data1.R2.G1 <- dat_sim[dat_sim$R2==1,c("Y","G1","Z","fZ")] # phase 2 data
    names(data1.R2.G1)[names(data1.R2.G1)=="G1"] <- "G"
    data0.R2.G1 <- dat_sim[dat_sim$R2==0,c("Y","Z","fZ")] # phase 2 data complement
    resHo.R2.G1 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R2.G1,data1.R2.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type , num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
    resHa.R2.G1 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R2.G1,data1.R2.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
    index_of_G1 <- which(names(resHa.R2.G1$theta) == "G")
    res.R2.G1 <-c( resHa.R2.G1$theta[index_of_G1], resHa.R2.G1$var_theta[index_of_G1], resHa.R2.G1$Wobs, resHo.R2.G1$Sobs, 2*(resHa.R2.G1$ll-resHo.R2.G1$ll) )

  }

  #############
  dat_nai.R1 <- dat_sim[dat_sim$R1==1,c("Y","G0", "G1", "Z","fZ")] # phase 2 data alone
  ordinalfit.nai.R1.G0 <- vglm(factor(Y,ordered = T) ~ G0+fZ, family = fam,  data = dat_nai.R1)
  ordinalfit.nai.Ho.R1 <- vglm(factor(Y,ordered = T) ~ fZ, family = fam,  data = dat_nai.R1)
  out_nai0.R1.G0 <- c(coef(ordinalfit.nai.R1.G0)[num_categories], diag(vcov(ordinalfit.nai.R1.G0))[num_categories])
  Wnai.R1.G0 <- out_nai0.R1.G0[1]^2/out_nai0.R1.G0[2]
  LRnai.R1.G0 <- as.numeric(2*(logLik(ordinalfit.nai.R1.G0) - logLik(ordinalfit.nai.Ho.R1)))
  Snai.R1.G0 <- (score.stat(ordinalfit.nai.R1.G0,orig.SE=T)^2)[1]
  out_nai.R1.G0 <- c( out_nai0.R1.G0, Wnai.R1.G0, Snai.R1.G0, LRnai.R1.G0)


  ordinalfit.nai.R1.G1 <- vglm(factor(Y,ordered = T) ~ G1+fZ, family = fam,  data = dat_nai.R1)
  out_nai0.R1.G1 <- c(coef(ordinalfit.nai.R1.G1)[num_categories], diag(vcov(ordinalfit.nai.R1.G1))[num_categories])
  Wnai.R1.G1 <- out_nai0.R1.G1[1]^2/out_nai0.R1.G1[2]
  LRnai.R1.G1 <- as.numeric(2*(logLik(ordinalfit.nai.R1.G1) - logLik(ordinalfit.nai.Ho.R1)))
  Snai.R1.G1 <- (score.stat(ordinalfit.nai.R1.G1,orig.SE=T)^2)[1]
  out_nai.R1.G1 <- c( out_nai0.R1.G1, Wnai.R1.G1, Snai.R1.G1, LRnai.R1.G1)


  dat_nai.R2 <- dat_sim[dat_sim$R2==1,c("Y","G0", "G1", "Z","fZ")] # phase 2 data alone
  ordinalfit.nai.R2.G0 <- vglm(factor(Y,ordered = T) ~ G0+fZ, family = fam,  data = dat_nai.R2)
  ordinalfit.nai.Ho.R2 <- vglm(factor(Y,ordered = T) ~ fZ, family = fam,  data = dat_nai.R2)
  out_nai0.R2.G0 <- c(coef(ordinalfit.nai.R2.G0)[num_categories], diag(vcov(ordinalfit.nai.R2.G0))[num_categories])
  Wnai.R2.G0 <- out_nai0.R2.G0[1]^2/out_nai0.R2.G0[2]
  LRnai.R2.G0 <- as.numeric(2*(logLik(ordinalfit.nai.R2.G0) - logLik(ordinalfit.nai.Ho.R2)))
  Snai.R2.G0 <- (score.stat(ordinalfit.nai.R2.G0,orig.SE=T)^2)[1]
  out_nai.R2.G0 <- c( out_nai0.R2.G0, Wnai.R2.G0, Snai.R2.G0, LRnai.R2.G0)


  ordinalfit.nai.R2.G1 <- vglm(factor(Y,ordered = T) ~ G1+fZ, family = fam,  data = dat_nai.R2)
  out_nai0.R2.G1 <- c(coef(ordinalfit.nai.R2.G1)[num_categories], diag(vcov(ordinalfit.nai.R2.G1))[num_categories])
  Wnai.R2.G1 <- out_nai0.R2.G1[1]^2/out_nai0.R2.G1[2]
  LRnai.R2.G1 <- as.numeric(2*(logLik(ordinalfit.nai.R2.G1) - logLik(ordinalfit.nai.Ho.R2)))
  Snai.R2.G1 <- (score.stat(ordinalfit.nai.R2.G1,orig.SE=T)^2)[1]
  out_nai.R2.G1 <- c( out_nai0.R2.G1, Wnai.R2.G1, Snai.R2.G1, LRnai.R2.G1)

  beta1_res <- data.frame()
  beta1_res <- data.frame(rbind(out_com.G0, res.R1.G0, out_nai.R1.G0, res.R2.G0, out_nai.R2.G0,
                                out_com.G1, res.R1.G1, out_nai.R1.G1, res.R2.G1, out_nai.R2.G1))
  G_vec <- c(rep("G0", 5), rep("G1", 5))
  Sample_vec <- rep(c("Complete", "SRS-EM", "SRS-Naive", "Balanced-EM", "Balanced-Naive"), 2)
  True_Beta1 <- rep(c(0,Beta1), each= 5)

  beta1_adj <- cbind(Sample_vec, G_vec, True_Beta1, beta1_res)
  row.names(beta1_adj) <- NULL







  adj_dat <- list(iteration= iteration,
                  seed = my_seed,
                  sample_freq=list(com_table=com_table, R1_table=R1_table, R2_table=R2_table),
                  beta1_table = beta1_adj)

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
  save(total.time, results, file = sprintf("adcat_data_arrayid_%d_n2_%d_cores_%d_iterations_%d.RData", array_id, n2, num_cores, m))
}

stopCluster(cl)
closeAllConnections()
