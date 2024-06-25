##################################
##  Genetic ##########
##################################


library(parallel)
library(TPOrd)

job_name <- Sys.getenv("JOB_NAME", "stopratanalysis")
data_name <- Sys.getenv("DATA_NAME")

Sys.setenv(OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           BLAS_NUM_THREADS = "1",
           LAPACK_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1",
           NUMEXPR_NUM_THREADS = "1")

Sys.getenv("MKL_NUM_THREADS","1")  # Should return '1'
Sys.getenv("MKL_DYNAMIC","FALSE")      # Should return 'FALSE'


array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset=4))



m <- as.numeric(Sys.getenv("ITERATIONS", "1")) # Number of iterations
mit <- as.numeric(Sys.getenv("MIT", "1"))
num_cores = as.numeric(Sys.getenv("NUM_CORES", "1"))



# Adjusting parameters using environment variables
target_cor <- as.numeric(Sys.getenv("TARGET_COR", "0.75"))
fractions_cat1 <- as.numeric(strsplit(Sys.getenv("FRACTIONS_CAT1", "0.35,0.25,0.2"), ",")[[1]])
fractions_cat2 <- as.numeric(strsplit(Sys.getenv("FRACTIONS_CAT2", "0.30,0.25,0.25,0.1"), ",")[[1]])
N <- as.numeric(Sys.getenv("N", "5000"))
cor_YZ_break <- as.numeric(Sys.getenv("COR_YZ_BREAK", "0.85"))
num_categories <- as.numeric(Sys.getenv("NUM_CATEGORIES", "5"))
n2_vector <- as.numeric(strsplit(Sys.getenv("N2_VECTOR", "2000"), ",")[[1]])
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
  
  
  my_seed=params$seed
  dat_sim=params$dat_sim
  sample_freq=params$sample_freq
  
  
  
  
  
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
  
  data1.R5.G0 <- dat_sim[dat_sim$R5==1,c("Y","G0","Z","fZ")] # phase 2 data
  names(data1.R5.G0)[names(data1.R5.G0)=="G0"] <- "G"
  data0.R5.G0 <- dat_sim[dat_sim$R5==0,c("Y","Z","fZ")] # phase 2 data complement
  resHo.R5.G0 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R5.G0,data1.R5.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )
  resHa.R5.G0 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R5.G0,data1.R5.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
  index_of_G0 <- which(names(resHa.R5.G0$theta) == "G")
  res.R5.G0<-  c( resHa.R5.G0$theta[index_of_G0], resHa.R5.G0$var_theta[index_of_G0], resHa.R5.G0$Wobs, resHo.R5.G0$Sobs, 2*(resHa.R5.G0$ll-resHo.R5.G0$ll) )
  
  
  data1.R5.G1 <- dat_sim[dat_sim$R5==1,c("Y","G1","Z","fZ")] # phase 2 data
  names(data1.R5.G1)[names(data1.R5.G1)=="G1"] <- "G"
  data0.R5.G1 <- dat_sim[dat_sim$R5==0,c("Y","Z","fZ")] # phase 2 data complement
  resHo.R5.G1 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R5.G1,data1.R5.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )
  resHa.R5.G1 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R5.G1,data1.R5.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
  index_of_G1 <- which(names(resHa.R5.G1$theta) == "G")
  res.R5.G1 <-c( resHa.R5.G1$theta[index_of_G1], resHa.R5.G1$var_theta[index_of_G1], resHa.R5.G1$Wobs, resHo.R5.G1$Sobs, 2*(resHa.R5.G1$ll-resHo.R5.G1$ll) )
  
  
  data1.R6.G0 <- dat_sim[dat_sim$R6==1,c("Y","G0","Z","fZ")] # phase 2 data
  names(data1.R6.G0)[names(data1.R6.G0)=="G0"] <- "G"
  data0.R6.G0 <- dat_sim[dat_sim$R6==0,c("Y","Z","fZ")] # phase 2 data complement
  resHo.R6.G0 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R6.G0,data1.R6.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )
  resHa.R6.G0 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R6.G0,data1.R6.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
  index_of_G0 <- which(names(resHa.R6.G0$theta) == "G")
  res.R6.G0<-  c( resHa.R6.G0$theta[index_of_G0], resHa.R6.G0$var_theta[index_of_G0], resHa.R6.G0$Wobs, resHo.R6.G0$Sobs, 2*(resHa.R6.G0$ll-resHo.R6.G0$ll) )
  
  
  data1.R6.G1 <- dat_sim[dat_sim$R6==1,c("Y","G1","Z","fZ")] # phase 2 data
  names(data1.R6.G1)[names(data1.R6.G1)=="G1"] <- "G"
  data0.R6.G1 <- dat_sim[dat_sim$R6==0,c("Y","Z","fZ")] # phase 2 data complement
  resHo.R6.G1 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R6.G1,data1.R6.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )
  resHa.R6.G1 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R6.G1,data1.R6.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
  index_of_G1 <- which(names(resHa.R6.G1$theta) == "G")
  res.R6.G1 <-c( resHa.R6.G1$theta[index_of_G1], resHa.R6.G1$var_theta[index_of_G1], resHa.R6.G1$Wobs, resHo.R6.G1$Sobs, 2*(resHa.R6.G1$ll-resHo.R6.G1$ll) )
  
  data1.R7.G0 <- dat_sim[dat_sim$R7==1,c("Y","G0","Z","fZ")] # phase 2 data
  names(data1.R7.G0)[names(data1.R7.G0)=="G0"] <- "G"
  data0.R7.G0 <- dat_sim[dat_sim$R7==0,c("Y","Z","fZ")] # phase 2 data complement
  resHo.R7.G0 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R7.G0,data1.R7.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )
  resHa.R7.G0 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R7.G0,data1.R7.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
  index_of_G0 <- which(names(resHa.R7.G0$theta) == "G")
  res.R7.G0<-  c( resHa.R7.G0$theta[index_of_G0], resHa.R7.G0$var_theta[index_of_G0], resHa.R7.G0$Wobs, resHo.R7.G0$Sobs, 2*(resHa.R7.G0$ll-resHo.R7.G0$ll) )
  
  data1.R7.G1 <- dat_sim[dat_sim$R7==1,c("Y","G1","Z","fZ")] # phase 2 data
  names(data1.R7.G1)[names(data1.R7.G1)=="G1"] <- "G"
  data0.R7.G1 <- dat_sim[dat_sim$R7==0,c("Y","Z","fZ")] # phase 2 data complement
  resHo.R7.G1 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R7.G1,data1.R7.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type , num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
  resHa.R7.G1 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R7.G1,data1.R7.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
  index_of_G1 <- which(names(resHa.R7.G1$theta) == "G")
  res.R7.G1 <-c( resHa.R7.G1$theta[index_of_G1], resHa.R7.G1$var_theta[index_of_G1], resHa.R7.G1$Wobs, resHo.R7.G1$Sobs, 2*(resHa.R7.G1$ll-resHo.R7.G1$ll) )
  
  data1.R8.G0 <- dat_sim[dat_sim$R8==1,c("Y","G0","Z","fZ")] # phase 2 data
  names(data1.R8.G0)[names(data1.R8.G0)=="G0"] <- "G"
  data0.R8.G0 <- dat_sim[dat_sim$R8==0,c("Y","Z","fZ")] # phase 2 data complement
  resHo.R8.G0 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R8.G0,data1.R8.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id )
  resHa.R8.G0 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R8.G0,data1.R8.G0,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
  index_of_G0 <- which(names(resHa.R8.G0$theta) == "G")
  res.R8.G0<-  c( resHa.R8.G0$theta[index_of_G0], resHa.R8.G0$var_theta[index_of_G0], resHa.R8.G0$Wobs, resHo.R8.G0$Sobs, 2*(resHa.R8.G0$ll-resHo.R8.G0$ll) )
  
  data1.R8.G1 <- dat_sim[dat_sim$R8==1,c("Y","G1","Z","fZ")] # phase 2 data
  names(data1.R8.G1)[names(data1.R8.G1)=="G1"] <- "G"
  data0.R8.G1 <- dat_sim[dat_sim$R8==0,c("Y","Z","fZ")] # phase 2 data complement
  resHo.R8.G1 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R8.G1,data1.R8.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type , num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
  resHa.R8.G1 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R8.G1,data1.R8.G1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type, num_categories = num_categories, N=N, iteration = iteration, array_id = array_id)
  index_of_G1 <- which(names(resHa.R8.G1$theta) == "G")
  res.R8.G1 <-c( resHa.R8.G1$theta[index_of_G1], resHa.R8.G1$var_theta[index_of_G1], resHa.R8.G1$Wobs, resHo.R8.G1$Sobs, 2*(resHa.R8.G1$ll-resHo.R8.G1$ll) )
  
  
  #####################
  
  
  
  #############
  dat_nai.R5 <- dat_sim[dat_sim$R5==1,c("Y","G0", "G1", "Z","fZ")] # phase 2 data alone
  ordinalfit.nai.R5.G0 <- vglm(factor(Y,ordered = T) ~ G0+fZ, family = fam,  data = dat_nai.R5)
  ordinalfit.nai.Ho.R5 <- vglm(factor(Y,ordered = T) ~ fZ, family = fam,  data = dat_nai.R5)
  out_nai0.R5.G0 <- c(coef(ordinalfit.nai.R5.G0)[num_categories], diag(vcov(ordinalfit.nai.R5.G0))[num_categories])
  Wnai.R5.G0 <- out_nai0.R5.G0[1]^2/out_nai0.R5.G0[2]
  LRnai.R5.G0 <- as.numeric(2*(logLik(ordinalfit.nai.R5.G0) - logLik(ordinalfit.nai.Ho.R5)))
  Snai.R5.G0 <- (score.stat(ordinalfit.nai.R5.G0,orig.SE=T)^2)[1]
  out_nai.R5.G0 <- c( out_nai0.R5.G0, Wnai.R5.G0, Snai.R5.G0, LRnai.R5.G0)
  
  
  ordinalfit.nai.R5.G1 <- vglm(factor(Y,ordered = T) ~ G1+fZ, family = fam,  data = dat_nai.R5)
  out_nai0.R5.G1 <- c(coef(ordinalfit.nai.R5.G1)[num_categories], diag(vcov(ordinalfit.nai.R5.G1))[num_categories])
  Wnai.R5.G1 <- out_nai0.R5.G1[1]^2/out_nai0.R5.G1[2]
  LRnai.R5.G1 <- as.numeric(2*(logLik(ordinalfit.nai.R5.G1) - logLik(ordinalfit.nai.Ho.R5)))
  Snai.R5.G1 <- (score.stat(ordinalfit.nai.R5.G1,orig.SE=T)^2)[1]
  out_nai.R5.G1 <- c( out_nai0.R5.G1, Wnai.R5.G1, Snai.R5.G1, LRnai.R5.G1)
  
  dat_nai.R6 <- dat_sim[dat_sim$R6==1,c("Y","G0", "G1", "Z","fZ")] # phase 2 data alone
  ordinalfit.nai.R6.G0 <- vglm(factor(Y,ordered = T) ~ G0+fZ, family = fam,  data = dat_nai.R6)
  ordinalfit.nai.Ho.R6 <- vglm(factor(Y,ordered = T) ~ fZ, family = fam,  data = dat_nai.R6)
  out_nai0.R6.G0 <- c(coef(ordinalfit.nai.R6.G0)[num_categories], diag(vcov(ordinalfit.nai.R6.G0))[num_categories])
  Wnai.R6.G0 <- out_nai0.R6.G0[1]^2/out_nai0.R6.G0[2]
  LRnai.R6.G0 <- as.numeric(2*(logLik(ordinalfit.nai.R6.G0) - logLik(ordinalfit.nai.Ho.R6)))
  Snai.R6.G0 <- (score.stat(ordinalfit.nai.R6.G0,orig.SE=T)^2)[1]
  out_nai.R6.G0 <- c( out_nai0.R6.G0, Wnai.R6.G0, Snai.R6.G0, LRnai.R6.G0)
  
  
  ordinalfit.nai.R6.G1 <- vglm(factor(Y,ordered = T) ~ G1+fZ, family = fam,  data = dat_nai.R6)
  out_nai0.R6.G1 <- c(coef(ordinalfit.nai.R6.G1)[num_categories], diag(vcov(ordinalfit.nai.R6.G1))[num_categories])
  Wnai.R6.G1 <- out_nai0.R6.G1[1]^2/out_nai0.R6.G1[2]
  LRnai.R6.G1 <- as.numeric(2*(logLik(ordinalfit.nai.R6.G1) - logLik(ordinalfit.nai.Ho.R6)))
  Snai.R6.G1 <- (score.stat(ordinalfit.nai.R6.G1,orig.SE=T)^2)[1]
  out_nai.R6.G1 <- c( out_nai0.R6.G1, Wnai.R6.G1, Snai.R6.G1, LRnai.R6.G1)
  
  
  dat_nai.R7 <- dat_sim[dat_sim$R7==1,c("Y","G0", "G1", "Z","fZ")] # phase 2 data alone
  ordinalfit.nai.R7.G0 <- vglm(factor(Y,ordered = T) ~ G0+fZ, family = fam,  data = dat_nai.R7)
  ordinalfit.nai.Ho.R7 <- vglm(factor(Y,ordered = T) ~ fZ, family = fam,  data = dat_nai.R7)
  out_nai0.R7.G0 <- c(coef(ordinalfit.nai.R7.G0)[num_categories], diag(vcov(ordinalfit.nai.R7.G0))[num_categories])
  Wnai.R7.G0 <- out_nai0.R7.G0[1]^2/out_nai0.R7.G0[2]
  LRnai.R7.G0 <- as.numeric(2*(logLik(ordinalfit.nai.R7.G0) - logLik(ordinalfit.nai.Ho.R7)))
  Snai.R7.G0 <- (score.stat(ordinalfit.nai.R7.G0,orig.SE=T)^2)[1]
  out_nai.R7.G0 <- c( out_nai0.R7.G0, Wnai.R7.G0, Snai.R7.G0, LRnai.R7.G0)
  
  
  ordinalfit.nai.R7.G1 <- vglm(factor(Y,ordered = T) ~ G1+fZ, family = fam,  data = dat_nai.R7)
  out_nai0.R7.G1 <- c(coef(ordinalfit.nai.R7.G1)[num_categories], diag(vcov(ordinalfit.nai.R7.G1))[num_categories])
  Wnai.R7.G1 <- out_nai0.R7.G1[1]^2/out_nai0.R7.G1[2]
  LRnai.R7.G1 <- as.numeric(2*(logLik(ordinalfit.nai.R7.G1) - logLik(ordinalfit.nai.Ho.R7)))
  Snai.R7.G1 <- (score.stat(ordinalfit.nai.R7.G1,orig.SE=T)^2)[1]
  out_nai.R7.G1 <- c( out_nai0.R7.G1, Wnai.R7.G1, Snai.R7.G1, LRnai.R7.G1)
  
  dat_nai.R8 <- dat_sim[dat_sim$R8==1,c("Y","G0", "G1", "Z","fZ")] # phase 2 data alone
  ordinalfit.nai.R8.G0 <- vglm(factor(Y,ordered = T) ~ G0+fZ, family = fam,  data = dat_nai.R8)
  ordinalfit.nai.Ho.R8 <- vglm(factor(Y,ordered = T) ~ fZ, family = fam,  data = dat_nai.R8)
  out_nai0.R8.G0 <- c(coef(ordinalfit.nai.R8.G0)[num_categories], diag(vcov(ordinalfit.nai.R8.G0))[num_categories])
  Wnai.R8.G0 <- out_nai0.R8.G0[1]^2/out_nai0.R8.G0[2]
  LRnai.R8.G0 <- as.numeric(2*(logLik(ordinalfit.nai.R8.G0) - logLik(ordinalfit.nai.Ho.R8)))
  Snai.R8.G0 <- (score.stat(ordinalfit.nai.R8.G0,orig.SE=T)^2)[1]
  out_nai.R8.G0 <- c( out_nai0.R8.G0, Wnai.R8.G0, Snai.R8.G0, LRnai.R8.G0)
  
  
  ordinalfit.nai.R8.G1 <- vglm(factor(Y,ordered = T) ~ G1+fZ, family = fam,  data = dat_nai.R8)
  out_nai0.R8.G1 <- c(coef(ordinalfit.nai.R8.G1)[num_categories], diag(vcov(ordinalfit.nai.R8.G1))[num_categories])
  Wnai.R8.G1 <- out_nai0.R8.G1[1]^2/out_nai0.R8.G1[2]
  LRnai.R8.G1 <- as.numeric(2*(logLik(ordinalfit.nai.R8.G1) - logLik(ordinalfit.nai.Ho.R8)))
  Snai.R8.G1 <- (score.stat(ordinalfit.nai.R8.G1,orig.SE=T)^2)[1]
  out_nai.R8.G1 <- c( out_nai0.R8.G1, Wnai.R8.G1, Snai.R8.G1, LRnai.R8.G1)
  
  beta1_res <- data.frame()
  beta1_res <- data.frame(rbind(out_com.G0, res.R5.G0, out_nai.R5.G0, res.R6.G0, out_nai.R6.G0, res.R7.G0, out_nai.R7.G0, res.R8.G0, out_nai.R8.G0,
                                out_com.G1, res.R5.G1, out_nai.R5.G1, res.R6.G1, out_nai.R6.G1, res.R7.G1, out_nai.R7.G1, res.R8.G1, out_nai.R8.G1))
  G_vec <- c(rep("G0", 9), rep("G1", 9))
  Sample_vec <- rep(c("Complete", "Fit_SRS-EM", "Fit-SRS-Naive", "GA-Y-EM", "GA-Y-Naive", "GA-YZ-EM", "GA-YZ-Naive", "GA-Z-EM", "GA-Z-Naive"), 2)
  True_Beta1 <- rep(c(0,Beta1), each=9)
  
  beta1_stop <- cbind(Sample_vec, G_vec, True_Beta1, beta1_res)
  row.names(beta1_stop) <- NULL
  
  
  
  
  
  adcat_dat <- list(
    iteration = iteration,
    seed = my_seed,
    sample_freq = sample_freq,
    beta1_table = beta1_stop
  )
  
  
  
  return(adcat_dat)
}

closeAllConnections()


cl <- makeCluster(num_cores)
clusterExport(cl, varlist = c(
  "array_id", "m", "num_cores",
  "target_cor", "fractions_cat1", "fractions_cat2",
  "N", "cor_YZ_break", "num_categories",
  "model_type", "Beta0", "Beta1",
  "fam", "run_simulation", "data_name"
))
clusterEvalQ(cl, {
  library(TPOrd)
})

# Results initialization
all_results <- list()

# Run simulations for each n2 value
for (n2 in n2_vector) {
  dataset_file_name <- sprintf("%s_arrayid_%d_n2_%d_cores_%d_iterations_%d.RData", data_name, array_id, n2, num_cores, m)
  
  # Load the corresponding dataset
  load(dataset_file_name)
  
  tasks <- lapply(1:m, function(iter) list(iter = iter, n2 = n2, dat_sim = results[[iter]]$dataset, sample_freq = results[[iter]][["sample_freq"]], seed = results[[iter]][["seed"]]))
  
  total.time <- system.time({
    results <- tryCatch({
      parLapply(cl, tasks, function(task) run_simulation(task))
    }, error = function(e) {
      print(paste("Error in parLapply:", e))
      NULL
    })
  })[3]
  
  
  # Save results
  save(total.time, results, file = sprintf("%s_arrayid_%d_n2_%d_cores_%d_iterations_%d.RData", job_name, array_id, n2, num_cores, m))
}

stopCluster(cl)
closeAllConnections()

