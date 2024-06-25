##################################
## Lagrange and Genetic ##########
##################################


library(parallel)
library(TPOrd)

job_name <- Sys.getenv("JOB_NAME", "propoddsdata_R5678")


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



m <- as.numeric(Sys.getenv("ITERATIONS", "128")) # Number of iterations
mit <- as.numeric(Sys.getenv("MIT", "1"))
num_cores = as.numeric(Sys.getenv("NUM_CORES", "32"))



# Adjusting parameters using environment variables
target_cor <- as.numeric(Sys.getenv("TARGET_COR", "0.75"))
fractions_cat1 <- as.numeric(strsplit(Sys.getenv("FRACTIONS_CAT1", "0.35,0.25,0.2"), ",")[[1]])
fractions_cat2 <- as.numeric(strsplit(Sys.getenv("FRACTIONS_CAT2", "0.30,0.25,0.25,0.1"), ",")[[1]])
N <- as.numeric(Sys.getenv("N", "5000"))
cor_YZ_break <- as.numeric(Sys.getenv("COR_YZ_BREAK", "0.85"))
num_categories <- as.numeric(Sys.getenv("NUM_CATEGORIES", "5"))
n2_vector <- as.numeric(strsplit(Sys.getenv("N2_VECTOR", "500,1000,1500,2000"), ",")[[1]])
model_type <- "Adjacent_Category"
Beta0 <- as.numeric(strsplit(Sys.getenv("BETA0", "-1.10580581,-0.510929214,2.260466119,2.33807856"), ",")[[1]])
Beta1 <- as.numeric(Sys.getenv("BETA1", "-0.1238386604"))
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






  dat_sim <- sim_categorical_data (Beta0 = Beta0 ,# intercept
                                   Beta1 = Beta1, # effect size of the sequencing variant
                                   fractions_cat1 =fractions_cat1,
                                   fractions_cat2 =fractions_cat2,
                                   target_cor = target_cor,
                                   N = N, # phase 1  sample size
                                   n2 = n2, # phase 2 sample size
                                   cor_YZ_break = cor_YZ_break,
                                   num_categories = num_categories,
                                   model_type = model_type)



  #######################################################
  #######################################################

  ### Data prep for GA ----
  dat_sim$YZ <- interaction(dat_sim$Y, dat_sim$Z)
  dat_sim$fZ <- factor(dat_sim$Z)
  dat_list <- TPOrd:::sample_dataframe_balanced_matrix(dat_sim, n2, 1000)

  dat = dat_sim[ , c("Y","Z","fZ","YZ")]
  dat$YZ <- as.numeric(dat$YZ)


  fit =  vglm(Y ~ fZ, family = fam, data = dat)
  fitcoef = c(fit@coefficients[1: (num_categories-1)], 0, fit@coefficients[num_categories: length(fit@coefficients)])





  pZ <- table(dat$Z)[1:(length(table(dat$Z))-1)]/N ## pZ
  ### summary(fit)


  Kind <- 1
  formula_Ha <- Y~G+fZ
  miscov_formula <- ~G


  pG <- fractions_cat1
  resultsgz <- simulateGZ(N, pG, pZ, target_cor, num_sim=3000, plot_results = T)

  p_gz0 <- convertMeanFrequenciesToTable(resultsgz$MeanFrequencies)


  ### The formula under the null is the same regardless
  formula_Ho <- Y~fZ
  auxvar = ~Z
  strataformula_YZ = ~YZ
  strataformula_Y = ~Y
  strataformula_Z = ~Z

  optMethod <- "Par-spec" # c("Par-spec","A-opt","D-opt")[1]

  ### Obtain GA design ----
  IM <- TPOrd:::obsIM_R_ord(formula=formula_Ha,
                            miscov=miscov_formula, auxvar=~Z,
                            model_type=model_type, dat,
                            beta=fitcoef,p_gz=p_gz0,num_categories=num_categories)



  indices_matrix_apply <- function(binary_matrix) {
    num_ones <- sum(binary_matrix[, 1])  # Assuming all columns have the same number of `1`s
    apply(binary_matrix, 2, function(col) which(col == 1))
  }

  pop1a_Y <- indices_matrix_apply(dat_list$R2)
  pop1a_YZ <- indices_matrix_apply(dat_list$R3)
  pop1a_Z <- indices_matrix_apply(dat_list$R4)


  fit1_Y <- apply(pop1a_Y,2,function(r){
    return(TPOrd:::fitnessTP_ord(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
  })

  fit1_YZ <- apply(pop1a_YZ,2,function(r){
    return(TPOrd:::fitnessTP_ord(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
  })

  fit1_Z <- apply(pop1a_Z,2,function(r){
    return(TPOrd:::fitnessTP_ord(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
  })





  pop2a <- sapply(1:1000,function(x){
    sa <- sample(N,n2)
  })
  fit2 <- apply(pop2a,2,function(r){
    return(TPOrd:::fitnessTP_ord(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
  })



  ### Put initial population (mix between LM and best candidates from SRS)
  pop1_Y <- pop1a_Y[,order(fit1_Y)[1:40]]
  pop1_YZ <- pop1a_YZ[,order(fit1_YZ)[1:40]]
  pop1_Z <- pop1a_Z[,order(fit1_Z)[1:40]]

  pop2 <- pop2a[,order(fit2)[1:20]]



  ### It may take a few minutes to run

  GA.sol1_Y <- TPOrd:::optimTP_GA_ord(ncores=1,
                                      formula=formula_Ha,
                                      miscov=miscov_formula,
                                      auxvar=~Z,
                                      model_type=model_type,
                                      n=n2,
                                      dat,
                                      beta=fitcoef,p_gz=p_gz0,
                                      ga.popsize=60,
                                      ga.propelit=0.9,
                                      ga.proptourney=0.9,
                                      ga.ngen=200, # can change to 500
                                      ga.mutrate=0.001,
                                      ga.initpop=t(cbind(pop1_Y,pop2)),
                                      optimMeasure=optMethod,K.idx=Kind,seed=1)


  GA.sol1_YZ <- TPOrd:::optimTP_GA_ord(ncores=1,
                                       formula=formula_Ha,
                                       miscov=miscov_formula,
                                       auxvar=~Z,
                                       model_type=model_type,
                                       n=n2,
                                       dat,
                                       beta=fitcoef,p_gz=p_gz0,
                                       ga.popsize=60,
                                       ga.propelit=0.9,
                                       ga.proptourney=0.9,
                                       ga.ngen=200, # can change to 500
                                       ga.mutrate=0.001,
                                       ga.initpop=t(cbind(pop1_YZ,pop2)),
                                       optimMeasure=optMethod,K.idx=Kind,seed=1)


  GA.sol1_Z <- TPOrd:::optimTP_GA_ord(ncores=1,
                                      formula=formula_Ha,
                                      miscov=miscov_formula,
                                      auxvar=~Z,
                                      model_type=model_type,
                                      n=n2,
                                      dat,
                                      beta=fitcoef,p_gz=p_gz0,
                                      ga.popsize=60,
                                      ga.propelit=0.9,
                                      ga.proptourney=0.9,
                                      ga.ngen=200, # can change to 500
                                      ga.mutrate=0.001,
                                      ga.initpop= t(cbind(pop1_Z,pop2)),
                                      optimMeasure=optMethod,K.idx=Kind,seed=1)


  ### Obtain the indicators
  dat_sim$R5 <- 1*(1:N %in% pop2a[,which(fit2 == min(fit2) )])
  dat_sim$R6 <- rep(0,N); dat_sim$R6[GA.sol1_Y$bestsol] <- 1
  dat_sim$R7 <- rep(0,N); dat_sim$R7[GA.sol1_YZ$bestsol] <- 1
  dat_sim$R8 <- rep(0,N); dat_sim$R8[GA.sol1_Z$bestsol] <- 1










  com_table <- table(dat_sim$Y, dat_sim$G1)
  R5_table <- table(dat_sim[dat_sim$R5==1,]$Y,dat_sim[dat_sim$R5==1,]$G1)
  R6_table <- table(dat_sim[dat_sim$R6==1,]$Y,dat_sim[dat_sim$R6==1,]$G1)
  R7_table <- table(dat_sim[dat_sim$R7==1,]$Y,dat_sim[dat_sim$R7==1,]$G1)
  R8_table <- table(dat_sim[dat_sim$R8==1,]$Y,dat_sim[dat_sim$R8==1,]$G1)




  ###############


  data_list <- list(iteration= iteration,
                    seed = my_seed,
                    sample_freq=list(com_table=com_table,  R5_table=R5_table, R6_table=R6_table, R7_table=R7_table, R8_table=R8_table),
                    dataset = dat_sim
  )



  return(data_list)

}

closeAllConnections()


cl <- makeCluster(num_cores)
clusterExport(cl, varlist = c(
  "array_id", "m", "num_cores",
  "target_cor", "fractions_cat1", "fractions_cat2",
  "N", "cor_YZ_break", "num_categories",
  "model_type", "Beta0", "Beta1",
  "fam", "run_simulation"
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
  save(total.time, results, file = sprintf("%s_arrayid_%d_n2_%d_cores_%d_iterations_%d.RData", job_name, array_id, n2, num_cores, m))
}

stopCluster(cl)
closeAllConnections()

