library(hglm)
DHGLM <- function(meanmodel_fixed=y~1, meanmodel_random=~id, dispersionmodel_fixed=~1, dispersionmodel_random=~id, data.in=NULL, meanmodel_RandC=NULL, dispersionmodel_RandC=NULL) {
  meanmodel_random <- as.formula(paste0("~0+",as.character(meanmodel_random)[2]))
  dispersionmodel_random <- as.formula(paste0("~0+",as.character(dispersionmodel_random)[2]))
  X <- model.matrix(meanmodel_fixed, data=data.in)
  Z <- model.matrix(meanmodel_random, data=data.in) #This might not work for multiple random effects, e.g. id+Box?? Might need to create separate Z for the different effects
  RandC1 <- ifelse(is.null(meanmodel_RandC), ncol(Z), meanmodel_RandC)
  X.disp <- model.matrix(dispersionmodel_fixed, data=data.in)
  Z.disp <- model.matrix(dispersionmodel_random, data=data.in)
  RandC2 <- ifelse(is.null(dispersionmodel_RandC), ncol(Z.disp), dispersionmodel_RandC)
  #This is the DHGLM fitting algorithm
  wghts1= rep(1,n)
  for (i in 1:10) {
    mean_model <- hglm(y=y, X=X, Z=Z, weights=wghts1, RandC=RandC1)
    y_disp <- (y-mean_model$fv)^2/(1-mean_model$hv[1:n])
    wghts2 <- (1-mean_model$hv[1:n])/2
    disp_model <- hglm(y=y_disp, X=X.disp, Z=Z.disp, weights=wghts2, family=Gamma(link=log), RandC=RandC2)
    wghts1 <- 1/disp_model$fv
  }
  
  cat("The final weighted residual variance is:",mean_model$varFix, " and should be 1.0 at convergence", "\n")
  print("Summary for the mean model")
  print(summary(mean_model))
  print("Summary for the dispersion model")
  print(summary(disp_model))
  list(mean_model=mean_model, disp_model=disp_model)
}


#simulate an example with random effects both in mean and dispersion
set.seed(1234)
k = 30 #number of individuals
m = 8 #number of observations per individual
n = m*k
mu = 10
c0 = 0
sigma_mean=sqrt(0.64)
sigma_disp=sqrt(0.25)
id <- factor(rep(1:k, each=m))
Z <- diag(k)%x%rep(1,m)
u <- rnorm(k, 0, sigma_mean)
v <- rnorm(k, 0, sigma_disp)
r <- rnorm(n, 0, sqrt(exp(c0+Z%*%v)))
y <- mu + Z%*%u + r
#This ends the simulation part
data1 <- data.frame(y=y, id=id)
fitted.model1 <- DHGLM(meanmodel_fixed=y~1, meanmodel_random=~id, dispersionmodel_fixed=~1, dispersionmodel_random=~id, data.in = data1)
