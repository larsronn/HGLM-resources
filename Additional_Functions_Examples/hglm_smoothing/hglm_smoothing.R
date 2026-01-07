library(ggplot2)
library(hglm)
#Simulated data
set.seed(12345)
n = 100
mu <- rep(1,n)
x <- 1:n
mu[(n/2+1):n] <- sin((1:(n/2))/pi*2)*(1:(n/2))

res.var <- exp(5-x/n)
err <- rnorm(n, 0, sqrt(res.var))
y <- mu + err

#Scale the data
y <- as.numeric(scale(y))
x <- as.numeric(scale(x))

#Code to plot the data and simple smoothing functions
dfsim <- data.frame(x=x, y=mu, residual.variance=res.var)
p1 <- ggplot(dfsim , aes(x=x, y=y))+
  geom_line()
p1 + ggtitle("Simulated function") +
  xlab("x") + ylab("y")
p2 <- ggplot(dfsim , aes(x=x, y=residual.variance))+
  geom_line()
p2 + ggtitle("Simulated Residual Variance") +
  xlab("x") + ylab("Residual variance")
sim.data <- data.frame(x=x, y=y)
p.data <- ggplot(sim.data, aes(x=x, y=y))+
  geom_point()

df <- data.frame(x=x, y=y)
p3 <- ggplot(df, aes(x=x, y=y)) +
  geom_point()
p4 <- ggplot(df, aes(x=x, y=y)) +
  geom_point()+
  geom_path()
p3b <- ggplot(df, aes(x=x, y=y)) +
  geom_point() + 
  geom_smooth() + 
  ggtitle("Default smoothing ggplot2")
###
#Here the computations start to get the smoothed curve assuming variable signal and noise
#The code should work on any data x and y
get.matern <-function(d, l, nu) {
  if (nu == 0.5) c <- exp(-d*l^(-1))
  if (nu == 1.5) c <- (1+sqrt(3)*d*l^(-1))*exp(-sqrt(3)*d*l^(-1))
  if (nu == 2.5) c <- (1+(sqrt(5)*d)*l^(-1)+(5*d^2)*3^(-1)*l^(-2))*exp(-sqrt(5)*d*l^(-1))
  return(c)
}
##
#see Ranjan et al (2011) Technometrics, 53:4, 366-378
bend <- function(C, OK.cond = 5) {
  eigVal <- eigen(C)$values
  #plot(-log10(eigVal))
  cond0 <- log10(max(eigVal)/min(eigVal)) 
  delta = max(eigVal)*(10^cond0-10^OK.cond)/(10^cond0*(10^OK.cond-1))
  C1 <- C
  diag(C1) = diag(C) + delta
  return(C1)
}
###
distMat <- as.matrix(dist(x)) #Compute the distance matrix

min.dist <- min(distMat+diag(dim(distMat)[1])*1e8)
if (min.dist < 1e-3) print("Consider clustering points by their means and use weights option")
n.rho = 40
like_results <- matrix(NA, n.rho, 4)
rho.seq <- seq(1e-2, 2*sd(x), length.out = n.rho)
maxiter=50
for(i in 1:n.rho) {
  rho = rho.seq[i]
  C1 <- get.matern(distMat, rho, 1.5) #Reasonable values for the Matern covariance function used
  BENDING = TRUE
  if (BENDING)     C1 <- bend(C1)
  ##### Fit the model in hglm ########
  library(hglm)
  n= length(y)
  X <- matrix(1, n, 1) #Design matrix for the intercept
  L <- t(chol(C1))    #Include the covariance in Z
  X.d <- model.matrix(~x)
  if (!is.null(mod2)) if(class(mod2) != "try-error") startval <- c(mod2$fixef, mod2$ranef, mod2$varRanef, mean(mod2$phi))
  if (is.null(mod2)) startval = NULL
  mod2 <- try( hglm(y = y, X = X, Z = L, calc.like = TRUE, maxit=maxiter), silent=TRUE)
  if (class(mod2) != "try-error") {
    if (mod2$Converge != "did not converge") {
      like_results[i,1:4] <- unlist(mod2$likelihood)
      curve2 <- data.frame(x=x, y = mod2$fv)
      p.mod2  <- p.data + geom_line(data=curve2)
      p2 <- p.mod2 + ggtitle("Model fitting variable signal and noise") +
        xlab("x") + ylab("y")
      plot(p2)
    }
  }
  print(i)
}

like_results
which.max(like_results[,3])


rho = rho.seq[which.max(like_results[,3])]
C1 <- get.matern(distMat, rho, 1.5) #Reasonable values for the Matern covariance function used
if (BENDING)     C1 <- bend(C1)
##### Fit the final hglm model and plot the results ########
n= length(y)
X <- matrix(1, n, 1) #Design matrix for the intercept
L <- t(chol(C1))    #Include the covariance in Z
X.d <- model.matrix(~x)
mod2.final <- hglm(y = y, X = X, Z = L, X.disp = X.d, X.rand.disp = list(X.d))

curve2 <- data.frame(x=x, y = mod2.final$fv)
p.mod2  <- p.data + geom_line(data=curve2)
p2 <- p.mod2 + ggtitle("Model fitting variable signal and noise") +
  xlab("x") + ylab("y")
plot(p2)


