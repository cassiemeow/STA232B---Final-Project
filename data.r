countn = c(15,1,
           6,1,2,
           6,6,
           7,2,3,2,
           16,9,3,3,1,
           57,38,17,2,2,
           119,81,45,6,1,1,
           173,118,57,16,3,1,
           136,103,50,13,6,1,1,
           54,51,32,5,1,1,
           13,15,12,3,1,
           4,3,1,
           1,1)

record = c(1,0,1,1,
           2,0,2,1,2,2,
           3,0,3,1,
           4,0,4,1,4,2,4,4,
           5,0,5,1,5,2,5,3,5,4,
           6,0,6,1,6,2,6,3,6,4,
           7,0,7,1,7,2,7,3,7,4,7,7,
           8,0,8,1,8,2,8,3,8,4,8,8,
           9,0,9,1,9,2,9,3,9,4,9,5,9,6,
           10,0,10,1,10,2,10,3,10,4,10,9,
           11,0,11,1,11,2,11,3,11,4,
           12,1,12,2,12,3,
           13,2,13,7)

record_m = matrix(record, ncol = 2, byrow = T)
data = matrix(ncol = 2)
for (i in 1:length(countn)) {
  data = rbind(data, matrix(rep(record_m[i,], countn[i]), ncol = 2, byrow = T))
}
data = cbind(data[-1,], Litter = 1:1328) %>% as.data.frame()
colnames(data)[1:2] = c("Nimplants","Ndead")

# MCEM
set.seed(232)
litter = 1328
mu = -2
sigma = 0.7
ite = 30
out = matrix(NA, ncol = 1, nrow = litter)
final_musig = matrix(NA, ncol = 3, nrow = ite)

for (i in 1:ite) {
  # E-step
  ## Metropolis-Hastings algorithm
  alpha_initial = rnorm(litter, 0, sigma)
  for (j in 1:1000) { # set a relatively large number for convergence
    alpha_update = rnorm(litter, 0, sigma)
    
    for (k in 1:litter) {
      a_new = exp(mu+alpha_update[k])^data[k,2]/(1+exp(mu+alpha_update[k]))^data[k,1]
      a_old = exp(mu+alpha_initial[k])^data[k,2]/(1+exp(mu+alpha_initial[k]))^data[k,1]
      alpha = pmin(1, a_new/a_old)
      
      u = runif(1)
      if (u < alpha) {alpha_initial[k] = alpha_update[k]}
    }
    out = cbind(out, alpha_initial)
  }
  out = out[,-(1:500)] # drop the non-converged values
  
  fr = function(x) {
    -mean(apply(out, 2, function(o) {sum((x+o)*data[,2]) -
        sum(log(1+exp(x+o))*data[,1]) # loss function
    }))}
  
  # M-step
  final_mu = NULL
  final_sigma = NULL
  final = optim(mu, fr, lower = -5, upper = 5, method = "L-BFGS-B")
  mu = final$par
  sigma = sqrt(mean(out^2))
  
  final_musig[i,] = cbind(i, mu, sigma)
  final_mu[i] = mu
  final_sigma[i] = sigma
}

final_musig = as.data.frame(final.musig)
colnames(final_musig) = c("iteration","mu","sigma")

ggplot(data=final_musig, aes(x = iteration, y = mu)) +
  geom_point() +
  xlab("Iteration") +
  ylab("MLE of mu") +
  theme_minimal()

ggplot(data=final_musig, aes(x = iteration, y = sigma)) +
  geom_point() +
  xlab("Iteration") +
  ylab("MLE of sigma") +
  theme_minimal()

#### SE
set.seed(232)
litter = 1328
mu = final_musig$mu[30]
sigma = final_musig$sigma[30]
out = matrix(NA, ncol = 1, nrow = litter)

  ## Metropolis-Hastings algorithm
  alpha_initial = rnorm(litter, 0, sigma)
  for (j in 1:1000) { # set a relatively large number for convergence
    alpha_update = rnorm(litter, 0, sigma)
    
    for (k in 1:litter) {
      a_new = exp(mu+alpha_update[k])^data[k,2]/(1+exp(mu+alpha_update[k]))^data[k,1]
      a_old = exp(mu+alpha_initial[k])^data[k,2]/(1+exp(mu+alpha_initial[k]))^data[k,1]
      alpha = pmin(1, a_new/a_old)
      
      u = runif(1)
      if (u < alpha) {alpha_initial[k] = alpha_update[k]}
    }
    out = cbind(out, alpha_initial)
  }
  out = out[,-(1:500)] # drop the non-converged values

I = matrix(0, 2, 2)  
for (i in 1:ncol(out)){
  alpha = out[,i]
  H = c(matrix(data[,1],nrow=1)%*%matrix(expit(mu+alpha)/(1+exp(mu+alpha)),ncol=1),
        3*alpha %*% alpha /(sigma^4)-1328/sigma^2)
  H = diag(H)
  S = c(sum(data[,2]) - data[,1]%*%expit(mu+alpha),
        alpha %*% alpha /sigma^3-1328/sigma)
  I = I + H - matrix(S, ncol = 1) %*% matrix(S, ncol = 2)
}
I  = I/ncol(out)
sqrt(solve(I)[1,1])
sqrt(solve(I)[2,2])

## Numerical Integration
s <- seq(1,25,1)
d <- length(s)
model1 <- glmer(cbind(Ndead, Nimplants-Ndead) ~ 1 + (1|litter), family = binomial, data = mice) #nAGQ = 1, laplace approximation
s_model1 <- summary(model1)
co <- s_model1$coefficients

mle <- matrix(0, nrow = d, ncol = 2)
a <- 0
mle_se <- matrix(0, nrow = d, ncol = 2)

for (i in s) {
  a <- a + 1
  model <- glmer(cbind(Ndead, Nimplants-Ndead) ~ 1 + (1|litter), family = binomial, nAGQ = i, data = mice) #take different values of nAGQ
  s_model <- summary(model)
  co <- s_model$coefficients
  mle[a,1] <- co[1]
  mle[a,2] <- sqrt(unlist(VarCorr(model)))
  asy = vcov(model, full=TRUE, ranpar="sd")
  mle_se[a,1] <- sqrt(asy[1,1])
  mle_se[a,2] <- sqrt(asy[2,2])
}

print(cbind(s,mle))
print(cbind(s,mle_se))
MLE_mu <- cbind(s,mle[,1],mle_se[,1])
colnames(MLE_mu) <- c('nAGQ', 'MLE_mu', 'StandardError_mu')
MLE_mu <- as.data.frame(MLE_mu)

MLE_mu %>%
ggplot(aes(nAGQ, MLE_mu)) +
geom_pointrange(aes(ymin = MLE_mu - StandardError_mu, ymax = MLE_mu + StandardError_mu))

MLE_sig <- cbind(s,mle[,2],mle_se[,2])
colnames(MLE_sig) <- c('nAGQ', 'MLE_sigma', 'StandardError_sigma')
MLE_sig <- as.data.frame(MLE_sig)

MLE_sig %>%
ggplot(aes(nAGQ, MLE_sigma)) +
geom_pointrange(aes(ymin = MLE_sigma - StandardError_sigma, ymax = MLE_sigma + StandardError_sigma))

## Bootstrap
mySumm <- function(mod) {
   c(sigma_e_square = getME(mod,"fixef"), sigma_s = sqrt(unlist(VarCorr(mod))))
}

booted <- bootMer(model, mySumm, nsim = 100, seed = 2047)
booted
