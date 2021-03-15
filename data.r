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
out = matrix(NA, ncol = 1, nrow = litter)

for (i in 1:30) {
  # E-step
  ## Metropolis-Hastings algorithm
  alpha_initial = rnorm(litter, mean=0, sd=sigma)
  for (k in 1:1000) { # set a relatively large number for convergence
    alpha_update = rnorm(litter, mean=0, sd=sigma)
    
    for (j in 1:litter) {
      a_new = exp(mu+alpha_update[j])^data[j,2]/(1+exp(mu+alpha_update[j]))^data[j,1]
      a_old = exp(mu+alpha_initial[j])^data[j,2]/(1+exp(mu+alpha_initial[j]))^data[j,1]
      alpha = min(1, a_new/a_old)
      
      u = runif(1)
      if (u < alpha) {alpha_initial[j] = alpha_update[j]}
    }
    out = cbind(out, alpha_initial)
  }
  out = out[,-(1:500)] # drop the non-converged values
  
  fr = function(x) {
    -mean(apply(out,2, function(i) {sum(x*data[,2]+data[,2]*i) -
        sum(data[,1]*log(1+exp(x+i)))
    }))}
  
  # M-step
  final_mu = NULL
  final_sigma = NULL
    final = optim(mu, fr, lower = -5, upper = 5, method = "L-BFGS-B")
    mu = final$par
    sigma = sqrt(mean(out^2))
    
    print(paste(i, mu, sigma))
    final_mu[i] = mu
    final_sigma[i] = sigma
}

final_result = c(final_mu[30], final_sigma[30])

mu_result = c(-2.201063438,-2.258719093,-2.268630974,-2.273614399,-2.274986966,-2.275624713,-2.275585099
,-2.275424387,-2.275470593,-2.275535516,-2.275852576,-2.274899202,-2.274859851,-2.274843385,-2.274928345
,-2.274776375,-2.274771227,-2.274588585,-2.274358027,-2.274551402,-2.274416098,-2.274216286,-2.274105974
,-2.274019603,-2.27403727,-2.273997129,-2.273955056,-2.273858039,-2.27388587,-2.273768701)

sigma_result = c(0.687109262,0.680981447,0.678560417,0.675848428,0.675100617,0.67434348,0.673529551
,0.673088305,0.672850585,0.672640734,0.672808548,0.671775368,0.67154516,0.671387591,0.671178017
,0.670688836,0.670425547,0.670056778,0.669656354,0.669726186,0.669533155,0.669422091,0.669232574,
0.668945662,0.668774407,0.668561676,0.668452382,0.668288421,0.66826806,0.668078129)

ok = cbind(seq(1:30), mu_result, sigma_result) %>% data.frame()

# save(final_result,file="/Users/xuchenghuiyun/Desktop/final.Rdata")
# load(file="/Users/xuchenghuiyun/Desktop/final.Rdata")

set.seed(232)
litter = 1328
######################
mu = -2.276878
sigma = 0.6742865
######################
out = matrix(NA, ncol = 1, nrow = litter)

  ## Metropolis-Hastings algorithm
  alpha_initial = rnorm(litter, 0, sigma)
  for (k in 1:1500) { # set a relatively large number for convergence
    alpha_update = rnorm(litter, 0, sigma)
    
    for (j in 1:litter) {
      p_new = exp(mu+alpha_update[j])^data[j,2]/(1+exp(mu+alpha_update[j]))^data[j,1]
      p_old = exp(mu+alpha_initial[j])^data[j,2]/(1+exp(mu+alpha_initial[j]))^data[j,1]
      alpha = min(1, p_new/p_old)
      
      u = runif(1)
      if (u < alpha) {alpha_initial[j] = alpha_update[j]}
    }
    out = cbind(out, alpha_initial)
  }
  out = out[,-(1:1000)] # drop the non-converged values

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
