library(rstan)

#initial parameter function
initf<-function(){
  list(N0=12.,mu=10.,sig=6.,pmr=8.,
       #gammaR=6.5,gammaT=6.5,gammaS=6.5,
       ec50R=20.5,ec50T=20.5,ec50S=20.5,
       emaxR=99.,emaxT=99.0,emaxS=99.99,SD=1.5)
}

#load fake data
source("bigcounts16.R")

#fit the model to the data
fit <- stan(file = 'SompobFixedRing.stan', data = bigcounts[[1]], iter = 100000, chains = 1, init = initf, algorithm = "NUTS");
saveRDS(fit, file = paste0("fit", ids[[1]], ".rds"))

#extract fitted data
source("extractdata4plotting.R")
FittedData("fitP0XX.rds")

#plot the output
source("plotgraphs.R")
csvfiles<-list.files(pattern = "fitP.*csv")
plot.counts(csvfiles[[1]],id=substr(csvfiles[[1]],4,7),split = TRUE)






