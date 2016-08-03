# Cluster shuffle: Determine best k for each value of n (number of agents) and generate plot

nvals = c(10000,50000,100000,250000,500000,750000,1000000)

bestk = c()
bestrate = c()
for (n in nvals) {
  filename = paste("output_var_C_", format(n, scientific=FALSE), ".csv", sep="")
  fr <- read.csv(filename)
  m <- max(fr$rate)
  L <- fr$rate == m
  fr <- fr[L,]
  k <- min(fr$k)
  bestk <- c(bestk, k)
  bestrate <- c(bestrate, m)
}
print("Best values of k for cluster shuffle")
print(data.frame(n=nvals,k=bestk,rate=bestrate))
bestkcluster = bestk # Save for comparison with Distribution match

options(scipen=2)
png('1_plotClusterShuffleBestK.png')
plot(nvals, bestk, type="b", main="Cluster shuffle",
      sub="Best number of neighbours for number of agents",
      ylab="Best value for k (neighbours)",  xlab="Number of agents (n)")
dev.off()

# Cluster shuffle: Determine best c for each value of n (number of agents) and generate plot

bestc = c()
bestrate = c()
for (n in nvals) {
  filename = paste("output_var_C_", format(n, scientific=FALSE), ".csv", sep="")
  fr <- read.csv(filename)
  m <- max(fr$rate)
  L <- fr$rate == m
  fr <- fr[L,]
  c <- min(fr$c)
  bestc <- c(bestc, c)
  bestrate <- c(bestrate, m)
}

print("Best values of c for cluster shuffle")
print(data.frame(n=nvals,c=bestc,rate=bestrate))
options(scipen=2)
png('2_plotClusterShuffleBestC.png')
plot(nvals, bestc, type="b", main="Cluster shuffle", sub="Best number of clusters for number of agents",
     ylab="Best value for c (clusters)",  xlab="Number of agents (n)")
dev.off()

# Cluster shuffle: Determine lowest c at best value of k and generate plot
l = length(nvals)
outfr = data.frame(k=1:l,c=1:l,rate=1:l,n=1:l)
i = 1
for (n in nvals) {
  filename = paste("output_var_C_", format(n, scientific=FALSE), ".csv", sep="")
  fr <- read.csv(filename)
  m <- max(fr$rate)
  L <- fr$rate == m
  fr <- fr[L,]
  k <- min(fr$k)
  L <- fr$k == k
  fr <- fr[L,]
  c = min(fr$c)
  row = c(k,c,n)
  outfr[i,] = c(k,c,rate=m,n)
  i <- i + 1
}
print("Best values of c for best values k for cluster shuffle")
print(outfr)

options(scipen=2)
png('3_plotClusterShuffleBestKC.png')
plot(nvals, outfr$k, type="b", main="Cluster shuffle", sub="Best number of clusters for number of agents",
     ylab="Best value for k (red) and c (green)",  xlab="Number of agents (n)", col="red")
lines(nvals,outfr$c,col="green", type="b")
dev.off()

# Distribution v Cluster: Determine best k for each value of n (number of agents) and generate plot

bestk = c()
bestrate = c()
for (n in nvals) {
  filename = paste("output_var_D_", format(n, scientific=FALSE), ".csv", sep="")
  fr <- read.csv(filename)
  m <- max(fr$rate)
  L <- fr$rate == m
  fr <- fr[L,]
  k <- min(fr$k)
  bestk <- c(bestk, k)
  bestrate <- c(bestrate, m)
}
print("Best values of k for distribution match")
print(data.frame(n=nvals,k=bestk,rate=bestrate))

options(scipen=2)
png('4_plotClusterShuffleVsDistributionBestK.png')
plot(nvals, bestk, type="b", main="Ideal number of neighbours",
     sub="Cluster shuffle (green) vs Distribution (red)",
     ylab="Best value for k (neighbours)",  xlab="Number of agents (n)", col="red")
lines(nvals,bestkcluster,col="green", type="b")
dev.off()

# Cluster shuffle: k,c, success

library(scatterplot3d)
fr <- read.csv("output_var_C_1000000.csv")
attach(fr)
k = fr$k
c = fr$c
rate = fr$rate * 100 # percent correct
png('5_plotClusterShuffleKCRate.png')
scatterplot3d(k,c,rate)
dev.off()
