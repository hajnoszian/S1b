library(tidyverse)

#myPaths <- .libPaths()
#myPaths <- c(myPaths, "C:/Users/hajno/Documents/R/win-library/4.1")
#.libPaths(myPaths)

#####################

ids <-sort(rep(101:105, 2, length.out = 10))
ids

x <- rep(
  c(1,1,2,2,3,3),
  2, length.out = 10)
time <- rep(
  c(1,2),
  5, length.out = 10)
y <- round(rnorm(10, 4, 1),2)
y2 <- -1*(y1 + rnorm(10, 0.5,0.5))
cor(y1,y2)#cool, we got negatively correlated Ys. About what we'd expect for our real data
m <- round(
  rep(rnorm(5, 3, 1), 1,each = 2),
         2)

df <- data.frame(
  id = ids,
  time = time,
  x = x,
  m = m,
  y = y
)
df

#df <- data.frame(
#  id = ids,
#  time = time,
#  x = x,
#  m = m,
#  y1 = y1,
#  y2 = y2
#)