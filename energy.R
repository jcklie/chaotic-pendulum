#!/bin/env R

data <- read.table("double_pendulum.csv", header=T,sep=";")

plot(x=data[ , "t"], y=data[ , "energy"], xlab="t [s]", ylab="E [J]")
