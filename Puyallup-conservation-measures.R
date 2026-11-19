library(tidyverse)
library(pegas)

loci <- 
  read.loci(file="../../Chinook/loci.txt") |>
  filter(population %in% c("White_R_H_DFW", "Puyallup_R_DFW")) |>
  as_tibble() |>
  mutate(population=as.character(population) |> gsub(pattern="\\w{4}$", replacement="") |> factor()) |>
  as.loci()

Fst(loci)[, 2] |> 
  as.double() |> 
  hist(breaks=30)
  
Fst(loci)[, 2] |> 
  as.double() |> 
  mean(na.rm=TRUE)

# loci$ots_isot <- NULL

# loci$ots_u07_25_325 <- NULL

set.seed(7)

B <- 
  slice_sample(loci, prop=0.5, by=population)

A <- 
  anti_join(x=loci, y=B)

Fst_A <- 
  Fst(A)[, 2] |> 
  as.double()

Fst_B <- 
  Fst(B)[, 2] |> 
  as.double()

mean(Fst_A, na.rm=TRUE)

mean(Fst_B, na.rm=TRUE)

hist(Fst_A, breaks=30, xlim=c(min(c(Fst_A, Fst_B), na.rm=TRUE)*1.2, max(c(Fst_A, Fst_B), na.rm=TRUE)*1.1), main="White R H vs. Puyallup R\n2021 (simulated)", sub=mean(Fst_A, na.rm=TRUE) |> round(digits=4), xlab=expression(italic(F)[ST]))

hist(Fst_B, breaks=30, xlim=c(min(c(Fst_A, Fst_B), na.rm=TRUE)*1.2, max(c(Fst_A, Fst_B), na.rm=TRUE)*1.1), main="White R H vs. Puyallup R\n2026 (simulated)", sub=mean(Fst_B, na.rm=TRUE) |> round(digits=4), xlab=expression(italic(F)[ST]))


mean(Fst_B, na.rm=TRUE) - 
  mean(Fst_A, na.rm=TRUE)


AB <- 
  c(Fst_A, Fst_B)


permute2samplet <- 
  function(x, n_loci){
    n <- n_loci*2
    idx_b <- sample(1:n, n_loci)
    idx_a <- setdiff(1:n, idx_b)
    mean_diff <- mean(x[idx_b], na.rm=TRUE) - mean(x[idx_a], na.rm=TRUE)
    return(mean_diff)}

permuteAB <- rep(NA, 1000)

set.seed(13)

for(i in 1:1000){permuteAB[i] <- permute2samplet(x=AB, n_loci=94)}

obsAB <-
  mean(Fst_B, na.rm=TRUE) - 
  mean(Fst_A, na.rm=TRUE)


# 2-tailed p-value
mean(abs(c(permuteAB, obsAB)) >= abs(obsAB))

hist(c(permuteAB, obsAB), xlab=expression(italic(F)[ST]~differences), breaks=30, main="Permutations")

abline(v=obsAB, lty=2, col="blue")

(c(permuteAB, obsAB) |> quantile(probs=c(0.025, 0.975), na.rm=TRUE))+mean(c(Fst_A, Fst_B), na.rm=TRUE)

set.seed(2)

C1 <- 
  filter(loci, population=="White_R_H") |> 
  slice_sample(prop=0.5) |> 
  as_tibble() |> 
  mutate(population=as.character(population))

C2 <- 
  anti_join(x=filter(loci, population=="White_R_H"), y=C1) |> 
  as_tibble() |> 
  mutate(population="Puyallup_R")

C <- 
  bind_rows(C1, C2) |> 
  mutate(population=factor(population)) |>
  as.loci()

Fst_C <- 
  Fst(C)[, 2] |> 
  as.double()

mean(Fst_A, na.rm=TRUE)

mean(Fst_B, na.rm=TRUE)

mean(Fst_C, na.rm=TRUE)

par(mfrow=c(3, 1))

hist(Fst_A, breaks=30, xlim=c(min(c(Fst_A, Fst_B, Fst_C), na.rm=TRUE)*1.2, max(c(Fst_A, Fst_B, Fst_C), na.rm=TRUE)*1.1), main="White R H vs. Puyallup R\n2021 (simulated)", sub=mean(Fst_A, na.rm=TRUE) |> round(digits=4), xlab=expression(italic(F)[ST]))

hist(Fst_B, breaks=20, xlim=c(min(c(Fst_A, Fst_B, Fst_C), na.rm=TRUE)*1.2, max(c(Fst_A, Fst_B, Fst_C), na.rm=TRUE)*1.1), main="White R H vs. Puyallup R\n2026 (simulated)", sub=mean(Fst_B, na.rm=TRUE) |> round(digits=4), xlab=expression(italic(F)[ST]))

hist(Fst_C, breaks=10, xlim=c(min(c(Fst_A, Fst_B, Fst_C), na.rm=TRUE)*1.2, max(c(Fst_A, Fst_B, Fst_C), na.rm=TRUE)*1.1), main="White R H vs. Puyallup R\n2031 (simulated)", sub=mean(Fst_C, na.rm=TRUE) |> round(digits=4), xlab=expression(italic(F)[ST]))


mean(Fst_C, na.rm=TRUE) - 
  mean(Fst_A, na.rm=TRUE)

BC <- 
  c(Fst_B, Fst_C)

permuteBC <- rep(NA, 1000)

set.seed(13)

for(i in 1:1000){permuteBC[i] <- permute2samplet(x=BC, n_loci=94)}

obsBC <-
  mean(Fst_C, na.rm=TRUE) - 
  mean(Fst_B, na.rm=TRUE)


# 2-tailed p-value
mean(abs(c(permuteBC, obsBC)) >= abs(obsBC))

par(mfrow=c(1, 1))

hist(c(permuteBC, obsBC), xlab=expression(italic(F)[ST]~differences), breaks=30, main=paste("p-value:", mean(abs(c(permuteBC, obsBC)) >= abs(obsBC)) |> round(digits=5)))

abline(v=obsBC, lty=2, col="blue")

