library(tidyverse)

gdatw <- 
  readRDS("../admixture/gdatw2.RDS")

# remove bad SNP
gdatw <- 
  select(gdatw, -ots_u07_07_161_a1, -ots_u07_07_161_a2)

write_delim(
  mutate(gdatw, ID=paste0("i", ID)), 
  col_names=FALSE, 
  file="../admixture/gdatw3.ped", 
  na="0")

length(names(gdatw)[7:508])/2

data.frame(Chromosome=1, SNP=paste0("rs", 10001:10251), GD=0, POS=seq(1, 251000000, 1000000)) |> 
  write_delim(col_names=FALSE, file="../admixture/gdatw3.map")

# bash
#
# conda activate admixture_Env
# cd /Users/tylermccraney/Google\ Drive/My\ Drive/Puyallup-White/admixture
# plink --file gdatw3 --make-bed --out gdatw3
# for K in 1 2 3 4; do admixture -s $RANDOM --cv=10 gdatw3.bed $K | tee gdatw3_log${K}.out; done
# grep -h CV gdatw3_log*  > gdatw3_log.txt



gdatw_cv <- 
  read_delim("../admixture/gdatw3_log.txt", col_names=FALSE) %>% 
  mutate(K=1:nrow(.)) |> 
  select(K, CV_error=X4)

plot(gdatw_cv)
lines(gdatw_cv)

tbl2 <- read.table("../admixture/gdatw3.2.Q")

tbl3 <- read.table("../admixture/gdatw3.3.Q")

tbl2_252 <- read.table("../admixture/gdatw2.2.Q")

unique(gdatw$pop)

gdatw_pop <- 
  gsub(pattern="River", replacement="R", gdatw$pop) |> 
  gsub(pattern="Creek", replacement="C") |> 
  gsub(pattern="Hatchery", replacement="H") |> 
  gsub(pattern="Chinook_", replacement="") |> 
  gsub(pattern="Buckley_Trap_", replacement="") |> 
  gsub(pattern="at_", replacement="_") |> 
  gsub(pattern="unclipped", replacement="Fall") |>
  gsub(pattern="_PBT", replacement="Fall") |>
  gsub(pattern="WDFW_FH", replacement="Fall") |>
  gsub(pattern="Spring_Fall_06KK", replacement="Spring_06KK_Baseline") |>
  gsub(pattern="R_24ND", replacement="R_Spring_24ND_Baseline") |>
  gsub(pattern="R_25NI", replacement="R_Spring_25NI_Baseline") |>
  gsub(pattern="fall", replacement="Fall") |>
  gsub(pattern="08HZ", replacement="08HZ_Baseline") |>
  gsub(pattern="98HB", replacement="98HB_Baseline") |>
  gsub(pattern="10JX", replacement="10JX_Baseline") |>
  gsub(pattern="12IP", replacement="12IP_Baseline") |>
  gsub(pattern="__", replacement="_") |>
  gsub(pattern="Hupp_Springs_23LQ", replacement="Hupp_S_Spring_23LQ_Baseline") |>
  gsub(pattern="Minter_C_23LR", replacement="Minter_C_Fall_23LR_Baseline")



unique(gdatw_pop)


rownames(tbl2) <- 
  gdatw$ID

rownames(tbl2_252) <- 
  gdatw$ID

rownames(tbl3) <- 
  gdatw$ID

tbl2$pop <- 
  gdatw_pop

tbl2_252$pop <- 
  gdatw_pop

tbl3$pop <- 
  gdatw_pop

tbl2 <- 
  bind_rows(tbl2, tibble(V1=rep(NA, 4), V2=rep(NA, 4), pop=rep(NA, 4)))

tbl2_252 <- 
  bind_rows(tbl2_252, tibble(V1=rep(NA, 4), V2=rep(NA, 4), pop=rep(NA, 4)))

tbl3 <- 
  bind_rows(tbl3, tibble(V1=rep(NA, 4), V2=rep(NA, 4), V3=rep(NA, 4), pop=rep(NA, 4)))

rownames(tbl2)[1317:1320] <- 
  paste0("spacer", 1:4)

rownames(tbl2_252)[1317:1320] <- 
  paste0("spacer", 1:4)

rownames(tbl3)[1317:1320] <- 
  paste0("spacer", 1:4)

rownames(tbl2)[1:1000]; rownames(tbl2)[1001:1320]
rownames(tbl3)[1:1000]; rownames(tbl3)[1001:1320]


tbl2 <- 
  tbl2[rownames(tbl2)[c(1:225, 1317, 226:966, 1318, 967:1002, 1319, 1003:1060, 1320, 1061:1316)], ]

tbl2_252 <- 
  tbl2_252[rownames(tbl2_252)[c(1:225, 1317, 226:966, 1318, 967:1002, 1319, 1003:1060, 1320, 1061:1316)], ]

tbl3 <- 
  tbl3[rownames(tbl3)[c(1:225, 1317, 226:966, 1318, 967:1002, 1319, 1003:1060, 1320, 1061:1316)], ]


pdf("../admixture/gdatw3_k2_3.pdf", width=13, height=6)
par(mar=c(7, 4, 4, 2) + 0.1, mfrow=c(2, 1))
barplot(t(as.matrix(tbl2[-3])), col=terrain.colors(4)[c(1, 3)], ylab="Ancestry", border=NA, names.arg=tbl2$pop, las=2, cex.names=0.3, space=0, main=expression(italic(K)==2*`*`), font.main=1)
abline(v=c(226, 968, 1005, 1064), col="white")
mtext(text="Baseline Spring", adj=0.08)
mtext(text="'Spring'")
mtext(text="?", adj=0.73)
mtext(text="'Fall'", adj=0.77)
mtext(text="Baseline Fall", adj=0.90)
barplot(t(as.matrix(tbl3[c(3, 1, 2)])), col=terrain.colors(4)[c(1, 2, 3)], ylab="Ancestry", border=NA, xaxt="n", las=2, cex.names=0.3, space=0, main=expression(italic(K)==3), font.main=1)
abline(v=c(226, 968, 1005, 1064), col="white")
mtext(text="Baseline Spring", adj=0.08)
mtext(text="'Spring'")
mtext(text="?", adj=0.73)
mtext(text="'Fall'", adj=0.77)
mtext(text="Baseline Fall", adj=0.90)
dev.off()


pdf("../admixture/gdatw2_k2_gdatw3_k2.pdf", width=13, height=6)
par(mar=c(5, 4, 4, 2) + 0.1, mfrow=c(2, 1))
barplot(t(as.matrix(tbl2_252[-3])), col=terrain.colors(4)[c(1, 3)], ylab="Ancestry", border=NA, xaxt="n", las=2, cex.names=0.3, space=0, main="252 loci (include Ots_u07-07.161)", font.main=1)
abline(v=c(226, 968, 1005, 1064), col="white")
mtext(text="Baseline Spring", adj=0.08)
mtext(text="'Spring'")
mtext(text="?", adj=0.73)
mtext(text="'Fall'", adj=0.77)
mtext(text="Baseline Fall", adj=0.90)
barplot(t(as.matrix(tbl2[-3])), col=terrain.colors(4)[c(1, 3)], ylab="Ancestry", border=NA, xaxt="n", las=2, cex.names=0.3, space=0, main="251 loci (remove Ots_u07-07.161)", font.main=1)
abline(v=c(226, 968, 1005, 1064), col="white")
mtext(text="Baseline Spring", adj=0.08)
mtext(text="'Spring'")
mtext(text="?", adj=0.73)
mtext(text="'Fall'", adj=0.77)
mtext(text="Baseline Fall", adj=0.90)
dev.off()

tbl2s <- 
  arrange(tbl2, desc(V1)) |>
  mutate(
    ASSIGNMENT=ifelse(grepl("_Spring_", pop), yes="SPRING", no=NA), 
    ASSIGNMENT=ifelse(grepl("_Unknown_", pop), yes="UNKNOWN", no=ASSIGNMENT), 
    ASSIGNMENT=ifelse(grepl("_Fall_", pop), yes="FALL", no=ASSIGNMENT), 
    BASELINE=ifelse(grepl("_Baseline", pop), yes=TRUE, no=FALSE),
    Q=cut(V1, breaks=c(0, 0.4, 0.6, 1)))

tbl2s$Q <- 
  factor(
    tbl2s$Q, 
    levels=c("(0.6,1]", "(0.4,0.6]", "(0,0.4]"), 
    labels=c("SPRING", "HYBRID", "FALL"), 
    ordered=TRUE)

tbl2s$ASSIGNMENT <- 
  factor(
    tbl2s$ASSIGNMENT, 
    levels=c("SPRING", "UNKNOWN", "FALL"), 
    ordered=TRUE)

tibl2s <- 
  as_tibble(tbl2s) |>
  filter(!is.na(pop)) |>
  mutate(collection=factor(pop), year=factor(pop))

tibl2s$collection <- 
  factor(
    tibl2s$collection, 
    levels=c(
      "Green_R_Fall_98HB_Baseline", "White_R_Spring_06KK_Baseline", 
      "Voights_C_H_Fall_08HZ_Baseline", "Soos_C_H_Fall_10JX_Baseline", 
      "Green_R_Soos_C_H_Fall_12IP_Baseline", 
      "White_R_Spring_21LC", "White_R_Unknown_21LC", "White_R_Fall_21LC", 
      "White_R_Spring_22LL", "White_R_Unknown_22LL", "White_R_Fall_22LL", 
      "White_R_Spring_23LI", "White_R_Unknown_23LI", "White_R_Fall_23LI", 
      "Hupp_S_Spring_23LQ_Baseline", "Minter_C_Fall_23LR_Baseline", 
      "White_R_Spring_24ND", "White_R_Unknown_24ND", "White_R_Fall_24ND", 
      "White_R_Spring_24ND_Baseline", 
      "White_R_Spring_25NI", "White_R_Unknown_25NI", "White_R_Fall_25NI", 
      "White_R_Spring_25NI_Baseline", 
      "White_R_Spring_26ON", "White_R_Unknown_26ON", "White_R_Fall_26ON"), 
    ordered=TRUE)
  
tibl2s$year <- 
  factor(
    tibl2s$collection, 
    levels=c(
      "Green_R_Fall_98HB_Baseline", "White_R_Spring_06KK_Baseline", 
      "Voights_C_H_Fall_08HZ_Baseline", "Soos_C_H_Fall_10JX_Baseline", 
      "Green_R_Soos_C_H_Fall_12IP_Baseline", 
      "White_R_Spring_21LC", "White_R_Unknown_21LC", "White_R_Fall_21LC", 
      "White_R_Spring_22LL", "White_R_Unknown_22LL", "White_R_Fall_22LL", 
      "White_R_Spring_23LI", "White_R_Unknown_23LI", "White_R_Fall_23LI", 
      "Hupp_S_Spring_23LQ_Baseline", "Minter_C_Fall_23LR_Baseline", 
      "White_R_Spring_24ND", "White_R_Unknown_24ND", "White_R_Fall_24ND", 
      "White_R_Spring_24ND_Baseline", 
      "White_R_Spring_25NI", "White_R_Unknown_25NI", "White_R_Fall_25NI", 
      "White_R_Spring_25NI_Baseline", 
      "White_R_Spring_26ON", "White_R_Unknown_26ON", "White_R_Fall_26ON"), 
    labels=c(
      1998, 2006, 2008, 2010, 2012, 
      2021, 2021, 2021, 
      2022, 2022, 2022, 
      2023, 2023, 2023, 2023, 2023, 
      2024, 2024, 2024, 2024, 
      2025, 2025, 2025, 2025, 
      2026, 2026, 2026)) |> 
  as.character() |> 
  as.numeric()


table(tibl2s$collection, tibl2s$BASELINE, useNA="ifany")

table(tibl2s$year, tibl2s$Q, useNA="ifany")

table(tibl2s$ASSIGNMENT, tibl2s$Q, useNA="ifany")
table(tibl2s$collection, tibl2s$Q, useNA="ifany")

table(filter(tibl2s, BASELINE==TRUE)$ASSIGNMENT, filter(tibl2s, BASELINE==TRUE)$Q, useNA="ifany")[-2, ]

table(filter(tibl2s, BASELINE==FALSE)$ASSIGNMENT, filter(tibl2s, BASELINE==FALSE)$Q, useNA="ifany")

table(filter(tibl2s, BASELINE==TRUE)$collection |> fct_drop(), filter(tibl2s, BASELINE==TRUE)$Q, useNA="ifany")[-2, ]

table(filter(tibl2s, BASELINE==FALSE)$collection |> fct_drop(), filter(tibl2s, BASELINE==FALSE)$Q, useNA="ifany")

saveRDS(tibl2s, file="../admixture/3-tibl2s.RDS")


pdf("../admixture/gdatw3_k2_sorted.pdf", width=13, height=3)
par(mar=c(7, 4, 4, 2) + 0.1)
barplot(t(as.matrix(tbl2s[1:2])), col=terrain.colors(4)[c(1, 3)], ylab="Ancestry", border=NA, names.arg=tbl2s$pop, las=2, cex.names=0.3, space=0, main=expression(italic(K)==2*`*`), font.main=1)
# abline(h=0.5, v=1026, col="white", lty=3)
abline(h=c(0.4, 0.6), v=c(952.5, 1062.5), col="white", lty=1)
mtext(text="'Spring'", adj=0.37)
mtext(text="'Hybrid'", adj=0.76)
mtext(text="'Fall'", adj=0.87)
dev.off()


saveRDS(gdatw, file="../admixture/gdatw3.RDS")

saveRDS(tbl2s, file="../admixture/3-tbl2s.RDS")

saveRDS(tbl2, file="../admixture/3-tbl2.RDS")

saveRDS(tbl3, file="../admixture/3-tbl3.RDS")


NOR <- 
  mutate(gdatw, population=gdatw_pop) |>
  filter(grepl(pattern="White_R", population)) |>
  filter(!grepl(pattern="Baseline", population))

table(NOR$population)

write_delim(mutate(NOR[-509], ID=paste0("i", ID)), col_names=FALSE, file="../admixture/NOR3.ped", na="0")


data.frame(Chromosome=1, SNP=paste0("rs", 10001:10251), GD=0, POS=seq(1, 251000000, 1000000)) |> 
  write_delim(col_names=FALSE, file="../admixture/NOR3.map")

# plink --file NOR3 --make-bed --out NOR3
# for K in 1 2 3 4; do admixture -s $RANDOM --cv=10 NOR3.bed $K | tee NOR3_log${K}.out; done
# grep -h CV NOR3_log*  > NOR3_log.txt


NOR_cv <- 
  read_delim("../admixture/NOR3_log.txt", col_names=FALSE) %>% 
  mutate(K=1:nrow(.)) |> 
  select(K, CV_error=X4)

plot(NOR_cv)
lines(NOR_cv)

NOR_1 <- read.table("../admixture/NOR3.1.Q")
NOR_2 <- read.table("../admixture/NOR3.2.Q")

NOR_1$id <- 
  NOR$ID
NOR_2$id <- 
  NOR$ID

NOR_1$pop <- 
  NOR$population
NOR_2$pop <- 
  NOR$population

NOR_2 <- 
  bind_rows(NOR_2, tibble(V1=rep(NA, 2), V2=rep(NA, 2)))

rownames(NOR_2)[836:837] <-
  paste0("spacer", 1:2)

rownames(NOR_2)


NOR_2 <-
  NOR_2[rownames(NOR_2)[c(1:741, 836, 742:777, 837, 778:835)], ]

pdf("../admixture/NOR3_k2_3.pdf", width=12, height=6)
par(mar=c(7, 4, 4, 2) + 0.1, mfrow=c(2, 1))
barplot(t(as.matrix(data.frame(V1=c(rep(1, 741), NA, rep(1, length(742:777)), NA, rep(1, length(778:835)))))), col=terrain.colors(4)[1], ylab="Ancestry", border=NA, names.arg=NOR_2$pop, las=2, cex.names=0.3, space=0, main=expression(italic(K)==1*`*`), font.main=1)
mtext(text="'Spring'", adj=0.37)
mtext(text="?", adj=0.80)
mtext(text="'Fall'", adj=0.97)
barplot(t(as.matrix(NOR_2[1:2])), col=terrain.colors(4)[c(1, 3)], ylab="Ancestry", border=NA, xaxt="n", las=2, cex.names=0.3, space=0, main=expression(italic(K)==2), font.main=1)
mtext(text="'Spring'", adj=0.37)
mtext(text="?", adj=0.89)
mtext(text="'Fall'", adj=0.97)
dev.off()
