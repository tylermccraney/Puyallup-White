library(readxl)
library(tidyverse)

name_repair <-
  function(names){
    names <- 
      tolower(names) |>
      gsub(pattern="[^A-Za-z0-9_%]+", replacement="_") |> 
      gsub(pattern="%", replacement="percent", fixed=TRUE)
    is_dup <- 
      names %in% names[duplicated(names)]
    if(any(is_dup)){
      suffixes <- 
        ave(names, names, FUN=function(x) paste0("_", seq_along(x)))
      names[is_dup] <- 
        paste0(names[is_dup], suffixes[is_dup])}
    return(names)}

r25 <- 
  read_excel(
    path="../data-from-Amelia/2025_WhiteRiverHatcherySpringChinook_GenotypingReport.xlsx", 
    .name_repair=name_repair, 
    na="NA")

r26 <- 
  read_excel(
    path="../data-from-Amelia/2026_WhiteRiverHatcherySpringChinook_GenotypingReport.xlsx", 
    .name_repair=name_repair, 
    na="NA")

r09_24 <- 
  read_excel(
    path="../White River NOR Chinook DNA Data 2009 to 2024.xlsx", 
    .name_repair=name_repair, 
    na="NA") |> 
  mutate(date=as.numeric(date) |> as.Date(origin="1899-12-30"))

rdat <- 
  bind_rows(r09_24, r25, r26)

table(rdat$sex)

rdat$lh <- 
  rdat$sex

rdat$sex <- 
  replace_values(
    rdat$sex,
    "f" ~ "F", 
    c("m", "J") ~ "M")

table(rdat$lh)

rdat$lh <- 
  replace_values(
    rdat$lh,
    c("f", "F", "m", "M") ~ "adult", 
    "J" ~ "jack")

table(rdat$assignment)

rdat$assignment <- 
  replace_values(
    rdat$assignment,
    "Unassigned" ~ "Unknown")



baseline <- 
  read_tsv(
    "../data-from-Amelia/Updated_WhiteRiverChinookBaseline_RubiasInput.txt", 
    col_types=str_flatten(rep("c", 594)), 
    na="NA",
    name_repair=name_repair)

glimpse(baseline)

baseline <- 
  reduce(
    names(baseline) |> 
      grep("^ots\\w+_[12]$", x=_, value=TRUE) |> 
      sub("_[12]$", "", x=_) |> 
      unique(),
    \(df, prefix){
      unite(
        df,
        col=!!prefix,
        matches(paste0("^", prefix, "_[12]$")),
        sep="",
        na.rm=TRUE)},
    .init=baseline) |>
  mutate(across(5:299, ~ na_if(., "")))



group_by(baseline, repunit, collection) |> count()



nor <- 
  read_csv(
    "../data-from-Amelia/WhiteRiver_NORGenotypeData_2021_2026.csv", 
    na="0", 
    name_repair=name_repair) |> 
  rename(id=sample_name)

names(nor)

rdat

dat <- 
  bind_rows(rename(baseline, id=indiv), nor) |> 
  select(-sample_type, -year)

names(dat)


gdat <- 
  select(dat, -299:-331, -335) |> 
  rowwise() |> 
  mutate(prop=sum(across(4:301, ~!is.na(.)))/298) |> 
  ungroup() |> 
  relocate(prop, .after=id) |> 
  filter(prop>=0.75)

gdat$prop |> summary()


gdat <- 
  left_join(x=gdat, y=rename(rdat, id=wdfw_id)) |> 
  relocate(303:317, .after=id) |> 
  filter(is.na(comment)) |>
  select(-genotypingsuccess, -comment) |> 
  mutate(year=ifelse(is.na(year), yes=(str_extract(string=id, pattern="^\\d{2}") |> as.numeric())+2000, no=year))

gdat$year <- 
  replace_values(gdat$year, 2098 ~ 1998)

gdat <- 
  arrange(gdat, year, id)

gdat <- 
  mutate(
    gdat, 
    collection=ifelse(
      test=is.na(collection), 
      yes=paste("White_River_Chinook_at_Buckley_Trap", assignment, str_extract(string=id, pattern="^\\w{4}"), sep="_"), 
      no=collection) |> 
      gsub(pattern="@", replacement="at", fixed=TRUE), 
    length=if_else(condition=length>100, true=length/10, false=length))


distinct(gdat, collection)

collnames <- 
  table(gdat$collection, useNA="ifany") |> 
  names()

gdat <- 
  filter(gdat, !collection %in% collnames[3:4]) |> 
  mutate(collection=factor(collection))

gdat$collection <- 
  factor(
    gdat$collection, 
    levels=collnames[c(23, 7:8, 17:22, 24:29, 9:14, 6, 2, 5, 1)], 
    ordered=TRUE)

gdat <- 
  arrange(gdat, collection, id)

gdat <- 
  distinct(gdat, id, .keep_all=TRUE) |> 
  filter(!is.na(collection))


names(gdat)


gsnps <- 
  enframe(colSums(!is.na(gdat[18:315]))/nrow(gdat)) |> 
  filter(value < 1) |> 
  filter(value >= 0.75)

gdatw <- 
  select(gdat, collection, id, gsnps$name) |> 
  separate_wider_position(cols=3:246, widths=c(a1=1, a2=1), names_sep="_")


gdatw <- 
  rename(gdatw, pop=collection, ID=id) |> 
  mutate(PaternalID=0, MaternalID=0, Sex=0, Phenotype=1) |> 
  relocate(491:494, .after=ID)

gdatw$pop <- 
  gsub(pattern="-", replacement="_", gdatw$pop)

write_delim(mutate(gdatw, ID=paste0("i", ID)), col_names=FALSE, file="../admixture/gdatw.ped", na="0")

length(names(gdatw)[7:494])/2

data.frame(Chromosome=1, SNP=paste0("rs", 10001:10244), GD=0, POS=seq(1, 244000000, 1000000)) |> 
  write_delim(col_names=FALSE, file="../admixture/gdatw.map")

str(gdat)

# bash
#
# conda activate admixture_Env
# cd /Users/tylermccraney/Google\ Drive/My\ Drive/Puyallup-White/admixture
# plink --file gdatw --make-bed --out gdatw
# for K in 1 2 3 4; do admixture -s $RANDOM --cv=10 gdatw.bed $K | tee gdatw_log${K}.out; done
# grep -h CV gdatw_log*  > gdatw_log.txt



gdatw_cv <- 
  read_delim("../admixture/gdatw_log.txt", col_names=FALSE) %>% 
  mutate(K=1:nrow(.)) |> 
  select(K, CV_error=X4)

plot(gdatw_cv)
lines(gdatw_cv)

tbl2 <- read.table("../admixture/gdatw.2.Q")

tbl3 <- read.table("../admixture/gdatw.3.Q")

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
  gsub(pattern="__", replacement="_")


unique(gdatw_pop)


rownames(tbl2) <- 
  gdatw$ID

rownames(tbl3) <- 
  gdatw$ID

tbl2$pop <- 
  gdatw_pop

tbl3$pop <- 
  gdatw_pop

tbl2 <- 
  bind_rows(tbl2, tibble(V1=rep(NA, 4), V2=rep(NA, 4), pop=rep(NA, 4)))

tbl3 <- 
  bind_rows(tbl3, tibble(V1=rep(NA, 4), V2=rep(NA, 4), V3=rep(NA, 4), pop=rep(NA, 4)))

rownames(tbl2)[1126:1129] <- 
  paste0("spacer", 1:4)

rownames(tbl3)[1126:1129] <- 
  paste0("spacer", 1:4)

rownames(tbl3)[1:1000]; rownames(tbl3)[1001:1129]


tbl2 <- 
  tbl2[rownames(tbl2)[c(1:129, 1126, 130:870, 1127, 871:906, 1128, 907:964, 1129, 965:1125)], ]

tbl3 <- 
  tbl3[rownames(tbl3)[c(1:129, 1126, 130:870, 1127, 871:906, 1128, 907:964, 1129, 965:1125)], ]


pdf("../admixture/gdatw_k2_3.pdf", width=12, height=6)
par(mar=c(7, 4, 4, 2) + 0.1, mfrow=c(2, 1))
barplot(t(as.matrix(tbl2[-3])), col=terrain.colors(4)[c(3, 1)], ylab="Ancestry", border=NA, names.arg=tbl2$pop, las=2, cex.names=0.3, space=0, main=expression(italic(K)==2*`*`), font.main=1)
barplot(t(as.matrix(tbl3[-4])), col=terrain.colors(4)[c(3, 1, 2)], ylab="Ancestry", border=NA, names.arg=tbl3$pop, las=2, cex.names=0.3, space=0, main=expression(italic(K)==3), font.main=1)
dev.off()

gdat

unique(gdatw_pop)

NOR <- 
  mutate(gdatw, population=gdatw_pop) |>
  filter(grepl(pattern="White_R", population)) |>
  filter(!grepl(pattern="Baseline", population))

table(NOR$population)

write_delim(mutate(NOR[-495], ID=paste0("i", ID)), col_names=FALSE, file="../admixture/NOR.ped", na="0")

length(names(gdatw)[7:494])/2

data.frame(Chromosome=1, SNP=paste0("rs", 10001:10244), GD=0, POS=seq(1, 244000000, 1000000)) |> 
  write_delim(col_names=FALSE, file="../admixture/NOR.map")

# plink --file NOR --make-bed --out NOR
# for K in 1 2 3 4; do admixture -s $RANDOM --cv=10 NOR.bed $K | tee NOR_log${K}.out; done
# grep -h CV NOR_log*  > NOR_log.txt


NOR_cv <- 
  read_delim("../admixture/NOR_log.txt", col_names=FALSE) %>% 
  mutate(K=1:nrow(.)) |> 
  select(K, CV_error=X4)

plot(NOR_cv)
lines(NOR_cv)

NOR_1 <- read.table("../admixture/NOR.1.Q")
NOR_2 <- read.table("../admixture/NOR.2.Q")

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

pdf("../admixture/NOR_k2_3.pdf", width=12, height=6)
par(mar=c(7, 4, 4, 2) + 0.1, mfrow=c(2, 1))
barplot(t(as.matrix(data.frame(V1=c(rep(1, 741), NA, rep(1, length(742:777)), NA, rep(1, length(778:835)))))), col=terrain.colors(4)[1], ylab="Ancestry", border=NA, names.arg=NOR_2$pop, las=2, cex.names=0.3, space=0, main=expression(italic(K)==1*`*`), font.main=1)
barplot(t(as.matrix(NOR_2[2:1])), col=terrain.colors(4)[c(3, 1)], ylab="Ancestry", border=NA, names.arg=NOR_2$pop, las=2, cex.names=0.3, space=0, main=expression(italic(K)==2), font.main=1)
dev.off()


saveRDS(dat, file="../data/dat.RDS")

saveRDS(rdat, file="../data/rdat.RDS")

saveRDS(select(gdat, 1:17, gsnps$name), file="../data/gdat.RDS")

saveRDS(gdatw, file="../admixture/gdatw.RDS")


