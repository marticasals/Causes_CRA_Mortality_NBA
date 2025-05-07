## R code of the analysis in
## Martínez et al. (2025): "Causes of Death in NBA Players"
## ========================================================
rm(list = ls(all = TRUE))

## R packages needed
## -----------------
library(survival)
library(dplyr)
library(descr)
library(etm)
library(ggplot2)

## Data set
## --------
load("NBA2019.RData")

## Complete data set
str(nba2019)
head(nba2019)
summary(nba2019)

## Subset of African American and White players
nba2019AAW <- droplevels(subset(nba2019, etni %in% c("Afr.American", "White")))
summary(nba2019AAW)


## Table 1
## =======
freq(subset(nba2019, status == 1)$dcause5, plot = FALSE)


## Preparation of multiple imputations
## ===================================
## Distribution of death causes (excluding idiopathic)
## ---------------------------------------------------
## All players
(tabAll <- with(subset(nba2019, status == 1 & dcauseCod != 3), table(dcauseCod)))
(ptabAll <- prop.table(tabAll))
rm(tabAll)
## Subset of African American and White players
(tabAAW <- with(subset(nba2019AAW, status == 1 & dcauseCod != 3),
                table(dcauseCod, etni)))
(ptabAAW <- prop.table(tabAAW, 2))
rm(tabAAW)

## Numbers of missing death causes
## -------------------------------
cod3All <- sum(nba2019$dcauseCod == 3)
cod3AA <- sum(nba2019AAW$etni == "Afr.American" & nba2019AAW$dcauseCod == 3)
cod3Wh <- sum(nba2019AAW$etni == "White" & nba2019AAW$dcauseCod == 3)

## Number of imputations
## ---------------------
M <- 15

## (Number of different) death times
## ---------------------------------
## All players
numtimes <- n_distinct(subset(nba2019, status == 1)$ageright)
Alltimes <- sort(unique(subset(nba2019, status == 1)$ageright))

## African American and White players
numtimesAAW <- nba2019AAW |>
  filter(status == 1) |>
  group_by(etni) |>
  summarise(n_distinct(ageright)) |>
  as.data.frame()
AAtimes <- sort(unique(subset(nba2019AAW, status == 1 &
                                          etni == "Afr.American")$ageright))
Whtimes <- sort(unique(subset(nba2019AAW, status == 1 &
                                          etni == "White")$ageright))

## Lists for cumulative incidence functions CIFs
## ---------------------------------------------
cifs <- vector("list", 4)
names(cifs) <- levels(nba2019AAW$dcause5)[c(1, 2, 4, 5)]
for (i in 1:4) {
  cifs[[i]] <- vector("list", 3)
  names(cifs[[i]]) <- c("All", levels(nba2019AAW$etni))
  cifs[[i]][[1]] <- cbind(times = Alltimes,
                          matrix(nrow = numtimes, ncol = M * 2,
                                 dimnames = list(1:numtimes,
                                                 paste0(rep(c("F", "V"), each = M),
                                                        1:M))))
  cifs[[i]][[2]] <- cbind(times = AAtimes,
                          matrix(nrow = numtimesAAW[1, 2], ncol = M * 2,
                                 dimnames = list(1:numtimesAAW[1, 2],
                                                 paste0(rep(c("F", "V"), each = M),
                                                        1:M))))
  cifs[[i]][[3]] <- cbind(times = Whtimes,
                          matrix(nrow = numtimesAAW[2, 2], ncol = M * 2,
                                 dimnames = list(1:numtimesAAW[2, 2],
                                                 paste0(rep(c("F", "V"), each = M),
                                                        1:M))))
}
rm(i)

## Lists for Cox models
## --------------------
betas <- vector("list", 4)
names(betas) <- levels(nba2019AAW$dcause5)[c(1, 2, 4, 5)]

## Lists for parameters (within betas)
for (i in 1:4){
  betas[[i]] <- vector("list", 2)
  names(betas[[i]]) <- c("Estimate", "Variance")
  betas[[i]][[1]] <- betas[[i]][[2]] <-
    matrix(ncol = M, nrow = 3, dimnames = list(c("Etni", "Height", "Year"),
                                               paste0("Imp.", 1:M)))
}
rm(i)


## Multiple imputations (with M iterations)
## ========================================
set.seed(1893)
for (m in 1:M) {
  nbaimpdatAll <- nba2019
  nbaimpdatAAW <- nba2019AAW |>
    mutate(etni = relevel(etni, ref = "White"))
  ## Imputation of death causes
  nbaimpdatAll[nbaimpdatAll$dcauseCod == 3,
               "dcauseCod"] <- sample(c(1, 2, 4, 5), cod3All, replace = TRUE,
                                      prob = ptabAll)
  nbaimpdatAAW[nbaimpdatAAW$dcauseCod == 3 & nbaimpdatAAW$etni == "White",
               "dcauseCod"] <- sample(c(1, 2, 4, 5), cod3Wh, replace = TRUE,
                                      prob = ptabAAW[, 2])
  nbaimpdatAAW[nbaimpdatAAW$dcauseCod == 3 & nbaimpdatAAW$etni == "Afr.American",
               "dcauseCod"] <- sample(c(1, 2, 4, 5), cod3AA, replace = TRUE,
                                      prob = ptabAAW[, 1])

  ## Cumulative incidence functions
  AllCIFsAll <- summary(etmCIF(Surv(ageleft, ageright, status) ~ 1, nbaimpdatAll,
                               etype = dcauseCod))
  AllCIFsAAW <- summary(etmCIF(Surv(ageleft, ageright, status) ~ etni, nbaimpdatAAW,
                               etype = dcauseCod))
  for (d in 1:4) {
    cifs[[d]]$All[, paste0(c("F", "V"), m)] <-
      as.matrix(AllCIFsAll[[1]][[d]][, c("P", "var")])
    cifs[[d]]$Afr.American[, paste0(c("F", "V"), m)] <-
      as.matrix(AllCIFsAAW$"etni=Afr.American"[[d]][, c("P", "var")])
    cifs[[d]]$White[, paste0(c("F", "V"), m)] <-
      as.matrix(AllCIFsAAW$"etni=White"[[d]][, c("P", "var")])
  }
  rm(d, AllCIFsAll, AllCIFsAAW, nbaimpdatAll)

  ## Cause-specific hazard models
  for (d in 1:4) {
    dc <- ifelse(d <= 2, d, d + 1)
    modCause <- coxph(Surv(ageleft, ageright,
                           dcauseCod == dc) ~ etni + height + debutyear,
                      nbaimpdatAAW)
    coefftable <- summary(modCause)$coef
    ## Parameters estimates
    betas[[d]]$Estimate[, m] <- coefftable[, 1]
    betas[[d]]$Variance[, m] <- coefftable[, 3]^2
    rm(modCause, coefftable, dc)
  }
  rm(d, nbaimpdatAAW)
}
rm(m, ptabAll, ptabAAW, cod3All, cod3AA, cod3Wh)


## Estimation of cumulative incidences
## ===================================
## Data frames with estimated cumulative incidence functions
## ---------------------------------------------------------
resultsAll <- data.frame(times = Alltimes, cif1 = NA, cif2 = NA, cif4 = NA, cif5 = NA)
resultsAA  <- data.frame(times = AAtimes, cif1 = NA, cif2 = NA, cif4 = NA, cif5 = NA)
resultsWh  <- data.frame(times = Whtimes, cif1 = NA, cif2 = NA, cif4 = NA, cif5 = NA)

for (d in 1:4) {
  ## CIFs estimates
  resultsAll[, d + 1] <- round(rowMeans(cifs[[d]]$All[, paste0("F", 1:M)]), 4)
  resultsAA[, d + 1]  <- round(rowMeans(cifs[[d]]$Afr.American[, paste0("F", 1:M)]), 4)
  resultsWh[, d + 1]  <- round(rowMeans(cifs[[d]]$White[, paste0("F", 1:M)]), 4)
}

## Table 2
## -------
## Look for values in times equal to or just below
## 50, 55, ..., 95 years.

## Figure 1
## --------
tiff("Figure1.tiff", width = 1200, height = 800, res = 90)
par(mfrow = c(2, 2), las = 1, font.lab = 4, font.axis = 2, font = 2,
    cex.lab = 1.3, cex.axis = 1.2)
for (d in 1:4) {
  dn <- ifelse(d <= 2, d, d + 1)
  with(resultsWh, plot(times, resultsWh[, d + 1], xlab = "Age",
                       ylab = "Cumulative incidence", type = "S", lwd = 3,
                       lty = 2, ylim = c(0, 0.35)))
  with(resultsAA, lines(times, resultsAA[, d + 1], type = "S", lwd = 3, col = 2,
                        lty = 1))
  axis(1, at = seq(20, 95, 5)[-c(5, 9, 13)])
  legend("topleft", c("African American", "White"), col = 2:1, lty = 1:2,
         bty = "n", lwd = 3, cex = 1.2)
  title(paste("Death cause:", levels(nba2019AAW$dcause5)[dn]), cex.main = 1.3)
}
dev.off()

## =============================
## Cause-specific hazards models
## =============================
## Data frame with summary statistics for each death cause
## -------------------------------------------------------
results <- data.frame(cause = factor(rep(levels(nba2019AAW$dcause5)[-3], each = 3)),
                      variables = c("African American", "Height [5 cm]",
                                    "NBA debut [5 years]"),
                      HR = NA, CIlow = NA, CIupp = NA)
results$cause <- relevel(results$cause, ref = "Neoplasm")

for (d in 1:4) {
  ## Parameter estimates
  betahat <- rowMeans(betas[[d]]$Estimate)

  ## Parameter variances (MI method)
  betahatVar <- rowMeans(betas[[d]]$Variance) + (1 + 1 / M) * apply(betas[[d]]$Estimate, 1, var)

  ## Hazard ratios
  hr1 <- round(exp(betahat[1]), 3)
  hr2 <- round(exp(5 * betahat[2]), 3)
  hr3 <- round(exp(5 * betahat[3]), 3)

  ## 95% confidence intervals of the hazard ratios
  ci1 <- round(exp(betahat[1] + qnorm(c(0.025, 0.975)) * sqrt(betahatVar[1])), 3)
  ci2 <- round(exp(5 * (betahat[2] + qnorm(c(0.025, 0.975)) * sqrt(betahatVar[2]))), 3)
  ci3 <- round(exp(5 * (betahat[3] + qnorm(c(0.025, 0.975)) * sqrt(betahatVar[3]))), 3)

  results[as.numeric(results$cause) == d, "HR"] <- c(hr1, hr2, hr3)
  results[as.numeric(results$cause) == d, "CIlow"] <- c(ci1[1], ci2[1], ci3[1])
  results[as.numeric(results$cause) == d, "CIupp"] <- c(ci1[2], ci2[2], ci3[2])
  rm(betahat, betahatVar, hr1, hr2, hr3, ci1, ci2, ci3)
}
rm(d, M)

results$variables <- factor(results$variables)
results$variables <- factor(results$variables,
                            levels = levels(results$variables)[3:1])

## Forest plots
## ------------
g <- ggplot(data = results, aes(x = variables, y = HR, ymin = CIlow,
                                ymax = CIupp, color = cause)) +
  geom_pointrange(aes(col = variables), lwd = 0.8) +
  geom_hline(aes(fill = variables), yintercept = 1, linetype = 2) +
  xlab("") +
  ylab("Hazard Ratio (95% Confidence Interval)") +
  geom_errorbar(aes(ymin = CIlow, ymax = CIupp, col = variables),
                width = 0.5, cex = 1) +
  facet_wrap(~cause, strip.position = "left", nrow = 9, scales = "free_y") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(hjust = 0,vjust = 1, angle = 180,
                                    face = "bold")) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_log10(n.breaks = 10, limits = c(min(results$CIlow), max(results$CIupp)))
g
ggsave("Figure2.tiff", height = 9)
rm(g)
