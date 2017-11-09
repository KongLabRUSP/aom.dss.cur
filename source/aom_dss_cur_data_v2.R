# |----------------------------------------------|
# | Project: Curcumin Treatment of AOM-DSS Model |
# | Study ID:                                    |
# | Scientist: Yue, Renyi Wu                     |
# | Data Analysis: Renyi Wu, Davit Sargsyan      |
# | Created: 11/06/2017                          |
# |----------------------------------------------|
# Header----
require(data.table)
require(ggplot2)

# Move up one directory
wd <- getwd()

# Part I: Load data----
setwd("data")
dList <- dir()
dList

# First batch
dt1 <- fread(dList[1])
dt1 <- dt1[, c(1:4,
               25:36)]
# dt1[, 5:16] <- lapply(dt1[, 5:16],
#                      function(a) {
#                        a[is.na(a)] <- 0
#                        a <- a + 1
#                      })
dt1

# "Bad" batch (Batch 2)
dt2 <- fread(dList[2])
# dt2[, 5:16] <- lapply(dt2[, 5:16],
#                       function(a) {
#                         a[is.na(a)] <- 0
#                         a <- a + 1
#                       })
dt2

# Rerun of the samples from Batch 2
dt3 <- fread(dList[3])
# Keep Cur samples only
dt3 <- dt3[, 1:16]
# dt3[, 5:16] <- lapply(dt3[, 5:16],
#                       function(a) {
#                         a[is.na(a)] <- 0
#                         a <- a + 1
#                       })
dt3

# Reset working directory
setwd(wd)

# Rename all samples using treatment names
sName <- c("Control_1",
           "Control_2",
           "AOM_DSS_1",
           "AOM_DSS_2",
           "Cur_1",
           "Cur_2")

t1 <- data.table(dt1[, 1:4], 
                 dt1[, seq(6, 16, 2), with = FALSE]/ 
                   dt1[, seq(5, 15, 2), with = FALSE])
names(t1)[5:10] <- sName
t1

t2 <- data.table(dt2[, 1:4], 
                 dt2[, seq(6, 16, 2), with = FALSE]/ 
                   dt2[, seq(5, 15, 2), with = FALSE])
names(t2)[5:10] <- sName
t2

t3 <- data.table(dt3[, 1:4], 
                 dt3[, seq(6, 16, 2), with = FALSE]/ 
                   dt3[, seq(5, 15, 2), with = FALSE])
names(t3)[5:10] <- sName
t3

rm(dt1, dt2, dt3, dList, sName)
gc()

# Part II: Compare Treatments and Batches----
# Combine first and third batches
dtc <- merge(t1, t3, by = names(t1)[1:4])

# Group averages: Batch 1----
dtc$ctrl_x <- (dtc$Control_1.x + dtc$Control_2.x)/2
dtc$aomdss_x <- (dtc$AOM_DSS_1.x + dtc$AOM_DSS_2.x)/2
dtc$cur_x <- (dtc$Cur_1.x + dtc$Cur_2.x)/2

# Fold-changes----
dtc$log2.aomdss.ctrl.x <- log2(dtc$aomdss_x/dtc$ctrl_x)
dtc$log2.aomdss.cur.x <- log2(dtc$aomdss_x/dtc$cur_x)

# Absolute differences----
dtc$diff.aomdss.ctrl.x <- dtc$aomdss_x - dtc$ctrl_x
dtc$diff.aomdss.cur.x <- dtc$aomdss_x - dtc$cur_x

# Plot proportions methylated----
hist(dtc$ctrl_x, 
     100, 
     xlab = "Proportion Methylated CpG",
     main = "Batch 1 Control")
hist(dtc$aomdss_x, 
     100, 
     xlab = "Proportion Methylated CpG",
     main = "Batch 1 AOM-DSS")
hist(dtc$cur_x, 
     100, 
     xlab = "Proportion Methylated CpG",
     main = "Batch 1 Curcumin")

# Differences in treatments within Batch 1-----
x <- dtc$ctrl_x
y <- dtc$aomdss_x
z <- dtc$cur_x

plot(y ~ x,
     xlab = "Control",
     ylab = "AOM-DSS",
     main = "Batch 1. RED:Diff <= -10% or >= 10%")
points(y[abs(dtc$diff.aomdss.ctrl.x) >= 0.1] ~ 
         x[abs(dtc$diff.aomdss.ctrl.x) >= 0.1],
       col = "red")

plot(z ~ y,
     xlab = "AOM-DSS",
     ylab = "Curcumin",
     main = "Batch 1. RED:Diff <= -10% or >= 10%")
points(z[abs(dtc$diff.aomdss.cur.x) >= 0.1] ~ 
         y[abs(dtc$diff.aomdss.cur.x) >= 0.1],
       col = "red")

# Fold-changes in treatments within Batch 1-----
x <- -log2(dtc$ctrl_x)
y <- -log2(dtc$aomdss_x)
z <- -log2(dtc$cur_x)

plot(y ~ x,
     xlab = "-log2(Control)",
     ylab = "-log2(AOM-DSS)",
     main = "Batch 1. RED:log2 <= -1 or >= 1")
points(y[abs(dtc$log2.aomdss.ctrl.x) >= 1] ~ 
         x[abs(dtc$log2.aomdss.ctrl.x) >= 1],
       col = "red")

plot(z ~ y,
     xlab = "-log2(AOM-DSS)",
     ylab = "-log2(Curcumin)",
     main = "Batch 1. RED:log2 <= -1 or >= 1")
points(z[abs(dtc$log2.aomdss.cur.x) >= 1] ~ 
         y[abs(dtc$log2.aomdss.cur.x) >= 1],
       col = "red")

# Group averages: Batch 2 Rerun----
dtc$ctrl_y <- (dtc$Control_1.y + dtc$Control_2.y)/2
dtc$aomdss_y <- (dtc$AOM_DSS_1.y + dtc$AOM_DSS_2.y)/2
dtc$cur_y <- (dtc$Cur_1.y + dtc$Cur_2.y)/2

# Differences in treatments within Batch 2 Rerun-----
# Fold-cahnges----
dtc$log2.aomdss.ctrl.y <- log2(dtc$aomdss_y/dtc$ctrl_y)
dtc$log2.aomdss.cur.y <- log2(dtc$aomdss_y/dtc$cur_y)

# Absolute differences----
dtc$diff.aomdss.ctrl.y <- dtc$aomdss_y - dtc$ctrl_y
dtc$diff.aomdss.cur.y <- dtc$aomdss_y - dtc$cur_y

# Plot proportions methylated
hist(dtc$ctrl_y, 
     100, 
     xlab = "Proportion Methylated CpG",
     main = "Batch 2 Rerun Control")
hist(dtc$aomdss_y, 
     100, 
     xlab = "Proportion Methylated CpG",
     main = "Batch 2 Rerun AOM-DSS")
hist(dtc$cur_y, 
     100, 
     xlab = "Proportion Methylated CpG",
     main = "Batch 2 Rerun Curcumin")

# Differences in treatments within Batch 2 Rerun-----
x <- dtc$ctrl_y
y <- dtc$aomdss_y
z <- dtc$cur_y

plot(y ~ x,
     xlab = "Control",
     ylab = "AOM-DSS",
     main = "Batch 2 Rerun. RED:Diff <= -10% or >= 10%")
points(y[abs(dtc$diff.aomdss.ctrl.y) >= 0.1] ~ 
         x[abs(dtc$diff.aomdss.ctrl.y) >= 0.1],
       col = "red")

plot(z ~ y,
     xlab = "AOM-DSS",
     ylab = "Curcumin",
     main = "Batch 2 Rerun. RED:Diff <= -10% or >= 10%")
points(z[abs(dtc$diff.aomdss.cur.y) >= 0.1] ~ 
         y[abs(dtc$diff.aomdss.cur.y) >= 0.1],
       col = "red")

# Fold-cahnges in treatments within Batch 2-----
x <- -log2(dtc$ctrl_y)
y <- -log2(dtc$aomdss_y)
z <- -log2(dtc$cur_y)

plot(y ~ x,
     xlab = "-log2(Control)",
     ylab = "-log2(AOM-DSS)",
     main = "Batch 2 Rerun. RED:log2 <= -1 or >= 1")
points(y[abs(dtc$log2.aomdss.ctrl.y) >= 1] ~ 
         x[abs(dtc$log2.aomdss.ctrl.y) >= 1],
       col = "red")

plot(z ~ y,
     xlab = "-log2(AOM-DSS)",
     ylab = "-log2(Curcumin)",
     main = "Batch 2 Rerun. RED:log2 <= -1 or >= 1")
points(z[abs(dtc$log2.aomdss.cur.y) >= 1] ~ 
         y[abs(dtc$log2.aomdss.cur.y) >= 1],
       col = "red")

# Fold-changes between batches----
plot(dtc$log2.aomdss.ctrl.x  ~ 
       dtc$log2.aomdss.ctrl.y,
     xlim = c(-6, 6),
     xlab = "Batch 2 Repeat",
     ylim = c(-6, 6),
     ylab = "Batch 1",
     main = "Log2 Differences in AOM-DSS vs. Control")
abline(0, 1, col = "red")

plot(dtc$log2.aomdss.cur.x ~ 
       dtc$log2.aomdss.cur.y,
     xlim = c(-6, 6),
     xlab = "Batch 2 Repeat",
     ylim = c(-6, 6),
     ylab = "Batch 1",
     main = "Log2 Differences in Curcumin vs. AOM-DSS")
abline(0, 1, col = "red")

# Absolute differences between batches----
plot(dtc$diff.aomdss.ctrl.x  ~ 
       dtc$diff.aomdss.ctrl.y,
     xlab = "Batch 2 Repeat",
     ylab = "Batch 1",
     main = "Absolute Differences in AOM-DSS vs. Control")
abline(0, 1, col = "red")

plot(dtc$diff.aomdss.cur.x ~ 
       dtc$diff.aomdss.cur.y,
     xlab = "Batch 2 Repeat",
     ylab = "Batch 1",
     main = "Absolute Differences in Curcumin vs. AOM-DSS")
abline(0, 1, col = "red")

# Plot Sample 1 vs. Sample 2 from Batch 1
x <- -log(dtc$Control_1.x)
y <- -log(dtc$Control_2.x)
  
plot(y ~ x,
     xlab = "-log2(Control) Batch 1",
     ylab = "-log2(Control) Batch 2 Rerun",
     main = "RED:log2(AOM-DSS vs. Control) <= -1 or >= 1")
points(y[abs(dtc$log2.aomdss.ctrl.x) >= 1] ~ 
         x[abs(dtc$log2.aomdss.ctrl.x) >= 1],
       col = "red")

x <- -log(dtc$AOM_DSS_1.x)
y <- -log(dtc$AOM_DSS_2.x)
plot(y ~ x,
     xlab = "-log2(AOM_DSS) Batch 1",
     ylab = "-log2(AOM_DSS) Batch 2 Rerun",
     main = "RED:log2(AOM-DSS vs. Control) <= -1 or >= 1")
points(y[abs(dtc$log2.aomdss.ctrl.x) >= 1] ~ 
         x[abs(dtc$log2.aomdss.ctrl.x) >= 1],
       col = "red")

x <- -log(dtc$Cur_1.x)
y <- -log(dtc$Cur_2.x)
plot(y ~ x,
     xlab = "-log2(Curcumin) Batch 1",
     ylab = "-log2(Curcumin) Batch 2 Rerun",
     main = "RED:log2(Curcumin vs. AOM-DSS) <= -1 or >= 1")
points(y[abs(dtc$log2.aomdss.cur.x) >= 1] ~ 
         x[abs(dtc$log2.aomdss.cur.x) >= 1],
       col = "red")

# MA Plot: X+Y vs. X-Y
M <- log2(dtc$Control_1.x) - log2(dtc$Control_2.x)
A <- (log2(dtc$Control_1.x) + log2(dtc$Control_2.x))/2
plot(M ~ A,
     main = "Batch 1 Controls")

M <- log2(dtc$AOM_DSS_1.x) - log2(dtc$AOM_DSS_2.x)
A <- (log2(dtc$AOM_DSS_1.x) + log2(dtc$AOM_DSS_2.x))/2
plot(M ~ A, 
     main = "Batch 1 AOM_DSS")

M <- log2(dtc$Cur_1.x) - log2(dtc$Cur_2.x)
A <- (log2(dtc$Cur_1.x) + log2(dtc$Cur_2.x))/2
plot(M ~ A, 
     main = "Batch 1 Curcumin")

M <- log2(dtc$aomdss_x) - log2(dtc$ctrl_x)
A <- (log2(dtc$aomdss_x) + log2(dtc$ctrl_x))/2
plot(M ~ A,
     main = "Batch 1 AOM-DSS vs. Control Averages")

M <- log2(dtc$cur_x) - log2(dtc$aomdss_x)
A <- (log2(dtc$cur_x) + log2(dtc$aomdss_x))/2
plot(M ~ A,
     main = "Batch 1 Curcumin vs. AOM-DSS Averages")

M <- log2(dtc$cur_x) - log2(dtc$ctrl_x)
A <- (log2(dtc$cur_x) + log2(dtc$ctrl_x))/2
plot(M ~ A,
     main = "Batch 1 Curcumin vs. Control Averages")

# Absolute vs. Fold-change differences----
plot(dtc$diff.aomdss.ctrl.x  ~ 
       dtc$log2.aomdss.ctrl.x,
     xlab = "log2(AOM_DSS/Control)",
     ylab = "AOM_DSS - Control",
     main = "Batch 1")

plot(dtc$diff.aomdss.ctrl.y  ~ 
       dtc$log2.aomdss.ctrl.y,
     xlab = "log2(AOM_DSS/Control)",
     ylab = "AOM_DSS - Control",
     main = "Batch 2 Rerun")


# CONTINUE HERE, DS 11/07/2017!
# DO PCA!
#   
# # Plot averages of same treatment in different batches----
# plot(dtc$ctrl_x ~ dtc$ctrl_y)
# plot(dtc$aomdss_x ~ dtc$aomdss_y)
# plot(dtc$cur_x ~ dtc$cur_y)