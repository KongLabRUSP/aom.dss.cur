# |----------------------------------------------|
# | Project: Curcumin Treatment of AOM-DSS Model |
# | Study ID:                                    |
# | Scientist: Yue, Renyi Wu                     |
# | Data Analysis: Renyi Wu, Davit Sargsyan      |
# | Created: 11/09/2017                          |
# |----------------------------------------------|
# Header----
# NOTE: run RStudio AS ADMINISDTRATOR!
# source("https://bioconductor.org/biocLite.R")
# biocLite("TxDb.Mmusculus.UCSC.mm9.knownGene",
#          suppressUpdates = TRUE)
# biocLite("ChIPseeker",
#          suppressUpdates = TRUE)
# biocLite("org.Mm.eg.db",
#          suppressUpdates = TRUE)

require(data.table)
require(ggplot2)
require(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm9.knownGene)

# Part I: Data----
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

# # "Bad" batch (Batch 2)
# dt2 <- fread(dList[2])
# # dt2[, 5:16] <- lapply(dt2[, 5:16],
# #                       function(a) {
# #                         a[is.na(a)] <- 0
# #                         a <- a + 1
# #                       })
# dt2

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

# t2 <- data.table(dt2[, 1:4], 
#                  dt2[, seq(6, 16, 2), with = FALSE]/ 
#                    dt2[, seq(5, 15, 2), with = FALSE])
# names(t2)[5:10] <- sName
# t2

t3 <- data.table(dt3[, 1:4], 
                 dt3[, seq(6, 16, 2), with = FALSE]/ 
                   dt3[, seq(5, 15, 2), with = FALSE])
names(t3)[5:10] <- sName
t3

# Combine first and third batches
dt13 <- merge(t1, 
              t3, 
              by = names(t1)[1:4],
              all = TRUE)
dt13

# Clean environment
rm(dt1, dt3, dList, sName, t1, t3)
gc()

# Part II: Annotate----
write.table(dt13, 
          file = "tmp/dt13.csv",
          row.names = FALSE,
          sep = "\t")

peakAnno1 <- annotatePeak(peak = "tmp/dt13.csv", 
                          tssRegion = c(-3000, 3000), 
                          TxDb = TxDb.Mmusculus.UCSC.mm9.knownGene,
                          annoDb = "org.Mm.eg.db")

# # CHECK: did I miss any region? All 7 region I found for TNF are in Promoter.
# peakAnno1 <- annotatePeak(peak = "data/comb2.csv", 
#                           tssRegion = c(-3000, 3000), 
#                           TxDb = TxDb.Mmusculus.UCSC.mm9.knownGene,
#                           annoDb = "org.Mm.eg.db")
# 
# peakAnno1 <- annotatePeak(peak = "data/Methyl_rep_2017.samdup.csv", 
#                           tssRegion = c(-3000, 3000), 
#                           TxDb = TxDb.Mmusculus.UCSC.mm9.knownGene,
#                           annoDb = "org.Mm.eg.db")
# # No, same results

dt1 <- data.table(as.data.frame(peakAnno1@anno@elementMetadata@listData))
dt1

# Remove unmapped regions
dt1 <- dt1[!is.na(dt1$SYMBOL == "NA"), ]
write.csv(dt1, 
            file = "tmp/dt13_annotated.csv",
            row.names = FALSE)

dt1[SYMBOL == "Tnf", ]

# Plot TNF
dtPlot <- dt1[SYMBOL == "Tnf", ]
dtPlot <- data.table(Location = paste(dtPlot$annotation,
                                   "\n#CpG =",
                                   dtPlot$CpG,
                                   ", ToTSS = ",
                                   dtPlot$distanceToTSS),
                     Control_1 = (dtPlot$Control_1.x + dtPlot$Control_2.x)/2,
                     Control_2 = (dtPlot$Control_1.y + dtPlot$Control_2.y)/2,
                     AOM_DSS_1 = (dtPlot$AOM_DSS_1.x + dtPlot$AOM_DSS_2.x)/2,
                     AOM_DSS_2 = (dtPlot$AOM_DSS_1.y + dtPlot$AOM_DSS_2.y)/2,
                     Curcumin_1 = (dtPlot$Cur_1.x + dtPlot$Cur_2.x)/2,
                     Curcumin_2 = (dtPlot$Cur_1.y + dtPlot$Cur_2.y)/2)
dtPlot <- melt.data.table(dtPlot,
                          id.vars = "Location",
                          variable.name = "Group",
                          value.name = "Methylation")

dtPlot$Location <- factor(dtPlot$Location,
                          levels = unique(dtPlot$Location))

dtPlot$Group <- as.character(dtPlot$Group)
dtPlot$Treatment <- gsub(x = dtPlot$Group, 
                         pattern = "_1", 
                         replacement = "")
dtPlot$Treatment <- gsub(x = dtPlot$Treatment, 
                         pattern = "_2", 
                         replacement = "")
dtPlot$Treatment <- factor(dtPlot$Treatment,
                          levels = unique(dtPlot$Treatment))

dtPlot$Batch <- factor(as.numeric(substr(dtPlot$Group, 
                                         nchar(dtPlot$Group),
                                         nchar(dtPlot$Group))))

p1 <- ggplot(dtPlot,
             aes(x = Treatment,
                 y = Methylation,
                 fill = Batch,
                 group = Batch)) +
  facet_wrap(~ Location,
             nrow = 2,
             scales = "free") +
  geom_bar(position = position_dodge(),
           stat="identity") +
  ggtitle("TNFa Methylation Level by DMR")
p1

tiff(filename = "tmp/tfn_methyl_barplot.tiff",
     height = 5,
     width = 12,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()