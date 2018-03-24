# |-----------------------------------------------------|
# | Study: Curcumin Treatment of AOM-DSS Model          |
# | Script: Editing figures for publication (grayscale) |
# | Scientist: Yue, Renyi Wu                            |
# | Data Analysis: Renyi Wu, Davit Sargsyan             |
# | Created: 03/24/2018                                 |
# |-----------------------------------------------------|
# Header----
require(data.table)
require(imager)
# Source: https://cran.r-project.org/web/packages/imager/vignettes/gettingstarted.html

dt1 <- fread("docs/renyi_figures_03222018/Fig 3D -methyl-seq-18wk-heatmap.csv")
dt1 <- dt1[!is.na(Control_18.mu), ]

# Make gene names unique----
dt1$N <- NULL
dt1[duplicated(dt1$gene), 
    N := 1:.N,
    by = gene]
dt1$N[is.na(dt1$N)] <- ""
dt1$reg <- paste(dt1$gene,
                 dt1$N,
                 sep = "")
dt1$reg

dt2 <- melt.data.table(dt1,
                       id.vars = "reg",
                       measure.vars = 1:3,
                       variable.name = "Treatment",
                       value.name = "Methylation Ratio")
summary(dt2)
hist(dt2$`Methylation (%)`, 100)

length(unique(dt2$reg))
dt2[substr(reg,
           1,
           3) == "Tnf",]

dt2$reg <- factor(dt2$reg,
                  levels = dt1$reg[order(dt1$Control_18.mu)])

# Figure 1: Hitmap
p1 <- ggplot(data = dt2) +
  geom_tile(aes(x =  Treatment,
                y = reg,
                fill = `Methylation Ratio`)) +
  scale_fill_gradient2(low = "white", 
                       high = "black") +
  scale_x_discrete("Treatment",
                   expand = c(0, 0),
                   breaks = levels(dt1$Treatment),
                   labels = c("Control",
                              "AOM+DSS",
                              "AOM+DSS\n+Curcumin")) + 
  scale_y_discrete("Gene",
                   expand = c(0, 0),
                   breaks = c("Tnf",
                              "Tnf1",
                              "Tnf2")) +
  ggtitle("") +
  # guides(fill = guide_legend(title.position = "top")) + 
  theme(legend.position = "top")
p1

tiff(filename = "tmp/aom_dss_cur_figure3d_heatmap.tiff",
     height = 7,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Figure 5A: convert heatmap to data----
p1 <- load.image("docs/renyi_figures_03222018/Fig 5A -heatmap 2.png")
plot(p1)
summary(p1)
class(p1)
p1
plot(grayscale(p1))

tiff(filename = "tmp/aom_dss_cur_figure5a_heatmap.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(grayscale(p1))
graphics.off()