library(haven)
datos <- read_dta("~/Downloads/Replication folder FE/programs/cc_xd_didtextbook_2025_9_30/Data sets/Moser and Voena 2012/moser_voena_didtextbook.dta")
install.packages("panelView")
library(panelView)
names(datos)
datos$D <- datos$treatmentgroup * datos$post
length(unique(datos$subclass))

panelview(data = datos, Y = "patents", D = "D",
          index = c("subclass", "year"),
          by.timing = TRUE, display.all = TRUE,
          type = "outcome", by.cohort = TRUE,
          xlab = "Year", ylab = "Patents",
          main = "Average Outcomes by Cohort",
          cex.main = 16, cex.axis = 10, cex.lab = 12, cex.legend = 12,
          axis.lab.gap = c(5, 0))

subs <- unique(datos$subclass[order(-datos$D, datos$subclass)])[1:400]
datos400 <- datos[datos$subclass %in% subs, ]

panelview(patents ~ D, data = datos400,
          index = c("subclass", "year"),
          xlab = "Year", ylab = "Subclass",
          display.all = TRUE, gridOff = TRUE, by.timing = TRUE,
          axis.lab.gap = c(5, 0),
          axis.lab = "both",
          axis.lab.angle = 45,
          background = "white",
          color = c("#D5D5D5", "#4A90D9"),
          legend.labs = c("Control", "Treatment"),
          main = "Treatment Status (Top 400 units, sorted)",
          cex.main = 16, cex.axis = 3, cex.lab = 12, cex.legend = 12)