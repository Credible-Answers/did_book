library(haven)
datos <- read_dta("~/Downloads/Replication folder FE/programs/cc_xd_didtextbook_2025_9_30/Data sets/Wolfers 2006/wolfers2006_didtextbook.dta")
install.packages("panelView")
library(panelView)
names(datos)
# We have to see the structure of treatment 
table(datos$udl)
table(datos$early_late_never)
table(datos$cohort)
length(unique(datos$state))
# Treatment Status
panelview(div_rate ~ udl, data = datos,
          index = c("state", "year"),
          xlab = "Year", ylab = "State",
          display.all = TRUE, gridOff = TRUE, by.timing = TRUE,
          axis.lab.gap = c(5, 0),
          axis.lab = "both",
          axis.lab.angle = 45,
          background = "white",
          color = c("#7BB3E0", "#2A6496", "#D5D5D5"),
          legend.labs = c("Treated (Pre)", "Treated (Post)", "Control"),
          main = "Treatment Status: Wolfers (2006)",
          cex.main = 16, cex.axis.x = 8, cex.axis.y = 6,
          cex.lab = 14, cex.legend = 12)
# Cohort outcome
panelview(data = datos, Y = "div_rate", D = "udl",
          index = c("state", "year"),
          by.timing = TRUE, display.all = TRUE,
          type = "outcome", by.cohort = TRUE,
          xlab = "Year", ylab = "Divorce Rate",
          main = "Average Outcomes by Cohort: Wolfers (2006)",
          cex.main = 16, cex.axis = 10, cex.lab = 14, cex.legend = 12,
          axis.lab.gap = c(5, 0))
