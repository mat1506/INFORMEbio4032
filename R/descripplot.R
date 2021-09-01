
library(here)
library(tidyverse)
library(ggplot2)
library(corrr)
library(faraway)
library("factoextra")
library("ggpubr")
library(rrcov)

y1 <- read_csv(here::here("analysis/data/raw_data/2021-08-31_informeMF2021_data_V01.csv"))


cory1 <- y1  %>%
  mutate(silt= silt_coarse + silt_fine,
         TOC.TN = (TOC/12)/(TN/14)) %>%
  dplyr::select(clay,silt,sand,
                TOC,TIC,TC,TN,TOC.TN,
                d13C, d18O) %>%
  rename(`Clay`=clay,
         `Silt`=silt,
         `Sand`=sand,
         `TOC/TN`=TOC.TN,
         `d13C`=d13C,
         `d18O`=d18O)

#variance inflation factors
#cuantifica la intensidad de la multicolinealidad en
#un análisis de regresión normal de mínimos cuadrados.
#Proporciona un índice que mide hasta qué punto
#la varianza de un coeficiente de regresión estimado
#se incrementa a causa de la colinealidad.

#1 = not correlated.
#Between 1 and 5 = moderately correlated.
#Greater than 5 = highly correlated.


sort(vif(cory1))

#####  Correlation plots

y1.cor <- corrr::correlate(cory1)
y1.cor %>% fashion()
p<- network_plot(y1.cor,min_cor = 0.4)

####  Principal component analysis
pcy1 <- y1  %>%
  mutate(silt= silt_coarse + silt_fine,
         TOC.TN = (TOC/12)/(TN/14)) %>%
  dplyr::select(clay,silt,
                TOC,TIC,TOC.TN,
                d13C, d18O) %>%
  rename(`Clay`=clay,
         `Silt`=silt,
         `TOC/TN`=TOC.TN,
         `d13C`=d13C,
         `d18O`=d18O)


y1.pca <- prcomp(x = pcy1, scale. = T)

p2 <- fviz_pca_biplot(y1.pca, labelsize = 6,
                     # Individuals
                     geom.ind = "point",col.ind = "black",
                     pointshape = 21, pointsize = 1.5,
                     palette = "Set3",
                     addEllipses = TRUE,
                     ellipse.type = "convex",
                     ellipse.level = 0.95,
                     ellipse.alpha = 0.25,
                     # Variables
                     col.var = "contrib",
                     legend.title = list(color = "Contrib", alpha = "Contrib")
) + scale_color_gradient2(low="blue",mid="black",
                          high="blue", midpoint=16)

p3 <- ggplot(data = get_eig(X = y1.pca),
             aes(x = as.factor(1:nrow(get_eig(X =y1.pca))),
                 y = cumulative.variance.percent,
                 group = 1)
) +
  geom_col(fill = "steelblue", color = "steelblue") +
  geom_line(color = "black", linetype = "solid") +
  geom_point(shape = 19, color = "black") +
  geom_text(aes(label = paste0(round(cumulative.variance.percent, 1), "%")),
            size = 3, vjust = -0.5, hjust = 0.7) +
  geom_hline(yintercept = 90, color = "firebrick", linetype = "dashed") +
  labs(title = "Accumulated variance (%)",
       x = "PC number",
       y = "Porcentage of explained variances") +
  theme_bw()


# Contribución a la primera componente.
p4 <- fviz_contrib(X = y1.pca, choice = "var", axes = 1, top = 10)
# Contribución a la segunda componente.
p5 <- fviz_contrib(X = y1.pca, choice = "var", axes = 2, top = 10)

Fig_3 <- ggpubr::ggarrange(ggarrange(p2, p3, ncol = 2, labels = c("a)", "b)"), align = "h",widths = c(2.8,2.6)),
                               ggarrange(p4, p5, ncol = 2, labels = c("c)", "d)"), align = "h",widths = c(3.5,3.5)),
                               nrow = 2,
                               heights = c(1, 1)
)
annotate_figure(Fig_3, top = text_grob("PCA Geoquimico sedimentary sequence", color = "black", face = "bold", size = 9))

ggsave(here::here("analysis/figures/Fig_3.pdf"))

#PCA robust XRF values

hpca<- PcaHubert(pcy1, k=7, scale = T)
rpca_eigen_var <- tibble(
  componente = seq_len(hpca$k),
  varianza   = rrcov::getSdev(hpca)^2,
  porcnt_varianza = 100 * varianza / sum(varianza),
  porcnt_varianza_acumulada = cumsum(porcnt_varianza)
)
#rpca_eigen_var
#Explained variances PCA robust plot
p6<-ggplot(data = rpca_eigen_var,
          aes(x = as.factor(componente),
              y = porcnt_varianza_acumulada,
              group = 1)
) +
  geom_col(fill = "steelblue", color = "steelblue") +
  geom_line(color = "black", linetype = "solid") +
  geom_point(shape = 19, color = "black") +
  geom_text(aes(label = paste0(round(porcnt_varianza_acumulada), "%")),
            size = 3, vjust = -0.5, hjust = 0.7) +
  geom_hline(yintercept = 90, color = "firebrick", linetype = "dashed") +
  lims(y = c(0, 100)) +
  labs(title = "Accumulated variance (%)",
       x = "PCs",
       y = "Explained variances (%)") +
  theme_bw()


ggpubr::ggarrange(p6, ncol = 1)

#PCA robusto
biplot(hpca, main="Figure 12S.Robust biplot", col=c("gray","blue"), scale =0, pc.biplot = FALSE)

#guardar

if(!file.exists(here::here("analysis/data/derived_data/datapca.csv"))){
  datapca<-as.data.frame(cbind(y1[,2:13],hpca@scores[,1:2]))
  path_out <- here::here("analysis/data/derived_data/","datapca.csv")
  write_csv(datapca,path_out)
}
