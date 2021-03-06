
library(here)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(tidypaleo)
theme_set(theme_bw(8))

y1 <- read_csv(here::here("analysis/data/raw_data/2021-08-31_informeMF2021_data_V01.csv"))
head(y1)

y1_longer <- y1 %>%
  mutate(silt= silt_coarse + silt_fine,
         TOC.TN = (TOC/12)/(TN/14)) %>%
  dplyr::select(depth,clay,silt,sand,
                TOC,TIC,TC,TN,TOC.TN,
                d13C, d18O) %>%
  rename(Depth = depth,`Clay`=clay,
         `Silt`=silt,
         `Sand`=sand,
         `TOC/TN`=TOC.TN,
         `d13C`=d13C,
         `d18O`=d18O) %>%
  pivot_longer(cols =Clay:d18O,names_to = "values", values_to = "count") %>%
  filter(values %in% c("Sand", "Silt","Clay",
                       "TOC","TIC","TC","TN","TOC/TN","d13C", "d18O")) %>%
  mutate(values = fct_relevel(values,"Sand", "Silt",
                              "Clay","TC","TIC","TOC","TN","TOC/TN","d13C", "d18O"))

zone_data <- tibble(ymin = c(38,140,240),
                    ymax = c(50,170,268),
                    xmin = -Inf, xmax = Inf)  ###### Stratiplot zone

coniss <- y1_longer %>% ###### Cluster plot
  nested_data(qualifiers = Depth, key = values, value = count, trans = scale) %>%
  nested_chclust_coniss()

alta_plot <- ggplot(y1_longer,aes(x = count, y = Depth)) +
  geom_lineh() +
  geom_point() +
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = zone_data,
            alpha = 0.4,
            fill = "blue",inherit.aes = FALSE ) +
  scale_y_reverse(breaks = c(40,60,80,100,120,140,
                             160,180,200,220,240,260,280)) +
  facet_geochem_gridh(
    vars(values),
    units = c("d13C" = "‰", "d18O" = "‰"),
    default_units = "%") +
  labs(x = NULL, y = "Depth (cm)") +
  scale_x_continuous(breaks = scales::breaks_extended(n = 3)) +
  theme_set(theme_paleo(8)) +
  theme(
    text = element_text(size=20),
    axis.ticks = element_line(colour = "grey70", size = 0.3),
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "grey50")
  )

wrap_plots(
  alta_plot +
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  ggplot() +
    layer_zone_boundaries(coniss, aes(y = Depth)) +
    layer_dendrogram(coniss, aes(y = Depth), param = "CONISS") +
    scale_x_continuous(breaks = scales::breaks_extended(n = 3)) +
    theme(axis.text.y.left = element_blank(),
          axis.ticks.y.left = element_blank(),
          text = element_text(size=15),
          panel.background = element_rect(fill = "white", colour = "grey50"))+
    labs(x = "coniss", y = NULL),
  nrow = 1,
  widths = c(8, 0.8)
)


png_out <- here::here("analysis/figures","Fig2.png")
ggsave(png_out,  width = 25,height = 12, units = 'cm',device = "png")


y2 <- read_csv(here::here("analysis/data/raw_data/2021-08-31_informeMF2021_polen_V01.csv"))
head(y2)

data.porcen <- (y2[,10:136]/(rowSums(y2[,10:136])))*100 # Calculate %
data.sum <- colSums(data.porcen) # Calculate column sums
data.pollen <- data.porcen[, which(data.sum > 10)] # Subset data 10%

data1 <- cbind(y2[,5],data.pollen)

## reordenamos los datos en un formato de tabla larga
y2data<- pivot_longer(data1,cols =Adiantum:Apiaceae_tipo_Bowlesia, names_to = "species", values_to = "rel_abund")

y2_plot <- ggplot(y2data, aes(x = rel_abund, y = median)) +
  geom_col_segsh() +
  geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2) +
  geom_point() +
  scale_y_reverse() +
  facet_abundanceh(vars(species)) +
  labs(x = "Abundence %", y = "Age")

keji_coniss <- y2data %>%
  nested_data(qualifiers = median, key = species, value = rel_abund, trans = scale) %>%
  nested_chclust_coniss()


library(patchwork)
wrap_plots(
  y2_plot +
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  ggplot() +
    layer_dendrogram(keji_coniss, aes(y = median)) +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank())+
    labs(y = NULL),
  nrow = 1,
  widths = c(8, 0.8)
)


