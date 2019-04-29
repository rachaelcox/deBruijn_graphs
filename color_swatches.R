library(tidyverse)
library(RColorBrewer)
library(colorspace)

colfunc <- colorRampPalette(c("dark blue", "light green", "white"))


palette <- colfunc(2000)
swatchplot(palette)

palette_LGL_format <- as_tibble(t(col2rgb(palette))) %>% rownames_to_column() %>%
  mutate(red = red/255) %>%
  mutate(green = green/255) %>%
  mutate(blue = blue/255)

