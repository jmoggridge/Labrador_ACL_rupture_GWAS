## Makes 2 plots from MDS data: 
# Fig1: first 2 coefficients
# Fig2: shows first 4 coefficients with 1d and 2d density

library(tidyverse)
library(janitor)
library(GGally)

# read mds data
mds <- read_delim("cr237_dryad_8.mds", delim = ' ') %>% 
  janitor::clean_names() %>% 
  mutate(across(c1:c10,  as.numeric)) %>% 
  select(id = iid, c1:c10)


# read phenotype data
dogs <- 
  read_delim("cr237_dryad_8.fam", delim = ' ', col_names = paste0('X', 1:6)) %>% 
  transmute(id = X1,
            sex = if_else(X5 == 1, 'M', 'F'),
            phenotype = case_when(
              X6 == 1 ~ 'Control', 
              X6 == 2 ~ 'ACL rupture', 
              TRUE ~ 'Missing')
  ) %>% 
  mutate_if(is_character, as_factor) %>% 
  # combine .fam and .mds data by id
 right_join(mds, by = 'id')


# plot mds
fig1 <- dogs %>% 
  ggplot(aes(c1, c2, color = phenotype, shape = sex)) +
  geom_point() +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(color = 'black', 
                                    size = 0.5, fill = NA)) +
  labs(x = 'MDS1', y = 'MDS2', 
       subtitle = paste("Metric multidimensional scaling of",
                        nrow(dogs),"dog genotypes"))
# fig1  


fig2 <- ggpairs(
  dogs, columns = (4:8),
  mapping = aes(color = phenotype),
  diag = list(continuous = wrap("densityDiag", alpha=0.5, size=0)),
  upper = list(continuous = wrap("density", alpha = 0.6),
               combo = "box_no_facet"),
  lower = list(continuous = wrap("points", alpha = 0.6, size=0.4), 
               combo = wrap("dot", alpha = 0.4, size=0.2)),
  legend = c(5,5)
  ) +
  theme_classic() +
  labs() + theme(legend.position = "bottom")


read_delim("result2.assoc.logistic")




write_rds(fig1, "MDS_figure1.rds")
write_rds(fig2, "MDS_figure2.rds")

# plot either of these into an rmarkdown with read_rds(file) 
read_rds("MDS_figure1.rds")
read_rds("MDS_figure2.rds")



