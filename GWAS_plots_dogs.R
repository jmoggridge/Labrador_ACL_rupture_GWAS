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

library(tidyverse)
# parse logistic regression results
logistic <- 
  read_delim("result2.assoc.logistic",delim = ' ') %>% 
  janitor::clean_names() %>% 
  mutate_all(trimws) %>% 
  transmute(chr, bp, p) %>% 
  mutate_all(as.numeric) 

glimpse(logistic)
library(patchwork)

man0 <- 
  logistic %>% 
  filter(chr <11) %>% 
  ggplot(aes(bp, -log(p, 10))) +
  geom_point(alpha = 0.5, size = 0.2) +
  facet_wrap(~chr, nrow = 1, scales = 'free_x') +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, 'lines')) +
  labs(x ='', y = expression(-log[10](p)))

man1 <-   logistic %>% 
  filter(chr > 10 & chr <21) %>% 
  ggplot(aes(bp, -log(p, 10))) +
  geom_point(alpha = 0.5, size = 0.2) +
  facet_wrap(~chr, nrow = 1, scales = 'free_x') +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, 'lines')) +
  labs(x = '', y = expression(-log[10](p)))

man2 <-   logistic %>% 
  filter(chr > 20 & chr <31) %>% 
  ggplot(aes(bp, -log(p, 10))) +
  geom_point(alpha = 0.5, size = 0.2) +
  facet_wrap(~chr, nrow = 1, scales = 'free_x') +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, 'lines')) +
  labs(x = '', y = expression(-log[10](p)))

man3 <-   logistic %>% 
  filter(chr > 30) %>% 
  ggplot(aes(bp, -log(p, 10))) +
  geom_point(alpha = 0.5, size = 0.2) +
  facet_wrap(~chr, nrow = 1, scales = 'free_x') +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, 'lines')) +
  labs(x = 'position', y = expression(-log[10](p)))

fig3_manhattan <- man0/man1/man2/man3 &
  plot_annotation(
    title = "Manhattan plot of logistic regression tests")


### write files

write_rds(fig1, "Figure1_mds1.rds")
write_rds(fig2, "Figure2_mds2.rds")
write_rds(fig3_manhattan, "Figure3_mahattan.rds")

# plot either of these into an rmarkdown with read_rds(file) 
read_rds("Figure1_mds1.rds")
read_rds("Figure2_mds2.rds")
read_rds("Figure3_mahattan.rds")



