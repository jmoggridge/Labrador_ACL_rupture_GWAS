## Makes 2 plots from MDS data: 
# Fig1: first 2 coefficients
# Fig2: shows first 4 coefficients with 1d and 2d density
# Fig3: manhattan plot
# Fig4: qq-plot

library(tidyverse)
library(janitor)
library(GGally)
library(ggrepel)
library(patchwork)
library(qqman)

# read mds data
mds <- read_delim("cr237_dryad_8.mds", delim = ' ') %>% 
  janitor::clean_names() %>% 
  mutate(across(c1:c10,  as.numeric)) %>% 
  select(id = iid, c1:c10)

# read phenotype data
dogs <- 
  read_delim("cr237_dryad_8.fam", delim = ' ', 
             col_names = paste0('X', 1:6)) %>% 
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



# parse logistic regression results
logistic <- 
  read_delim("result2.assoc.logistic",delim = ' ') %>% 
  janitor::clean_names() %>% 
  mutate_all(trimws) %>% 
  transmute(chr, bp, p, snp) %>% 
  mutate(across(c(chr, bp, p), as.numeric))

glimpse(logistic)



manhattan(logistic, chr = 'chr', bp = 'bp' , p = 'p', snp = 'snp', 
          annotatePval = 0.0001,
          main = "Manhattan plot for logistic regression tests")

qq(logistic$p, main = "Q-Q plot of GWAS p-values from logisitic regression")


### write files

write_rds(fig1, "Figure1_mds1.rds")
write_rds(fig2, "Figure2_mds2.rds")
write_rds(fig3_manhattan, "Figure3_mahattan.rds")

# plot either of these into an rmarkdown with read_rds(file) 
read_rds("Figure1_mds1.rds")
read_rds("Figure2_mds2.rds")
read_rds("Figure3_mahattan.rds")



## ggplot manhattan 
# (adapted from https://www.r-graph-gallery.com/101_Manhattan_plot.html)


gg_manhattan_plot <- function(assoc_data, snpsOfInterest = NULL){

  assoc_data <- assoc_data %>% 
    group_by(chr) %>% 
    summarise(chr_len = max(bp)) %>% 
    # cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    # join cumulative position to initial data
    left_join(logistic, by = 'chr') %>%
    # get the cumulative position of each SNP
    arrange(chr, bp) %>%
    mutate(BPcum = bp + tot) %>%
    # Add highlight and annotation information
    mutate(is_highlight = ifelse(snp %in% snpsOfInterest, "yes", "no")) %>%
    mutate(is_annotate = ifelse(-log10(p) > 4, "yes", "no"))
  
  # Prepare X axis breaks + labels
  axisdf <- assoc_data %>% 
    group_by(chr) %>%
    summarize(center=(max(BPcum) + min(BPcum)) / 2)
  
  # Make the manhattan plot showing -logP ~ position:
  ggplot(assoc_data, aes(BPcum, -log10(p))) +
    geom_point(aes(color = as.factor(chr)), alpha = 0.8, size = 1) +
    scale_color_manual(values = rep(c("grey35", "skyblue"), 100)) +
    # Add highlighted points
    geom_point(
      data = assoc_data %>% filter(is_highlight == "yes"),
      color = "orange",
      size = 2) +
    # Add labels using ggrepel to avoid overlap
    geom_text_repel(
      data = assoc_data %>% filter(is_annotate == "yes"),
      aes(label = snp),
      size = 2) +
    # Set x-axis breaks + labels
    scale_x_continuous(
      label = axisdf$chr, 
      breaks = axisdf$center,
      expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.05, 0.05)) +
    # Custom the theme:
    labs(x = 'position') +
    theme_bw() +
    theme(legend.position = "none", 
          axis.text.x = element_text(size = 7),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
    )
}


# create a vector of snps to we want to highlight
snpsOfInterest <- logistic %>% filter(p < 10e-6) %>% pull(snp)

gg_manhattan_plot(logistic, snpsOfInterest) +
  scale_y_continuous(expand = c(0.05, 0.05), limits = c(0, 5.5),
                     breaks = seq(0, 5.5, 0.5)) +
  labs(subtitle = "Logistic regression test for SNP association")


# QQ plot

ci <- 0.95
nSNPs <- nrow(logistic)

# create a sorted vector of p-values,
# use the ppoints function to get the expected distribution
plotdata <- tibble(
  observed = sort(logistic$p),
  expected = ppoints(nSNPs),
  ci_lower   = qbeta(p = (1 - ci) / 2, 
                   shape1 = seq(nSNPs), 
                   shape2 = rev(seq(nSNPs))),
  ci_upper   = qbeta(p = (1 + ci) / 2, 
                   shape1 = seq(nSNPs), 
                   shape2 = rev(seq(nSNPs)))
) %>% 
  mutate_all(~ -log10(.x))

plotdata %>% 
  ggplot(aes(expected, observed)) +
  geom_abline(intercept = 0, slope = 1, lty = 2, color = 'gray') +
  geom_point(shape = 1, ) +
  coord_equal() +
  theme_classic() +
  labs(
    x = expression(Expected ~ -log[10](p)),
    y = expression(Observed ~ -log[10](p)),
    title = "Q-Q plot for logistic regression association model"
    )

