
# Human CRISPR Knockout Pooled Library (GeCKO v2)
# addgene catalog# (Pooled Library #1000000048, #1000000049) 
# Improved vectors and genome-wide libraries for CRISPR screening . 
# Sanjana NE, Shalem O, Zhang F. Nat Methods. 2014 Aug;11(8):783-4. 
# doi: 10.1038/nmeth.3047. 10.1038/nmeth.3047 PubMed 25075903

# This R script was used after the NGS reads extracted by mageck software.
# Li W , Xu H , Xiao T , et al. 
# MAGeCK enables robust identification of essential genes from genome-scale 
# CRISPR/Cas9 knockout screens[J]. Genome Biology, 2014, 15(12):554.

# Detecting the required packages -----------------------------------------
# To detect the required packages if they are not installed, 
# the packages will installed then.

if(!require('tidyverse')){
  install.packages('tidyverse')
  library(tidyverse)
}
if(!require('compiler')){
  install.packages('compiler')
  library(compiler)
}

# tidy and mark the NGS reads count data ----------------------------------



# plot functions ----------------------------------------------------------

NGS_plot <- function(nm = '', readcount = '', i = ''){
  obj <- filter(readcount, Gene == nm)
  
  NonTargetingCon <- filter(readcount, str_detect(Gene,'^Non'))
  NonTargetingCon$NonTargetingCon <- 'NonTargettingControlGuideRNA'
  
  label_for_lines <- tribble(
    ~x,~y,~z,
    1, 100,'100-fold',
    1, 50, '50-fold',
    1, 20, '20-fold',
    1, 10, '10-fold',
    1, 5, '5-fold',
    1, 3, '3-fold',
    1, 1, 'equal')
  gg <- ggplot(readcount, aes(x = con, y = GFP)) +
    geom_point(aes(color = readcount,
                   shape = readcount),
               size =0.5,
               alpha = (1/5)) +
    geom_abline(intercept = 0, slope = 1, colour = 'white', size = 0.5) +
    geom_abline(intercept = log10(3), slope = 1, colour = 'white', size = 0.5) +
    geom_abline(intercept = log10(5), slope = 1, colour = 'red', size = 0.5) +
    geom_abline(intercept = log10(10), slope = 1, colour = 'white', size = 0.5) +
    geom_abline(intercept = log10(20), slope = 1, colour = 'white', size = 0.5) +
    geom_abline(intercept = log10(50), slope = 1, colour = 'white', size = 0.5) +
    geom_abline(intercept = log10(100), slope = 1, colour = 'white', size = 0.5) +
    geom_text(data = label_for_lines, aes(x,y,label = z), angle = 27.5, size = 3) +
    geom_point(data = NonTargetingCon, 
               aes(con, 
                   GFP, 
                   color = NonTargetingCon,
                   shape = NonTargetingCon),
               size = 1,
               alpha = (1/2)) +
    geom_point(data = obj, 
               aes(con, 
                   GFP, 
                   color = Gene,
                   shape = Gene),
               size = 3) +
    labs(y = ' Sorted GFP- sgRNA reads', 
         x = ' sgRNA reads', 
         title = 'Enriched gene by FACs sorted'
    ) +
    scale_color_manual(name = 'sgRNAs', 
                       limits = c('Total','NonTargettingControlGuideRNA',nm),
                       label = c('Total','NT Control',nm), 
                       values = c('grey50','black','blue')) +
    scale_shape_manual(name = 'sgRNAs',
                       limits = c('Total','NonTargettingControlGuideRNA',nm),
                       label = c('Total','NT Control',nm),
                       values = c(10,17,15)) +
    scale_x_log10(breaks = c(1, 10, 40,60,80,100, 200,400,
                             600,1000, 10000, 100000),
                  labels = c('1', '10', '40','60','80','100','200','400',
                             '600','1000', '10000', '100000')) +
    scale_y_log10(breaks = c(1, 10, 100, 200,400,600,1000, 
                             2000,4000,6000,
                             10000, 100000),
                  labels = c('1', '10', '100','200','400','600','1000',
                             '2000','4000','6000',
                             '10000', '100000')) + 
    theme(axis.title = element_text(family = "Helvetica"), 
          axis.text = element_text(family = "Helvetica"), 
          axis.text.x = element_text(family = "Helvetica"), 
          axis.text.y = element_text(family = "Helvetica"), 
          plot.title = element_text(family = "Helvetica", hjust = 0.5), 
          legend.text = element_text(family = "Helvetica"), 
          legend.title = element_text(family = "Helvetica")) +
    guides(color = guide_legend(override.aes = list(alpha = 1,size = 4)))
  ggsave(str_c(i, '.', nm, ".GFPminus Sorted sgRNA reads ",".tiff"),width = 10,height = 6.18)
  
}

