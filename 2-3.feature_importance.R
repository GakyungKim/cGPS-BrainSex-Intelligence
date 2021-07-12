#######################################################################
##################### 2-3. Feature Importance #####################
#######################################################################

### install packages
#install.packages("remotes")
remotes::install_github("LCBC-UiO/ggseg", build_vignettes = F)
##https://github.com/LCBC-UiO/ggseg
remotes::install_github("LCBC-UiO/ggseg3d", build_vignettes = T)
remotes::install_github("LCBC-UiO/ggsegDesterieux", build_vignettes = T)
remotes::install_github("LCBC-UiO/ggsegExtra", build_vignettes = T)

### load packages
library(ggseg)
library(ggseg3d)
library(ggsegDesterieux)
library(ggsegExtra)
library(ggplot2)

# check the names of atlas
atlas_name = ggseg_atlas_repos()

# make the repository for chosen atlas
dk_rep = ggseg_atlas_repos('desterieux', ignore.case = T)
dk_region <- ggsegDesterieux::desterieux_3d$ggseg_3d
dk_name <- data.frame(region=dk_region[[1]]["region"],label=dk_region[[1]]["label"])

## default atlas of dk
ggseg(atlas = dk)
dk$region[grep("insula",dk$region)]
aseg$region[grep("cingulate",aseg$region)]

### collect important features to classify sex
brain = data.frame(
  #          3        5       6        7      9       10
  hemi = c("left", "right","right", "right", "left", "left",
  #           2      4       9   
           "left", "left", "left"),
  #                    3            5             6            7                    9                  10
  region = c("postcentral", "middle temporal","insula", "parahippocampal", "superior temporal","transverse temporal",
  #                     2                4                 9   
             "lateral occipital", "precuneus", "rostral middle frontal"),
  weight = c(0.031, 0.027, 0.028, 0.027, 0.023, 0.023,
             -0.034, -0.027, -0.023),
  stringsAsFactors = FALSE)

brain2 = data.frame(
  #           5,      8      10
  hemi = c("right","left", "left"),
  
  #                  5,            8         10                 
  region = c("thalamus proper", "putamen", "amygdala"),
             #      3,5                 7
             #"cerebellum cortex", "4th ventricle"),
  #            5         8        10       
  weight = c("0.028","-0.024", "-0.022"),
             #3 ,5          7      
             #"-0.0275",  "-0.025"),
  stringAsFactors = FALSE
)


brain3 = data.frame(
  region = c("cerebellum cortex", "4th ventricle"),
  weight = c("-0.0275", "-0.025"),
  stringAsFactors = FALSE
)

### visualize brain
ggseg(.data=brain, mapping=aes(fill=weight), color = 'white', position = "stacked")+
  scale_fill_gradient2(limits = c(-0.04, 0.04), breaks = seq(-0.04, 0.04, 0.01), low = "blue", high = "red", mid = "white") +
  theme(legend.position = 'bottom', legend.key.width = unit(2.1, 'cm'),
        legend.title = element_blank())

ggseg(.data=brain2, atlas = "aseg",
      mapping=aes(fill=as.numeric(weight)), color = 'white')+
  scale_fill_gradient2(limits = c(-0.04, 0.00), breaks = seq(-0.04, 0.04, 0.01), low = "blue", high = "red", mid = "white") +
  theme(legend.position = 'bottom', legend.key.width = unit(2.1, 'cm'),
        legend.title = element_blank())

ggseg(.data=brain3, atlas = "aseg" ,mapping=aes(fill=as.numeric(weight)), color = 'white')+
  scale_fill_gradient2(limits = c(-0.04, 0.00), breaks = seq(-0.04, 0.00, 0.01), low = "blue", high = "red", mid = "white") +
  theme(legend.position = 'bottom', legend.key.width = unit(2.1, 'cm'),
        legend.title = element_blank())
