.First.sys()

#https://gist.github.com/stevenworthington/3178163
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("MASS", "gplots", "infotheo", "lme4", "car", "rgl", "lattice", "ggplot2", "mlr")
ipak(packages)
