# Mitochondrial ubiquitin landscape in neurons

Software to accompany manuscript [Antico et al., "Global ubiquitylation analysis of mitochondria in primary neurons identifies endogenous Parkin targets following activation of PINK1", Sci. Adv. 7, eabj0722 (2021)](https://www.science.org/doi/10.1126/sciadv.abj0722).

## Usage

We suggest using RStudio. Start in the top project directory. The first step is to create environment using 'renv':

```
install.packages("renv")
renv::restore()
```

This will install all necessary packages. Run the `targets` pipeline.

```
targets::tar_make()
```

This will analyse data, create figures (in directory `fig`) and a report (in directory `doc`). 
