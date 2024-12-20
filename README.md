# SLDivMaxR v. 0.1.0
### A R package to calculate between and within single allele divergence

## Package Summary
`SLDivMaxR` is a simple package that allows calculation of divergence (or distance) on a single allele at each sequence level, by defined groups (which can be species, subspecies, lineagesm, haplotypes, etc.) and within each defined groups. The package was designed to analize large datasets, increase replicability, and produce outputs in a `data.frame` format. `SLDivMaxR` is built upon the `ape` funtion `dist.dna` to calculate genetic divergence and the `vegan` funtion `meandist` to calculate mean divergence between and within groups. This approach shortcots manual designation of groups in other commonly used softwares, and has potential to be used sistematically in large datasets and as part of more complex population genetics analises. 

## Instalation
Make sure to have `devtools` package installed. 
```{r}
install.packages("devtools") 
```
To install `SLDivMaxR` simply run: 
```{r}
devtools::install_github("rcarrmart/SLDivMaxR") 
```
Packages `ape` and `vegan` will be needed. 

## Functions Overlook
The package currently has six funtions. 

The funtion `DNA.GroupDist.SE` applies a modified version of `meandist` to calculate the Standard Error of the divergences calculated by `DNA.Dist.df`. 

## Citation
Please cite this package as follows: 
Carrera-Martínez, R (2025) SLDivMaxR: A R package to calculate between and within single allele divergence. Version 0.1.0 url<<https://github.com/rcarrmart/SLDivMaxR/>>. 
