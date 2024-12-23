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

Carrera-Martínez, R (2025) SLDivMaxR: A R package to calculate between and within single allele divergence. Version 0.1.0 URL:<<https://github.com/rcarrmart/SLDivMaxR/>>. 

To cite `vegan`: 

Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P, O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, Barbour M, Bedward M, Bolker B, Borcard D, Borman T, Carvalho G, Chirico M, De Caceres M, Durand S, Evangelista H, FitzJohn R, Friendly M, Furneaux B, Hannigan G, Hill M, Lahti L, McGlinn D, Ouellette M, Ribeiro Cunha E, Smith T, Stier A, Ter Braak C, Weedon J (2024). vegan: Community Ecology Package. R package version 2.7-0, https://github.com/vegandevs/vegan, <<https://vegandevs.github.io/vegan/>>.

To cite `ape`: 

Paradis E, Schliep K (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. _Bioinformatics_, 35, 526-528. doi:10.1093/bioinformatics/bty633. 
