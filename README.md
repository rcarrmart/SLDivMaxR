# SLDivMaxR v. 0.1.0
### A R package to calculate single allele divergence matrix in between and within taxa

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
### `Dist2DF.Long`: `dist` to `data.frame` long format conversion
A function to change the data class from `dist` to `data.frame` in long format. 
```{r}
Dist2DF.Long (inDist,
              Convert100 = FALSE)
```
The input of the function (`inDist`) is a matrix of class `dist`, as in the output of `ape::dist.dna`. The output can be converted to percentage if the input file is in proportions and `Convert100` is set to `TRUE`. Otherwise, when `Convert100 = FALSE` (the default) the output will remain the same.  

Example: 
```{r}

```

### `Dist2DF.Wide`: `dist` to `data.frame` wide format conversion
A function to change the data class from `dist` to `data.frame` in wide format. 
```{r}
Dist2DF.Wide (inDist,
              Convert100 = FALSE)
```
The input of the function (`inDist`) is a matrix of class `dist`, as in the output of `ape::dist.dna`. The output can be converted to percentage if the input file is in proportions and `Convert100` is set to `TRUE`. Otherwise, when `Convert100 = FALSE` (the default) the output will remain the same.  

Example: 
```{r}

```
### `Dist2DF.Groups`: Mean Distances between Group
Applies `hclust` and `cutree` to a matrix of class `"dist"` to cluster samples by a  maximum distance threshold or to a set number of groups.
```
Dist2DF.Groups(inDist,
               Convert100 = FALSE,
               get.group = NULL,
               n.groups = NULL)
```
`Dist2DF.Groups` takes a distance matrix (`inDist`) and calculates the mean distances between descreet grouping using maximum divergence threshold between sequences (argument `get.group`) or by creating a _N_ number of groups regardless their within group divergences (argument `n.groups`). 

Example: 
```{r}

```
### `DNA.Dist.df`: Single locus divergence data.frame
Applies `ape::dist.dna()` but changes `"dist"` class output from `ape::dist.dna()` to `data.frame`.

```{r}
DNA.Dist.df(inDNA,
            Dist.type = "percent",
            Out.Format = "long",
            Model = "K80",
            GAMMA = FALSE,
            PW.deletion = TRUE,
            Var = FALSE,
            by.group = NULL,
            group.summary = FALSE,
            within.group = FALSE)
```
Arguments: 

`inDNA` the input aligned sequences in `.fasta` format.

`Dist.type` specifies if the distance output should be in proportions (`"prop"`), percentages (`"percentage"`), or number of nucleotides that are different (`"nucleotide"`).

`Model` specifies the evolutionary model to use. Same as in `ape::dist.dna`. Options are: `"RAW"`, `"JC69"`, `"K80"` (default), `"F81"`, `"K81"`, `"F84"`, `"T92"`, `"TN93"`, `"GG95"`, `"LOGDET"`, `"BH87"`, `"PARALIN"`, `"N"`, `"TS"`, `"TV"`, `"INDEL"`, and `"INDELBLOCK"`. Note: if `Model` is `"N"` and `Dist.type` is set for `"prop"` or `"percentage"`, then an error will be return.

`GAMMA` same as `gamma` in `ape::dist.dna()`.

`PW.deletion` specifies if Pairwise deletion should be considered (`TURE`, default) or not (`FALSE`).

`Var` if `TRUE` output is the calculated variation (`FALSE` is default)

`by.group` a vector to group sequences on different categories (e.g., populations, species, lineages, etc.). The vector should follow the same order as the sequence alignment, by applying `vegan::meandist()`. The vector can be the `.csv` output of species delimitation analyses such as `ASAP` and `ABGD`.

`group.summary` applies `vegan::meandist()` to formulate a data.frame with within groups mean distance, between groups mean distance and overall distance. A vector with discreet groupings should be imported in `by.group`.

`within.group` applies `vegan::meandist()` to create a `data.frame` with each within group mean divergence, and number of samples ("N") by group. A vector with discreet groupings should be imported in `by.group`. Groupls with N = 1 will return `NA`'s.


Example:
```{r}
```


### `SE`: Standard Error 
Calculates standard error on a vector `x`. 
```{r}
SE(x)
```
`SE` automatically ignores missing data (`NA`'s). `x` must be a numeric vector. 

Example:
```{r}
a <- c(1:10)
SE(a)

a <- c(1:10, NA, NA)
SE(a)

a <- c(NA, NA, NA)
SE(a)
```

### `DNA.GroupDist.SE`: Divergence Standard Error between and within groups
Calculates standard error of specified between and within group genetic divergences. The funtion applies a modified version of `vegan::meandist` to calculate the Standard Error of the divergences calculated by `DNA.Dist.df`. 

```{r}
DNA.GroupDist.SE(inDNA,
                 Dist.type = "percent",
                 Out.Format = "long",
                 Model = "K80",
                 GAMMA = FALSE,
                 PW.deletion = TRUE,
                 by.group,
                 within.group = FALSE)
```
`inDNA` the imput aligned sequences in `.fasta` format.

`Dist.type` specifies if the distance output should be in proportions (`"prop"`), percentages (`"percentage"`), or number of nucleotides that are different (`"nucleotide"`).

`Model` specifies the evolutionary model to use. Same as in `ape::dist.dna`. Options are: `"RAW"`, `"JC69"`, `"K80"` (default), `"F81"`, `"K81"`, `"F84"`, `"T92"`, `"TN93"`, `"GG95"`, `"LOGDET"`, `"BH87"`, `"PARALIN"`, `"N"`, `"TS"`, `"TV"`, `"INDEL"`, and `"INDELBLOCK"`. Note: if `Model` is `"N"` and `Dist.type` is set for `"prop"` or `"percentage"`, then an error will be return.

`GAMMA` same as `gamma` in `ape::dist.dna()`.

`PW.deletion` specifies if Pairwise deletion should be considered (`TURE`, default) or not (`FALSE`).

`by.group` a vector to group sequences on different categories (e.g., populations, species, lineages, etc.). The vector should follow the same order as the sequence alignment. The vector can be the `.csv` output of species delimitation analyses such as `ASAP` and `ABGD`.

`within.group` creates a `data.frame` with each within group standard error, and number of saples ("N") by group. A vector with discreet groupings should be imported in `by.group`. Groups with N = 1 will return `NA`'s.


Example:
```{r}
```

## Citation
Please cite this package as follows: 

Carrera-Martínez, R (2025) SLDivMaxR: A R package to calculate between and within single allele divergence. Version 0.1.0 URL:<<https://github.com/rcarrmart/SLDivMaxR/>>. 

To cite `vegan`: 

Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P, O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, Barbour M, Bedward M, Bolker B, Borcard D, Borman T, Carvalho G, Chirico M, De Caceres M, Durand S, Evangelista H, FitzJohn R, Friendly M, Furneaux B, Hannigan G, Hill M, Lahti L, McGlinn D, Ouellette M, Ribeiro Cunha E, Smith T, Stier A, Ter Braak C, Weedon J (2024). vegan: Community Ecology Package. R package version 2.7-0, https: <<https://vegandevs.github.io/vegan/>>.

To cite `ape`: 

Paradis E, Schliep K (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. _Bioinformatics_, 35, 526-528. doi:10.1093/bioinformatics/bty633. 
