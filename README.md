# SLDivMaxR v. 0.1.0
### A R Package for Calculating Single Allele Divergence Matrices Between and Within Taxa

## Package Summary
`SLDivMaxR` is a straightforward R package designed to calculate genetic divergence (or distance) for single alleles at each sequence level, defined by user-specified groups (such as species, subspecies, lineages, haplotypes, etc.), as well as within each defined group. The package is optimized for analyzing large datasets, enhancing replicability, and producing outputs in a `data.frame` format. Built upon the `ape` function `dist.dna` to compute genetic divergence and the `vegan` function `meandist` for calculating mean divergence between and within groups, SLDivMaxR simplifies the manual designation of groups required by other commonly used software. This makes it a useful tool for large datasets and can be seamlessly integrated into more complex population genetics analyses.

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
Dist2DF.Long(inDist,
             Convert100 = FALSE)
```
The input for the function (`inDist`) should be a matrix of class `dist`, such as the output from `ape::dist.dna`. The output can be converted to percentages if the input is in proportions and `Convert100 = TRUE`. By default, when `Convert100 = FALSE`, the output will remain unchanged.


Example: 
```{r}
library(ape)
test.dist <- dist.dna(test.DNA)

Dist.df <- Dist2DF.Long(test.dist)

Dist.df.100 <- Dist2DF.Long(test.dist, Convert100 = TRUE)

Dist.df.100x <- Dist2DF.Long(test.dist*100) # equivalent to `Convert100 = TRUE
```

### `Dist2DF.Wide`: `dist` to `data.frame` wide format conversion
A function to change the data class from `dist` to `data.frame` in wide format. 
```{r}
Dist2DF.Wide(inDist,
             Convert100 = FALSE)
```
The input for the function (`inDist`) should be a matrix of class `dist`, such as the output from `ape::dist.dna`. The output can be converted to percentages if the input is in proportions and `Convert100 = TRUE`. By default, when `Convert100 = FALSE`, the output will remain unchanged.

Example: 
```{r}
library(ape)
test.dist <- dist.dna(test.DNA)

Dist.df <- Dist2DF.Wide(test.dist)

Dist.df.100 <- Dist2DF.Wide(test.dist, Convert100 = TRUE)

Dist.df.100x <- Dist2DF.Wide(test.dist*100) # equivalent to `Convert100 = TRUE
```
### `Dist2DF.Groups`: Mean Distances between Group
Applies `hclust` and `cutree` to a matrix of class `"dist"` to cluster samples by a  maximum distance threshold or to a set number of groups.
```
Dist2DF.Groups(inDist,
               Convert100 = FALSE,
               get.group = NULL,
               n.groups = NULL,
               method = "PW",
               G.Name = "group")
````
`Dist2DF.Groups` takes a distance matrix (`inDist`) and calculates the mean distances between distinct groups. It does this either by using the maximum divergence threshold between sequences (via the `get.group` argument, with either `method` options) or by creating _N_ groups, regardless of within-group divergences (via the `n.groups` argument, and `method = "CT"`). 

Example: 
```{r}
test.dist <- dist.dna(test.DNA)

# threshold of 0.7 
dist.g.prop0.7 <- Dist2DF.Groups(test.dist, get.group = 0.7) 

# threshold of 7%
dist.g.per7 <- Dist2DF.Groups(test.dist, get.group = 7, Convert100 = TRUE)  

# threshold of 10% using the `hclust` and `cutree` approach
# Setting group names to "Species"
bySpecies <- Dist2DF.Groups(test.dist, get.group = 10, Convert100 = T,
                             method = "CT", G.Name = "Species")

# Cluster to 3 groups
Dist2DF.Groups(test.dist, n.groups = 3, Convert100 = T, method = "CT")
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
**Arguments: **
**`inDNA`**: The input aligned sequences in `.fasta` format.

**`Dist.type`**: Specifies the distance output format. Options are: `"prop"`: Proportions; `"percentage"`: Percentages; and `"nucleotide"`: Number of nucleotides that are different.

**`Model`**: Specifies the evolutionary model to use, as in `ape::dist.dna`. Options include: `"RAW"`, `"JC69"`, `"K80"` (default), `"F81"`, `"K81"`, `"F84"`, `"T92"`, `"TN93"`, `"GG95"`, `"LOGDET"`, `"BH87"`, `"PARALIN"`, `"N"`, `"TS"`, `"TV"`, `"INDEL"`, and `"INDELBLOCK"`. Note: If `Model = "N"` and `Dist.type` is set to `"prop"` or `"percentage"`, an error will be returned.

**`GAMMA`**: Same as the `gamma` parameter in `ape::dist.dna()`.

**`PW.deletion`**: Specifies whether pairwise deletion should be considered (`TRUE`, default) or not (`FALSE`).

**`Var`**: If `TRUE`, the output will include calculated variation. Default is `FALSE`.

**`by.group`**: A vector for grouping sequences into different categories (e.g., populations, species, lineages). The order should match the sequence alignment and be used with `vegan::meandist()`.

**`group.summary`**: Applies `vegan::meandist()` to generate a `data.frame` with within-group mean distance, between-group mean distance, and overall distance. A vector with discrete groupings should be provided in `by.group`.

**`within.group`**: Applies `vegan::meandist()` to create a `data.frame` with each group’s mean divergence and the number of samples ("N") by group. A vector with discrete groupings should be provided in `by.group`. Groups with `N = 1` will return `NA`s.

Examples:
```{r}
# Distance matrix class == data.frame of all specimens
DNA.Dist.df(test.DNA, Out.Format = "wide")

# same as above but with JC96 model
DNA.Dist.df(test.DNA, Out.Format = "wide", Model = "JC69")

# as above, but long pair-wise comparison format
DNA.Dist.df(test.DNA, Out.Format = "long", Model = "JC69")

####
# Defining groupings
####
#create data.frame with specimen and grouping (e.g., species) IDs
Species.info <- data.frame(Specimen = paste0("Ind",LETTERS[1:9]),
                           Species = paste0("Sp",c(rep(1,3),rep(2,3),rep(3,3))),
                            Seq.names = names(test.DNA))

Species.info  # Check that they match.
# The funtion uses the order of the `by.group` vector to group the sequences

DNA.Dist.df(test.DNA, Out.Format = "long", Model = "JC69", by.group = Species.info$Species)

# Reorganize the output 
Comparison <- DNA.Dist.df(test.DNA, Out.Format = "long", Model = "JC69", by.group = Species.info$Species)
Comparison <- data.frame(Pairs = paste(Comparison$x1, Comparison$x2, sep = "_vs_"),
                         Distance = Comparison$Distance)
Comparison

#wide format
WDist <- DNA.Dist.df(test.DNA, Out.Format = "wide", Model = "JC69", by.group = Species.info$Species)
WDist
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
**Arguments:**
**`inDNA`**: The input aligned sequences in `.fasta` format.

**`Dist.type`**: Specifies the distance output format. Options are `"prop"`: Proportions; `"percentage"`: Percentages; and `"nucleotide"`: Number of nucleotides that differ.

**`Model`**: Specifies the evolutionary model to use, as in `ape::dist.dna`. Options include: `"RAW"`, `"JC69"`, `"K80"` (default), `"F81"`, `"K81"`, `"F84"`, `"T92"`, `"TN93"`, `"GG95"`, `"LOGDET"`, `"BH87"`, `"PARALIN"`, `"N"`, `"TS"`, `"TV"`, `"INDEL"`, and `"INDELBLOCK"`. Note: If `Model = "N"` and `Dist.type` is set to `"prop"` or `"percentage"`, an error will be returned.

**`GAMMA`**: Same as the `gamma` parameter in `ape::dist.dna()`.

**`PW.deletion`**: Specifies whether pairwise deletion should be considered (`TRUE`, default) or not (`FALSE`).

**`by.group`**: A vector for grouping sequences into different categories (e.g., populations, species, lineages). The vector should follow the same order as the sequence alignment.

**`within.group`**: Creates a `data.frame` with the within-group standard error and the number of samples ("N") by group. A vector with discrete groupings should be provided in `by.group`. Groups with `N = 1` will return `NA`s.


Example:
```{r}
# Defining groupings
####
#create data.frame with specimen and grouping (e.g., species) IDs
Species.info <- data.frame(Specimen = paste0("Ind",LETTERS[1:9]),
                           Species = paste0("Sp",c(rep(1,3),rep(2,3),rep(3,3))),
                           Seq.names = names(test.DNA))

Species.info  # Check that they match.
# The funtion uses the order of the `by.group` vector to group the sequences

DNA.GroupDist.SE(test.DNA, Out.Format = "long", Model = "JC69", by.group = Species.info$Species)

# Reorganize the output 
CompSE<- DNA.GroupDist.SE(test.DNA, Out.Format = "long", Model = "JC69", by.group = Species.info$Species)
CompSE <- data.frame(Pairs = paste(CompSE$x1, CompSE$x2, sep = "_vs_"),
                         SE = CompSE$Standard_Error)
CompSE

#wide format
WDistSE <- DNA.GroupDist.SE(test.DNA, Out.Format = "wide", Model = "JC69", by.group = Species.info$Species)
WDistSE
```

## Citation
Please cite this package as follows: 

Carrera-Martínez, R (2025) SLDivMaxR: A R Package for Calculating Single Allele Divergence Matrices Between and Within Taxa. Version 0.1.0 URL:<<https://github.com/rcarrmart/SLDivMaxR/>>. 

To cite `vegan`: 

Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P, O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, Barbour M, Bedward M, Bolker B, Borcard D, Borman T, Carvalho G, Chirico M, De Caceres M, Durand S, Evangelista H, FitzJohn R, Friendly M, Furneaux B, Hannigan G, Hill M, Lahti L, McGlinn D, Ouellette M, Ribeiro Cunha E, Smith T, Stier A, Ter Braak C, Weedon J (2024). vegan: Community Ecology Package. R package version 2.7-0, https: <<https://vegandevs.github.io/vegan/>>.

To cite `ape`: 

Paradis E, Schliep K (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. _Bioinformatics_, 35, 526-528. doi:10.1093/bioinformatics/bty633. 
