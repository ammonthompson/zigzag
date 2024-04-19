# <strong><em>Z<sup>i</sup><sub>g</sub> Z<sup>a</sup><sub>g</sub></em></strong> v1.0.0

# Quick Installation Guide and Simple Analysis

The `zigzag` R package is designed to compute the posterior probability that genes are actively expressed, given a set of RNA-seq relative expression estimates. Requires at least 2 replicates (expression estimates from RNA-seq libraries).

## Installation

To install the `zigzag` package, you can use the following commands in your R environment:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("ammonthompson/zigzag")
```
or you can clone the repository and install it locally:

``` r
install.packages("path/to/zigzag", type = "source", repos = NULL)
```

## Simple Analysis Example

You can perform a basic analysis to get a feel for `zigzag`'s functionality. Here's a simple example to get you started:

``` r
library(zigzag)

# Load your RNA-seq expression data with at least 2 replicate libraries.
expression_data <- read.csv('path/to/your/data.csv', header = T, row.names = 1)
gene_length_data <- read.csv('/path/to/your/gene_length_data.csv', header = T, row.names = 1)

# Load data into a zigzag object
# Compute the posterior probability of active expression by
# first running burnin until the chain is stationary,
# then run mcmc to sample from the posterior distribution.
my_zigzag <- zigzag$new(expression_data, 
                        gene_lengths_data, 
                        output_dir = "my_output")
my_zigzag$burnin()
my_zigzag$mcmc()

```
View the results located in "my_output/" and evaluate quality of the MCMC samples and priors in "my_output/*mcmc_output/mcmc_report". 

You will likely not want to use all default settings. Use ? to view documentation on functions. For example

``` r
?zigzag
?zigzag::burnin 
?zigzag::mcmc
```

## Documentation and Tutorial

For a more comprehensive guide and detailed tutorials on how to use `zigzag`, please visit our [tutorial page](https://ammonthompson.github.io/zigzag_user_guide/).

## Citation

If you use `zigzag` in your research, please cite our work using the following reference:

[Thompson et al., 2020, PNAS](https://doi.org/10.1073/pnas.1919748117)

