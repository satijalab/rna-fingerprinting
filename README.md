# RNA Fingerprinting 

RNA fingerprinting is a statistical framework for scRNA-seq data that maps transcriptional responses observed in new experiments to perturbation dictionaries. This takes place in two main steps. First, fingerprints (i.e. denoised representations of perturbation effects) are estimated for each reference perturbation in the dictionary using a multi-condition latent factor model. Second, queries are mapped to these fingerprints using a Bayesian regression framework. Queries can be mapped either at the level of individual cells or as samples (i.e. groups of cells with the same label). Further details on the statistical methodology and biological applications are in our [preprint](https://www.biorxiv.org/content/10.1101/2025.09.19.676866v1). 

## Installation

This package can be installed as:

```
# install.packages("devtools")
devtools::install_github("satijalab/rna-fingerprinting")
```

If you would like to use any of our precomputed fingerprints (from Genome-Wide Perturb-seq or the Immune Dictionary), these can be accessed by installing the suggested dependency [`RNAFingerprintingData`](https://github.com/satijalab/rna-fingerprinting-data) as follows:

```
# install.packages("remotes")
remotes::install_url("https://github.com/satijalab/rna-fingerprinting-data/raw/refs/heads/main/RNAFingerprintingData_0.1.0.tar.gz")
```

## Vignettes

Knitted vignettes (HTML files) showing example usage can be found in the `docs` folder of this repository. See [example usage](https://satijalab.github.io/rna-fingerprinting/jost-example-usage.html) (recommended start) and [example usage with precomputed fingerprints](https://satijalab.github.io/rna-fingerprinting/jost-gwps-example.html).
