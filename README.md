# RNA Fingerprinting 

RNA fingerprinting is a statistical framework for scRNA-seq data that maps transcriptional responses observed in new experiments to perturbation dictionaries. This takes place in two main steps. First, fingerprints (i.e. denoised representations of perturbation effects) are estimated for each reference perturbation in the dictionary using a multi-condition latent factor model. Second, queries are mapped to these fingerprints using a Bayesian regression framework. Queries can be mapped either at the level of individual cells or as samples (i.e. groups of cells with the same label). Further details on the statistical methodology and biological applications are in our preprint. 

## Installation

If you would like to use any of our precomputed fingerprints (in addition or instead of estimating fingerprints from your own perturbation dictionaries), these can be accessed by installing the suggested dependency `RNAFingerprintingData` as follows:

```
# install.packages("remotes"")
remotes::install_url("https://github.com/satijalab/rna-fingerprinting-data/raw/refs/heads/main/RNAFingerprintingData_0.1.0.tar.gz")
```

Then you can install this package as:

```
# install.packages("devtools")
devtools::install_github("satijalab/rna-fingerprinting")
```

## Vignettes

Knitted vignettes (HTML files) showing example usage can be found in the `docs` folder of this repository. See [example usage](https://github.com/satijalab/rna-fingerprinting/blob/main/docs/jost-example-usage.html) (recommended start) and [example usage with precomputed fingerprints](https://github.com/satijalab/rna-fingerprinting/blob/main/docs/jost-gwps-example.html).
