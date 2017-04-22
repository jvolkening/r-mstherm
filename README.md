MSTherm :: Analyze MS/MS protein melting data
=======

[![Build Status](https://travis-ci.org/jvolkening/r-mstherm.svg?branch=master)](https://travis-ci.org/jvolkening/r-mstherm)
[![Coverage Status](https://coveralls.io/repos/github/jvolkening/r-mstherm/badge.svg?branch=master)](https://coveralls.io/github/jvolkening/r-mstherm?branch=master)


**MSTherm** is an R package to aid in modeling and analyzing mass
spectrometry-based protein melting data. Quantitative data is imported and
normalized and thermal behavior is modeled at the protein level. Methods exist
for normalization, visualization, exploratory analysis, and data export.

Installation
------------

In Linux, do:

    git clone https://github.com/jvolkening/r-MSTherm.git
    cd r-MSTherm
    R CMD INSTALL

You may need to install some prerequisites first depending on your environment
-- the installer will notify you of this.

Usage
-----

Once you have produced spectral quantification tables and created the
necessary metadata files describing the experimental setup (as described in
the vignette), a modeling session can be a simple as:

```
library(mstherm)
expt <- MSThermExperiment("control.tsv", "annotations.tsv")
expt <- normalize_to_std(expt, "BOVINE_SERUM_ALBUMIN")
res  <- model_experiment(expt)
tbl  <- as.data.frame(res)
pdf("plots.pdf", 5, 5, pointsize=10)
plot(res)
dev.off()
```

(Of couse, you will likely want to experiment with the various options
available to the modeling and plotting methods).

Further details of the required input files and available methods are
available in the package vignette.
