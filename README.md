MSTherm :: Analyze MS/MS protein melting data
=======

[![Build Status](https://travis-ci.org/jvolkening/r-mstherm.svg?branch=master)](https://travis-ci.org/jvolkening/r-mstherm)
[![Coverage Status](https://coveralls.io/repos/github/jvolkening/r-mstherm/badge.svg?branch=master)](https://coveralls.io/github/jvolkening/r-mstherm?branch=master)


**MSTherm** is an R package to aid in modeling and analyzing
mass-spectrometry-based protein melting data. Quantitative data is imported
and normalized and thermal behavior is modeled at the protein level. Functions
exist for visualization, exploratory analysis, and analsysis of melting
temperature shifts.

Installation
------------

In Linux, do:

    git clone https://github.com/jvolkening/r-MSTherm.git
    cd r-MSTherm
    R CMD INSTALL

You may need to install some prerequisites first depending on your environment
-- the installer will notify you of this.
