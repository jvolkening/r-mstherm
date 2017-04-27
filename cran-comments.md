## Resubmission
This is a resubmission. In this version I have:

* Made all examples in documentation executable

* Added a reference to a relevant paper as background, to be replaced with a
  citation to work from our own lab on next release.

* Moved an experimental block of code into a development branch


## Test environments
* local Debian 8 install, R 3.3.3
* Ubuntu 12.04.5 LTS, R 3.3.3, R-devel (on Travis-CI)
* Windows (win-builder), R-devel and R-release
* Windows 10, R 3.4.1 i386 and x64


## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Possibly mis-spelled words in DESCRIPTION:
    Savitski (12:61)
    al (12:73)
    et (12:70)
    proteome (9:5)
    spectrometry (8:61)

  These words are spelled correctly


## Downstream dependencies
There are currently no downstream dependencies for this package
