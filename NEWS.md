# mstherm 0.4.7

* Added additional method parameters related to filtering of protein results

# mstherm 0.4.6

* Fixed conversion of MSThermResultSet to data frame

# mstherm 0.4.5

* Added export to SQLite database
* Added RMSD as measure of fit
* Various modifications to satisfy R CMD check
* Renaming of some methods for the sake of clarity
* Package name to lowercase

# MSTherm 0.4.3

* Additional normalization method (by Tm distribution)
* Fixed warnings due to missing 'score' column (should be optional)
* Fixed bug due to spectra with all missing quant values (from some quant
  software)

# MSTherm 0.4.1

* The plot.MSThermResult method has two additional arguments (CI.points and
  CI.Tm) which control whether to plot the confidence intervals if present.
  Previously these were always plotted if present. The default of both is still
  TRUE.
