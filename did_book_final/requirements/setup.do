* ============================================================
* did_book - setup.do
* Installs the Stata packages used in the book
* ============================================================

* ------------------------------------------------------------
* 1. SHARED DEPENDENCIES
* ------------------------------------------------------------

* reghdfe (also needed by eventstudyinteract and did_imputation)
ssc install ftools, replace
ssc install require, replace
ssc install reghdfe, replace

* Fast data tools
ssc install gtools, replace

* panelview
net install grc1leg, from("http://www.stata.com/users/vwiggins") replace
net install gr0075, from("http://www.stata-journal.com/software/sj18-4") replace
ssc install labutil, replace
ssc install sencode, replace
ssc install palettes, replace
ssc install colrspace, replace
ssc install panelview, replace

* ------------------------------------------------------------
* 2. PACKAGES BY CHAPTER
* (a package used in several chapters is listed once)
* ------------------------------------------------------------

* Chapter 3
ssc install pretrends, replace

* Chapter 4
ssc install did_multiplegt_dyn, replace
net install fect, from("https://raw.githubusercontent.com/xuyiqing/fect_stata/master/") replace
ssc install honestdid, replace
ssc install sdid_event, replace

* Chapter 5
ssc install twowayfeweights, replace

* Chapter 6
ssc install avar, replace
ssc install eventstudyinteract, replace
ssc install drdid, replace
ssc install csdid, replace
ssc install did_imputation, replace

* Chapter 7
net install nprobust, from("https://raw.githubusercontent.com/nppackages/nprobust/master/stata") replace
ssc install did_had, replace

* Chapter 8
ssc install did_multiplegt_stat, replace

display as result "did_book: all Stata packages installed."
