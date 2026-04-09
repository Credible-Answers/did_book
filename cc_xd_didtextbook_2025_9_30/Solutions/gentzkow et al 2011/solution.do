use "C:\Users\134476.SCIENCESPO\C DE CHAISEMARTIN Dropbox\clément de chaisemartin\A Mini course DID\Applications\Data sets\Gentzkow et al 2011\gentzkowetal_didtextbook.dta", clear
use "C:\Users\134476\C DE CHAISEMARTIN Dropbox\clément de chaisemartin\A Mini course DID\Applications\Data sets\Gentzkow et al 2011\gentzkowetal_didtextbook.dta", clear

************************** Chapter 5

*** 1. TWFE Regression

*Estimation
areg prestout i.year numdailies, absorb(cnty90) cluster(cnty90)

*Decomposition
twowayfeweights prestout cnty90 year numdailies, type(feTR)

*** 2. TWFE with state-specific trends

qui tab styr, gen(styr)

*Estimation
qui areg prestout i.year i.styr numdailies, absorb(cnty90) cluster(cnty90)
di _b[numdailies], _se[numdailies]

*Decomposition
twowayfeweights prestout cnty90 year numdailies, type(feTR) controls(styr1-styr683)

*** 3. FD Regression with state-specific trends

*Estimation
areg changeprestout changedailies, absorb(styr) cluster(cnty90)

*Decomposition
twowayfeweights changeprestout cnty90 year changedailies numdailies, type(fdTR) controls(styr1-styr683) 

*Assessing if weights correlated with year variable
twowayfeweights changeprestout cnty90 year changedailies numdailies, type(fdTR) controls(styr1-styr683) test_random_weights(year)

************************** Chapter 8

*** 1. Testing whether change in daily newspapers as good as random

reg changedailies lag_numdailies, cluster(cnty90)
reg changedailies lag_ishare_urb, cluster(cnty90)

*** 2. Distributed lage TWFE Regression

*Estimation
areg prestout i.year numdailies lag_numdailies, absorb(cnty90) cluster(cnty90)

*Decomposition
twowayfeweights prestout cnty90 year numdailies, other_treatments(lag_numdailies) type(feTR)
twowayfeweights prestout cnty90 year lag_numdailies, other_treatments(numdailies) type(feTR)

*** 3. Non-normalized event-study effects

did_multiplegt_dyn prestout cnty90 year numdailies, effects(4) placebo(4) effects_equal(all)

// With graph
*did_multiplegt_dyn prestout cnty90 year numdailies, effects(4) placebo(4) graphoptions(ylabel(-0.04(0.01)0.05) xlabel(-4(1)4) yscale(range(-0.04 0.05)) legend(off) xtitle(Relative time to change in newspapers) title(Non-normalized DID estimates) ytitle(Effect))
*graph export "C:\Users\134476\C DE CHAISEMARTIN Dropbox\clément de chaisemartin\A Mini course DID\Textbook\graphnewspapers_dCDH.pdf", replace

*** 4. Analyzing the paths whose effect is averaged in the non-normalized event-study effects

did_multiplegt_dyn prestout cnty90 year numdailies, effects(1) design(0.8,console) graph_off

did_multiplegt_dyn prestout cnty90 year numdailies, effects(2) design(0.8,console) graph_off

did_multiplegt_dyn prestout cnty90 year numdailies, effects(4) design(0.8,console) graph_off

*** 5. Normalized event-study effects

did_multiplegt_dyn prestout cnty90 year numdailies, effects(4) placebo(4) normalized normalized_weights effects_equal(all)

// With graph
*did_multiplegt_dyn prestout cnty90 year numdailies, effects(4) placebo(4) normalized normalized_weights effects_equal(all) graphoptions(ylabel(-0.01(0.01)0.02) xlabel(-4(1)4) yscale(range(-0.01 0.02)) legend(off) xtitle(Relative time to change in newspapers) title(Normalized DID estimates) ytitle(Effect))
*graph export "C:\Users\134476\C DE CHAISEMARTIN Dropbox\clément de chaisemartin\A Mini course DID\Textbook\graphnewspapers_dCDH_normalized.pdf", replace

*** 6. Testing if the lagged number of newspapers affects turnout

did_multiplegt_dyn prestout cnty90 year numdailies if year<=first_change|same_treat_after_first_change==1, effects(2) effects_equal(all) same_switchers graph_off

*** 7. Estimators assuming away effects of lagged treatments on the outcome

egen election_number=group(year)
did_multiplegt_stat prestout cnty90 election_number numdailies, placebo(1) exact_match
tab lag_numdailies if year==first_change
tab lag_numdailies if changedailies!=0&changedailies!=.&year!=1868


