capture use "C:\Users\134476.SCIENCESPO\C DE CHAISEMARTIN Dropbox\clément de chaisemartin\A Mini course DID\Applications\Data sets\Wolfers 2006\wolfers2006_didtextbook.dta", clear

if _rc==601{
use "C:\Users\134476\C DE CHAISEMARTIN Dropbox\clément de chaisemartin\A Mini course DID\Applications\Data sets\Wolfers 2006\wolfers2006_didtextbook.dta", clear
}

*** 1) Static TWFE regression
 
reg div_rate udl i.state i.year [w=stpop], vce(cluster state)

*** 2) Decomposing the static TWFE regression
 
*ssc install twowayfeweights
*help twowayfeweights
twowayfeweights div_rate state year udl, type(feTR) test_random_weights(exposurelength) weight(stpop)

*** 3) Testing randomized treatment timing

reg div_rate i.early_late_never if cohort!=1956&year<=1968 [w=stpop], vce(cluster state)

*** 4) Event-study TWFE regression
 
reg div_rate rel_time* i.state i.year [w=stpop], vce(cluster state)
test rel_timeminus1 rel_timeminus2 rel_timeminus3 rel_timeminus4 rel_timeminus5 rel_timeminus6 rel_timeminus7 rel_timeminus8 rel_timeminus9

*** 5) Decomposing the first estimated effect in the event-study TWFE regression

twowayfeweights div_rate state year rel_time1, type(feTR) test_random_weights(year) weight(stpop) other_treatments(rel_time2-rel_time16) controls(rel_timeminus1-rel_timeminus9)

*** 6) Sun and Abraham event-study estimators

*ssc install eventstudyinteract
*help eventstudyinteract
replace cohort=. if cohort==0
eventstudyinteract div_rate rel_time* [aweight=stpop], ///
	absorb(i.state i.year) cohort(cohort) control_cohort(controlgroup) ///
	vce(cluster state)
	
*** 7) Estimators of Callaway and Sant'Anna	

*ssc install csdid
*help csdid
replace cohort=0 if cohort==.
csdid div_rate [weight=stpop], ivar(state) time(year) gvar(cohort) notyet agg(event)

*** 8) Estimators of de Chaisemartin and D'Haultfoeuille

*ssc install did_multiplegt_dyn
*help did_imputation_dyn
did_multiplegt_dyn div_rate state year udl, effects(16) placebo(9) weight(stpop)

*** 9) Estimators of Borusyak et. al.

*ssc install did_imputation
*help did_imputation
replace cohort=. if cohort==0
did_imputation div_rate state year cohort [aweight=stpop], horizons(0/15) autosample minn(0) pre(9)


