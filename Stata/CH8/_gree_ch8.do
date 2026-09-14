
cd "/Users/anzony.quisperojas/Documents/GitHub/did_book"
use "_data/gentzkowetal_didtextbook.dta", replace


/* Using gentzkowetal_didtextbook, regress changedailies on
lag_numdailies, clustering standard errors at the county level. Interpret the results */

reg changedailies lag_numdailies, vce(cluster cnty90 )


/* Using gentzkowetal_didtextbook, regress changedailies on lag_ishare_urb, counties' lagged
urbanization rate, clustering standard errors at the county level. Interpret the results. */

reg changedailies lag_ishare_urb, vce(cluster cnty90 )

/* Using gentzkowetal_didtextbook, regress turnout
on the number of newspapers and its lag and on county and year FEs, clustering standard errors
at the county level. Interpret the results. 317 */


areg prestout i.year numdailies lag_numdailies, absorb(cnty90) cluster(cnty90)

/* Use the twowayfeweights Stata package and the other_treatments option
to decompose bβ dl0 and bβ dl 1 , and interpret the results. */




did_multiplegt_dyn prestout cnty90 year numdailies, effects(4) placebo(4) effects_equal(all)


did_multiplegt_dyn prestout cnty90 year numdailies, effects(1) design(0.8,console) graph_off



did_multiplegt_dyn prestout cnty90 year numdailies, effects(4) placebo(4) normalized effects_equal(all)


did_multiplegt_dyn prestout cnty90 year numdailies if year<=first_change|same_treat_after_first_change==1, effects(2) effects_equal(all) same_switchers graph_off

did_multiplegt_dyn prestout cnty90 year numdailies, effects(4) placebo(4) normalized normalized_weights effects_equal(all)

