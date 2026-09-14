
/****************************************************************************************
  DID Textbook – Applications
  ***************************************************************/

version 17.0
clear all

* Change this path to a folder on your computer
cd "/Users/karlavega/Documents/GitHub/did_book"

* Download textbook datasets and do-files
ssc describe cc_xd_didtextbook
net get cc_xd_didtextbook


use "/Users/karlavega/Documents/GitHub/did_book/cc_xd_didtextbook_2025_9_30/Data sets/Moser and Voena 2012/moser_voena_didtextbook.dta", clear


******************* Chapter 3

*** 1) Static TWFE regression
 
* Two-way fixed effects regression with subclass and year FE.
* Standard errors clustered at the subclass level.

xtreg patents twea i.year, fe i(subclass) cluster(subclass)
reghdfe patents twea, absorb(subclass year) cluster(subclass)


*** 2) Equivalence between static TWFE regression and DID
 
* Canonical DiD regression.
* The coefficient on twea equals β_fe.

reg patents treatmentgroup post twea, cluster(subclass)


*** 3) Testing randomized treatment

* Regress patents on the treatment group indicator
* using only pre-treatment years (before 1919).

reg patents treatmentgroup if year<=1918, cluster(subclass)



*** 4) Event-study TWFE regression
 
* Event-study with year FE, treatment group FE,
* leads (reltimeminus*) and lags (reltimeplus*).
* Standard errors clustered at the subclass level.
reg patents i.year treatmentgroup reltimeminus* reltimeplus*, cluster(subclass)


* Test of pre-trends (NA / Parallel Trends)
* Joint test: all pre-treatment coefficients equal zero.

test reltimeminus1 reltimeminus2 reltimeminus3 reltimeminus4 reltimeminus5 reltimeminus6 reltimeminus7 reltimeminus8 reltimeminus9 reltimeminus10 reltimeminus11 reltimeminus12 reltimeminus13 reltimeminus14 reltimeminus15 reltimeminus16 reltimeminus17 reltimeminus18

* Generating corresponding ES plot
{
/* We produce the E-S graph by creating a matrix (res) gathering the time to the 
event, the point estimates and the CI */

reg patents reltimeminus* reltimeplus* i.year treatmentgroup, cluster(subclass)

matrix temp=r(table)'
matrix res=J(40,4,0)
matrix res[19,1]=0
forvalues x = 1/18 {
matrix res[19-`x',1]=-`x'
matrix res[19-`x',2]=temp[`x',1]
matrix res[19-`x',3]=temp[`x',5]
matrix res[19-`x',4]=temp[`x',6]
}
forvalues x = 1/21 {
matrix res[`x'+19,1]=`x'
matrix res[`x'+19,2]=temp[`x'+18,1]
matrix res[`x'+19,3]=temp[`x'+18,5]
matrix res[`x'+19,4]=temp[`x'+18,6]
}

// Store sub matrix for the figure with the Borusyak et al estimator
matrix res_post=res["r19".."r40","c1".."c4"]

preserve
drop _all
svmat res
twoway (scatter res2 res1, msize(medlarge) msymbol(o) mcolor(navy) legend(off)) ///
	(line res2 res1, lcolor(navy)) (rcap res4 res3 res1, lcolor(maroon)), ///
	 title("TWFE Event-study estimates") xtitle("Relative time to year before TWEA") ///
	 ytitle("Effect") xlabel(-18(3)21) yscale(range(-0.25 1)) ylabel(-0.25(.25)1)
graph export "/Users/karlavega/Documents/GitHub/did_book/graphES_moser1.pdf", replace
restore
}

* Verifying "by hand" that event-study coefficients are simple DIDs

sum patents if year==1919&treatmentgroup==1
scalar m1=r(mean)
sum patents if year==1918&treatmentgroup==1
scalar m2=r(mean)
sum patents if year==1919&treatmentgroup==0
scalar m3=r(mean)
sum patents if year==1918&treatmentgroup==0
scalar m4=r(mean)
di m1-m2-(m3-m4)


*** 5) Event-study TWFE regression, without pre-trends estimates

* Event-study specification excluding pre-treatment leads.
* Used to verify equation (3.7) for ℓ = 1 by computing the DID on the RHS.
reg patents i.yearpost treatmentgroup reltimeplus*, cluster(subclass)

* Generating corresponding ES plot
{
reg patents reltimeplus* i.yearpost treatmentgroup, cluster(subclass)

matrix temp=r(table)'
matrix res=J(22,4,0)
matrix res[1,1]=0
forvalues x = 1/21 {
matrix res[`x'+1,1]=`x'
matrix res[`x'+1,2]=temp[`x',1]
matrix res[`x'+1,3]=temp[`x',5]
matrix res[`x'+1,4]=temp[`x',6]
}

preserve
drop _all
svmat res
svmat res_post
twoway (scatter res2 res1, msize(small) msymbol(o) mcolor(midblue) legend(order(2 /*"  &" 3*/ "Without Pre-Periods" 4 /*"  &" 5*/ "With Pre-Periods") pos(6) col(2))) ///
	(line res2 res1, lcolor(midblue)) (rcap res4 res3 res1, lcolor(midblue)) ///
	(line res_post2 res1, lcolor(red)) (rcap res_post4 res_post3 res1, lcolor(red)) ///
	(scatter res_post2 res1, msize(small) msymbol(o) mcolor(red)), ///
	 title("TWFE Event-study estimates") xtitle("Relative time to year before TWEA") ///
	 ytitle("Effect") xlabel(0(3)21) yscale(range(-0.25 1)) ylabel(-0.25(.25)1)
graph export "C:\Users\134476\C DE CHAISEMARTIN Dropbox\clément de chaisemartin\A Mini course DID\Textbook\graphES_moser2.pdf", replace
restore
}

*Numerically equivalent to Borusyak et al, Gardner, Liu et al

gen cohort=1919 if treatmentgroup==1

did_imputation patents subclass year cohort, horizons(0/20) autosample minn(0)


*Numerically equivalent to having all year FEs in the regression

reg patents i.year treatmentgroup reltimeplus*, cluster(subclass)



*** 6) Linear pre-trends we could fail to detect

* This exercise evaluates whether the estimated effects of compulsory licensing
* could be driven by differential linear trends between treated and control subclasses
* that we do not have enough power to detect.

*Estimation based on our first 6 pre-trends estimates only, because otherwise command takes too long to run.
local github https://raw.githubusercontent.com
net install pretrends, from(`github'/mcaceresb/stata-pretrends/main) replace
reghdfe patents reltimeminus* reltimeplus*, absorb(treatmentgroup year) cluster(subclass)
pretrends power 0.5, numpre(6)

* Generating corresponding ES plot
{
local slope=r(slope)	
	
preserve
	
keep if year>=1912&year<=1939
drop reltimeminus7-reltimeminus18
	
reg patents reltimeminus* reltimeplus* i.year treatmentgroup, cluster(subclass)

restore 

matrix temp=r(table)'
matrix res=J(28,4,0)
matrix res[7,1]=0
forvalues x = 1/6 {
matrix res[7-`x',1]=-`x'
matrix res[7-`x',2]=temp[`x',1]
matrix res[7-`x',3]=temp[`x',5]
matrix res[7-`x',4]=temp[`x',6]
}
forvalues x = 1/21 {
matrix res[`x'+7,1]=`x'
matrix res[`x'+7,2]=temp[`x'+6,1]
matrix res[`x'+7,3]=temp[`x'+6,5]
matrix res[`x'+7,4]=temp[`x'+6,6]
}

preserve
drop _all
svmat res
twoway (scatter res2 res1, msize(medlarge) msymbol(o) mcolor(navy) legend(off)) ///
	(line res2 res1, lcolor(navy)) (rcap res4 res3 res1, lcolor(maroon)) (function y=x*`slope', range(-6 21) lcolor(gray) lpattern(dash)), ///
	 title("TWFE Event-study estimates") xtitle("Relative time to year before TWEA") ///
	 ytitle("Effect") xlabel(-6(3)21) yscale(range(-0.25 1)) ylabel(-0.25(.25)1)
graph export "C:\Users\fe-kn\C DE CHAISEMARTIN Dropbox\RAs De Chaisemartin\Mini course DID\Applications\Solutions\Moser and Voena 2012\graphs\graphES_moser3.pdf", replace
restore
}




*** 7) The variance of the effect of having been exposed to treatment for 14 years

* Compare the variance of outcomes between treated and control subclasses
* 14 years after treatment (year = 1932).
sdtest diffpatentswrt1918 if year==1932, by(treatmentgroup)

* Difference in standard deviations between treated and control groups
di r(sd_2)-r(sd_1)
scalar sd_effects=r(sd_2)-r(sd_1)

* Point estimate of the treatment effect
reg diffpatentswrt1918 treatmentgroup if year==1932

* Confidence interval using the estimated variance difference
di _b[treatmentgroup]-1.96*sd_effects,_b[treatmentgroup]+1.96*sd_effects






*** 8) Placebo test of the assumptions underlying the estimation of the variance of treatment effects

*** Treatment effect heterogeneity (14 years)
* Variance test 14 years after treatment.
* H0: v14 = 0 (no heterogeneity).
sdtest diffpatentswrt1918 if year==1904, by(treatmentgroup)

*** Placebo tests (pre-treatment)
forvalue i=1900/1939{
sdtest diffpatentswrt1918 if year==`i', by(treatmentgroup)
}



******************* Chapter 4

*** 1) Estimators with controls

*Controlling for patents in 1900, TWFE (in blue on figure, legend="TWFE controlling for baseline patents")
reghdfe patents reltimeminus* reltimeplus*, absorb(year#patents1900 treatmentgroup) cluster(subclass)
test reltimeminus1 reltimeminus2 reltimeminus3 reltimeminus4 reltimeminus5 reltimeminus6 reltimeminus7 reltimeminus8 reltimeminus9 reltimeminus10 reltimeminus11 reltimeminus12 reltimeminus13 reltimeminus14 reltimeminus15 reltimeminus16 reltimeminus17 reltimeminus18

*** Event-study plot with controls
*Producing the graph
reghdfe patents reltimeminus* reltimeplus*, absorb(year#patents1900 treatmentgroup) cluster(subclass)
matrix temp=r(table)'
matrix res=J(40,4,0)
matrix res[19,1]=0
forvalues x = 1/18 {
matrix res[19-`x',1]=-`x'
matrix res[19-`x',2]=temp[`x',1]
matrix res[19-`x',3]=temp[`x',5]
matrix res[19-`x',4]=temp[`x',6]
}
forvalues x = 1/21 {
matrix res[`x'+19,1]=`x'
matrix res[`x'+19,2]=temp[`x'+18,1]
matrix res[`x'+19,3]=temp[`x'+18,5]
matrix res[`x'+19,4]=temp[`x'+18,6]
}

preserve
drop _all
svmat res
twoway (scatter res2 res1, msize(medlarge) msymbol(o) mcolor(navy) legend(off)) ///
	(line res2 res1, lcolor(navy)) (rcap res4 res3 res1, lcolor(maroon)), ///
	 title("TWFE Event-study estimates") xtitle("Relative time to year before TWEA") ///
	 ytitle("Effect") xlabel(-18(3)21) yscale(range(-0.25 1)) ylabel(-0.25(.25)1)
	 
graph export "/Users/karlavega/Documents/GitHub/did_book/graphES_moser_controls.pdf", replace
restore

*Testing that treatment group indicator and covariate are correlated

reg patents1900 treatmentgroup if year==1900

*Controlling for patents in 1900, DID (not on figure)
did_multiplegt_dyn patents subclass year twea, effects(21) placebo(18) trends_nonparam(patents1900)



*** 2) Interactive fixed effects

net install fect, from(https://raw.githubusercontent.com/xuyiqing/fect_stata/master/) replace
ssc install _gwtmean, replace

// Optimal number of factors = 2
fect patents, treat(twea) unit(subclass) time(year) method("ife") r(4) tol(1e-4) cv

// If cross-validation on the treated only, optimal number of factors =1
fect patents, treat(twea) unit(subclass) time(year) method("ife") r(4) cv tol(1e-4) cvtreat  

timer clear
timer on 1
set seed 1 
fect patents, treat(twea) unit(subclass) time(year) method("ife") r(2) tol(1e-4) se
timer off 1
timer list
matrix list e(ATT)

// Run time on Dell desktop computer, processor 11th Gen Intel(R) Core(TM) i7-11700T @ 1.40GHz 1.39 GHz, Stata MP 18: 554 seconds. 

matrix res_ife=J(21,4,0)
forvalues x = 1/21 {
matrix res_ife[`x',1]=`x'
matrix res_ife[`x',2]=e(ATTs)[`x'+19,3]
matrix res_ife[`x',3]=e(ATTs)[`x'+19,6]
matrix res_ife[`x',4]=e(ATTs)[`x'+19,7]
}

matrix res_post=res_post[2..22,1..4]

preserve
drop _all
svmat res_ife
svmat res_post
twoway (scatter res_ife2 res_ife1, msize(small) msymbol(o) mcolor(midblue) legend(order(2 /*"  &" 3*/ "IFE" 4 /*"  &" 5*/ "TWFE") pos(6) col(2))) ///
	(line res_ife2 res_ife1, lcolor(midblue)) (rcap res_ife4 res_ife3 res_ife1, lcolor(midblue)) ///
	(line res_post2 res_post1, lcolor(red)) (rcap res_post4 res_post3 res_post1, lcolor(red)) ///
	(scatter res_post2 res_post1, msize(small) msymbol(o) mcolor(red)), ///
	 title("Event-study estimates") xtitle("Years after TWEA") ///
	 ytitle("Effect") xscale(range(1 21)) xlabel(1(5)21) yscale(range(-0.5 1.5)) ylabel(-0.5(.25)1.5)
graph export "/Users/karlavega/Documents/GitHub/did_book/graphES_moser_ife_did.pdf", replace
restore


*** 3) Synthetic control

* Estimate dynamic treatment effects using synthetic control.
* Standard errors obtained via bootstrap.

ssc install sdid_event, replace

timer clear
timer on 1
set seed 1 
sdid_event patents subclass year twea, method("sc") brep(200)
timer off 1
timer list

matrix res_sc=J(21,4,0)
forvalues x = 1/21 {
matrix res_sc[`x',1]=`x'
matrix res_sc[`x',2]=e(H)[`x'+1,1]
matrix res_sc[`x',3]=e(H)[`x'+1,3]
matrix res_sc[`x',4]=e(H)[`x'+1,4]
}

preserve
drop _all
svmat res_sc
svmat res_post
twoway (scatter res_sc2 res_sc1, msize(small) msymbol(o) mcolor(midblue) legend(order(2 /*"  &" 3*/ "SC" 4 /*"  &" 5*/ "TWFE") pos(6) col(2))) ///
	(line res_sc2 res_sc1, lcolor(midblue)) (rcap res_sc4 res_sc3 res_sc1, lcolor(midblue)) ///
	(line res_post2 res_post1, lcolor(red)) (rcap res_post4 res_post3 res_post1, lcolor(red)) ///
	(scatter res_post2 res_post1, msize(small) msymbol(o) mcolor(red)), ///
	 title("Event-study estimates") xtitle("Years after TWEA") ///
	 ytitle("Effect") xscale(range(1 21)) xlabel(1(5)21) yscale(range(-2 2)) ylabel(-2(.5)2)
graph export "/Users/karlavega/Documents/GitHub/did_book/graphES_moser_sc_did.pdf", replace
restore

*Demeaned

bys subclass: egen pre_mean_temp=mean(patents) if year<=1918
bys subclass: egen pre_mean=mean(pre_mean_temp)
gen patents_demeaned=patents-pre_mean
timer clear
timer on 1
set seed 1 
sdid_event patents_demeaned subclass year twea, method("sc") brep(200)
timer off 1
timer list

matrix res_sc=J(21,4,0)
forvalues x = 1/21 {
matrix res_sc[`x',1]=`x'
matrix res_sc[`x',2]=e(H)[`x'+1,1]
matrix res_sc[`x',3]=e(H)[`x'+1,3]
matrix res_sc[`x',4]=e(H)[`x'+1,4]
}

preserve
drop _all
svmat res_sc
svmat res_post
twoway (scatter res_sc2 res_sc1, msize(small) msymbol(o) mcolor(midblue) legend(order(2 /*"  &" 3*/ "Demeaned SC" 4 /*"  &" 5*/ "TWFE") pos(6) col(2))) ///
	(line res_sc2 res_sc1, lcolor(midblue)) (rcap res_sc4 res_sc3 res_sc1, lcolor(midblue)) ///
	(line res_post2 res_post1, lcolor(red)) (rcap res_post4 res_post3 res_post1, lcolor(red)) ///
	(scatter res_post2 res_post1, msize(small) msymbol(o) mcolor(red)), ///
	 title("Event-study estimates") xtitle("Years after TWEA") ///
	 ytitle("Effect") xscale(range(1 21)) xlabel(1(5)21) yscale(range(-2 2)) ylabel(-2(.5)2)
graph export "/Users/karlavega/Documents/GitHub/did_book/graphES_moser_demeaned_sc_did.pdf", replace
restore




*** 4) Synthetic did

timer clear
timer on 1
set seed 1 
sdid_event patents subclass year twea, brep(200)
timer off 1
timer list

matrix res_sd=J(21,4,0)
forvalues x = 1/21 {
matrix res_sd[`x',1]=`x'
matrix res_sd[`x',2]=e(H)[`x'+1,1]
matrix res_sd[`x',3]=e(H)[`x'+1,3]
matrix res_sd[`x',4]=e(H)[`x'+1,4]
}

preserve
drop _all
svmat res_sd
svmat res_post
twoway (scatter res_sd2 res_sd1, msize(small) msymbol(o) mcolor(midblue) legend(order(2 /*"  &" 3*/ "SD" 4 /*"  &" 5*/ "TWFE") pos(6) col(2))) ///
	(line res_sd2 res_sd1, lcolor(midblue)) (rcap res_sd4 res_sd3 res_sd1, lcolor(midblue)) ///
	(line res_post2 res_post1, lcolor(red)) (rcap res_post4 res_post3 res_post1, lcolor(red)) ///
	(scatter res_post2 res_post1, msize(small) msymbol(o) mcolor(red)), ///
	 title("Event-study estimates") xtitle("Years after TWEA") ///
	 ytitle("Effect") xscale(range(1 21)) xlabel(1(5)21) yscale(range(-2 2)) ylabel(-2(.5)2)
graph export "/Users/karlavega/Documents/GitHub/did_book/graphES_moser_sd_did.pdf", replace
restore



*** 5) Sensitivity analysis of Rambachan and Roth 

reghdfe patents reltimeminus* reltimeplus*, absorb(year treatmentgroup) cluster(subclass)
matrix l_vec = J(21,1,0)
matrix l_vec[14,1]=1
honestdid, pre(1/18) post(19/39) mvec(0.5(0.5)2) l_vec(l_vec) coefplot
honestdid, pre(1/3) post(19) mvec(0.5(0.5)2)
honestdid, pre(1/7) post(19) mvec(0.5(0.5)2)

*** Rambachan & Roth (2023) sensitivity (honestdid) using β_fe,-1 and β_fe,1

preserve
keep if year==1918|year==1904|year==1932
reghdfe patents reltimeminus14 reltimeplus14 /// 
, absorb(year treatmentgroup) cluster(subclass)
honestdid, pre(1) post(2) mvec(2(2)10) coefplot xtitle(M, size(large)) ytitle(95% Robust CI, size(large))
graph export "/Users/karlavega/Documents/GitHub/did_book/coefplot.pdf", replace
restore


use "/Users/karlavega/Documents/GitHub/did_book/cc_xd_didtextbook_2025_9_30/Data sets/gentzkow et al 2011/gentzkowetal_didtextbook.dta", clear

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


use "/Users/karlavega/Documents/GitHub/did_book/cc_xd_didtextbook_2025_9_30/Data sets/Wolfers 2006/wolfers2006_didtextbook.dta", clear

************************** Chapter 6


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




use "/Users/karlavega/Documents/GitHub/did_book/cc_xd_didtextbook_2025_9_30/Data sets/Pierce and Schott 2016/pierce_schott_didtextbook.dta", clear

************************** Chapter 7

*1) TWFE regressions
reg delta2001 ntrgap, vce(hc2, dfadjust)
reg delta2002 ntrgap, vce(hc2, dfadjust)
reg delta2004 ntrgap, vce(hc2, dfadjust)
reg delta2005 ntrgap, vce(hc2, dfadjust)

*2) Weights analysis
twowayfeweights delta2001 indusid cons ntrgap ntrgap, type(fdTR)

*3) Test that the NTR-gap treatment is as good as randomly assigned
reg ntrgap lemp1997 lemp1998 lemp1999 lemp2000, vce(hc2, dfadjust)

*4) Stute test


stute_test delta2001 ntrgap, seed(1)
stute_test delta2002 ntrgap, seed(1)
stute_test delta2004 ntrgap, seed(1)
stute_test delta2005 ntrgap, seed(1)

*Joint test
preserve
reshape long delta deltalintrend, i(indusid) j(year)
stute_test delta ntrgap indusid year if year>=2001, seed(1)
restore

*Test that quasi stayers
sort ntrgap
scalar stat_test_qs=ntrgap[1]/(ntrgap[2]-ntrgap[1])
di stat_test_qs

*5) Pre-trends test: linear
reg delta1999 ntrgap, vce(hc2, dfadjust)
reg delta1998 ntrgap, vce(hc2, dfadjust)
reg delta1997 ntrgap, vce(hc2, dfadjust)

*6) Pre-trends test, with industry-specific linear trends: linear
reg deltalintrend1998 ntrgap, vce(hc2, dfadjust)
reg deltalintrend1997 ntrgap, vce(hc2, dfadjust)

*Pre-trends test, with industry-specific linear trends: non-parametric
stute_test deltalintrend1998 ntrgap, order(0) seed(1)
stute_test deltalintrend1997 ntrgap, order(0) seed(1)
*Joint test
preserve
reshape long delta deltalintrend, i(indusid) j(year)
stute_test deltalintrend ntrgap indusid year if year<=1998, order(0) seed(1)
restore

*7) Stute test, linear trends
stute_test deltalintrend2001 ntrgap, seed(1)
stute_test deltalintrend2002 ntrgap, seed(1)
stute_test deltalintrend2004 ntrgap, seed(1)
stute_test deltalintrend2005 ntrgap, seed(1)
*Joint Stute tests
preserve
reshape long delta deltalintrend, i(indusid) j(year)
stute_test deltalintrend ntrgap indusid year if year>=2001, seed(1)
restore

*9) Estimators with linear trends
reg deltalintrend2001 ntrgap, vce(hc2, dfadjust)
reg deltalintrend2002 ntrgap, vce(hc2, dfadjust)
reg deltalintrend2004 ntrgap, vce(hc2, dfadjust)
reg deltalintrend2005 ntrgap, vce(hc2, dfadjust)




use "/Users/karlavega/Documents/GitHub/did_book/cc_xd_didtextbook_2025_9_30/Data sets/gentzkow et al 2011/gentzkowetal_didtextbook.dta", clear

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







