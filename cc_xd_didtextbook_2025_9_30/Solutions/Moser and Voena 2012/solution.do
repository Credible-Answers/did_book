use "~/Dropbox/Mini course DID/Applications/Data sets/Moser and Voena 2012/moser_voena_didtextbook.dta", clear
use "~/Dropbox/Mini course DID/Applications/Data sets/Moser and Voena 2012/moser_voena_didtextbook.dta", clear

******************* Chapter 3

*** 1) Static TWFE regression
 
xtreg patents twea i.year, fe i(subclass) cluster(subclass)
reghdfe patents twea, absorb(subclass year) cluster(subclass)

*** 2) Equivalence between static TWFE regression and DID
 
reg patents treatmentgroup post twea, cluster(subclass)

*** 3) Testing randomized treatment

reg patents treatmentgroup if year<=1918, cluster(subclass)

*** 4) Event-study TWFE regression
 
reg patents i.year treatmentgroup reltimeminus* reltimeplus*, cluster(subclass)

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
graph export "C:\Users\fe-kn\C DE CHAISEMARTIN Dropbox\RAs De Chaisemartin\Mini course DID\Applications\Solutions\Moser and Voena 2012\graphs\graphES_moser1.pdf", replace
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

sdtest diffpatentswrt1918 if year==1932, by(treatmentgroup)
di r(sd_2)-r(sd_1)
scalar sd_effects=r(sd_2)-r(sd_1)
reg diffpatentswrt1918 treatmentgroup if year==1932
di _b[treatmentgroup]-1.96*sd_effects,_b[treatmentgroup]+1.96*sd_effects

*** 8) Placebo test of the assumptions underlying the estimation of the variance of treatment effects

sdtest diffpatentswrt1918 if year==1904, by(treatmentgroup)

forvalue i=1900/1939{
sdtest diffpatentswrt1918 if year==`i', by(treatmentgroup)
}

******************* Chapter 4

*** 1) Estimators with controls

*Controlling for patents in 1900, TWFE (in blue on figure, legend="TWFE controlling for baseline patents")
reghdfe patents reltimeminus* reltimeplus*, absorb(year#patents1900 treatmentgroup) cluster(subclass)
test reltimeminus1 reltimeminus2 reltimeminus3 reltimeminus4 reltimeminus5 reltimeminus6 reltimeminus7 reltimeminus8 reltimeminus9 reltimeminus10 reltimeminus11 reltimeminus12 reltimeminus13 reltimeminus14 reltimeminus15 reltimeminus16 reltimeminus17 reltimeminus18

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
// twoway (scatter res2 res1, msize(medlarge) msymbol(o) mcolor(navy) legend(off)) ///
// 	(line res2 res1, lcolor(navy)) (rcap res4 res3 res1, lcolor(maroon)), ///
// 	 title("TWFE Event-study estimates") xtitle("Relative time to year before TWEA") ///
// 	 ytitle("Effect") xlabel(-18(3)21) yscale(range(-0.25 1)) ylabel(-0.25(.25)1) 
twoway ///
    (scatter res2 res1, msize(medlarge) msymbol(o) ///
        mcolor(gs4) legend(off)) ///
    (line res2 res1, lcolor(gs4) lpattern(solid) lwidth(medthick)) ///
    (rcap res4 res3 res1, lcolor(gs10) lpattern(solid) lwidth(medthin)), ///
    scheme(s1mono) ///
    title("TWFE Event-study estimates") ///
    xtitle("Relative time to year before TWEA") ///
    ytitle("Effect") ///
    xlabel(-18(3)21) ///
    yscale(range(-0.25 1)) ///
    ylabel(-0.25(.25)1)
graph export "${figures}/graphES_moser_controls.pdf", replace
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
// twoway (scatter res_ife2 res_ife1, msize(small) msymbol(o) mcolor(midblue) legend(order(2 /*"  &" 3*/ "IFE" 4 /*"  &" 5*/ "TWFE") pos(6) col(2))) ///
// 	(line res_ife2 res_ife1, lcolor(midblue)) (rcap res_ife4 res_ife3 res_ife1, lcolor(midblue)) ///
// 	(line res_post2 res_post1, lcolor(red)) (rcap res_post4 res_post3 res_post1, lcolor(red)) ///
// 	(scatter res_post2 res_post1, msize(small) msymbol(o) mcolor(red)), ///
// 	 title("Event-study estimates") xtitle("Years after TWEA") ///
// 	 ytitle("Effect") xscale(range(1 21)) xlabel(1(5)21) yscale(range(-0.5 1.5)) ylabel(-0.5(.25)1.5)
twoway ///
    (scatter res_ife2 res_ife1, msize(small) msymbol(o) mcolor(gs4) ///
        legend(order(2 "IFE" 4 "TWFE") pos(6) col(2))) ///
    (line  res_ife2 res_ife1, lcolor(gs4)) ///
    (rcap  res_ife4 res_ife3 res_ife1, lcolor(gs4)) ///
    (line  res_post2 res_post1, lcolor(gs10)) ///
    (rcap  res_post4 res_post3 res_post1, lcolor(gs10)) ///
    (scatter res_post2 res_post1, msize(small) msymbol(o) mcolor(gs10)), ///
    scheme(s1mono) ///
    title("Event-study estimates") xtitle("Years after TWEA") ///
    ytitle("Effect") xscale(range(1 21)) xlabel(1(5)21) ///
    yscale(range(-0.5 1.5)) ylabel(-0.5(.25)1.5)
graph export "${figures}/graphES_moser_ife_did.pdf", replace
restore


*** 3) Synthetic control

*Standard

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
graph export "C:\Users\134476\C DE CHAISEMARTIN Dropbox\clément de chaisemartin\A Mini course DID\Textbook\graphES_moser_sc_did.pdf", replace
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
// twoway (scatter res_sc2 res_sc1, msize(small) msymbol(o) mcolor(midblue) legend(order(2 /*"  &" 3*/ "Demeaned SC" 4 /*"  &" 5*/ "TWFE") pos(6) col(2))) ///
// 	(line res_sc2 res_sc1, lcolor(midblue)) (rcap res_sc4 res_sc3 res_sc1, lcolor(midblue)) ///
// 	(line res_post2 res_post1, lcolor(red)) (rcap res_post4 res_post3 res_post1, lcolor(red)) ///
// 	(scatter res_post2 res_post1, msize(small) msymbol(o) mcolor(red)), ///
// 	 title("Event-study estimates") xtitle("Years after TWEA") ///
// 	 ytitle("Effect") xscale(range(1 21)) xlabel(1(5)21) yscale(range(-2 2)) ylabel(-2(.5)2)
twoway ///
    (scatter res_sc2 res_sc1, msize(small) msymbol(o) mcolor(gs4) ///
        legend(order(2 "Demeaned SC" 4 "TWFE") pos(6) col(2))) ///
    (line  res_sc2   res_sc1, lcolor(gs4)) ///
    (rcap  res_sc4   res_sc3 res_sc1, lcolor(gs4)) ///
    (line  res_post2 res_post1, lcolor(gs10)) ///
    (rcap  res_post4 res_post3 res_post1, lcolor(gs10)) ///
    (scatter res_post2 res_post1, msize(small) msymbol(o) mcolor(gs10)), ///
    scheme(s1mono) ///
    title("Event-study estimates") xtitle("Years after TWEA") ///
    ytitle("Effect") xscale(range(1 21)) xlabel(1(5)21) ///
    yscale(range(-2 2)) ylabel(-2(.5)2)
graph export "${figures}/graphES_moser_demeaned_sc_did.pdf", replace
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
graph export "C:\Users\134476\C DE CHAISEMARTIN Dropbox\clément de chaisemartin\A Mini course DID\Textbook\graphES_moser_sd_did.pdf", replace
restore

*** 5) Sensitivity analysis of Rambachan and Roth 
// All black and white
clear all 
clear matrix
use "~/Dropbox/Mini course DID/Applications/Data sets/Moser and Voena 2012/moser_voena_didtextbook.dta", clear

set scheme plotplain
reghdfe patents reltimeminus* reltimeplus*, absorb(year treatmentgroup) cluster(subclass)
matrix l_vec = J(21,1,0)
matrix l_vec[14,1]=1
honestdid, pre(1) post(2) mvec(2(2)10) coefplot ///
    scheme(s1mono) ///
    plotopts(lcolor(black) mcolor(black)) ///
    ciopts(lcolor(black))
	
honestdid, pre(1) post(2) mvec(2(2)10) coefplot
	
	
	
graph display, scheme(plotplain)

* Edit graph to make all colors black
gr_edit plotregion1.plot1.style.editstyle line(color(black)) editcopy
gr_edit plotregion1.plot2.style.editstyle line(color(black)) editcopy

honestdid, pre(1) post(2) mvec(2(2)10) coefplot xtitle(M, size(large)) ytitle(95% Robust CI, size(large))
graph display, scheme(plotplain)
graph export "coefplot_bw.pdf", replace


honestdid, pre(1/18) post(19/39) mvec(0.5(0.5)2) l_vec(l_vec) coefplot
honestdid, pre(1/3) post(19) mvec(0.5(0.5)2)
honestdid, pre(1/7) post(19) mvec(0.5(0.5)2)

global figures "/Users/anzony.quisperojas/Documents/GitHub/did_book/figures"

preserve
keep if year==1918|year==1904|year==1932
reghdfe patents reltimeminus14 reltimeplus14 /// 
, absorb(year treatmentgroup) cluster(subclass)
honestdid, pre(1) post(2) mvec(2(2)10) coefplot xtitle(M, size(large)) ytitle(95% Robust CI, size(large)) ///
preserve
keep if year==1918|year==1904|year==1932
reghdfe patents reltimeminus14 reltimeplus14 /// 
, absorb(year treatmentgroup) cluster(subclass)
honestdid, pre(1) post(2) mvec(2(2)10) coefplot xtitle(M, size(large)) ytitle(95% Robust CI, size(large)) ///
	scheme(plotplain) ///
    plotopts(lcolor(black) mcolor(black)) ///
    ciopts(lcolor(black))
graph export "${figures}/coefplot.pdf", replace
restore
// graph export "C:\Users\134476.SCIENCESPO\C DE CHAISEMARTIN Dropbox\clément de chaisemartin\A Mini course DID\Textbook\coefplot.pdf", replace
restore
