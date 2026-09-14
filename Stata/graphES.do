cd "~/Documents/GitHub/did_book"

global figures "figures"
use "data/Divorce-Wolfers-AER.dta", clear 
set matsize 1000


/* Matrix saving the average effects from 0 to 7 years after the law, from the 
different methods */
matrix res_avg = J(5,2,0)

* Part 1 Wolfers' results in table 2, w/o linear time trends
xi i.years_unilateral i.st i.year 
quietly: reg div_rate _I* if year>1955 & year<1989 [w=stpop]

matrix temp = e(V)
matrix V_b = temp[1..4,1..4]
matrix e2 = J(4,1,1/4)

matrix temp=r(table)'
matrix res = ((0\1\2\3), temp[1..4,1], temp[1..4,5], temp[1..4,6])
matrix list res

matrix res_avg[1,1] = res[1..4,2]'*e2

matrix temp = e2'*(V_b*e2)
matrix res_avg[1,2] = sqrt(temp[1,1])

/* Part 2: TWFE reg., without gathering by groups of two years (as opposed to
 Wolfers) and with leads of the treatment */
gen Dur=year-lfdivlaw  if lfdivlaw<2000
replace Dur=min(15,max(Dur,-10)) if lfdivlaw<2000
replace Dur=-10 if lfdivlaw==2000

forvalues x = 0/15 {
	g Dt`x'=(Dur==`x')
}
* For placebo estimates
forvalues x = 2/10 {
	g Dt_`x'=(Dur==-`x')
}

drop _I*
xi i.st i.year 
quiet: reg div_rate Dt* _I* if year>1955 & year<1989 [w=stpop], vce(cluster st)

/* We produce the E-S graph by creating a matrix (res) gathering the time to the 
event, the point estimates and the CI */

matrix temp = e(V)
matrix V_b = temp[1..8,1..8]
matrix e = J(8,1,1/8)
matrix temp=r(table)'
matrix res=J(26,4,0)
matrix res[10,1]=0
forvalues x = 2/10 {
matrix res[11-`x',1]=-`x'+1
matrix res[11-`x',2]=temp[`x'+15,1]
matrix res[11-`x',3]=temp[`x'+15,5]
matrix res[11-`x',4]=temp[`x'+15,6]
}
forvalues x = 0/15 {
matrix res[`x'+11,1]=`x'+1
matrix res[`x'+11,2]=temp[`x'+1,1]
matrix res[`x'+11,3]=temp[`x'+1,5]
matrix res[`x'+11,4]=temp[`x'+1,6]
}

matrix res_avg[2,1] = res[11..18,2]'*e

matrix temp = e'*(V_b*e)
matrix res_avg[2,2] = sqrt(temp[1,1])

preserve
drop _all
svmat res
// twoway (scatter res2 res1, msize(medlarge) msymbol(o) mcolor(navy) legend(off)) ///
// 	(line res2 res1, lcolor(navy)) (rcap res4 res3 res1, lcolor(maroon)), ///
// 	 title("TWFE estimates") xtitle("Relative time to year before law") ///
// 	 ytitle("Effect") ylabel(-1(.5)0.5) xlabel(-9(3)16) yscale(range(-1.1 0.8)) name(g1)
twoway ///
    (scatter res2 res1, msize(medlarge) msymbol(o) mcolor(gs4) legend(off)) ///
    (line    res2 res1, lcolor(gs4)) ///
    (rcap    res4 res3 res1, lcolor(gs10)), ///
    scheme(s1mono) ///
    title("TWFE estimates") ///
    xtitle("Relative time to year before law") ///
    ytitle("Effect") ///
    xlabel(-9(3)16, grid glcolor(gs14) glwidth(vthin)) ///
    yscale(range(-1.1 0.8)) ///
    ylabel(-1(.5)0.5, grid glcolor(gs14) glwidth(vthin)) ///
    name(g1)
graph export "${figures}/graphES.pdf", replace
restore
