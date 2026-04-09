capture use "C:\Users\134476\C DE CHAISEMARTIN Dropbox\clément de chaisemartin\A Mini course DID\Applications\Data sets\Pierce and Schott 2016\pierce_schott_didtextbook.dta", clear

if _rc==601{
use "C:\Users\134476.SCIENCESPO\C DE CHAISEMARTIN Dropbox\clément de chaisemartin\A Mini course DID\Applications\Data sets\Pierce and Schott 2016\pierce_schott_didtextbook.dta", clear
}

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



