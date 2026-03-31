
*********************************************************************************************************************************
cap program drop Sustained_IPW_stab
cap program define Sustained_IPW_stab, rclass
args A Y den newid

preserve
cap drop p_num*
cap drop p_den*
cap drop sw*
cap drop cumA*
cap drop A_lag
cap drop cumY
cap drop cumden


sort `newid' t
qui  by `newid'(t): gen A_lag = `A'[_n-1]

	
* -------------------------------------------------------------
* 1. PS model
* -------------------------------------------------------------

*t=1
foreach t of numlist 1{
	
	logit `A'                                                      if t==`t', nolog
	predict double p_num`t'                                        if t==`t'
	
	logit `A'   $base_conf                                         if t==`t', nolog
	predict double p_den`t'                                        if t==`t'
	
*	twoway (kdensity p_den`t' if SEN==1 & t==`t', lcol(red)) (kdensity p_den`t' if SEN==0 & t==`t',lcol(blue)) , name(PS`t', replace) xtitle(PS) title("PS for time `t'")
	replace p_num`t'=1-p_num`t'                                    if t==`t' & `A'==0
	replace p_den`t'=1-p_den`t'                                    if t==`t' & `A'==0

}


*t>1
foreach t of numlist 2/3{
	
	logit `A' A_lag i.t                                            if t==`t', nolog
	predict double p_num`t'                                        if t==`t'
	
	logit `A' A_lag  $tv_conf $base_conf                        if t==`t', nolog
	predict double p_den`t'                                        if t==`t'
	
*	twoway (kdensity p_den`t' if SEN==1 & t==`t', lcol(red)) (kdensity p_den`t' if SEN==0 & t==`t',lcol(blue)) , name(PS`t', replace) xtitle(PS) title("PS for time `t'")
	replace p_num`t'=1-p_num`t'                                    if t==`t' & `A'==0
	replace p_den`t'=1-p_den`t'                                    if t==`t' & `A'==0

}

gen p_num=.
gen p_den=.
foreach t of numlist 1/3{
	replace p_num=p_num`t' if t==`t' 
	replace p_den=p_den`t' if t==`t' 
}

* -------------------------------------------------------------
* 2. Build stabilized weight per row, then cumulative product per newid
* -------------------------------------------------------------

* Per-time stabilized weight factor
gen double lsw_t = log(p_num) - log(p_den)

tabstat  p_num p_den lsw_t , by(t) s(count mean min max) c(s)

* Cumulative product across time for each subject (product of sw_t up to t on log scale to improve stability)
bys `newid' (t): gen double lsw = lsw_t
bys `newid': replace lsw = lsw[_n-1] + lsw_t if _n>1
gen double sw = exp(lsw)

* -------------------------------------------------------------
* 3. Diagnostics: examine sw_final distribution, truncate if needed
* -------------------------------------------------------------
/* Inspect percentiles. If extreme weights, truncate at, e.g., 1st/99th or 0.01/99.99 */
centile sw, centile(1 99)
local p1 = r(c_1)
local p99 = r(c_2)
display "1st pctile = `p1', 99th pctile = `p99'"
*/

* Truncate example
gen double sw_trim = sw if sw<.
replace sw_trim = `p99' if sw_trim > `p99' & sw_trim<.
replace sw_trim = `p1' if sw_trim < `p1' 
centile sw_trim, centile(1 99)

* ---------------------------------------------------------------------
* 4. Create cumulative exposure and cumulative outcome and denominator
* ---------------------------------------------------------------------
ta `A'
bys `newid' (t): gen cumA = `A'
bys `newid': replace cumA = cumA[_n-1] + `A' if _n>1
bys `newid' (t): gen cumA_lag = cumA[_n-1]
ta t cumA_lag

bys `newid' (t): gen sw_lag = sw[_n-1]
bys `newid' (t): gen sw_trim_lag = sw_trim[_n-1]
tabstat sw_trim sw_trim_lag, by(t) s(count min max)


bys `newid' (t): gen cumY = `Y' if _n>1
bys `newid'  (t): replace cumY = cumY[_n-1] + `Y' if _n>2

bys `newid' (t): gen cumden = `den' if _n>1
bys `newid'  (t): replace cumden = cumden[_n-1] + `den' if _n>2

*first record not needed anymore- note Y starts from Year 2 
drop if t==1

* -------------------------------------------------------------
* 6. Fit the Marginal Structural Model (weighted)
* -------------------------------------------------------------
* Simple MSM: 
poisson cumY cumA_lag i.t  [pw=sw_trim_lag], e(cumden) vce(cluster `newid') nolog


*more general MSM
poisson cumY c.cumA_lag##c.cumA_lag A_lag i.t  [pw=sw_trim_lag], e(cumden) vce(cluster `newid') nolog

*results for Year 4
	di "base rate at t=4"
	nlcom exp(_b[_cons]+_b[4.t])
	scalar rate_000=exp(_b[_cons]+_b[4.t])
	return scalar  rate_000=rate_000
	
	di "exposed rates at t=`t'"
	nlcom exp(_b[_cons]+_b[4.t]+_b[cumA_lag]*3+_b[c.cumA_lag#c.cumA_lag]*9+_b[A_lag])
	scalar rate_111=exp(_b[_cons]+_b[4.t]+_b[cumA_lag]*3+_b[c.cumA_lag#c.cumA_lag]*9+_b[A_lag])
	return scalar  rate_111=rate_111
	
	di "RR at t=4"
	nlcom exp(_b[cumA_lag]*3+_b[c.cumA_lag#c.cumA_lag]*9+_b[A_lag])
    scalar  ATE_RR=exp(_b[cumA_lag]*3+_b[c.cumA_lag#c.cumA_lag]*9+_b[A_lag])
    return scalar ATE_RR=ATE_RR
	
	di "RD at t=4"
	nlcom exp(_b[_cons]+_b[4.t]+_b[cumA_lag]*3+_b[c.cumA_lag#c.cumA_lag]*9+_b[A_lag])-exp(_b[_cons]+_b[4.t])
	scalar ATE_RD=exp(_b[_cons]+_b[4.t]+_b[cumA_lag]*3+_b[c.cumA_lag#c.cumA_lag]*9+_b[A_lag])-exp(_b[_cons]+_b[4.t])
	return scalar ATE_RD=ATE_RD
	
	scalar list
	

restore

end

