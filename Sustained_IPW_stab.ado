
*********************************************************************************************************************************
cap program drop Sustained_IPW_stab
cap program define Sustained_IPW_stab, rclass
args A1 A2 A3 A4 Y1 Y2 Y3 Y4 den1 den2 den3 den4


preserve

reshape long `Y' inter  Hosp `A' `den' SEN_lag1 Hosp_lag1 Y_lag1 SEN_lag2 Hosp_lag2 SEN_lag3 Hosp_lag3 cumden cumY
ta `A'
sort id t
qui  by id: gen A_lag = `A'[_n-1]


	
* -------------------------------------------------------------
* 1. PS model
* -------------------------------------------------------------

*t=1
foreach t of numlist 1{
	
	logit `A'                                                      if t==`t', nolog
	predict double p_num`t'                                        if t==`t'
	
	logit `A'   $base_conf                                         if t==`t', nolog
	predict double p_den`t'                                        if t==`t'
	
	replace p_num`t'=1-p_num`t'                                    if t==`t' & `A'==0
	replace p_den`t'=1-p_den`t'                                    if t==`t' & `A'==0

}


*t>1
foreach t of numlist 2/3{
	
	logit `A' A_lag                                             if t==`t', nolog
	predict double p_num`t'                                        if t==`t'
	
	logit `A' A_lag  $tv_conf $base_conf                        if t==`t', nolog
	predict double p_den`t'                                        if t==`t'
	

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
* 2. Build stabilized weight per row, then cumulative product per id
* -------------------------------------------------------------

* Per-time stabilized weight factor
gen double sw_t = p_num / p_den
*tabstat  p_num p_den sw_t , by(t) s(count mean min max) c(s)

* Cumulative product across time for each subject (product of sw_t up to t)
bys id (t): gen double sw = sw_t
bys id: replace sw = sw[_n-1] * sw_t if _n>1
*tabstat   sw , by(`A') s(count mean min max) c(s)


* ---------------------------------------------------------------------
* 3. Create cumulative exposure and cumulative outcome and denominator
* ---------------------------------------------------------------------
cap drop cumA
cap drop cumY
cap drop cumden

bys id (t): gen cumA = `A'
bys id: replace cumA = cumA[_n-1] + `A' if _n>1
bys id: gen cumA_lag = cumA[_n-1]
bys id: gen sw_lag = sw[_n-1]



bys id (t): gen cumY = `Y' if _n>1
bys id: replace cumY = cumY[_n-1] + `Y' if _n>2

bys id (t): gen cumden = `den' if _n>1
bys id: replace cumden = cumden[_n-1] + `den' if _n>2


drop if t==1

* -------------------------------------------------------------
* 5. Fit the Marginal Structural Model (weighted)
* -------------------------------------------------------------
* Simple MSM: 
poisson cumY cumA_lag i.t  [pw=sw_lag], e(cumden) vce(cluster id) nolog


*more general MSM
poisson cumY c.cumA_lag##c.cumA_lag A_lag i.t [pw=sw_lag], e(cumden) vce(cluster id) nolog

*results
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

