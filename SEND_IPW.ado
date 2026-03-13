
************************************************************************************************************************
cap program drop SEND_IPW
cap program define SEND_IPW, rclass

args  Y A D H
	

cap drop num 
cap drop den 
cap drop wt

*PS for receiving SEND 
*numerator
logit `A' , nolog
predict num
*denominator
logit `A' $base_conf `H', nolog
predict den
*twoway (kdensity den if `A'==1, lcol(red)) (kdensity den if `A'==0,lcol(blue)), name(byhand, replace) xtitle(PS)

*ATE
gen wt=cond(`A', num/den, (1-num)/(1-den))
*Poisson 
poisson `Y' `A' [pw=wt], exp(`D') nolog
*ATE,ratio, rate
return scalar ATE_rr=exp(_b[`A'])
*ATE,dif,  rate
return scalar ATE_rd=exp(_b[_cons]+_b[`A']) -exp(_b[_cons])

*ATT
cap drop wt
gen wt=cond(`A', 1, den/(1-den))
poisson `Y' `A' [pw=wt], exp(`D') nolog
*ATE,ratio, rate
return scalar ATT_rr=exp(_b[`A'])
*ATE,dif,  rate
return scalar ATT_rd=exp(_b[_cons]+_b[`A']) -exp(_b[_cons])

end



