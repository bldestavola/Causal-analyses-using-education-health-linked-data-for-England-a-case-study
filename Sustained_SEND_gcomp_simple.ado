
*********************************************************************************************************************************
cap program drop Sustained_SEND|_gcomp_simple
cap program define Sustained_SEND_gcomp_simple, rclass

	
*****************************************************************************************************
*MONTE CARLO STEP
******************************************************************************************************


preserve
cap drop Y2_* H2_* Y3_* Y4_*

expand 1000

sort id t
by id t:gen original=_n==1
by id t:gen counter=_n

**************************************************************
*ESTIMATED Potential Outcomes at (end of time 1) time 2
**************************************************************

*unexposed
poisson Y  male white idacibin eyfsp_bin Hosp_lag1  if SEN_lag1==0 & t==2 & original , e(den) nolog
predict Y2_0 
*exposed
  poisson Y  male white idacibin eyfsp_bin Hosp_lag1  if SEN_lag1==1 & t==2 & original,e(den) nolog
predict Y2_1 

order counter original id t
sort  counter id t 
qui by counter  id: replace Y2_0=Y2_0[2]
qui by counter  id: replace Y2_1=Y2_1[2]

**************************************************************
*potential confounder at time 2
**************************************************************
logit Hosp male white idacibin eyfsp Hosp_lag1  if t==2  & original, nolog

*unexposed
qui gen p2_0= _b[_cons] + _b[male] * male+ _b[white] * white + _b[idacibin] * (idacibin==1) + _b[eyfsp]*eyfsp+ _b[Hosp_lag1] * Hosp_lag1 
qui gen H2_0= runiform() < 1/(1+ exp(-p2_0)) 

*exposed
qui gen p2_1= _b[_cons] + _b[male] * male+ _b[white] * white + _b[idacibin] * (idacibin==1) + _b[eyfsp]*eyfsp+ _b[Hosp_lag1] * Hosp_lag1 
qui gen H2_1= runiform() < 1/(1+ exp(-p2_1)) 

sort counter  id t
qui by counter  id: replace H2_0=H2_0[2]
qui by counter  id: replace H2_1=H2_1[2]
drop p2_*


**************************************************************
*potential outcomes at time 3
**************************************************************
 poisson Y  male white idacibin eyfsp_bin Hosp_lag1 Hosp_lag2 SEN_lag1  if  t==3  & original, e(den)nolog

cap drop m3_*
*unexposed
qui gen m3_00= exp(_b[_cons] + _b[male] * male + _b[white] * white + _b[idacibin] * (idacibin==1) + _b[eyfsp_bin]*eyfsp_bin +  _b[Hosp_lag1] * H2_0  + _b[Hosp_lag2] * Hosp_lag2 +_b[SEN_lag1] * 0  )   *den

*exposed
qui gen m3_11= exp(_b[_cons] + _b[male] * male + _b[white] * white + _b[idacibin] * (idacibin==1) + _b[eyfsp_bin]*eyfsp_bin +  _b[Hosp_lag1] * H2_1  + _b[Hosp_lag2] * Hosp_lag2 +_b[SEN_lag1] * 1 )  *den

sort counter  id t
qui by counter  id: replace m3_00=m3_00[3]
qui by counter  id: replace m3_11=m3_11[3]

qui gen Y3_00=rpoisson(m3_00)
qui gen Y3_11=rpoisson(m3_11)
drop m3_*

sort counter  id t
qui by counter  id: replace Y3_00=Y3_00[3]
qui by counter  id: replace Y3_11=Y3_11[3]


**************************************************************
*potential confounder at time 3
**************************************************************
logit Hosp male white idacibin Hosp_lag1 if t==3  & original, nolog

*unexposed
qui gen p3_00= _b[_cons] + _b[male] * male+ _b[white] * white + _b[idacibin] * (idacibin==1) + _b[Hosp_lag1] * H2_0 
qui gen H3_00= runiform() < 1/(1+ exp(-p3_00)) 

*exposed
qui gen p3_11= _b[_cons] + _b[male] * male+ _b[white] * white + _b[idacibin] * (idacibin==1) + _b[Hosp_lag1] * H2_1 
qui gen H3_11= runiform() < 1/(1+ exp(-p3_11)) 

sort counter  id t
qui by counter  id: replace H3_00=H3_00[3]
qui by counter  id: replace H3_11=H3_11[3]
drop p3_*


**************************************************************
*potential outcomes at time 4
**************************************************************
poisson Y  male  idacibin eyfsp_bin Hosp_lag1  SEN_lag1  if  t==4  & original, e(den) nolog

cap drop m4_*
*unexposed
qui gen m4_000= exp(_b[_cons] + _b[male] * male  + _b[idacibin] * (idacibin==1) + _b[eyfsp_bin]*eyfsp_bin +  _b[Hosp_lag1] * H3_00   +_b[SEN_lag1] * 0 )  *den

*exposed
qui gen m4_111= exp(_b[_cons] + _b[male] * male  + _b[idacibin] * (idacibin==1) + _b[eyfsp_bin]*eyfsp_bin +  _b[Hosp_lag1] * H3_11   +_b[SEN_lag1] * 1  )  *den

sort counter  id t
qui by counter  id: replace m4_000=m4_000[4]
qui by counter  id: replace m4_111=m4_111[4]
su m4_* if t==4

qui gen Y4_000=rpoisson(m4_000)
qui gen Y4_111=rpoisson(m4_111)
drop m4_*

sort counter  id t
qui by counter  id: replace Y4_000=Y4_000[4]
qui by counter  id: replace Y4_111=Y4_111[4]


**************************************************************
*Estimands
**************************************************************
*difference in mean rates and rate ratios: 
cap drop rate*
qui gen rate_111=(Y2_1+Y3_11+Y4_111)/cumden4
qui gen rate_000=(Y2_0+Y3_00+Y4_000)/cumden4

*ATE
gen ATE_RD=rate_111-rate_000

qui su rate_111 if t==4
scalar mean1=r(mean)
qui su rate_000 if t==4
scalar mean0=r(mean)  
gen ATE_RR=mean1/mean0


*POST THE RESULTS
	su rate_000 if t==4
   return scalar rate_000=r(mean)
   
	su rate_111 if t==4
   return scalar rate_111=r(mean)
   
    su ATE_RD if t==4
	return scalar ate_RD=r(mean)
	
    su ATE_RR if t==4
	return scalar ate_RR=r(mean)



restore

end


